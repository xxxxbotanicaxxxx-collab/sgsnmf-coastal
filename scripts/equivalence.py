"""Run upstream Octave on the original full demo; replay exact inputs in Python.

Only upstream plotting is suppressed in a temporary copy. No optimized
equations or NNLS implementation are replaced. Optional reference MAT reuse.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import time
import numpy as np
from scipy.io import loadmat
from sgsnmf_py.core.sgsnmf import sgsnmf
from sgsnmf_py.core.slic_hsi import Segmentation
from sgsnmf_py.core.metrics import compare


UPSTREAM_COMMIT = "e124bc76499bbbeebe10001b410a11b1d4945966"


def quote(path):
    return str(Path(path).resolve()).replace("\\", "/").replace("'", "''")


def run(upstream, output, octave="octave-cli", max_iter=30, reuse=False):
    upstream, output = Path(upstream).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    fixture = output / "octave_reference.mat"
    if not reuse:
        revision = subprocess.check_output(["git", "-C", str(upstream), "rev-parse", "HEAD"], text=True).strip()
        if revision != UPSTREAM_COMMIT:
            raise ValueError(f"Expected upstream {UPSTREAM_COMMIT}; got {revision}")
        source = (upstream / "code/Func/slic_HSI.m").read_text(errors="replace")
        begin = source.index("   % Plot results")
        end = source.index("\nend", begin)
        source = source[:begin] + "   % Plotting suppressed by validation harness.\n" + source[end:]
        with tempfile.TemporaryDirectory() as directory:
            temp = Path(directory)
            (temp / "slic_HSI.m").write_text(source)
            script = f"""
warning('off','Octave:possible-matlab-short-circuit-operator');
pkg load image;
addpath(genpath('{quote(upstream / 'code')}'));
addpath('{quote(temp)}','-begin');
rand('seed',17); randn('seed',17);
disp('STAGE load-demo'); fflush(stdout);
load('{quote(upstream / 'data/F1_A_9.mat')}');
load('{quote(upstream / 'data/F1_S_9.mat')}');
[rows,cols,M]=size(S); L=size(A,1); N=rows*cols;
Htrue=reshape(S,N,M)';
Y=reshape((A*Htrue)',rows,cols,L);
[X,noise,Cn]=addNoise(Y,'additive',20,0,0);
X=max(X,eps); Y=reshape(X',rows,cols,L);
disp('STAGE original-SLIC'); fflush(stdout);
seg=slic_HSI(Y,round(N/36),0.5);
disp('STAGE original-VCA-FCLS'); fflush(stdout);
[W0,unused]=hyperVca(seg.X_c,M);
H0=fcls(W0,X);
para=struct('X',X,'W',W0,'H',H0,'M',M,'lambda',0.3,...
 'tol',0.05,'maxiter',{max_iter},'timelimit',1e9,'verbose',1,'print_iter',5);
save('-mat7-binary','{quote(output / 'octave_inputs.mat')}','X','Y','A','Htrue','W0','H0','seg','para');
disp('STAGE original-optimizer'); fflush(stdout);
[W,H]=sgsnmf(para,seg);
save('-mat7-binary','{quote(fixture)}','X','Y','A','Htrue','W0','H0','seg','W','H','para');
"""
            (temp / "reference.m").write_text(script)
            started = time.monotonic()
            with (output / "octave.log").open("w") as logfile:
                process = subprocess.run([octave, "--no-gui", "--quiet", str(temp / "reference.m")],
                                         stdout=logfile, stderr=subprocess.STDOUT, text=True, timeout=7200)
            process.check_returncode()
            print(f"Octave completed in {time.monotonic() - started:.1f} s", flush=True)
    data = loadmat(fixture, simplify_cells=True)
    s = data["seg"]
    labels = np.asarray(s["labels"], int).ravel() - 1
    seg = Segmentation(labels, s["X_c"], np.asarray(s["Cj"]).ravel(), float(s["Sw"]),
                       labels.reshape(data["Y"].shape[:2], order="F"))
    result = sgsnmf(data["X"], data["W0"], data["H0"], seg,
                    max_iter=int(data["para"]["maxiter"]), time_limit=1e9)
    metrics = compare(data["W"], data["H"], result.endmembers, result.abundances)
    discrepancies = {k: v for k, v in metrics.items() if k != "permutation" and np.any(np.asarray(v) > 1e-3)}
    report = {"upstream_commit": UPSTREAM_COMMIT, "fixture_sha256": hashlib.sha256(fixture.read_bytes()).hexdigest(),
              "scope": "Full original 20 dB demo; shared Octave SLIC, noise, VCA and FCLS initialization",
              "iterations_requested": int(data["para"]["maxiter"]),
              "python_iterations": result.iterations, "python_stop": result.stop_reason,
              "python_vs_octave": metrics, "discrepancies_over_1e-3": discrepancies,
              "octave_vs_truth": compare(data["A"], data["Htrue"], data["W"], data["H"]),
              "python_vs_truth": compare(data["A"], data["Htrue"], result.endmembers, result.abundances),
              "pass": not discrepancies}
    (output / "equivalence.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--upstream", required=True)
    parser.add_argument("--output", default="results/equivalence")
    parser.add_argument("--octave", default="octave-cli")
    parser.add_argument("--max-iter", type=int, default=30)
    parser.add_argument("--reuse", action="store_true")
    args = parser.parse_args()
    report = run(args.upstream, args.output, args.octave, args.max_iter, args.reuse)
    raise SystemExit(0 if report["pass"] else 1)
