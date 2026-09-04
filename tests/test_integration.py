import importlib.util
import os
from pathlib import Path
import numpy as np
import pytest
from sgsnmf_py.io.loaders import public_file, read, prepare
from sgsnmf_py.core.sgsnmf import fit


@pytest.mark.integration
def test_fixed_public_pace(tmp_path):
    if os.environ.get("RUN_NETWORK_TESTS") != "1":
        pytest.skip("Set RUN_NETWORK_TESTS=1 to run real PACE public smoke test")
    path = Path(os.environ["PACE_FIXTURE"]) if "PACE_FIXTURE" in os.environ else public_file(tmp_path)
    # File-like interface also exercises the form returned by earthaccess.open.
    with open(path, "rb") as handle:
        ds = read(handle)
        try:
            scene = prepare(ds, bbox=(-91, 27, -88, 29), max_side=48)
            result = fit(scene.cube, mask=scene.mask, p=4, n_segments=10, max_iter=2, seed=17)
            assert scene.mask.sum() >= 50
            assert result.endmembers.shape == (len(scene.wavelength), 4)
            assert np.isfinite(result.abundances).all()
        finally:
            ds.close()


@pytest.mark.octave
def test_octave_original_demo(tmp_path):
    if "SGSNMF_UPSTREAM" not in os.environ or "OCTAVE" not in os.environ:
        pytest.skip("Set SGSNMF_UPSTREAM and OCTAVE to run original MATLAB code in Octave")
    script = Path(__file__).parents[1] / "scripts/equivalence.py"
    spec = importlib.util.spec_from_file_location("equivalence", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    report = module.run(os.environ["SGSNMF_UPSTREAM"], tmp_path, os.environ["OCTAVE"], max_iter=5)
    assert report["pass"], report["discrepancies_over_1e-3"]
