"""Docker/Figshare replay entry point. Never persist or print credentials."""
import os
from pathlib import Path
from sgsnmf_py.cli import main


def run():
    output = Path(os.environ.get("RESULTS_DIR", "/results"))
    scene = os.environ.get("SGSNMF_SCENE")
    arguments = ["run", "--output", str(output), "--p", "4", "--seed", "17",
                 "--segments", "25", "--max-side", "96", "--max-iter", "100"]
    snapshot = os.environ.get("SGSNMF_SNAPSHOT")
    bundled = Path(__file__).resolve().parents[1] / "data/coastal_snapshot.npz"
    if snapshot or (not scene and bundled.exists()):
        arguments += ["--snapshot", snapshot or str(bundled)]
    elif scene:
        arguments += ["--scene", scene, "--collection", os.environ.get("SGSNMF_COLLECTION", "PACE_OCI_L2_AOP")]
        bbox = os.environ.get("SGSNMF_BBOX", "-77.5,7.5,-74,12").split(",")
        if len(bbox) != 4:
            raise ValueError("SGSNMF_BBOX must be west,south,east,north")
        arguments += ["--bbox", *bbox]
        if os.environ.get("SGSNMF_CHL_SCENE"):
            arguments += ["--chlorophyll-scene", os.environ["SGSNMF_CHL_SCENE"]]
    elif os.environ.get("SGSNMF_PUBLIC_FILE"):
        arguments += ["--file", os.environ["SGSNMF_PUBLIC_FILE"], "--bbox", "-91", "27", "-88", "29"]
    else:
        arguments += ["--public-demo", "--cache", os.environ.get("SGSNMF_CACHE", "/data/cache")]
    main(arguments)


if __name__ == "__main__":
    run()
