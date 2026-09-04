"""Interactive credential entry; secrets live only in this process.

Run from a user-controlled terminal. Production/capsule uses environment
secrets and the noninteractive CLI instead.
"""
import argparse
import getpass
import os
from pathlib import Path
from sgsnmf_py.cli import main


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/caribbean")
    parser.add_argument("--scene", default="PACE_OCI.20240427T174513.L2.OC_AOP.V3_2.nc")
    parser.add_argument("--bbox", nargs=4, default=["-74.5", "11.5", "-71", "13.5"])
    args = parser.parse_args()
    print("NASA Earthdata: introduce tus credenciales en esta terminal, no en el chat.", flush=True)
    username = input("EARTHDATA_USERNAME: ").strip()
    password = getpass.getpass("EARTHDATA_PASSWORD (oculta): ")
    if not username or not password:
        raise ValueError("Se requieren usuario y contraseña")
    os.environ["EARTHDATA_USERNAME"] = username
    os.environ["EARTHDATA_PASSWORD"] = password
    try:
        main(["run", "--scene", args.scene,
              "--collection", "PACE_OCI_L2_AOP", "--bbox", *args.bbox,
              "--p", "4", "--segments", "25", "--max-side", "96",
              "--chlorophyll-scene", args.scene.replace(".OC_AOP.", ".OC_BGC."),
              "--output", str(Path(args.output).resolve())])
    finally:
        os.environ.pop("EARTHDATA_USERNAME", None)
        os.environ.pop("EARTHDATA_PASSWORD", None)
        password = None


if __name__ == "__main__":
    run()
