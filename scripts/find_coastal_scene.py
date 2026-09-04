"""Bounded, quality-driven search; no selection against NMF/chlorophyll fit."""
import argparse
import getpass
import json
import os
from pathlib import Path
import numpy as np
from sgsnmf_py.io.loaders import search, open_granule, prepare, attach_chlorophyll
from sgsnmf_py.core.sgsnmf import fit
from sgsnmf_py.products import save_products


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/caribbean")
    parser.add_argument("--bbox", nargs=4, type=float, default=[-74.5, 11.5, -71, 13.5])
    parser.add_argument("--dates", nargs=2, default=["2024-04-20", "2024-05-10"])
    parser.add_argument("--count", type=int, default=12)
    args = parser.parse_args()
    os.environ["EARTHDATA_USERNAME"] = input("EARTHDATA_USERNAME: ").strip()
    os.environ["EARTHDATA_PASSWORD"] = getpass.getpass("EARTHDATA_PASSWORD (oculta): ")
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    records = []
    try:
        found = search(args.bbox, args.dates, count=args.count)
        names = sorted({g["umm"]["GranuleUR"] for g in found})
        for name in names:
            print("Inspecting", name, flush=True)
            try:
                with open_granule(name) as (aop, metadata):
                    try:
                        initial = prepare(aop, bbox=args.bbox, max_side=96, metadata=metadata)
                        print("AOP valid pixels:", int(initial.mask.sum()), flush=True)
                        bgc_name = name.replace(".OC_AOP.", ".OC_BGC.")
                        with open_granule(bgc_name, "PACE_OCI_L2_BGC") as (bgc, _):
                            try:
                                combined = attach_chlorophyll(aop, bgc)
                                scene = prepare(combined, bbox=args.bbox, max_side=96,
                                                metadata=dict(metadata, chlorophyll_scene=bgc_name))
                            finally:
                                bgc.close()
                    finally:
                        aop.close()
                pairs = int(np.sum(scene.mask & np.isfinite(scene.chlorophyll) & (scene.chlorophyll > 0)))
                print("Valid chlorophyll pairs:", pairs, flush=True)
                if pairs < 20:
                    records.append({"scene": name, "status": "too_few_chlorophyll_pairs", "pairs": pairs})
                    continue
                result = fit(scene.cube, mask=scene.mask, p=4, n_segments=25, max_iter=100, seed=17, time_limit=1e9)
                save_products(result, scene, output, parameters={"p": 4, "segments": 25, "seed": 17, "max_iter": 100,
                                                               "selection": "first CMR granule in chronological order with >=16 valid spectra and >=20 chlor_a pairs"})
                # Small, explicit analysis snapshot enables credential-free replay.
                np.savez_compressed(output / "input_snapshot.npz", cube=scene.cube, wavelength=scene.wavelength,
                                    latitude=scene.latitude, longitude=scene.longitude, mask=scene.mask,
                                    scale=scene.scale, chlorophyll=scene.chlorophyll,
                                    metadata=json.dumps(scene.metadata))
                records.append({"scene": name, "status": "selected", "pairs": pairs})
                (output / "scene_selection.json").write_text(json.dumps(records, indent=2))
                print("SUCCESS", output.resolve(), flush=True)
                return
            except (ValueError, RuntimeError) as error:
                record = {"scene": name, "status": "rejected", "reason": str(error)}
                records.append(record)
                print(record, flush=True)
            (output / "scene_selection.json").write_text(json.dumps(records, indent=2))
        raise RuntimeError("No qualifying scene in the bounded search; see scene_selection.json")
    finally:
        os.environ.pop("EARTHDATA_USERNAME", None)
        os.environ.pop("EARTHDATA_PASSWORD", None)


if __name__ == "__main__":
    run()
