"""Command line entry point."""
import argparse
import json
from pathlib import Path
from .core.sgsnmf import fit
from .io.loaders import read, prepare, open_granule, search, public_file, sha256, attach_chlorophyll, load_snapshot, PUBLIC_SCENE, PACE_COLLECTION
from .domain import library_initialization
from .products import save_products


def main(argv=None):
    parser = argparse.ArgumentParser(description="SGSNMF coastal: exploratory hyperspectral unmixing")
    subs = parser.add_subparsers(dest="command", required=True)
    find = subs.add_parser("search")
    find.add_argument("--bbox", nargs=4, type=float, required=True)
    find.add_argument("--dates", nargs=2, required=True)
    find.add_argument("--collection", default=PACE_COLLECTION)
    find.add_argument("--count", type=int, default=5)
    run = subs.add_parser("run")
    src = run.add_mutually_exclusive_group(required=True)
    src.add_argument("--scene", help="Exact NASA granule id")
    src.add_argument("--file", type=Path)
    src.add_argument("--snapshot", type=Path, help="Replay archived analysis subset without NASA credentials")
    src.add_argument("--public-demo", action="store_true", help="Fixed 2024 public HyperCoast PACE fixture")
    run.add_argument("--mission", choices=["pace", "emit"], default="pace")
    run.add_argument("--collection", default=PACE_COLLECTION)
    run.add_argument("--bbox", nargs=4, type=float)
    run.add_argument("--p", type=int, choices=[3, 4], default=4)
    run.add_argument("--segments", type=int, default=40)
    run.add_argument("--compactness", type=float, default=0.5)
    run.add_argument("--lambda", dest="lam", type=float, default=0.3)
    run.add_argument("--max-iter", type=int, default=100)
    run.add_argument("--max-side", type=int, default=96)
    run.add_argument("--seed", type=int, default=17)
    run.add_argument("--library", type=Path, help="Measured Rrs domain library CSV; no extrapolation")
    chl = run.add_mutually_exclusive_group()
    chl.add_argument("--chlorophyll-file", type=Path, help="Standard PACE OC_BGC product, same acquisition and processing")
    chl.add_argument("--chlorophyll-scene", help="Exact PACE OC_BGC NASA granule; geolocation must match")
    run.add_argument("--chlorophyll-collection", default="PACE_OCI_L2_BGC")
    run.add_argument("--allow-missing-flags", action="store_true")
    run.add_argument("--output", type=Path, default=Path("results/coastal"))
    run.add_argument("--cache", type=Path, default=Path(".cache"))
    args = parser.parse_args(argv)
    if args.command == "search":
        results = search(args.bbox, args.dates, args.collection, args.count)
        print(json.dumps([{"scene": g["umm"]["GranuleUR"], "links": g.data_links()} for g in results], indent=2))
        return
    def process_with_chlorophyll(ds, metadata):
        if args.chlorophyll_scene:
            with open_granule(args.chlorophyll_scene, args.chlorophyll_collection, "pace") as (bgc, _):
                try:
                    merged = attach_chlorophyll(ds, bgc)
                    process(merged, dict(metadata, chlorophyll_scene=args.chlorophyll_scene))
                finally:
                    bgc.close()
        elif args.chlorophyll_file:
            bgc = read(args.chlorophyll_file)
            try:
                process(attach_chlorophyll(ds, bgc), dict(metadata, chlorophyll_file=args.chlorophyll_file.name,
                                                       chlorophyll_sha256=sha256(args.chlorophyll_file)))
            finally:
                bgc.close()
        else:
            process(ds, metadata)
    def process(ds, metadata):
        bbox = args.bbox
        if args.public_demo and bbox is None:
            bbox = [-91, 27, -88, 29]
        scene = prepare(ds, bbox=bbox, max_side=args.max_side, metadata=metadata,
                        allow_missing_flags=args.allow_missing_flags)
        process_scene(scene)
    def process_scene(scene):
        library = library_initialization(args.library, scene.wavelength, args.p) if args.library else None
        result = fit(scene.cube, p=args.p, mask=scene.mask, n_segments=args.segments,
                     compactness=args.compactness, lam=args.lam, max_iter=args.max_iter,
                     seed=args.seed, time_limit=1e9, w0=library / scene.scale if library is not None else None)
        parameters = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}
        report = save_products(result, scene, args.output, library, parameters)
        print(json.dumps({"output": str(args.output.resolve()), "valid_pixels": scene.metadata["valid_pixels"],
                          "iterations": report["iterations"], "chlorophyll_check": report["chlorophyll_sanity_check"]}, indent=2))
    if args.snapshot:
        if args.bbox or args.chlorophyll_file or args.chlorophyll_scene:
            parser.error("Snapshot already contains a fixed crop and chlorophyll reference")
        process_scene(load_snapshot(args.snapshot))
    elif args.scene:
        with open_granule(args.scene, args.collection, args.mission) as (ds, metadata):
            try:
                process_with_chlorophyll(ds, metadata)
            finally:
                ds.close()
    else:
        path = public_file(args.cache) if args.public_demo else args.file
        ds = read(path, args.mission)
        try:
            process_with_chlorophyll(ds, {"scene": PUBLIC_SCENE if args.public_demo else path.name,
                         "sha256": sha256(path), "access": "automatic public cache" if args.public_demo else "local"})
        finally:
            ds.close()


if __name__ == "__main__":
    main()
