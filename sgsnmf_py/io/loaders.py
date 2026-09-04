"""Mission loaders, geolocation-aware cropping and conservative quality masks."""
from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import urllib.request
import numpy as np

PUBLIC_SCENE = "PACE_OCI.20240423T184658.L2.OC_AOP.V1_0_0.NRT.nc"
PUBLIC_URL = "https://github.com/opengeos/datasets/releases/download/netcdf/" + PUBLIC_SCENE
PUBLIC_SHA256 = "1c7aca086fa4681c961c59f0ad7b28a1be1909d9ece793f258b5aa6bb6c83311"
PACE_COLLECTION = "PACE_OCI_L2_AOP"
REJECT_FLAGS = {"ATMFAIL", "LAND", "HIGLINT", "HILT", "HISATZEN", "STRAYLIGHT", "CLDICE", "COCCOLITH", "HISOLZEN", "LOWLW", "CHLFAIL", "NAVWARN", "MAXAERITER", "ATMWARN", "NAVFAIL", "FILTER", "PRODFAIL", "SEAICE", "BOWTIEDEL", "HIPOL"}


@dataclass
class Scene:
    cube: np.ndarray
    wavelength: np.ndarray
    latitude: np.ndarray
    longitude: np.ndarray
    mask: np.ndarray
    scale: float
    metadata: dict
    chlorophyll: object = None


def login():
    import earthaccess
    if not (os.environ.get("EARTHDATA_USERNAME") and os.environ.get("EARTHDATA_PASSWORD")):
        raise RuntimeError("Set EARTHDATA_USERNAME and EARTHDATA_PASSWORD as environment secrets; interactive login is disabled")
    auth = earthaccess.login(strategy="environment", persist=False)
    if not auth.authenticated:
        raise RuntimeError("Earthdata authentication failed")
    return auth


def search(bbox, temporal, collection=PACE_COLLECTION, count=5):
    import earthaccess
    validate_bbox(bbox)
    return earthaccess.search_data(short_name=collection, bounding_box=tuple(bbox),
                                   temporal=tuple(temporal), count=count)


def validate_bbox(bbox):
    if bbox is None:
        return
    west, south, east, north = bbox
    if not (-180 <= west < east <= 180 and -90 <= south < north <= 90):
        raise ValueError("bbox must be west south east north, without dateline crossing")


@contextmanager
def open_granule(scene, collection=PACE_COLLECTION, mission="pace"):
    """Authenticated HTTPS byte-range access; handles stay alive until crop loads."""
    import earthaccess
    login()
    found = earthaccess.search_data(short_name=collection, granule_name=scene, count=2)
    if len(found) != 1:
        raise ValueError(f"Expected one exact granule, found {len(found)}; specify collection/versioned name")
    handles = earthaccess.open(found)
    if len(handles) != 1:
        for handle in handles:
            handle.close()
        raise ValueError("Granule resolves to multiple files; select a single spectral product")
    try:
        yield read(handles[0], mission), {"scene": scene, "collection": collection,
                                         "access": "earthaccess HTTPS streaming"}
    finally:
        handles[0].close()


def read(path_or_file, mission="pace"):
    import hypercoast
    if mission == "pace":
        return hypercoast.read_pace(path_or_file, engine="h5netcdf")
    if mission == "emit":
        # Native geometry retains per-pixel lat/lon; no interpolation of spectra.
        return hypercoast.read_emit(path_or_file, ortho=False, engine="h5netcdf")
    raise ValueError("mission must be pace or emit")


def public_file(cache):
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    target = cache / PUBLIC_SCENE
    if not target.exists():
        temporary = target.with_suffix(".part")
        urllib.request.urlretrieve(PUBLIC_URL, temporary)
        if sha256(temporary) != PUBLIC_SHA256:
            raise ValueError("Public fixture checksum mismatch; refusing changed/corrupted input")
        temporary.replace(target)
    if sha256(target) != PUBLIC_SHA256:
        raise ValueError("Cached public fixture checksum mismatch")
    return target


def prepare(ds, bbox=None, max_side=96, wavelength_range=(400, 700), metadata=None,
            allow_missing_flags=False, pixel_window=None):
    """Load only an enclosing swath crop; invalid pixels remain excluded from fit.

    Pixel-window is (row_start, row_stop, col_start, col_stop), useful for the
    fixed public fixture. A single global scale preserves spectral shape.
    """
    import xarray as xr
    validate_bbox(bbox)
    if max_side < 2:
        raise ValueError("max_side must be >= 2")
    variable = "Rrs" if "Rrs" in ds else "reflectance"
    if variable not in ds:
        raise ValueError("No Rrs/reflectance spectral variable")
    spectral = "wavelength" if "wavelength" in ds[variable].dims else "wavelengths"
    spatial = [d for d in ds[variable].dims if d != spectral]
    if len(spatial) != 2:
        raise ValueError("Expected two spatial dimensions")
    lat_name = "latitude" if "latitude" in ds else "lat"
    lon_name = "longitude" if "longitude" in ds else "lon"
    lat, lon = xr.broadcast(ds[lat_name], ds[lon_name])
    lat, lon = lat.transpose(*spatial), lon.transpose(*spatial)
    latv, lonv = lat.values, lon.values
    in_box = np.isfinite(latv) & np.isfinite(lonv)
    if bbox is not None:
        w, s, e, n = bbox
        in_box &= (lonv >= w) & (lonv <= e) & (latv >= s) & (latv <= n)
    if pixel_window is not None:
        window = np.zeros(in_box.shape, bool)
        r0, r1, c0, c1 = pixel_window
        if not (0 <= r0 < r1 <= window.shape[0] and 0 <= c0 < c1 <= window.shape[1]):
            raise ValueError("Pixel window outside swath")
        window[r0:r1, c0:c1] = True
        in_box &= window
    rows, cols = np.where(in_box)
    if not len(rows):
        raise ValueError("No scene pixels intersect bbox/window")
    stride = max(1, int(np.ceil(max(np.ptp(rows) + 1, np.ptp(cols) + 1) / max_side)))
    rs, cs = slice(rows.min(), rows.max() + 1, stride), slice(cols.min(), cols.max() + 1, stride)
    subset = ds.isel({spatial[0]: rs, spatial[1]: cs})
    wl = np.asarray(subset[spectral].values)
    bands = (wl >= wavelength_range[0]) & (wl <= wavelength_range[1])
    if "good_wavelengths" in subset and subset.good_wavelengths.ndim == 1:
        bands &= subset.good_wavelengths.values.astype(bool)
    if bands.sum() < 4:
        raise ValueError("Too few usable visible bands")
    cube = subset[variable].isel({spectral: np.where(bands)[0]}).transpose(*spatial, spectral).values.astype(float)
    diagnostics = {"in_bbox": int(in_box[rs, cs].sum()),
                   "finite_spectra": int(np.isfinite(cube).all(axis=2).sum()),
                   "nonnegative_spectra": int((np.isfinite(cube).all(axis=2) & (cube >= 0).all(axis=2)).sum()),
                   "entirely_missing_wavelengths_nm": wl[bands][~np.isfinite(cube).any(axis=(0, 1))].tolist()}
    mask = in_box[rs, cs] & np.isfinite(cube).all(axis=2) & (cube >= 0).all(axis=2) & (np.linalg.norm(cube, axis=2) > 0)
    quality = "unavailable"
    if "l2_flags" in subset:
        flags = subset.l2_flags
        meanings = str(flags.attrs.get("flag_meanings", "")).split()
        bits = np.asarray(flags.attrs.get("flag_masks", []))
        if len(meanings) != len(bits) or not len(bits):
            raise ValueError("l2_flags lacks valid flag_meanings/flag_masks metadata")
        flag_values = flags.transpose(*spatial).values
        mask &= np.isfinite(flag_values)
        values = np.nan_to_num(flag_values).astype(np.uint32)
        rejected = [name for name in meanings if name.upper() in REJECT_FLAGS]
        diagnostics["flag_counts_in_bbox"] = {name: int(np.count_nonzero((values & np.uint32(int(bit) & 0xffffffff))[in_box[rs, cs]])) for name, bit in zip(meanings, bits) if name.upper() in REJECT_FLAGS}
        reject = np.uint32(0)
        for name, bit in zip(meanings, bits):
            if name.upper() in REJECT_FLAGS:
                reject |= np.uint32(int(bit) & 0xffffffff)
        mask &= (values & reject) == 0
        quality = {"rejected_flags": rejected}
    elif not allow_missing_flags:
        raise ValueError("No l2_flags; provide a quality-controlled product or explicitly allow missing flags for exploratory work")
    if mask.sum() < 16:
        raise ValueError(f"Only {mask.sum()} valid pixels after quality/band filtering. Diagnostics: {diagnostics}")
    scale = float(cube[mask].max())
    cube /= scale
    cube[~mask] = np.nan
    chlorophyll = None
    if "chlor_a" in subset:
        chlorophyll = subset.chlor_a.transpose(*spatial).values
    info = dict(metadata or {})
    info.update({"variable": variable, "units": ds[variable].attrs.get("units", "unknown"),
                 "quality_diagnostics": diagnostics,
                 "bbox": bbox, "stride": stride, "valid_pixels": int(mask.sum()),
                 "quality": quality, "wavelength_range_nm": list(wavelength_range),
                 "global_scale": scale, "negative_policy": "exclude entire pixel; no offset or clipping"})
    return Scene(cube, wl[bands], latv[rs, cs], lonv[rs, cs], mask, scale, info, chlorophyll)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def attach_chlorophyll(aop, bgc):
    """Attach standard chlor_a only after confirming exact swath geolocation.

    Refuse nearest-neighbor/time-mismatched products. Reference quality masks
    chlorophyll only; BGC product failures must not erase valid AOP spectra.
    """
    if "chlor_a" not in bgc:
        raise ValueError("BGC product has no chlor_a")
    for coordinate in ("latitude", "longitude"):
        a, b = aop[coordinate].values, bgc[coordinate].values
        if a.shape != b.shape or not np.allclose(a, b, atol=1e-5, rtol=0, equal_nan=True):
            raise ValueError("AOP/BGC geolocation mismatch; provide the same acquisition and processing version")
    if "l2_flags" not in bgc:
        raise ValueError("Standard BGC reference needs l2_flags")
    flags = bgc.l2_flags
    meanings = str(flags.attrs.get("flag_meanings", "")).split()
    bits = np.asarray(flags.attrs.get("flag_masks", []))
    if len(meanings) != len(bits) or not len(bits):
        raise ValueError("BGC flag definitions missing")
    reject = np.uint32(0)
    for name, bit in zip(meanings, bits):
        if name.upper() in REJECT_FLAGS:
            reject |= np.uint32(int(bit) & 0xffffffff)
    valid = (flags.astype("uint32") & reject) == 0
    return aop.assign(chlor_a=bgc.chlor_a.where(valid & (bgc.chlor_a > 0)))


def load_snapshot(path):
    """Read a small analysis snapshot, with pickle disabled and schema checks."""
    with np.load(path, allow_pickle=False) as data:
        required = {"cube", "wavelength", "latitude", "longitude", "mask", "scale", "chlorophyll", "metadata"}
        if not required.issubset(data.files):
            raise ValueError("Incomplete analysis snapshot")
        cube = np.array(data["cube"], float)
        mask = np.array(data["mask"], bool)
        wl = np.array(data["wavelength"], float)
        lat, lon, chl = [np.array(data[k], float) for k in ("latitude", "longitude", "chlorophyll")]
        scale = float(data["scale"])
        metadata = json.loads(str(data["metadata"]))
    if cube.ndim != 3 or mask.shape != cube.shape[:2] or wl.shape != (cube.shape[2],):
        raise ValueError("Incompatible snapshot dimensions")
    if any(a.shape != mask.shape for a in (lat, lon, chl)) or not mask.any():
        raise ValueError("Invalid snapshot coordinates/mask")
    if not np.isfinite(cube[mask]).all() or np.any(cube[mask] < 0) or not np.isfinite(scale) or scale <= 0:
        raise ValueError("Snapshot requires finite, nonnegative valid spectra and positive scale")
    metadata.update(snapshot_sha256=sha256(path), access="archived analysis snapshot")
    return Scene(cube, wl, lat, lon, mask, scale, metadata, chl)
