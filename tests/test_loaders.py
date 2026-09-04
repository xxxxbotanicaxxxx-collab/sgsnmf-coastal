import numpy as np
import pytest
import xarray as xr
from sgsnmf_py.io.loaders import prepare, login, attach_chlorophyll


def dataset():
    lat, lon = np.meshgrid(np.linspace(8, 9, 8), np.linspace(-77, -76, 8), indexing="ij")
    ds = xr.Dataset({"Rrs": (("y", "x", "wavelength"), np.full((8, 8, 5), 0.01)),
                     "l2_flags": (("y", "x"), np.zeros((8, 8), dtype=np.int32))},
                    coords={"latitude": (("y", "x"), lat), "longitude": (("y", "x"), lon),
                            "wavelength": [400, 450, 500, 600, 700]})
    ds.l2_flags.attrs.update(flag_meanings="ATMFAIL LAND", flag_masks=[1, 2])
    return ds


def test_bbox_quality_and_negative_policy():
    ds = dataset()
    ds.l2_flags.values[1, 1] = 2
    ds.Rrs.values[2, 2, 2] = -0.01
    scene = prepare(ds, bbox=(-77, 8, -76, 9))
    assert scene.mask.sum() == 62
    assert not scene.mask[1, 1] and not scene.mask[2, 2]
    assert np.isnan(scene.cube[1, 1]).all()
    assert scene.scale == 0.01
    with pytest.raises(ValueError, match="No scene"):
        prepare(ds, bbox=(0, 0, 1, 1))


def test_missing_flags_and_credentials_fail_explicitly(monkeypatch):
    with pytest.raises(ValueError, match="No l2_flags"):
        prepare(dataset().drop_vars("l2_flags"))
    monkeypatch.delenv("EARTHDATA_USERNAME", raising=False)
    monkeypatch.delenv("EARTHDATA_PASSWORD", raising=False)
    with pytest.raises(RuntimeError, match="environment secrets"):
        login()


def test_chlorophyll_requires_matching_geolocation_and_separate_quality():
    aop = dataset()
    bgc = aop.drop_vars("Rrs").copy(deep=True)
    bgc["chlor_a"] = (("y", "x"), np.ones((8, 8)))
    bgc.l2_flags.values[0, 0] = 1
    merged = attach_chlorophyll(aop, bgc)
    assert merged.l2_flags.values[0, 0] == 0
    assert np.isnan(merged.chlor_a.values[0, 0])
    assert "chlor_a" in merged
    bgc.latitude.values[0, 0] += 0.01
    with pytest.raises(ValueError, match="geolocation mismatch"):
        attach_chlorophyll(aop, bgc)
