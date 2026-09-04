import numpy as np
import pytest
from sgsnmf_py.domain import library_initialization, chlorophyll_check
from sgsnmf_py.io.loaders import load_snapshot


def test_domain_library_refuses_extrapolation(tmp_path):
    path = tmp_path / "library.csv"
    path.write_text("wavelength_nm,water,phyto,cdom\n400,0.01,0.02,0.03\n700,0.001,0.002,0.003\n")
    result = library_initialization(path, np.array([400., 550., 700.]), 3)
    assert result.shape == (3, 3)
    with pytest.raises(ValueError, match="extrapolation"):
        library_initialization(path, np.array([350., 700.]), 3)


def test_no_false_chlorophyll_validation():
    maps = np.ones((5, 5, 3))
    mask = np.ones((5, 5), bool)
    assert chlorophyll_check(maps, None, mask)["status"] == "not_available"
    assert chlorophyll_check(maps, np.full((5, 5), np.nan), mask)["status"] == "insufficient_pairs"


def test_snapshot_rejects_incomplete_schema(tmp_path):
    path = tmp_path / "bad.npz"
    np.savez(path, cube=np.ones((2, 2, 4)))
    with pytest.raises(ValueError, match="Incomplete"):
        load_snapshot(path)
