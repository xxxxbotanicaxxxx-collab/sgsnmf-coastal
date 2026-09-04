import numpy as np
import pytest
from sgsnmf_py.core.sgsnmf import sgsnmf, _subproblem
from sgsnmf_py.core.slic_hsi import from_labels, slic_hsi
from sgsnmf_py.core.metrics import compare, sad


def test_permutation_aware_metrics():
    w = np.array([[1., 0.1], [0.2, 1.], [0.4, 0.3]])
    h = np.array([[0.2, 0.7], [0.8, 0.3]])
    metrics = compare(w, h, w[:, ::-1], h[::-1])
    assert metrics["permutation"] == [1, 0]
    assert metrics["abundance_rmse"] == 0
    assert max(metrics["sad_radians"]) < 1e-7
    with pytest.raises(ValueError):
        sad(w, np.zeros_like(w))


def test_fortran_group_order_and_mask():
    cube = np.arange(1, 25, dtype=float).reshape(3, 4, 2)
    labels = np.array([[7, 7, 9, -1], [7, 9, 9, 9], [7, 7, 9, 9]])
    seg = from_labels(cube, labels)
    np.testing.assert_array_equal(seg.labels, np.array([0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1]))
    np.testing.assert_allclose(seg.centers[:, 0], cube[labels == 7].mean(axis=0))
    assert np.isfinite(seg.confidence).all()


def test_zero_abundance_subgradient_is_finite():
    h, g, _ = _subproblem(np.ones((3, 5)), np.ones((3, 2)), np.zeros((2, 5)),
                           tol=0.001, wp=np.ones(2), confidence=np.ones(5), lam=0.3)
    assert np.isfinite(g).all() and np.isfinite(h).all()
    assert (h >= 0).all()


def test_optimizer_improves_controlled_reconstruction():
    rng = np.random.default_rng(7)
    w = rng.uniform(0.1, 0.8, (8, 3))
    h = rng.dirichlet(np.ones(3), 100).T
    x = w @ h
    cube = x.T.reshape(10, 10, 8, order="F")
    seg = from_labels(cube, np.zeros((10, 10), int))
    result = sgsnmf(x, w * 0.8, h * 0.8, seg, lam=0, max_iter=15)
    assert result.reconstruction_sse[-1] < 0.01 * result.reconstruction_sse[0]
    assert np.max(np.abs(result.abundances.sum(axis=0) - 1)) < 0.01
    assert (result.endmembers >= 0).all() and (result.abundances >= 0).all()
    with pytest.raises(ValueError):
        sgsnmf(-x, w, h, seg)


def test_slic_excludes_masked_pixels():
    rng = np.random.default_rng(8)
    cube = rng.uniform(0.1, 1, (12, 12, 5))
    mask = np.ones((12, 12), bool)
    mask[:3, :] = False
    cube[~mask] = np.nan
    seg = slic_hsi(cube, 6, mask=mask)
    assert len(seg.labels) == mask.sum()
    assert np.all(seg.label_image[~mask] == -1)
