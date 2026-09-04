"""Permutation-aware metrics; spectral angles are in radians."""
import numpy as np
from scipy.optimize import linear_sum_assignment


def sad(reference, estimate):
    a, b = np.asarray(reference, float), np.asarray(estimate, float)
    if a.shape != b.shape or a.ndim != 2:
        raise ValueError("Expected matching bands x components matrices")
    norms = np.linalg.norm(a, axis=0) * np.linalg.norm(b, axis=0)
    if np.any(norms == 0):
        raise ValueError("SAD is undefined for zero spectra")
    return np.arccos(np.clip(np.sum(a * b, axis=0) / norms, -1, 1))


def rmse(reference, estimate, axis=None):
    a, b = np.asarray(reference), np.asarray(estimate)
    if a.shape != b.shape:
        raise ValueError("RMSE shapes must match")
    return np.sqrt(np.mean((a - b) ** 2, axis=axis))


def compare(reference_w, reference_h, w, h):
    if reference_w.shape != w.shape or reference_h.shape != h.shape:
        raise ValueError("Comparison shapes must match")
    na = np.linalg.norm(reference_w, axis=0)
    nb = np.linalg.norm(w, axis=0)
    if np.any(na == 0) or np.any(nb == 0):
        raise ValueError("Cannot match zero endmembers")
    cost = np.arccos(np.clip((reference_w.T @ w) / (na[:, None] * nb), -1, 1))
    _, order = linear_sum_assignment(cost)
    return {"permutation": order.tolist(),
            "sad_radians": sad(reference_w, w[:, order]).tolist(),
            "endmember_rmse": float(rmse(reference_w, w[:, order])),
            "abundance_rmse": float(rmse(reference_h, h[order])),
            "abundance_rmse_per_component": rmse(reference_h, h[order], axis=1).tolist()}
