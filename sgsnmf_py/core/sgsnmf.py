"""Projected-gradient port of Wang et al. SGSNMF, revised 2018.

Not multiplicative NMF. Attribution and redistribution scope: THIRD_PARTY.md.
"""
from dataclasses import dataclass
from time import process_time
import numpy as np
from scipy.optimize import nnls
from .slic_hsi import slic_hsi


@dataclass
class Result:
    endmembers: np.ndarray
    abundances: np.ndarray
    reconstruction_sse: list
    iterations: int
    stop_reason: str
    segmentation: object = None


def fcls(a, x):
    """Original augmented NNLS: [1e-5*A; ones], not exact equality FCLS."""
    a, x = np.asarray(a, float), np.asarray(x, float)
    augmented = np.vstack((1e-5 * a, np.ones(a.shape[1])))
    targets = np.vstack((1e-5 * x, np.ones(x.shape[1])))
    h = np.column_stack([nnls(augmented, t, maxiter=100 * a.shape[1])[0]
                         for t in targets.T])
    return np.clip(h, np.finfo(float).eps, 1)


def _subproblem(v, w, h, tol, wp=None, confidence=None, lam=0, maxiter=100):
    h = h.copy()
    wt_v, wt_w = w.T @ v, w.T @ w
    alpha, beta = 1., 0.1
    for iteration in range(1, maxiter + 1):
        grad = wt_w @ h - wt_v
        if wp is not None:
            weighted = wp[:, None] * h
            norm = np.linalg.norm(weighted, axis=0)
            # Defined zero subgradient at zero columns (upstream yields NaN).
            unit = np.divide(weighted, norm, out=np.zeros_like(h), where=norm > 0)
            grad += lam * wp[:, None] * confidence[None, :] * unit
        if np.linalg.norm(grad[(grad < 0) | (h > 0)]) < tol:
            break
        for inner in range(20):
            trial = np.maximum(h - alpha * grad, 0)
            delta = trial - h
            sufficient = 0.99 * np.sum(grad * delta) + 0.5 * np.sum((wt_w @ delta) * delta) < 0
            if inner == 0:
                decrease = not sufficient
                previous = h
            if decrease:
                if sufficient:
                    h = trial
                    break
                alpha *= beta
            else:
                if not sufficient or np.array_equal(previous, trial):
                    h = previous
                    break
                alpha /= beta
                previous = trial
    return h, grad, iteration


def sgsnmf(x, w0, h0, segmentation, lam=0.3, tol=0.05,
           max_iter=100, time_limit=600.):
    """x: bands x pixels; w0: bands x p; h0: p x pixels.

    Shared initial matrices and segmentation are necessary for equivalence.
    ASC appends 15, matching upstream. No post-hoc abundance normalization.
    """
    x, w, h = [np.array(v, dtype=float, copy=True) for v in (x, w0, h0)]
    if x.ndim != 2 or w.ndim != 2 or h.ndim != 2 or (w.shape[0], h.shape[1]) != x.shape or w.shape[1] != h.shape[0]:
        raise ValueError("Incompatible factorization shapes")
    if any(not np.isfinite(v).all() or np.any(v < 0) for v in (x, w, h)):
        raise ValueError("NMF requires finite nonnegative inputs")
    if lam < 0 or tol < 0 or max_iter < 1 or time_limit <= 0:
        raise ValueError("Invalid optimizer parameters")
    seg = segmentation
    if len(seg.labels) != x.shape[1] or seg.centers.shape[0] != x.shape[0] or len(seg.confidence) != x.shape[1]:
        raise ValueError("Segmentation does not match x")
    if not np.array_equal(np.unique(seg.labels), np.arange(seg.centers.shape[1])):
        raise ValueError("Groups must be contiguous and nonempty")
    if not np.isfinite(seg.confidence).all() or np.any(seg.confidence < 0):
        raise ValueError("Invalid confidence")
    grad_w, grad_h = w @ (h @ h.T) - x @ h.T, (w.T @ w) @ h - w.T @ x
    initial = np.sqrt(np.sum(grad_w**2) + np.sum(grad_h**2))
    tol_w = max(0.001, tol) * initial
    tol_h = np.full(seg.centers.shape[1], tol_w)
    history = [float(np.sum((x - w @ h)**2)), 0.]
    consecutive = cumulative = completed = 0
    start = process_time()
    reason = "max_iter"
    for iteration in range(1, max_iter + 1):
        if process_time() - start > time_limit:
            reason = "time_limit"
            break
        if history[-2] - history[-1] > 0.0001:
            consecutive = 0
        else:
            consecutive += 1
            cumulative += 1
        if iteration < 5:
            consecutive = 0
        if consecutive >= 5 and cumulative >= 20:
            reason = "stagnation"
            break
        weights = 1 / (w.shape[1]**2 * fcls(w, seg.centers) + 1)
        tw = np.vstack((w, 15 * np.ones(w.shape[1])))
        tx = np.vstack((x, 15 * np.ones(x.shape[1])))
        for group in range(seg.centers.shape[1]):
            mask = seg.labels == group
            h[:, mask], _, count = _subproblem(tx[:, mask], tw, h[:, mask], tol_h[group],
                                               weights[:, group], seg.confidence[mask], lam)
            if count == 1:
                tol_h[group] *= 0.1
        wt, _, count = _subproblem(x.T, h.T, w.T, tol_w)
        w = wt.T
        if count == 1:
            tol_w *= 0.1
        history.append(float(np.sum((x - w @ h)**2)))
        completed = iteration
    if not np.isfinite(w).all() or not np.isfinite(h).all():
        raise FloatingPointError("Nonfinite optimizer result")
    return Result(w, h, history[:1] + history[2:], completed, reason, seg)


def initialize(x, p, seed=0):
    """Deterministic farthest-point initialization; distinct from original VCA."""
    if p < 2 or p > min(x.shape):
        raise ValueError("p must be in [2, min(bands, pixels)]")
    rng = np.random.default_rng(seed)
    chosen = [int(rng.integers(x.shape[1]))]
    distance = np.full(x.shape[1], np.inf)
    for _ in range(p - 1):
        distance = np.minimum(distance, np.sum((x - x[:, chosen[-1], None])**2, axis=0))
        distance[chosen] = -1
        chosen.append(int(np.argmax(distance)))
    w = np.maximum(x[:, chosen], np.finfo(float).eps)
    return w, fcls(w, x)


def fit(cube, p=4, n_segments=100, compactness=0.5, mask=None, seed=0, w0=None, **kwargs):
    cube = np.asarray(cube, float)
    seg = slic_hsi(cube, n_segments, compactness, mask=mask)
    x = cube.reshape(-1, cube.shape[-1], order="F")[seg.label_image.ravel(order="F") >= 0].T
    if w0 is None:
        w, h = initialize(x, p, seed)
    else:
        w = np.asarray(w0, float)
        if w.shape != (x.shape[0], p):
            raise ValueError("Initial domain spectra must be bands x p")
        h = fcls(w, x)
    return sgsnmf(x, w, h, seg, **kwargs)
