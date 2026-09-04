"""SLIC adapter; differs from Wang's angular/hexagonal SLIC (see docs)."""
from dataclasses import dataclass
import numpy as np
from skimage.segmentation import slic


@dataclass
class Segmentation:
    labels: np.ndarray  # zero-based, valid pixels only, Fortran spatial order
    centers: np.ndarray  # bands x groups
    confidence: np.ndarray
    spacing: float
    label_image: np.ndarray  # -1 at invalid pixels


def from_labels(cube, labels, compactness=0.5, spacing=None):
    cube = np.asarray(cube, float)
    labels = np.asarray(labels, int)
    if cube.ndim != 3 or labels.shape != cube.shape[:2]:
        raise ValueError("Expected rows x cols x bands and rows x cols labels")
    valid = labels >= 0
    if not valid.any() or not np.isfinite(cube[valid]).all():
        raise ValueError("Segmentation needs finite valid pixels")
    ids = np.unique(labels[valid])
    image = np.full(labels.shape, -1, int)
    image[valid] = np.searchsorted(ids, labels[valid])
    flat = image.ravel(order="F")
    x = cube.reshape(-1, cube.shape[-1], order="F")[flat >= 0]
    labs = flat[flat >= 0]
    centers = np.array([x[labs == i].mean(axis=0) for i in range(len(ids))]).T
    rr, cc = np.indices(labels.shape)
    r, c = rr.ravel(order="F")[flat >= 0], cc.ravel(order="F")[flat >= 0]
    spacing = float(spacing or np.sqrt(valid.sum() / len(ids)))
    confidence = np.empty(len(x))
    for i in range(len(ids)):
        sel = labs == i
        # MATLAB round for positive one-based positions, then return to zero-based.
        cr = np.floor(np.mean(r[sel]) + 1.5) - 1
        cl = np.floor(np.mean(c[sel]) + 1.5) - 1
        norms = np.linalg.norm(x[sel], axis=1) * np.linalg.norm(centers[:, i])
        cosine = np.divide(x[sel] @ centers[:, i], norms,
                           out=np.ones(sel.sum()), where=norms > 0)
        dc = np.arccos(np.clip(cosine, 0, 1))
        ds = (r[sel] - cr) ** 2 + (c[sel] - cl) ** 2
        confidence[sel] = 1 / (np.sqrt(dc + ds / spacing**2 * compactness) + np.finfo(float).eps)
    return Segmentation(labs, centers, confidence, spacing, image)


def slic_hsi(cube, n_segments=100, compactness=0.5, max_iter=10, mask=None):
    cube = np.asarray(cube, float)
    if cube.ndim != 3 or n_segments < 1 or compactness <= 0 or max_iter < 1:
        raise ValueError("Invalid cube or SLIC parameters")
    if mask is None:
        mask = np.isfinite(cube).all(axis=2)
    mask = np.asarray(mask, bool)
    if mask.shape != cube.shape[:2] or not mask.any():
        raise ValueError("Empty or incompatible SLIC mask")
    if not np.isfinite(cube[mask]).all():
        raise ValueError("Valid spectra must be finite")
    image = np.where(mask[:, :, None], cube, 0)
    labels = slic(image, n_segments=min(n_segments, int(mask.sum())),
                  compactness=compactness, max_num_iter=max_iter, sigma=0,
                  convert2lab=False, channel_axis=-1, start_label=1,
                  enforce_connectivity=True, mask=mask) - 1
    return from_labels(cube, labels, compactness)
