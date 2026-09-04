"""Explicit domain hypotheses; NMF factors are not constituent concentrations."""
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr

HYPOTHESES = ["agua clara (proxy; no agua pura aislada)",
              "fitoplancton (hipótesis espectral)",
              "CDOM (hipótesis espectral)",
              "sedimentos (hipótesis espectral)"]


def library_initialization(csv_file, wavelengths, p):
    """CSV columns: wavelength_nm followed by p measured, nonnegative Rrs spectra.

    Order: clear water, phytoplankton-rich, CDOM-rich, sediment-rich. The
    library must use the same radiometric quantity and units as the scene.
    """
    library = np.loadtxt(csv_file, delimiter=",", skiprows=1)
    if library.ndim != 2 or library.shape[1] != p + 1 or not np.isfinite(library).all():
        raise ValueError("Domain library needs wavelength_nm plus p finite spectra")
    if np.any(np.diff(library[:, 0]) <= 0) or np.any(library[:, 1:] < 0):
        raise ValueError("Library wavelengths must increase and spectra be nonnegative")
    if wavelengths.min() < library[0, 0] or wavelengths.max() > library[-1, 0]:
        raise ValueError("Library does not cover the scene; extrapolation forbidden")
    return np.column_stack([np.interp(wavelengths, library[:, 0], library[:, i + 1]) for i in range(p)])


def label_factors(w, initial_library=None):
    if initial_library is None:
        return [f"Componente {i + 1} (sin asignación bioóptica)" for i in range(w.shape[1])]
    a = initial_library / np.maximum(np.linalg.norm(initial_library, axis=0), 1e-300)
    b = w / np.maximum(np.linalg.norm(w, axis=0), 1e-300)
    rows, cols = linear_sum_assignment(np.arccos(np.clip(a.T @ b, -1, 1)))
    names = [""] * w.shape[1]
    for i, j in zip(rows, cols):
        names[j] = HYPOTHESES[i]
    return names


def chlorophyll_check(abundance_maps, chlorophyll, mask):
    """Descriptive rank association only; no tuning against L2 chlorophyll."""
    if chlorophyll is None:
        return {"status": "not_available", "reason": "No colocated standard chlor_a product supplied"}
    if chlorophyll.shape != mask.shape:
        raise ValueError("Chlorophyll must be geolocated on the abundance grid")
    valid = mask & np.isfinite(chlorophyll) & (chlorophyll > 0)
    if valid.sum() < 20:
        return {"status": "insufficient_pairs", "n": int(valid.sum())}
    rho = []
    for i in range(abundance_maps.shape[-1]):
        values = abundance_maps[:, :, i][valid]
        result = spearmanr(values, np.log10(chlorophyll[valid])).statistic
        rho.append(float(result) if np.isfinite(result) else None)
    return {"status": "exploratory_only", "n": int(valid.sum()), "spearman_vs_log10_chlor_a": rho,
            "interpretation": "Spatially dependent samples; correlations are not constituent identification or ground truth accuracy"}
