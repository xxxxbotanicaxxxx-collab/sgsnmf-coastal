"""Portable NetCDF results and noninteractive scientific figures."""
from pathlib import Path
import json
import importlib.metadata
import textwrap
import numpy as np
import xarray as xr
from .domain import label_factors, chlorophyll_check


def save_products(result, scene, output, library=None, parameters=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    p = result.endmembers.shape[1]
    valid = result.segmentation.label_image.ravel(order="F") >= 0
    flat = np.full((valid.size, p), np.nan)
    flat[valid] = result.abundances.T
    maps = flat.reshape((*scene.mask.shape, p), order="F")
    x = scene.cube.reshape(-1, scene.cube.shape[-1], order="F")[valid].T
    residual = np.full(valid.size, np.nan)
    residual[valid] = np.sqrt(np.mean((x - result.endmembers @ result.abundances)**2, axis=0)) * scene.scale
    residual = residual.reshape(scene.mask.shape, order="F")
    labels = label_factors(result.endmembers * scene.scale, library)
    ds = xr.Dataset({
        "abundance": (("row", "col", "component"), maps),
        "endmember": (("wavelength", "component"), result.endmembers * scene.scale),
        "reconstruction_rmse": (("row", "col"), residual),
        "valid_mask": (("row", "col"), scene.mask.astype(np.int8)),
        "superpixel": (("row", "col"), result.segmentation.label_image),
    }, coords={"latitude": (("row", "col"), scene.latitude),
               "longitude": (("row", "col"), scene.longitude),
               "wavelength": scene.wavelength, "component": np.arange(1, p + 1)})
    ds.abundance.attrs.update(units="1", description="Exploratory spectral fractions; soft sum-to-one constraint; not concentrations")
    ds.endmember.attrs.update(units=scene.metadata.get("units", "unknown"))
    ds.wavelength.attrs["units"] = "nm"
    ds.attrs["component_labels"] = json.dumps(labels, ensure_ascii=False)
    if scene.chlorophyll is not None:
        ds["chlor_a_reference"] = (("row", "col"), scene.chlorophyll)
        ds.chlor_a_reference.attrs.update(units="mg m-3", description="Standard PACE L2 chlor_a, independently quality masked; not ground truth")
    ds.to_netcdf(output / "abundances.nc", engine="h5netcdf")
    fig, axes = plt.subplots(1, p, figsize=(4 * p, 4), constrained_layout=True, squeeze=False)
    for i, ax in enumerate(axes[0]):
        art = ax.scatter(scene.longitude[scene.mask], scene.latitude[scene.mask], c=maps[:, :, i][scene.mask],
                         s=12, cmap="viridis", vmin=0, vmax=1, marker="s")
        ax.set(title=textwrap.fill(labels[i], 28), xlabel="Longitud", ylabel="Latitud")
        ax.title.set_fontsize(10)
        ax.set_aspect(1 / max(np.cos(np.deg2rad(np.nanmean(scene.latitude))), 0.1))
        fig.colorbar(art, ax=ax, label="Fracción espectral", shrink=0.7)
    fig.savefig(output / "abundance_maps.png", dpi=170)
    plt.close(fig)
    sanity = chlorophyll_check(maps, scene.chlorophyll, scene.mask)
    if sanity["status"] == "exploratory_only":
        paired = scene.mask & np.isfinite(scene.chlorophyll) & (scene.chlorophyll > 0)
        fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
        for i, ax in enumerate(axes.flat):
            if i >= p:
                ax.axis("off")
                continue
            ax.scatter(np.log10(scene.chlorophyll[paired]), maps[:, :, i][paired], s=8, alpha=0.35)
            rho = sanity["spearman_vs_log10_chlor_a"][i]
            rho_text = f"{rho:.3f}" if rho is not None else "indefinido"
            ax.set(xlabel="log10(chlor_a L2 / mg m⁻³)", ylabel="Fracción espectral",
                   title=f"Componente {i + 1} · Spearman ρ = {rho_text}")
        fig.suptitle(f"Sanity-check exploratorio · {sanity['n']} pares · no ground truth", fontsize=12)
        fig.savefig(output / "chlorophyll_sanity.png", dpi=160)
        plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 4.5), constrained_layout=True)
    for i in range(p):
        ax.plot(scene.wavelength, result.endmembers[:, i] * scene.scale, label=labels[i])
    ax.set(xlabel="Longitud de onda (nm)", ylabel=scene.metadata.get("units", "Reflectancia"), title="Endmembers estimados")
    ax.legend(fontsize=8)
    fig.savefig(output / "endmembers.png", dpi=170)
    plt.close(fig)
    report = {"scene": scene.metadata, "parameters": parameters or {}, "labels": labels,
              "iterations": result.iterations, "stop_reason": result.stop_reason,
              "reconstruction_sse_scaled": result.reconstruction_sse,
              "mean_reconstruction_rmse_physical_units": float(np.nanmean(residual)),
              "max_sum_to_one_error": float(np.max(np.abs(result.abundances.sum(axis=0) - 1))),
              "chlorophyll_sanity_check": sanity,
              "versions": {key: importlib.metadata.version(key) for key in
                           ["numpy", "scipy", "scikit-image", "xarray", "hypercoast", "earthaccess"]},
              "limitations": ["NMF components require external bio-optical validation",
                              "Water constituents interact nonlinearly in Rrs; fractions are not mg/m3 or absorption coefficients",
                              "skimage SLIC is an approximate preprocessing replacement"]}
    (output / "run.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return report
