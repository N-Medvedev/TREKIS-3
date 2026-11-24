"""
! This file is part of TREKIS-3
Plots of material and ion data: electron / hole / photon mean free paths, elastic mean
free paths, density of states, and ion mean free path, energy loss, range and charge.

Plot functions (each takes the path of a data file and returns a matplotlib Figure):

    plot_mfp        OUTPUT_Electron_IMFPs*, OUTPUT_Hole_IMFPs*, OUTPUT_Photon_IMFPs*
                    inelastic MFP (photons: attenuation length) vs. energy, per shell + total
    plot_emfp       OUTPUT_Electron*EMFPs*, OUTPUT_Hole*EMFPs*      elastic MFP vs. energy
    plot_dos        OUTPUT_<material>_DOS_analysis.dat              DOS and effective mass
    plot_ion_imfp   OUTPUT_<ion>*_IMFP.dat                          ion inelastic MFP vs. ion energy
    plot_ion_dedx   OUTPUT_<ion>*_dEdx.dat                          electronic energy loss
    plot_ion_zeff   OUTPUT_<ion>*_effective_charges.dat             effective ion charge
    plot_ion_range  OUTPUT_<ion>*_Range.dat                         ion range (micrometres)

The MFP-type files have no header, their columns are
    energy | one column per shell (atom by atom) | total
Shell names for the legends are taken from the header of a Total_energies.txt found near
the file (run folders contain it) or passed as `shell_names=[...]`; if neither matches
the number of columns, generic names ('Shell 1', ...) are used. Values of 1e30 (TREKIS-3's
"infinite") and non-positive values are not drawn.

Common keyword arguments: xlog, ylog, xlim, ylim, ydecades (dynamic range of a log y axis,
counted from the smallest total MFP), figsize and title (True, False or a string).
"""
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from . import labels
from .common import add_title, require_data, set_ylabel_fit
from .readers import detect_kind, find_shell_names, read_dos, read_table
from .style import LINESTYLES, PALETTE, SINGLE, TOTAL_COLOUR, apply_style

ANG = labels.ANGSTROM
HUGE = 1e29            # TREKIS-3 writes 1e30 for "infinite" mean free paths


def _clean(y):
    """Copy of y with infinite (>= 1e29), non-positive and non-finite values set to NaN."""
    y = np.array(y, float)
    y[~np.isfinite(y) | (y >= HUGE) | (y <= 0)] = np.nan
    return y


def _require(Y):
    require_data(_clean(Y), "mean free path / energy-loss data")


def _shell_labels(path, n_shells, shell_names):
    """Legend names of the shell columns (see the module docstring)."""
    names = shell_names if shell_names is not None else find_shell_names(path)
    if names is not None and len(names) == n_shells:
        return [labels.series_label(n) for n in names]
    if names is not None:
        print(f"Note: {Path(path).name} has {n_shells} shell columns but {len(names)} shell names "
              f"were found; using generic names (pass shell_names=[...] to set them).")
    return [f"Shell {i + 1}" for i in range(n_shells)]


def _new_figure(figsize):
    return plt.subplots(figsize=figsize or SINGLE, layout="constrained")


def _finish(fig, ax, xlabel, ylabel, xlog, ylog, xlim=None, ylim=None, nlines=0):
    """Axis labels, scales, limits and legend (if there is more than one curve)."""
    ax.set_xlabel(xlabel)
    set_ylabel_fit(ax, ylabel, fig=fig)
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    if xlim:
        ax.set_xlim(*xlim)
    else:
        ax.set_xlim(ax.dataLim.x0, ax.dataLim.x1)
    if ylim:
        ax.set_ylim(*ylim)
    if nlines > 1:
        ax.legend(loc="best", ncol=2 if nlines > 6 else 1)


def _shells_and_total(ax, x, Y, names, total_label="Total"):
    """Shell curves (coloured, thin) and the total (dark, thick, last column of Y)."""
    n = Y.shape[1] - 1
    for i in range(n):
        ax.plot(x, _clean(Y[:, i]), color=PALETTE[i % len(PALETTE)],
                ls=LINESTYLES[(i // len(PALETTE)) % len(LINESTYLES)], lw=1.0, label=names[i])
    ax.plot(x, _clean(Y[:, -1]), color=TOTAL_COLOUR, lw=1.5, label=total_label)


# --------------------------------------------------------------------------- #
#  Electron / hole / photon inelastic mean free paths
# --------------------------------------------------------------------------- #
_MFP = {      # kind: (x label, title, y label)
    "mfp_electron": ("Electron energy (eV)", "Electron inelastic mean free path", "Mean free path"),
    "mfp_hole": ("Hole energy (eV)", "Hole inelastic mean free path", "Mean free path"),
    "mfp_photon": ("Photon energy (eV)", "Photon attenuation length", "Attenuation length"),
}


def plot_mfp(path, kind=None, shell_names=None, xlog=True, ylog=True, xlim=None, ylim=None,
             ydecades=None, figsize=None, title=True):
    """Inelastic mean free path versus energy, for each shell and in total.
    Default y range: from just below the smallest total MFP up to `ydecades`
    decades above it (5 for electrons/holes, 8 for photons)."""
    apply_style()
    kind = kind or detect_kind(path)
    xlabel, main, ylabel = _MFP[kind]
    data, _ = read_table(path)
    x, Y = data[:, 0], data[:, 1:]
    _require(Y)
    fig, ax = _new_figure(figsize)
    names = _shell_labels(path, Y.shape[1] - 1, shell_names)
    _shells_and_total(ax, x, Y, names)
    if ylim is None and ylog:
        tot = _clean(Y[:, -1])
        if np.isfinite(tot).any():
            lo = np.nanmin(tot)
            ylim = (lo / 2, lo * 10 ** (ydecades or (8 if kind == "mfp_photon" else 5)))
    _finish(fig, ax, xlabel, f"{ylabel} ({ANG})", xlog, ylog, xlim, ylim, Y.shape[1])
    add_title(ax, path, main, title)
    return fig


def plot_emfp(path, xlog=True, ylog=True, xlim=None, ylim=None, ydecades=5, figsize=None,
              title=True):
    """Elastic mean free path versus energy (total, and emission / absorption if present)."""
    apply_style()
    data, _ = read_table(path)
    x, Y = data[:, 0], data[:, 1:]
    _require(Y)
    who = "Hole" if Path(path).name.startswith("OUTPUT_Hole") else "Electron"
    fig, ax = _new_figure(figsize)
    for i, nm in enumerate(["Total", "Emission", "Absorption"][:Y.shape[1]]):
        ax.plot(x, _clean(Y[:, i]), color=TOTAL_COLOUR if i == 0 else PALETTE[i - 1],
                lw=1.5 if i == 0 else 1.0, label=nm)
    if ylim is None and ylog and np.isfinite(_clean(Y[:, 0])).any():
        lo = np.nanmin(_clean(Y[:, 0]))
        ylim = (lo / 2, lo * 10 ** ydecades)
    _finish(fig, ax, f"{who} energy (eV)", f"Mean free path ({ANG})", xlog, ylog, xlim, ylim,
            Y.shape[1])
    add_title(ax, path, f"{who} elastic mean free path", title)
    return fig


# --------------------------------------------------------------------------- #
#  Density of states
# --------------------------------------------------------------------------- #
def plot_dos(path, xlim=None, figsize=None, title=True):
    """DOS (left axis) and, if present in the file, the effective mass (right axis,
    dashed) versus energy."""
    apply_style()
    cols, units = read_dos(path)
    E = cols[next(k for k in cols if k.lower().startswith("e"))]
    dos = next(k for k in cols if k.upper().startswith("DOS"))
    mass = next((k for k in cols if "eff" in k.lower()), None)
    fig, ax = _new_figure(figsize)
    ax.plot(E, cols[dos], color=PALETTE[0], label="DOS")
    ax.set_xlabel("Energy (eV)")
    set_ylabel_fit(ax, f"Density of states ({units.get(dos) or 'arb. units'})", fig=fig)
    ax.set_ylim(bottom=0)
    ax.set_xlim(*(xlim or (E.min(), E.max())))
    if mass:
        ax2 = ax.twinx()
        ax2.spines["right"].set_visible(True)
        ax2.grid(False)
        ax2.plot(E, cols[mass], color=PALETTE[1], ls="--", label="Effective mass")
        unit = units.get(mass, "")
        unit = "$m_e$" if unit.lower() in ("", "me", "m_e") else unit
        ax2.set_ylabel(f"Effective mass ({unit})")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax2.legend(h1 + h2, l1 + l2, loc="best")
    add_title(ax, path, "Density of states", title)
    return fig


# --------------------------------------------------------------------------- #
#  Ion data: IMFP, dE/dx, effective charge, range (energies in MeV)
# --------------------------------------------------------------------------- #
def _ion_shell_plot(path, shell_names, figsize):
    """Figure with the per-shell curves and the total of an ion file (energy in MeV)."""
    apply_style()
    data, _ = read_table(path)
    x, Y = data[:, 0], data[:, 1:]
    _require(Y)
    fig, ax = _new_figure(figsize)
    names = _shell_labels(path, Y.shape[1] - 1, shell_names)
    _shells_and_total(ax, x, Y, names)
    return fig, ax, Y


def plot_ion_imfp(path, shell_names=None, xlog=True, ylog=True, xlim=None, ylim=None,
                  ydecades=5, figsize=None, title=True):
    """Ion inelastic mean free path versus ion energy, per shell and in total."""
    fig, ax, Y = _ion_shell_plot(path, shell_names, figsize)
    if ylim is None and ylog and np.isfinite(_clean(Y[:, -1])).any():
        lo = np.nanmin(_clean(Y[:, -1]))
        ylim = (lo / 2, lo * 10 ** ydecades)
    _finish(fig, ax, "Ion energy (MeV)", f"Mean free path ({ANG})", xlog, ylog, xlim, ylim,
            Y.shape[1])
    add_title(ax, path, "Ion inelastic mean free path", title)
    return fig


def plot_ion_dedx(path, shell_names=None, xlog=True, ylog=False, xlim=None, ylim=None,
                  ydecades=5, figsize=None, title=True):
    """Ion electronic energy loss dE/dx versus ion energy, per shell and in total."""
    fig, ax, Y = _ion_shell_plot(path, shell_names, figsize)
    if not ylog:
        ax.set_ylim(bottom=0)
    elif ylim is None and np.isfinite(_clean(Y[:, -1])).any():
        hi = np.nanmax(_clean(Y[:, -1]))
        ylim = (hi / 10 ** ydecades, hi * 2)
    _finish(fig, ax, "Ion energy (MeV)", f"Energy loss $\\mathit{{dE/dx}}$ (eV {ANG}$^{{-1}}$)",
            xlog, ylog, xlim, ylim, Y.shape[1])
    add_title(ax, path, "Electronic energy loss", title)
    return fig


def plot_ion_zeff(path, xlog=True, xlim=None, ylim=None, figsize=None, title=True):
    """Ion effective charge versus ion energy."""
    apply_style()
    data, _ = read_table(path)
    fig, ax = _new_figure(figsize)
    ax.plot(data[:, 0], data[:, 1], color=PALETTE[0])
    ax.set_xlabel("Ion energy (MeV)")
    set_ylabel_fit(ax, "Effective charge $\\mathit{Z}_{\\mathrm{eff}}$", fig=fig)
    if xlog:
        ax.set_xscale("log")
    ax.set_xlim(*(xlim or (data[:, 0].min(), data[:, 0].max())))
    ax.set_ylim(*(ylim or (0, data[:, 1].max() * 1.05)))
    add_title(ax, path, "Ion effective charge", title)
    return fig


def plot_ion_range(path, xlog=True, ylog=True, xlim=None, ylim=None, figsize=None, title=True):
    """Ion range (micrometres) versus ion energy. Input energies are read as eV
    and the range as Angstrom (units are taken from the file header when given)."""
    apply_style()
    data, header = read_table(path)
    units = re.findall(r"\[([^\]]*)\]", " ".join(header[1:2])) or ["eV", "eV/A", "A"]
    e_scale = {"ev": 1e-6, "kev": 1e-3, "mev": 1.0}.get(units[0].lower(), 1e-6)
    x = data[:, 0] * e_scale                       # MeV
    rng = data[:, 2] * 1e-4                        # Angstrom -> micrometre
    _require(rng)
    fig, ax = _new_figure(figsize)
    ax.plot(x, np.where(rng > 0, rng, np.nan), color=PALETTE[0])
    ax.set_xlabel("Ion energy (MeV)")
    set_ylabel_fit(ax, "Range ($\\mu$m)", fig=fig)
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    ax.set_xlim(*(xlim or (ax.dataLim.x0, ax.dataLim.x1)))
    if ylim:
        ax.set_ylim(*ylim)
    add_title(ax, path, "Ion range", title)
    return fig
