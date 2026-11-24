"""
! This file is part of TREKIS-3
Plots of material and ion data: electron / hole / photon mean free paths, elastic MFPs,
DOS, and ion MFP, energy loss, range and effective charge.

The MFP-type files have no header, so the columns are
    energy | one column per shell (atom by atom) | total
Shell names are taken from the header of a Total_energies.txt found near the file
(run folders contain it) or passed as `shell_names`; if neither matches the number of
columns, generic names are used.
"""
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from . import labels
from .readers import find_shell_names, read_dos, read_table
from .style import LINESTYLES, PALETTE, SINGLE, apply_style

ANG = labels.ANGSTROM
HUGE = 1e29            # TREKIS-3 writes 1e30 for "infinite" mean free paths


def _clean(y):
    y = np.array(y, float)
    y[~np.isfinite(y) | (y >= HUGE) | (y <= 0)] = np.nan
    return y


def _require(Y):
    from .plots import require_data
    require_data(_clean(Y), "mean free path / energy-loss data")


def _shell_labels(path, n_shells, shell_names):
    names = shell_names if shell_names is not None else find_shell_names(path)
    if names is not None and len(names) == n_shells:
        return [labels.series_label(n) for n in names]
    if names is not None:
        print(f"Note: {Path(path).name} has {n_shells} shell columns but {len(names)} shell names "
              f"were found; using generic names (pass shell_names=[...] to set them).")
    return [f"Shell {i + 1}" for i in range(n_shells)]


def _finish(ax, xlabel, ylabel, xlog, ylog, xlim=None, ylim=None, nlines=0):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
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
    """Shell curves (coloured, thin) and the total (black, thick, last column of Y)."""
    n = Y.shape[1] - 1
    for i in range(n):
        ax.plot(x, _clean(Y[:, i]), color=PALETTE[i % len(PALETTE)],
                ls=LINESTYLES[(i // len(PALETTE)) % len(LINESTYLES)], lw=1.1, label=names[i])
    ax.plot(x, _clean(Y[:, -1]), color="black", lw=1.8, label=total_label)


# --------------------------------------------------------------------------- #
#  Electron / hole / photon inelastic mean free paths
# --------------------------------------------------------------------------- #
_MFP = {
    "mfp_electron": ("Electron energy (eV)", "Electron inelastic MFP"),
    "mfp_hole": ("Hole energy (eV)", "Hole inelastic MFP"),
    "mfp_photon": ("Photon energy (eV)", "Photon attenuation length"),
}


def plot_mfp(path, kind=None, shell_names=None, xlog=True, ylog=True, xlim=None, ylim=None,
             ydecades=None, figsize=None):
    """Inelastic mean free path versus energy, for each shell and in total.
    Default y range: from just below the smallest total MFP up to `ydecades`
    decades above it (5 for electrons/holes, 8 for photons)."""
    from .readers import detect_kind
    apply_style()
    kind = kind or detect_kind(path)
    xlabel, what = _MFP[kind]
    data, _ = read_table(path)
    x, Y = data[:, 0], data[:, 1:]
    _require(Y)
    fig, ax = plt.subplots(figsize=figsize or SINGLE)
    names = _shell_labels(path, Y.shape[1] - 1, shell_names)
    _shells_and_total(ax, x, Y, names)
    if ylim is None and ylog:
        tot = _clean(Y[:, -1])
        if np.isfinite(tot).any():
            lo = np.nanmin(tot)
            ylim = (lo / 2, lo * 10 ** (ydecades or (8 if kind == "mfp_photon" else 5)))
    _finish(ax, xlabel, f"{what} ({ANG})", xlog, ylog, xlim, ylim, Y.shape[1])
    return fig


def plot_emfp(path, xlog=True, ylog=True, xlim=None, ylim=None, ydecades=5, figsize=None):
    """Elastic mean free path versus energy (total, and emitted / absorbed if present)."""
    apply_style()
    data, _ = read_table(path)
    x, Y = data[:, 0], data[:, 1:]
    _require(Y)
    who = "Hole" if Path(path).name.startswith("OUTPUT_Hole") else "Electron"
    fig, ax = plt.subplots(figsize=figsize or SINGLE)
    for i, nm in enumerate(["Total", "Emission", "Absorption"][:Y.shape[1]]):
        ax.plot(x, _clean(Y[:, i]), color="black" if i == 0 else PALETTE[i - 1],
                lw=1.8 if i == 0 else 1.1, label=nm)
    if ylim is None and ylog and np.isfinite(_clean(Y[:, 0])).any():
        lo = np.nanmin(_clean(Y[:, 0]))
        ylim = (lo / 2, lo * 10 ** ydecades)
    _finish(ax, f"{who} energy (eV)", f"{who} elastic MFP ({ANG})", xlog, ylog, xlim, ylim, 3)
    return fig


# --------------------------------------------------------------------------- #
#  Density of states
# --------------------------------------------------------------------------- #
def plot_dos(path, xlim=None, figsize=None):
    """DOS (left axis) and its integral (right axis) versus energy."""
    apply_style()
    cols, units = read_dos(path)
    E = cols[next(k for k in cols if k.lower().startswith("e"))]
    dos = next(k for k in cols if k.upper().startswith("DOS"))
    #integ = next((k for k in cols if "int" in k.lower()), None)
    integ = next((k for k in cols if "eff" in k.lower()), None)
    fig, ax = plt.subplots(figsize=figsize or SINGLE)
    ax.plot(E, cols[dos], color=PALETTE[0], label="DOS")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel(f"Density of states ({units.get(dos, 'arb. units')})")
    ax.set_ylim(bottom=0)
    if xlim:
        ax.set_xlim(*xlim)
    if integ:
        ax2 = ax.twinx()
        ax2.plot(E, cols[integ], color=PALETTE[1], ls="--", label="Effective mass")
        ax2.set_ylabel(f"Effective mass ({units.get(integ, '$m_e$')})")
        ax2.minorticks_on()
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)
    return fig


# --------------------------------------------------------------------------- #
#  Ion data: IMFP, dE/dx, effective charge, range (energies in MeV)
# --------------------------------------------------------------------------- #
def _ion_shell_plot(path, ylabel, ylog, ydecades, shell_names, xlog, xlim, ylim, figsize):
    apply_style()
    data, _ = read_table(path)
    x, Y = data[:, 0], data[:, 1:]                 # energy in MeV
    _require(Y)
    fig, ax = plt.subplots(figsize=figsize or SINGLE)
    names = _shell_labels(path, Y.shape[1] - 1, shell_names)
    _shells_and_total(ax, x, Y, names)
    return fig, ax, Y


def plot_ion_imfp(path, shell_names=None, xlog=True, ylog=True, xlim=None, ylim=None,
                  ydecades=5, figsize=None):
    """Ion inelastic mean free path versus ion energy, per shell and in total."""
    fig, ax, Y = _ion_shell_plot(path, "", ylog, ydecades, shell_names, xlog, xlim, ylim, figsize)
    if ylim is None and ylog and np.isfinite(_clean(Y[:, -1])).any():
        lo = np.nanmin(_clean(Y[:, -1]))
        ylim = (lo / 2, lo * 10 ** ydecades)
    _finish(ax, "Ion energy (MeV)", f"Ion inelastic MFP ({ANG})", xlog, ylog, xlim, ylim, Y.shape[1])
    return fig


def plot_ion_dedx(path, shell_names=None, xlog=True, ylog=False, xlim=None, ylim=None,
                  ydecades=5, figsize=None):
    """Ion electronic energy loss dE/dx versus ion energy, per shell and in total."""
    fig, ax, Y = _ion_shell_plot(path, "", ylog, ydecades, shell_names, xlog, xlim, ylim, figsize)
    if not ylog:
        ax.set_ylim(bottom=0)
    elif ylim is None and np.isfinite(_clean(Y[:, -1])).any():
        hi = np.nanmax(_clean(Y[:, -1]))
        ylim = (hi / 10 ** ydecades, hi * 2)
    _finish(ax, "Ion energy (MeV)", f"Energy loss $\\mathit{{dE/dx}}$ (eV {ANG}$^{{-1}}$)",
            xlog, ylog, xlim, ylim, Y.shape[1])
    return fig


def plot_ion_zeff(path, xlog=True, xlim=None, ylim=None, figsize=None):
    """Ion effective charge versus ion energy."""
    apply_style()
    data, _ = read_table(path)
    fig, ax = plt.subplots(figsize=figsize or SINGLE)
    ax.plot(data[:, 0], data[:, 1], color=PALETTE[0])
    ax.set_xlabel("Ion energy (MeV)")
    ax.set_ylabel("Effective charge $\\mathit{Z}_{\\mathrm{eff}}$")
    if xlog:
        ax.set_xscale("log")
    ax.set_xlim(*(xlim or (data[:, 0].min(), data[:, 0].max())))
    ax.set_ylim(*(ylim or (0, data[:, 1].max() * 1.05)))
    return fig


def plot_ion_range(path, xlog=True, ylog=True, xlim=None, ylim=None, figsize=None):
    """Ion range (micrometres) versus ion energy. Input energies are read as eV
    and the range as Angstrom (units are taken from the file header when given)."""
    apply_style()
    data, header = read_table(path)
    units = re.findall(r"\[([^\]]*)\]", " ".join(header[1:2])) or ["eV", "eV/A", "A"]
    e_scale = {"ev": 1e-6, "kev": 1e-3, "mev": 1.0}.get(units[0].lower(), 1e-6)
    x = data[:, 0] * e_scale                       # MeV
    rng = data[:, 2] * 1e-4                        # Angstrom -> micrometre
    _require(rng)
    fig, ax = plt.subplots(figsize=figsize or SINGLE)
    ax.plot(x, np.where(rng > 0, rng, np.nan), color=PALETTE[0])
    _finish_simple(ax, "Ion energy (MeV)", "Range ($\\mu$m)", xlog, ylog, xlim, ylim)
    return fig


def _finish_simple(ax, xlabel, ylabel, xlog, ylog, xlim, ylim):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    ax.set_xlim(*(xlim or (ax.dataLim.x0, ax.dataLim.x1)))
    if ylim:
        ax.set_ylim(*ylim)
