"""Plot functions for TREKIS-3 output files (publication-quality matplotlib figures)."""
import os
import re
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable

from . import labels
from . import material_plots as mp
from .readers import (_SKIP_DIRS, detect_kind, read_multi_time, read_time_series)
from .style import (DOUBLE, LINESTYLES, PALETTE, SINGLE, TIME_CMAP, apply_style,
                    time_colours, time_label, time_norm)

MAX_LEGEND_TIMES = 8      # more time columns than this -> colour bar instead of a legend


# --------------------------------------------------------------------------- #
#  Helpers
# --------------------------------------------------------------------------- #
def _select_times(mt, times):
    """Keep only the columns closest to the requested times [fs]."""
    if times is None:
        return mt.times, mt.y
    idx = []
    for t in times:
        i = int(np.argmin(np.abs(mt.times - t)))
        if i not in idx:
            idx.append(i)
        if not np.isclose(mt.times[i], t, rtol=0.05):
            print(f"Note: requested time {t:g} fs, using the closest available: {mt.times[i]:g} fs")
    idx.sort()
    return mt.times[idx], mt.y[:, idx]


def require_data(y, what="data"):
    """Raise ValueError if there is nothing to draw (all values zero, negative or missing)."""
    y = np.asarray(y, float)
    if not np.any(np.isfinite(y) & (y > 0)):
        raise ValueError(f"no non-zero {what} (all values are zero or missing), nothing to plot")


def _positive(y):
    y = np.array(y, float)
    y[y <= 0] = np.nan
    return y


def _finish_multi_time(fig, ax, times, colours, norm, curves_labels, legend_loc="best"):
    """Legend (few times) or colour bar (many times)."""
    if len(times) <= MAX_LEGEND_TIMES:
        ax.legend(loc=legend_loc, title="Time", title_fontsize=7.5)
    else:
        sm = ScalarMappable(norm=norm, cmap=TIME_CMAP)
        cb = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.06)
        cb.set_label("Time (fs)")
        cb.ax.tick_params(direction="out")


def _plot_curves(ax, x, y, times, ylog, marker="o"):
    """One curve per time. Small markers keep isolated non-zero points visible
    (e.g. holes that exist in one radial bin only)."""
    colours, norm = time_colours(times)
    for j, t in enumerate(times):
        yy = _positive(y[:, j]) if ylog else y[:, j]
        ax.plot(x, yy, color=colours[j], label=time_label(t),
                marker=marker, ms=2.2, mew=0)
    return colours, norm


def _set_log_ylim(ax, y, decades):
    pos = y[y > 0]
    if pos.size == 0:
        return
    ymax = pos.max()
    ymin = max(pos.min(), ymax * 10.0 ** (-decades))
    ymin = min(ymin, ymax / 100.0)             # always show at least two decades
    ax.set_ylim(ymin / 1.5, ymax * 3)


def _xlim_nonzero(x, y, pad=1.25):
    """x range over which at least one curve is non-zero."""
    nz = np.where((y > 0).any(axis=1))[0]
    if nz.size == 0:
        return x[0], x[-1]
    lo = x[0] / pad if x[0] > 0 else x[0]
    hi = max(x[nz[-1]] * pad, x[0] * 10)       # always show at least one decade
    return lo, hi


# --------------------------------------------------------------------------- #
#  Time series: Total_numbers.txt, Total_energies.txt
# --------------------------------------------------------------------------- #
def plot_time_series(path, hide_empty=True, xlog=True, columns=None, figsize=None):
    """Quantities versus time. Energy columns (unit eV) and particle numbers are
    put in separate panels. All-zero columns are hidden unless hide_empty=False;
    `columns` restricts the plot to the given column names."""
    apply_style()
    ts = read_time_series(path)
    energy, number = [], []
    for i in range(1, len(ts.names)):
        name = ts.names[i]
        if columns is not None and name not in columns:
            continue
        if hide_empty and not np.any(ts.data[:, i] != 0):
            continue
        (energy if ts.units[i].lower() == "ev" else number).append(i)
    groups = [(g, lab) for g, lab in ((number, "Number of particles"),
                                      (energy, "Energy (eV)")) if g]
    if not groups:
        raise ValueError(f"Nothing to plot in {path}")

    fs = figsize or (SINGLE if len(groups) == 1 else (SINGLE[0], SINGLE[1] * 1.55))
    fig, axes = plt.subplots(len(groups), 1, sharex=True, figsize=fs, squeeze=False)
    t = ts.time
    pos = t > 0 if xlog else np.ones_like(t, bool)

    k = 0
    for ax, (idx, ylabel) in zip(axes[:, 0], groups):
        for n, i in enumerate(idx):
            ax.plot(t[pos], ts.data[pos, i], color=PALETTE[n % len(PALETTE)],
                    ls=LINESTYLES[(n // len(PALETTE)) % len(LINESTYLES)],
                    label=labels.series_label(ts.names[i]))
        ax.set_ylabel(ylabel)
        if all(ts.data[pos, i].min() >= 0 for i in idx):
            ax.set_ylim(bottom=0)
        top = max(ts.data[pos, i].max() for i in idx)
        if top > 0:
            ax.set_ylim(top=top * (1.12 + 0.12 * ((len(idx) + 2) // 3)))
        ax.legend(loc="upper center", ncol=min(3, len(idx)), columnspacing=1.0,
                  handletextpad=0.5)
        k += 1
    ax = axes[-1, 0]
    if xlog:
        ax.set_xscale("log")
    ax.set_xlim(t[pos].min(), t.max())
    ax.set_xlabel("Time (fs)")
    fig.align_ylabels(axes[:, 0])
    return fig


# --------------------------------------------------------------------------- #
#  Radial profiles: Radial_*.txt
# --------------------------------------------------------------------------- #
def plot_radial(path, times=None, xlog=True, ylog=True, ydecades=8, xlim=None, ylim=None,
                figsize=None):
    """Radial profile versus distance from the ion trajectory, one curve per time."""
    apply_style()
    mt = read_multi_time(path)
    tt, y = _select_times(mt, times)
    require_data(y)
    fig, ax = plt.subplots(figsize=figsize or SINGLE)
    colours, norm = _plot_curves(ax, mt.x, y, tt, ylog)
    ax.set_xlabel(f"Radius $\\mathit{{r}}$ ({labels.ANGSTROM})")
    ax.set_ylabel(labels.radial_ylabel(path))
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
        _set_log_ylim(ax, y, ydecades)
    ax.set_xlim(*(xlim or _xlim_nonzero(mt.x, y)))
    if ylim:
        ax.set_ylim(*ylim)
    _finish_multi_time(fig, ax, tt, colours, norm, None)
    return fig


# --------------------------------------------------------------------------- #
#  Energy spectra: *_distribution_vs_E[1_eV].txt
# --------------------------------------------------------------------------- #
def plot_spectrum(path, times=None, xlog=True, ylog=True, ydecades=8, xlim=None, ylim=None,
                  figsize=None):
    """Energy distribution of particles, one curve per time."""
    apply_style()
    mt = read_multi_time(path)
    tt, y = _select_times(mt, times)
    require_data(y)
    fig, ax = plt.subplots(figsize=figsize or SINGLE)
    colours, norm = _plot_curves(ax, mt.x, y, tt, ylog)
    ax.set_xlabel("Energy $\\mathit{E}$ (eV)")
    ax.set_ylabel(f"{labels.multi_body(path)} spectrum (arb. units)")
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
        _set_log_ylim(ax, y, ydecades)
    ax.set_xlim(*(xlim or _xlim_nonzero(mt.x, y)))
    if ylim:
        ax.set_ylim(*ylim)
    _finish_multi_time(fig, ax, tt, colours, norm, None)
    return fig


# --------------------------------------------------------------------------- #
#  Angular distributions: *_theta_distribution.txt
# --------------------------------------------------------------------------- #
def plot_theta(path, times=None, ylog=False, xlim=(0, 180), ylim=None, figsize=None, ydecades=6):
    """Angular distribution of particles, one curve per time."""
    apply_style()
    mt = read_multi_time(path)
    tt, y = _select_times(mt, times)
    require_data(y)
    fig, ax = plt.subplots(figsize=figsize or SINGLE)
    colours, norm = _plot_curves(ax, mt.x, y, tt, ylog)
    ax.set_xlabel(r"Polar angle $\theta$ (deg)")
    ax.set_ylabel(f"{labels.multi_body(path)} angular distribution (arb. units)")
    ax.set_xlim(*xlim)
    ax.set_xticks(np.arange(0, 181, 30))
    if ylog:
        ax.set_yscale("log")
        _set_log_ylim(ax, y, ydecades)
    else:
        ax.set_ylim(bottom=0)
    if ylim:
        ax.set_ylim(*ylim)
    _finish_multi_time(fig, ax, tt, colours, norm, None)
    return fig


# --------------------------------------------------------------------------- #
#  Dispatcher, saving, folders
# --------------------------------------------------------------------------- #
_PLOTTERS = {
    "timeseries": plot_time_series, "radial": plot_radial,
    "spectrum": plot_spectrum, "theta": plot_theta,
    "mfp_electron": mp.plot_mfp, "mfp_hole": mp.plot_mfp, "mfp_photon": mp.plot_mfp,
    "emfp": mp.plot_emfp, "dos": mp.plot_dos,
    "ion_imfp": mp.plot_ion_imfp, "ion_dedx": mp.plot_ion_dedx,
    "ion_zeff": mp.plot_ion_zeff, "ion_range": mp.plot_ion_range,
}

# keyword arguments understood by each kind (others are silently ignored in folder mode)
_ACCEPTED = {
    "timeseries": {"hide_empty", "xlog", "columns", "figsize"},
    "radial": {"times", "xlog", "ylog", "ydecades", "xlim", "ylim", "figsize"},
    "spectrum": {"times", "xlog", "ylog", "ydecades", "xlim", "ylim", "figsize"},
    "theta": {"times", "ylog", "xlim", "ylim", "figsize", "ydecades"},
    "mfp_electron": {"shell_names", "xlog", "ylog", "xlim", "ylim", "ydecades", "figsize"},
    "emfp": {"xlog", "ylog", "xlim", "ylim", "ydecades", "figsize"},
    "dos": {"xlim", "figsize"},
    "ion_imfp": {"shell_names", "xlog", "ylog", "xlim", "ylim", "ydecades", "figsize"},
    "ion_dedx": {"shell_names", "xlog", "ylog", "xlim", "ylim", "ydecades", "figsize"},
    "ion_zeff": {"xlog", "xlim", "ylim", "figsize"},
    "ion_range": {"xlog", "ylog", "xlim", "ylim", "figsize"},
}
_ACCEPTED["mfp_hole"] = _ACCEPTED["mfp_photon"] = _ACCEPTED["mfp_electron"]


def accepted_kwargs(kind, kwargs):
    return {k: v for k, v in kwargs.items() if k in _ACCEPTED[kind]}


def plot_file(path, **kwargs):
    """Plot any known TREKIS-3 output file; returns the matplotlib Figure.
    Keyword arguments are passed to the plot function of that file type."""
    kind = detect_kind(path)
    if kind is None:
        raise ValueError(f"Unknown TREKIS-3 file type: {Path(path).name}")
    return _PLOTTERS[kind](path, **kwargs)


def save_figure(fig, path, out_dir=None, formats=("png", "pdf"), dpi=300):
    """Save fig next to the data file (or in out_dir) with the data file's base name."""
    path = Path(path)
    out_dir = Path(out_dir) if out_dir else path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    base = re.sub(r"[^\w\-]+", "_", path.stem).strip("_")
    saved = []
    for fmt in formats:
        target = out_dir / f"{base}.{fmt}"
        fig.savefig(target, dpi=dpi)
        saved.append(target)
    return saved


def _data_files(folder):
    return sorted(f for f in Path(folder).iterdir()
                  if f.is_file() and f.suffix.lower() in (".txt", ".dat") and detect_kind(f))


def plot_folder(folder, out_dir=None, formats=("png", "pdf"), only=None, recursive=False,
                **kwargs):
    """Plot every recognised TREKIS-3 file in a folder (and, with recursive=True, in all
    sub-folders; 'diff_CS' and 'figures' folders are skipped).

    out_dir: None -> a 'figures' sub-folder next to the data of each folder;
             otherwise this folder (recursive: the sub-folder structure is mirrored in it).
    only: list of words; only files whose name contains one of them are plotted.
    """
    folder = Path(folder)
    saved = []
    if recursive:
        dirs = []
        for root, subdirs, _ in os.walk(folder):
            subdirs[:] = sorted(d for d in subdirs if d not in _SKIP_DIRS)
            dirs.append(Path(root))
    else:
        dirs = [folder]
    for d in dirs:
        if out_dir is None:
            target = d / "figures"
        else:
            target = Path(out_dir) / d.relative_to(folder) if recursive else Path(out_dir)
        for f in _data_files(d):
            if only and not any(s.lower() in f.name.lower() for s in only):
                continue
            kind = detect_kind(f)
            fig = None
            try:
                fig = _PLOTTERS[kind](f, **accepted_kwargs(kind, kwargs))
                saved += save_figure(fig, f, target, formats)
            except Exception as exc:                   # report and keep going with the other files
                print(f"WARNING: skipped {f}: {exc}")
            finally:
                if fig is not None:
                    plt.close(fig)
    return saved
