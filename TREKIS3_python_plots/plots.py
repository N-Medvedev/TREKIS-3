"""
! This file is part of TREKIS-3
Plots of the dynamics of a TREKIS-3 run, the plot dispatcher and the folder walker.

Plot functions (each takes the path of a data file and returns a matplotlib Figure):

    plot_time_series   Total_energies.txt, Total_numbers.txt    quantities vs. time
    plot_radial        Radial_*.txt                              radial profiles, one curve per time
    plot_spectrum      *_distribution_vs_E[1_eV].txt             energy spectra, one curve per time
    plot_theta         *_theta_distribution.txt                  angular distributions, one curve per time

Common keyword arguments (where they make sense):
    times      only the time steps closest to these times [fs]
    xlog, ylog logarithmic axes
    xlim, ylim axis limits
    ydecades   dynamic range shown on a logarithmic y axis
    figsize    figure size in inches (default: style.SINGLE)
    title      True (automatic title and sub-title), False (none) or a string

Other functions:
    plot_file(path, **kw)    plot any recognised file (kind detected from the file name)
    save_figure(fig, path)   save as png / pdf / ... next to the data file or in a folder
    plot_folder(folder)      plot all recognised files of a folder, optionally recursively;
                             empty files are skipped with a warning
    accepted_kwargs(kind, kw) the keyword arguments understood by a file kind

The mean-free-path, DOS and ion plots are in material_plots.py.
"""
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable

from . import labels
from . import material_plots as mp
from .common import add_title, require_data, set_ylabel_fit
from .readers import _SKIP_DIRS, detect_kind, read_multi_time, read_time_series
from .style import (LINESTYLES, PALETTE, SINGLE, TIME_CMAP, apply_style, time_colours,
                    time_label)

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


def _positive(y):
    """Copy of y with non-positive values replaced by NaN (for log axes)."""
    y = np.array(y, float)
    y[y <= 0] = np.nan
    return y


def _finish_multi_time(fig, ax, times, norm):
    """Legend (few times) or colour bar (many times)."""
    if len(times) <= MAX_LEGEND_TIMES:
        ax.legend(loc="best", title="Time")
    else:
        sm = ScalarMappable(norm=norm, cmap=TIME_CMAP)
        cb = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.06)
        cb.set_label("Time (fs)")
        cb.outline.set_visible(False)
        cb.ax.tick_params(direction="out")
        cb.ax.grid(False)


def _plot_curves(ax, x, y, times, ylog, marker="o"):
    """One curve per time. Small markers keep isolated non-zero points visible
    (e.g. holes that exist in one radial bin only)."""
    colours, norm = time_colours(times)
    for j, t in enumerate(times):
        yy = _positive(y[:, j]) if ylog else y[:, j]
        ax.plot(x, yy, color=colours[j], label=time_label(t), marker=marker, ms=2.2, mew=0)
    return colours, norm


def _set_log_ylim(ax, y, decades):
    """y limits for a log axis: from the maximum down `decades` decades (at least two)."""
    pos = y[y > 0]
    if pos.size == 0:
        return
    ymax = pos.max()
    ymin = max(pos.min(), ymax * 10.0 ** (-decades))
    ymin = min(ymin, ymax / 100.0)             # always show at least two decades
    ax.set_ylim(ymin / 1.5, ymax * 3)


def _xlim_nonzero(x, y, pad=1.25):
    """x range over which at least one curve is non-zero (at least one decade wide)."""
    nz = np.where((y > 0).any(axis=1))[0]
    if nz.size == 0:
        return x[0], x[-1]
    lo = x[0] / pad if x[0] > 0 else x[0]
    hi = max(x[nz[-1]] * pad, x[0] * 10)
    return lo, hi


# --------------------------------------------------------------------------- #
#  Time series: Total_numbers.txt, Total_energies.txt
# --------------------------------------------------------------------------- #
def _load_extra_column(base_ts, companion_path, col_name="Energy[eV]"):
    """Load `col_name` from `companion_path` aligned to `base_ts.time`.
    Returns (name, unit, values) or None if unavailable."""
    companion_path = Path(companion_path)
    if not companion_path.is_file():
        return None
    extra_ts = read_time_series(companion_path)

    target_base = col_name.split("[", 1)[0].strip()
    idx = None
    for i, name in enumerate(extra_ts.names):
        if name == col_name or name.strip() == target_base:
            idx = i
            break
    if idx is None:
        return None

    values = np.asarray(extra_ts.data[:, idx], float)
    if not np.array_equal(extra_ts.time, base_ts.time):
        values = np.interp(base_ts.time, extra_ts.time, values)

    unit = extra_ts.units[idx] if idx < len(extra_ts.units) and extra_ts.units[idx] else ""
    if not unit and "[" in col_name and "]" in col_name:
        unit = col_name.split("[", 1)[1].split("]", 1)[0].strip()

    return extra_ts.names[idx], unit, values


def plot_time_series(path, hide_empty=True, xlog=True, columns=None, figsize=None, title=True):
    """Quantities versus time. Energy columns (unit eV) and particle numbers are
    put in separate panels. All-zero columns are hidden unless hide_empty=False;
    `columns` restricts the plot to the given column names. For Total_energies.txt the
    total energy of a Total_numbers.txt in the same folder is added."""
    apply_style()
    path = Path(path)
    ts = read_time_series(path)

    series_list = [(ts.names[i], ts.units[i], ts.data[:, i]) for i in range(1, len(ts.names))]

    if path.name.lower() == "total_energies.txt":
        extra = _load_extra_column(ts, path.with_name("Total_numbers.txt"), "Energy[eV]")
        if extra is not None and all(s[0] != extra[0] for s in series_list):
            series_list.append(extra)

    energy, number = [], []
    for name, unit, data in series_list:
        if columns is not None and name not in columns:
            continue
        if hide_empty and not np.any(data != 0):
            continue
        (energy if unit.lower() == "ev" else number).append((name, data))

    groups = [(g, lab) for g, lab in ((number, "Number of particles"),
                                      (energy, "Energy (eV)")) if g]
    if not groups:
        raise ValueError(f"Nothing to plot in {path}")

    fig, axes = plt.subplots(len(groups), 1, sharex=True, figsize=figsize or SINGLE,
                             squeeze=False, layout="constrained")
    t = ts.time
    pos = t > 0 if xlog else np.ones_like(t, bool)

    for ax, (items, ylabel) in zip(axes[:, 0], groups):
        for n, (name, data) in enumerate(items):
            ax.plot(t[pos], data[pos], color=PALETTE[n % len(PALETTE)],
                    ls=LINESTYLES[(n // len(PALETTE)) % len(LINESTYLES)],
                    label=labels.series_label(name))
        set_ylabel_fit(ax, ylabel, fig=fig)
        if all(data[pos].min() >= 0 for _, data in items):
            ax.set_ylim(bottom=0)
        top = max(data[pos].max() for _, data in items)
        if top > 0:
            # leave room above the curves for the boxed legend (rows of 3 entries)
            rows = (len(items) + 2) // 3
            axes_h = (figsize or SINGLE)[1] * 0.62 / len(groups)       # approx. axes height, in
            frac = min(0.6, (0.17 * rows + 0.14) / axes_h)
            ax.set_ylim(top=top / (1.0 - frac))
        ax.legend(loc="upper center", ncol=min(3, len(items)), columnspacing=1.0,
                  handletextpad=0.5)

    ax = axes[-1, 0]
    if xlog:
        ax.set_xscale("log")
    ax.set_xlim(t[pos].min(), t.max())
    ax.set_xlabel("Time (fs)")
    main = ("Energy balance" if path.name.lower().startswith("total_energies")
            else "Particle numbers and energy")
    add_title(axes[0, 0], path, main, title)
    return fig


# --------------------------------------------------------------------------- #
#  Radial profiles: Radial_*.txt
# --------------------------------------------------------------------------- #
def plot_radial(path, times=None, xlog=True, ylog=True, ydecades=8, xlim=None, ylim=None,
                figsize=None, title=True):
    """Radial profile versus distance from the ion trajectory, one curve per time."""
    apply_style()
    mt = read_multi_time(path)
    tt, y = _select_times(mt, times)
    require_data(y, allow_negative=not ylog)
    fig, ax = plt.subplots(figsize=figsize or SINGLE, layout="constrained")
    colours, norm = _plot_curves(ax, mt.x, y, tt, ylog)
    ax.set_xlabel(f"Radius $\\mathit{{r}}$ ({labels.ANGSTROM})")
    set_ylabel_fit(ax, labels.radial_ylabel(path), fig=fig)
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
        _set_log_ylim(ax, y, ydecades)
    ax.set_xlim(*(xlim or _xlim_nonzero(mt.x, y)))
    if ylim:
        ax.set_ylim(*ylim)
    _finish_multi_time(fig, ax, tt, norm)
    add_title(ax, path, f"Radial profile: {labels.radial_body(path)}", title)
    return fig


# --------------------------------------------------------------------------- #
#  Energy spectra: *_distribution_vs_E[1_eV].txt
# --------------------------------------------------------------------------- #
def plot_spectrum(path, times=None, xlog=True, ylog=True, ydecades=8, xlim=None, ylim=None,
                  figsize=None, title=True):
    """Energy distribution of particles, one curve per time."""
    apply_style()
    mt = read_multi_time(path)
    tt, y = _select_times(mt, times)
    require_data(y, allow_negative=not ylog)
    fig, ax = plt.subplots(figsize=figsize or SINGLE, layout="constrained")
    colours, norm = _plot_curves(ax, mt.x, y, tt, ylog)
    ax.set_xlabel("Energy $\\mathit{E}$ (eV)")
    set_ylabel_fit(ax, "Spectrum (arb. units)", fig=fig)
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
        _set_log_ylim(ax, y, ydecades)
    ax.set_xlim(*(xlim or _xlim_nonzero(mt.x, y)))
    if ylim:
        ax.set_ylim(*ylim)
    _finish_multi_time(fig, ax, tt, norm)
    add_title(ax, path, f"{labels.multi_body(path)} energy spectrum", title)
    return fig


# --------------------------------------------------------------------------- #
#  Angular distributions: *_theta_distribution.txt
# --------------------------------------------------------------------------- #
def plot_theta(path, times=None, ylog=True, xlim=(0, 180), ylim=None, figsize=None,
               ydecades=6, title=True):
    """Angular distribution of particles, one curve per time."""
    apply_style()
    mt = read_multi_time(path)
    tt, y = _select_times(mt, times)
    require_data(y, allow_negative=not ylog)
    fig, ax = plt.subplots(figsize=figsize or SINGLE, layout="constrained")
    colours, norm = _plot_curves(ax, mt.x, y, tt, ylog)
    ax.set_xlabel(r"Polar angle $\theta$ (deg)")
    set_ylabel_fit(ax, "Angular distribution (arb. units)", fig=fig)
    ax.set_xlim(*xlim)
    ax.set_xticks(np.arange(0, 181, 30))
    if ylog:
        ax.set_yscale("log")
        _set_log_ylim(ax, y, ydecades)
    else:
        ax.set_ylim(bottom=0)
    if ylim:
        ax.set_ylim(*ylim)
    _finish_multi_time(fig, ax, tt, norm)
    add_title(ax, path, f"{labels.multi_body(path)} angular distribution", title)
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
    "timeseries": {"hide_empty", "xlog", "columns", "figsize", "title"},
    "radial": {"times", "xlog", "ylog", "ydecades", "xlim", "ylim", "figsize", "title"},
    "spectrum": {"times", "xlog", "ylog", "ydecades", "xlim", "ylim", "figsize", "title"},
    "theta": {"times", "ylog", "xlim", "ylim", "figsize", "ydecades", "title"},
    "mfp_electron": {"shell_names", "xlog", "ylog", "xlim", "ylim", "ydecades", "figsize", "title"},
    "emfp": {"xlog", "ylog", "xlim", "ylim", "ydecades", "figsize", "title"},
    "dos": {"xlim", "figsize", "title"},
    "ion_imfp": {"shell_names", "xlog", "ylog", "xlim", "ylim", "ydecades", "figsize", "title"},
    "ion_dedx": {"shell_names", "xlog", "ylog", "xlim", "ylim", "ydecades", "figsize", "title"},
    "ion_zeff": {"xlog", "xlim", "ylim", "figsize", "title"},
    "ion_range": {"xlog", "ylog", "xlim", "ylim", "figsize", "title"},
}
_ACCEPTED["mfp_hole"] = _ACCEPTED["mfp_photon"] = _ACCEPTED["mfp_electron"]


def accepted_kwargs(kind, kwargs):
    """The subset of `kwargs` that the plot function of this file kind understands."""
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
        fig.savefig(target, dpi=dpi, bbox_inches="tight", pad_inches=0.05)
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
    Files without non-zero data, or that fail for another reason, are skipped with a
    warning; the other files are still plotted.
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
