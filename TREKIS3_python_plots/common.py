"""
! This file is part of TREKIS-3
Helpers shared by the plot modules (plots.py and material_plots.py).

    require_data(y, what, allow_negative)
        raise ValueError if there is nothing to draw (all zero / missing), so that empty
        output files are skipped with a clear message instead of failing in matplotlib
    add_title(ax, path, main, title)
        bold title above the axes plus a small grey sub-title with the context of the
        run (ion, energy, material), taken from the folder names
    set_ylabel_fit(ax, text, fig)
        y label at full size, reduced in font size only if it would run off the figure
"""
import matplotlib.pyplot as plt
import numpy as np

from . import labels
from .style import CONTEXT_COLOUR


def require_data(y, what="data", allow_negative=False):
    """Raise ValueError if there is nothing finite and plottable to draw."""
    y = np.asarray(y, float)
    mask = np.isfinite(y) & (y != 0 if allow_negative else y > 0)
    if not np.any(mask):
        raise ValueError(f"no non-zero {what} (all values are zero, negative or missing), nothing to plot")


def add_title(ax, path, main, title=True):
    """Title above the axes: `main` in bold, and below it, in grey, the context of the run
    (e.g. 'Ar 1 MeV in Al2O3') found from the folder names of `path`.

    title: True -> automatic text; False/None -> no title; str -> this text instead of `main`
    (the context line is kept)."""
    if title is False or title is None:
        return
    if isinstance(title, str):
        main = title
    context = labels.run_context(path)
    ax.set_title(main, loc="left", pad=17 if context else 6)
    if context:
        ax.annotate(context, xy=(0, 1), xycoords="axes fraction", xytext=(0, 3.5),
                    textcoords="offset points", ha="left", va="bottom", fontsize=7.5,
                    color=CONTEXT_COLOUR, annotation_clip=False)


def _renderer(fig):
    try:
        return fig.canvas.get_renderer()
    except AttributeError:                      # backends without get_renderer()
        return fig._get_renderer()


def _measure_top_margin_px(fig, ax):
    """Pixels by which the top edge of the y label is inside (> 0) or above (< 0) the figure."""
    fig.canvas.draw()
    renderer = _renderer(fig)
    fig_bbox = fig.get_window_extent(renderer=renderer)
    ylabel_bbox = ax.yaxis.label.get_window_extent(renderer=renderer)
    return float(fig_bbox.y1 - ylabel_bbox.y1)


def set_ylabel_fit(ax, text, fig=None, min_pad_px=2.0, min_scale=0.6, step=0.05):
    """Keep the y-axis label on a single line at full font size unless its top edge
    would go past the upper side of the figure; then shrink it step by step (down to
    `min_scale` of the normal size)."""
    if fig is None:
        fig = ax.figure
    base_fs = float(plt.rcParams.get("axes.labelsize", 9))

    ax.set_ylabel(text, fontsize=base_fs)
    margin_px = _measure_top_margin_px(fig, ax)
    if margin_px >= min_pad_px:
        return

    best_fs, best_margin = base_fs, margin_px
    for scale in np.arange(1.0 - step, min_scale - 1e-9, -step):
        fs = base_fs * float(scale)
        ax.set_ylabel(text, fontsize=fs)
        cur_margin = _measure_top_margin_px(fig, ax)
        if cur_margin > best_margin:
            best_fs, best_margin = fs, cur_margin
        if cur_margin >= min_pad_px:
            return
    ax.set_ylabel(text, fontsize=best_fs)
