"""Publication-quality matplotlib style and colour handling for TREKIS-3 plots."""
import matplotlib as mpl
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize

# Figure widths in inches (typical journal single / double column)
SINGLE = (3.4, 2.7)
DOUBLE = (7.0, 3.0)

# Okabe-Ito colour-blind-safe palette (for different quantities)
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
           "#56B4E9", "#F0E442", "#000000"]
LINESTYLES = ["-", "--", "-.", ":"]

# Sequential colour map for ordering by time (dark -> light, avoids the washed-out ends)
_base = mpl.colormaps["viridis"]
TIME_CMAP = LinearSegmentedColormap.from_list("trekis_time", _base(np.linspace(0.0, 0.88, 256)))


def apply_style():
    """Set global rcParams: clean look, inward ticks on all sides, no grid."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.labelsize": 9,
        "axes.linewidth": 0.8,
        "axes.grid": False,
        "axes.prop_cycle": mpl.cycler(color=PALETTE),
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.major.size": 3.5, "ytick.major.size": 3.5,
        "xtick.minor.size": 2.0, "ytick.minor.size": 2.0,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.minor.width": 0.6, "ytick.minor.width": 0.6,
        "lines.linewidth": 1.3,
        "lines.markersize": 3.5,
        "legend.frameon": False,
        "legend.fontsize": 7.5,
        "legend.handlelength": 1.8,
        "legend.borderaxespad": 0.4,
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
        "pdf.fonttype": 42,       # editable text in PDF
        "svg.fonttype": "none",
        "mathtext.default": "regular",
    })


def time_norm(times):
    """Normalisation for colouring by time: logarithmic if times span > 1 decade."""
    t = np.asarray(times, float)
    if len(t) > 1 and t.min() > 0 and t.max() / t.min() > 10:
        return LogNorm(t.min(), t.max())
    return Normalize(t.min(), t.max() if t.max() > t.min() else t.min() + 1)


def time_colours(times):
    norm = time_norm(times)
    return [TIME_CMAP(norm(t)) for t in times], norm


def time_label(t):
    """0.1 -> '0.1 fs', 1500 -> '1.5 ps'."""
    if t >= 1000:
        return f"{t / 1000:g} ps"
    return f"{t:g} fs"
