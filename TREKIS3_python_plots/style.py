"""
! This file is part of TREKIS-3
Look of the TREKIS-3 figures: fonts, colours, sizes, axes, grid and legend boxes.

Everything that decides how the figures look is collected here, so that the look can be
changed in one place:

    SINGLE, DOUBLE   figure sizes in inches (journal single / double column width)
    PALETTE          colours for different quantities (colour-blind-safe, high contrast)
    TOTAL_COLOUR     colour of "total" curves (shells are drawn in PALETTE colours)
    TIME_CMAP        sequential colour map for curves ordered by time (early -> late)
    apply_style()    sets the matplotlib rcParams (called by every plot function)

The style is a clean, modern one: no top/right frame lines, outward ticks, a light grid,
semi-transparent boxed legends and a bold left-aligned title above each figure.

Helper functions for the time colouring (used for radial profiles, spectra and angular
distributions, where one curve is drawn per time step):

    time_norm(times)     matplotlib normalisation (logarithmic if the times span > 1 decade)
    time_colours(times)  list of colours (one per time) and the normalisation used
    time_label(t)        legend text for a time in fs, e.g. '0.1 fs' or '1.5 ps'
"""
import matplotlib as mpl
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize

# --------------------------------------------------------------------------- #
#  Sizes (inches)
# --------------------------------------------------------------------------- #
SINGLE = (3.4, 3.0)      # one journal column
DOUBLE = (7.0, 3.3)      # full page width

# --------------------------------------------------------------------------- #
#  Colours
# --------------------------------------------------------------------------- #
# High-contrast, colour-blind-safe (Paul Tol "vibrant" + gold and grey)
PALETTE = ["#0077BB", "#EE7733", "#009988", "#CC3311",
           "#33BBEE", "#EE3377", "#DDAA33", "#555555"]
LINESTYLES = ["-", "--", "-.", ":"]       # used after the palette has been used up

TOTAL_COLOUR = "#1A1A1A"                  # total over all shells
TEXT_COLOUR = "#222222"
CONTEXT_COLOUR = "#666666"                # sub-title (ion, energy, material)
GRID_COLOUR = "#DADADA"
EDGE_COLOUR = "#444444"

# Sequential map for time (dark purple -> orange); the very light end is cut off
TIME_CMAP = LinearSegmentedColormap.from_list(
    "trekis_time", mpl.colormaps["plasma"](np.linspace(0.02, 0.86, 256)))


def apply_style():
    """Set the global matplotlib rcParams of the TREKIS-3 look."""
    mpl.rcParams.update({
        # fonts and text
        "font.family": "sans-serif",
        "font.sans-serif": ["Inter", "Roboto", "Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "text.color": TEXT_COLOUR,
        "mathtext.default": "regular",
        # axes: only left and bottom lines, light grid below the data
        "axes.labelsize": 9,
        "axes.labelcolor": TEXT_COLOUR,
        "axes.edgecolor": EDGE_COLOUR,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "axes.axisbelow": True,
        "grid.color": GRID_COLOUR,
        "grid.linewidth": 0.6,
        "grid.linestyle": "-",
        "axes.prop_cycle": mpl.cycler(color=PALETTE),
        # title: bold, left-aligned
        "axes.titlesize": 9,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.titlecolor": TEXT_COLOUR,
        # ticks: outward, minor ticks only where matplotlib adds them (log axes)
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.color": EDGE_COLOUR, "ytick.color": EDGE_COLOUR,
        "xtick.labelcolor": TEXT_COLOUR, "ytick.labelcolor": TEXT_COLOUR,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0,
        "xtick.minor.size": 1.8, "ytick.minor.size": 1.8,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.minor.width": 0.6, "ytick.minor.width": 0.6,
        # lines
        "lines.linewidth": 1.0,
        "lines.markersize": 1.5,
        "lines.solid_capstyle": "round",
        # legend: semi-transparent rounded box
        "legend.frameon": True,
        "legend.fancybox": True,
        "legend.framealpha": 0.8,
        "legend.facecolor": "white",
        "legend.edgecolor": "#C8C8C8",
        "legend.fontsize": 7.5,
        "legend.title_fontsize": 7.5,
        "legend.borderpad": 0.2,
        "legend.labelspacing": 0.35,
        "legend.handlelength": 1.8,
        "legend.borderaxespad": 0.2,
        # figure and output
        "figure.dpi": 120,
        "figure.facecolor": "white",
        "savefig.dpi": 300,
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "pdf.fonttype": 42,       # editable text in PDF
        "svg.fonttype": "none",
    })


# --------------------------------------------------------------------------- #
#  Colouring by time
# --------------------------------------------------------------------------- #
def time_norm(times):
    """Normalisation for colouring by time: logarithmic if times span > 1 decade."""
    t = np.asarray(times, float)
    if len(t) > 1 and t.min() > 0 and t.max() / t.min() > 10:
        return LogNorm(t.min(), t.max())
    return Normalize(t.min(), t.max() if t.max() > t.min() else t.min() + 1)


def time_colours(times):
    """Colours (one per time, from TIME_CMAP) and the normalisation that was used."""
    norm = time_norm(times)
    return [TIME_CMAP(norm(t)) for t in times], norm


def time_label(t):
    """0.1 -> '0.1 fs', 1500 -> '1.5 ps'."""
    if t >= 1000:
        return f"{t / 1000:g} ps"
    return f"{t:g} fs"
