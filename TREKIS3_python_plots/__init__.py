"""
! This file is part of TREKIS-3
TREKIS3_python_plots
====================
Modern, publication-quality plots of TREKIS-3 output files with matplotlib
(numpy and matplotlib are the only requirements).

Quick use (from Python):
    from TREKIS3_python_plots import plot_file, save_figure
    fig = plot_file("Radial_electron_energy[eV_A^-3].txt", times=[0.1, 1, 10, 100])
    save_figure(fig, "Radial_electron_energy[eV_A^-3].txt", out_dir="figures")

Quick use (command line, from the folder that contains TREKIS3_python_plots):
    python -m TREKIS3_python_plots <file or run folder> [more files] [--out DIR] ...

Modules:
    readers          read the output files (by column names from the headers) and detect their kind
    plots            dynamics plots (time series, radial profiles, spectra, angles), dispatcher,
                     saving and walking through folders
    material_plots   mean free paths, density of states and ion data
    style            look of the figures (palette, fonts, sizes, legend boxes)
    labels           axis labels and titles from file and folder names
    common           small helpers shared by the plot modules
    __main__, run.py command line (python -m ..., and a launcher that works from any folder)
"""
__version__ = "0.2"

from .material_plots import (plot_dos, plot_emfp, plot_ion_dedx, plot_ion_imfp, plot_ion_range,
                             plot_ion_zeff, plot_mfp)
from .plots import (plot_file, plot_folder, plot_radial, plot_spectrum, plot_theta,
                    plot_time_series, save_figure)
from .readers import detect_kind, read_multi_time, read_time_series
from .style import apply_style

__all__ = ["plot_file", "plot_folder", "plot_radial", "plot_spectrum", "plot_theta",
           "plot_time_series", "plot_mfp", "plot_emfp", "plot_dos", "plot_ion_imfp",
           "plot_ion_dedx", "plot_ion_zeff", "plot_ion_range", "save_figure", "detect_kind",
           "read_multi_time", "read_time_series", "apply_style"]
