#!/usr/bin/env python3
"""
! This file is part of TREKIS-3
Example: plot a TREKIS-3 run folder, or single files, with custom settings.
Edit the CONFIG block and run:  python example_usage.py
"""
import sys
from pathlib import Path

# make the package importable when this script lives inside the package folder
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib.pyplot as plt
from TREKIS3_python_plots import plot_file, plot_folder, save_figure

# ----------------------------- CONFIG ------------------------------------- #
RUN_FOLDER = Path(".")                 # folder with the TREKIS-3 output files
OUT_DIR = RUN_FOLDER / "figures"       # where the figures are saved
FORMATS = ("png", "pdf")               # any matplotlib format
TIMES = [0.1, 1, 10, 100]              # fs; None = all time steps in the files
# -------------------------------------------------------------------------- #

# 1) everything that is recognised in the folder
plot_folder(RUN_FOLDER, OUT_DIR, FORMATS, times=TIMES)

# 2) or a single file with its own settings
# fig = plot_file(RUN_FOLDER / "Radial_electron_energy[eV_A^-3].txt",
#                 times=[1, 10, 100], xlim=(1, 300), ylim=(1e-5, 20))
# save_figure(fig, RUN_FOLDER / "Radial_electron_energy[eV_A^-3].txt", OUT_DIR, FORMATS)
# plt.show()
