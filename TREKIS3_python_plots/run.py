#!/usr/bin/env python3
"""
! This file is part of TREKIS-3
Launcher for TREKIS3_python_plots: works from any working directory, without installing
the package. This is what TREKIS-3 calls (see FORTRAN_INTEGRATION.md):

    python TREKIS3_python_plots/run.py <run folder / files> [options]

Same arguments as  python -m TREKIS3_python_plots  (see README.md).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))   # folder containing the package

from TREKIS3_python_plots.__main__ import main

if __name__ == "__main__":
    main()
