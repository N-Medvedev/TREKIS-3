"""
! This file is part of TREKIS-3

Readers for TREKIS-3 output files.

All files are self-describing through their '#' header lines, so columns are
looked up by name and never by position:

  * time series  (Total_numbers.txt, Total_energies.txt)
        #Time[fs]  Electrons[eV]  Atoms[eV] ...
  * "x versus time" tables (Radial_*.txt, *_theta_distribution.txt,
    *_distribution_vs_E[1_eV].txt)
        #Radius[A]  0.01[fs]  0.10[fs] ...                  (time in the header)
        #Energy[eV]  Spectrum[arb.units]
        #time[fs]:  0.01  0.10 ...                          (spectra: 2 header lines)

  * headerless tables (OUTPUT_*_IMFPs*.dat, ion *_IMFP / _dEdx / _effective_charges /
    _Range files): plain numeric columns (energy | shells | total); the Range file has two
    header lines '# Energy dEdx Range' / '# [eV] [eV/A] [A]'
  * density of states (OUTPUT_<material>_DOS_analysis.dat): one header line
        #E[eV]  k[1/m]  DOS[a.u.]  Int_DOS[a.u.]  Eff_mass[me]

Main functions and classes:
    read_time_series(path) -> TimeSeries     columns by name (Total_numbers, Total_energies)
    read_multi_time(path)  -> MultiTime      x and y(x) at several times (radial, spectra, angles)
    read_table(path)       -> (data, header) plain numeric table
    read_dos(path)         -> (columns, units) dictionaries
    find_shell_names(path)                   shell names from a nearby Total_energies.txt
    detect_kind(path)                        kind of a file from its name (None if unknown)

Numbers written by Fortran without the exponent letter (1.0-100) are understood.

File names may contain either brackets (as written by TREKIS-3) or underscores
(as produced when files are copied or downloaded), both are handled.
"""
import os
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np


# --------------------------------------------------------------------------- #
#  Helpers
# --------------------------------------------------------------------------- #
def split_name_unit(token):
    """'Energy[eV]' -> ('Energy', 'eV');  'Time' -> ('Time', '')."""
    m = re.match(r"^([^\[\]]*)(?:\[([^\]]*)\])?:?$", token.strip())
    if not m:
        return token.strip(), ""
    return m.group(1).strip(), (m.group(2) or "").strip()


def _header_lines(path):
    """Header ('#') lines of a file, without the leading '#' and line ends."""
    lines = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                lines.append(s.lstrip("#").strip())
            else:
                break
    return lines


_FORTRAN_EXP = re.compile(r"(?<=\d)([+-]\d{3})(?=\s|$)")   # 1.0-100 -> 1.0E-100


def _load_numeric(path):
    """Numeric table; tolerates Fortran numbers without 'E' (1.0-100) and NaN/Infinity."""
    try:
        data = np.loadtxt(path, comments="#", ndmin=2)
    except ValueError:
        rows = []
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                rows.append([float(v) for v in _FORTRAN_EXP.sub(r"E\1", line).split()])
        n = max(len(r) for r in rows)
        data = np.array([r + [np.nan] * (n - len(r)) for r in rows])
    if data.size == 0:
        raise ValueError(f"No numeric data in {path}")
    return data


_FS_TOKEN = re.compile(r"(?<![A-Za-z])([-+]?\d[\d.]*(?:[eE][-+]?\d+)?)\s*\[fs\]")
_NUMBER = re.compile(r"[-+]?\d[\d.]*(?:[eE][-+]?\d+)?")


def _parse_times(header, ncols):
    """Times [fs] from the header; None if they cannot be found."""
    for line in header:                                  # '#time[fs]:  0.01  0.10 ...'
        if line.lower().startswith("time[fs]"):
            vals = [float(v) for v in _NUMBER.findall(line.split(":", 1)[-1])]
            if len(vals) == ncols:
                return np.array(vals)
    for line in header:                                  # '#Radius[A]  0.01[fs]  0.10[fs] ...'
        vals = [float(v) for v in _FS_TOKEN.findall(line)]
        if len(vals) == ncols:
            return np.array(vals)
    return None


# --------------------------------------------------------------------------- #
#  Data containers
# --------------------------------------------------------------------------- #
@dataclass
class TimeSeries:
    """Columns versus time: Total_numbers.txt, Total_energies.txt."""
    names: list          # column names, names[0] is the time column
    units: list
    data: np.ndarray     # shape (n_rows, n_columns)

    @property
    def time(self):
        return self.data[:, 0]

    def column(self, name):
        return self.data[:, self.names.index(name)]


@dataclass
class MultiTime:
    """A quantity y(x) given at several times."""
    x: np.ndarray        # (nx,)
    y: np.ndarray        # (nx, nt)
    times: np.ndarray    # (nt,) [fs]
    xname: str
    xunit: str
    yname: str = ""
    yunit: str = ""


# --------------------------------------------------------------------------- #
#  Readers
# --------------------------------------------------------------------------- #
def read_time_series(path):
    header = _header_lines(path)
    data = _load_numeric(path)
    ncol = data.shape[1]
    tokens = header[0].split() if header else []
    if len(tokens) != ncol:
        tokens = ["Time[fs]"] + [f"col{i}" for i in range(1, ncol)]
    names, units = zip(*(split_name_unit(t) for t in tokens))
    return TimeSeries(list(names), list(units), data)


def read_multi_time(path):
    """Radial profiles, spectra and angular distributions (x + one column per time)."""
    header = _header_lines(path)
    data = _load_numeric(path)
    x, y = data[:, 0], data[:, 1:]
    times = _parse_times(header, y.shape[1])
    if times is None:
        times = np.arange(y.shape[1], dtype=float)
    first = header[0].split() if header else ["x"]
    xname, xunit = split_name_unit(first[0])
    yname = yunit = ""
    if len(first) > 1 and "[fs]" not in first[1]:
        yname, yunit = split_name_unit(first[1])
    return MultiTime(x, y, times, xname, xunit, yname, yunit)


def read_table(path):
    """Plain numeric table: (data, header_lines). Used for the headerless MFP files."""
    return _load_numeric(path), _header_lines(path)


def read_dos(path):
    """OUTPUT_<material>_DOS_analysis.dat: dict of column name -> array."""
    header = _header_lines(path)
    data = _load_numeric(path)
    tokens = header[0].split() if header else []
    if len(tokens) != data.shape[1]:
        tokens = [f"col{i}" for i in range(data.shape[1])]
    return {split_name_unit(t)[0]: data[:, i] for i, t in enumerate(tokens)}, \
           {split_name_unit(t)[0]: split_name_unit(t)[1] for t in tokens}


# --------------------------------------------------------------------------- #
#  Shell names (the MFP files themselves contain no column names)
# --------------------------------------------------------------------------- #
_SKIP_DIRS = {"diff_CS", "figures"}
_NOT_SHELLS = {"time", "electrons", "atoms", "field", "photons"}


def _find_total_energies(base, max_depth=3):
    base = Path(base)
    base_depth = len(base.parts)
    for root, dirs, files in os.walk(base):
        dirs[:] = sorted(d for d in dirs if d not in _SKIP_DIRS)
        if len(Path(root).parts) - base_depth >= max_depth:
            dirs[:] = []
        if "Total_energies.txt" in files:
            return Path(root) / "Total_energies.txt"
    return None


def find_shell_names(start):
    """Shell names (in file column order) taken from the header of a Total_energies.txt
    of a run in the neighbourhood of `start` (a file or folder); None if not found.
    Assumes the shell order of Total_energies.txt equals that of the MFP files."""
    start = Path(start).resolve()
    d = start if start.is_dir() else start.parent
    for base in (d, d.parent, d.parent.parent):
        f = _find_total_energies(base)
        if f is not None:
            ts = read_time_series(f)
            return [n for n in ts.names if n.lower() not in _NOT_SHELLS]
    return None


# --------------------------------------------------------------------------- #
#  File kind detection
# --------------------------------------------------------------------------- #
def detect_kind(path):
    """Kind of a TREKIS-3 output file, or None if it is not recognised.

    dynamics:   'timeseries', 'radial', 'spectrum', 'theta'
    material:   'mfp_electron', 'mfp_hole', 'mfp_photon', 'emfp', 'dos'
    ion:        'ion_imfp', 'ion_dedx', 'ion_zeff', 'ion_range'
    """
    name = Path(path).name
    if name.startswith(("Total_numbers", "Total_energies")):
        return "timeseries"
    if name.startswith("Radial_"):
        return "radial"
    if "_theta_distribution" in name:
        return "theta"
    if "_distribution_vs_E" in name:
        return "spectrum"
    if name.startswith("OUTPUT_"):
        if name.endswith("_DOS_analysis.dat"):
            return "dos"
        if name.endswith("_effective_charges.dat"):
            return "ion_zeff"
        if name.endswith("_Range.dat"):
            return "ion_range"
        if name.endswith("_dEdx.dat"):
            return "ion_dedx"
        if name.endswith("_IMFP.dat"):
            return "ion_imfp"
        if name.startswith(("OUTPUT_Electron", "OUTPUT_Hole")) and \
                ("EMFPs" in name or "_elas" in name):
            return "emfp"
        if name.startswith("OUTPUT_Electron_IMFPs"):
            return "mfp_electron"
        if name.startswith("OUTPUT_Hole") and "IMFPs" in name:
            return "mfp_hole"
        if name.startswith("OUTPUT_Photon_IMFPs"):
            return "mfp_photon"
    return None
