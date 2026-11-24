"""
! This file is part of TREKIS-3
Axis labels and legend texts derived from TREKIS-3 file names and headers.
"""
import re
from pathlib import Path

ANGSTROM = "\u00c5"

_QTY = {
    "energy": "energy density",
    "kin_energy": "kinetic energy density",
    "pot_energy": "potential energy density",
    "density": "density",
    "temperature": "temperature",
}
_QTY_DEFAULT_UNIT = {
    "energy": f"eV {ANGSTROM}$^{{-3}}$",
    "kin_energy": f"eV {ANGSTROM}$^{{-3}}$",
    "pot_energy": f"eV {ANGSTROM}$^{{-3}}$",
    "density": "cm$^{-3}$",
    "temperature": "K",
}
_RADIAL_RE = re.compile(
    r"^Radial_(?P<body>.+?)_(?P<qty>kin_energy|pot_energy|energy|density|temperature)(?P<unit>.*)$")

_SERIES_NAMES = {            # column names of Total_numbers.txt
    "Ne": "Electrons",
    "Ne_Emitted": "Emitted electrons",
    "N_photons": "Photons",
    "Energy": "Total energy",
    "Energy_Emitted": "Emitted energy",
    "Valence": "Valence band",
}


def _capitalise(s):
    return s[:1].upper() + s[1:]


def pretty_unit(unit_text, qty=None):
    """Unit from a file-name fragment ('[eV_A^-3]', '_eV_A_-3_', '[1_cm^-3]', '[K]')."""
    u = unit_text.lower()
    if "ev_a" in u:
        return f"eV {ANGSTROM}$^{{-3}}$"
    if "cm" in u:
        return "cm$^{-3}$"
    if re.fullmatch(r"[\[_]*k[\]_]*", u):
        return "K"
    return _QTY_DEFAULT_UNIT.get(qty, "")


def radial_ylabel(path):
    """'Radial_Si_L23-shell_holes_energy[eV_A^-3].txt' ->
       'Si L23-shell holes energy density (eV Å$^{-3}$)'."""
    stem = Path(path).stem
    m = _RADIAL_RE.match(stem)
    if not m:
        return "Value"
    body = _capitalise(m.group("body").replace("_", " "))
    qty = m.group("qty")
    unit = pretty_unit(m.group("unit"), qty)
    return f"{body} {_QTY[qty]} ({unit})"


def series_label(name):
    """Legend text of a time-series column ('Si_L23-shell' -> 'Si L23-shell')."""
    return _SERIES_NAMES.get(name, name.replace("_", " "))


def multi_body(path):
    """Subject of a spectrum / angular distribution file name."""
    stem = Path(path).stem
    m = re.match(r"^(?P<body>.+?)_(?:theta_distribution|distribution_vs_E)", stem)
    body = m.group("body") if m else stem
    body = body.replace("_", " ")
    body = {"VB holes": "Valence-band holes", "Electron": "Electron", "Electrons": "Electron"}.get(body, body)
    return _capitalise(body)
