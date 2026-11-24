"""
! This file is part of TREKIS-3
Axis labels, legend texts and figure titles derived from TREKIS-3 file and folder names.

File names carry the information about what is plotted, e.g.
    Radial_Si_L23-shell_holes_energy[eV_A^-3].txt  ->  y label 'Si L23-shell holes energy
                                                       density (eV Å^-3)'
and the folder names carry the context of the run, e.g.
    OUTPUT_Al2O3/OUTPUT_Ar_in_Al2O3/Ar_E_1.00_MeV_150.00_fs/  ->  'Ar 1 MeV in Al2O3'
which is shown as the sub-title of the figures.

Functions:
    radial_body(path), radial_ylabel(path)   subject and y label of a Radial_*.txt file
    pretty_unit(text, qty)                   unit text from a file-name fragment
    series_label(name)                       legend text of a time-series column
    multi_body(path)                         subject of a spectrum / angular-distribution file
    pretty_formula(text)                     'Al2O3' -> 'Al$_2$O$_3$' (mathtext subscripts)
    run_context(path)                        ion / energy / material from the folder names
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

_SERIES_NAMES = {            # column names of Total_numbers.txt and Total_energies.txt
    "Ne": "Electrons",
    "Ne_Emitted": "Emitted electrons",
    "N_photons": "Photons",
    "Energy": "Total energy",
    "Energy_Emitted": "Emitted energy",
    "Valence": "Valence band",
}

# folder names: Ar_E_1.00_MeV_150.00_fs  /  OUTPUT_Ar_in_Al2O3  /  OUTPUT_Al2O3
_RUN_RE = re.compile(r"^(?P<ion>.+?)_E_(?P<E>[0-9.]+)_MeV_(?P<t>[0-9.]+)_fs")
_ION_RE = re.compile(r"^OUTPUT_(?P<ion>.+?)_in_(?P<mat>.+)$")
_MAT_RE = re.compile(r"^OUTPUT_(?P<mat>.+)$")


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


def radial_body(path):
    """Subject of a radial file: 'Radial_Si_L23-shell_holes_energy...' -> 'Si L23-shell holes'."""
    m = _RADIAL_RE.match(Path(path).stem)
    return _capitalise(m.group("body").replace("_", " ")) if m else Path(path).stem


def radial_ylabel(path):
    """'Radial_Si_L23-shell_holes_energy[eV_A^-3].txt' ->
       'Si L23-shell holes energy density (eV Å$^{-3}$)'."""
    m = _RADIAL_RE.match(Path(path).stem)
    if not m:
        return "Value"
    qty = m.group("qty")
    unit = pretty_unit(m.group("unit"), qty)
    return f"{radial_body(path)} {_QTY[qty]} ({unit})"


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


def pretty_formula(text):
    """'Al2O3' -> 'Al$_{2}$O$_{3}$' (digits after a letter become subscripts)."""
    return re.sub(r"(?<=[A-Za-z)])(\d+)", r"$_{\1}$", text)


def run_context(path):
    """Context of a file from the names of the folders it lies in, for the figure sub-title:
        .../OUTPUT_Al2O3/OUTPUT_Ar_in_Al2O3/Ar_E_1.00_MeV_150.00_fs/<file>  ->  'Ar 1 MeV in Al2O3'
        .../OUTPUT_Al2O3/OUTPUT_Ar_in_Al2O3/<file>                          ->  'Ar in Al2O3'
        .../OUTPUT_Al2O3/<file>                                             ->  'Al2O3'
    Returns '' if the folder names do not follow this scheme."""
    p = Path(path).resolve()
    folder = p if p.is_dir() else p.parent
    ion = energy = mat = None
    for d in [folder] + list(folder.parents)[:3]:
        m = _RUN_RE.match(d.name)
        if m and ion is None:
            ion, energy = m.group("ion"), float(m.group("E"))
            continue
        m = _ION_RE.match(d.name)
        if m:
            ion = ion or m.group("ion")
            mat = m.group("mat")
            break
        m = _MAT_RE.match(d.name)
        if m:
            mat = m.group("mat")
            break
    mat = pretty_formula(mat) if mat else None
    if ion and energy is not None:
        return f"{ion} {energy:g} MeV" + (f" in {mat}" if mat else "")
    if ion and mat:
        return f"{ion} in {mat}"
    return mat or ""
