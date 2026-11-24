# TREKIS3_python_plots

Modern, publication-quality figures for TREKIS-3 output, made with matplotlib. It is an alternative
(or addition) to the gnuplot scripts written by `Gnuplotting_subs.f90`: the files are read by their
header names, the plot type is recognised from the file name, and titles, axis labels and units are
set automatically.

**Requirements:** Python 3.9+ with `numpy` and `matplotlib` (3.5 or newer). No installation: keep the
folder `TREKIS3_python_plots` where TREKIS-3 is run (next to `INPUT_PARAMETERS.txt`), or anywhere you like.

## Quick start

From the folder that contains `TREKIS3_python_plots/`:

    python -m TREKIS3_python_plots OUTPUT_Al2O3                 # the whole tree below this folder
    python -m TREKIS3_python_plots path\to\run_folder --no-recursive
    python -m TREKIS3_python_plots "Radial_electron_energy[eV_A^-3].txt" -t 0.1 1 10 100 -o figures

From any folder (this is also what TREKIS-3 calls, see *Calling it from TREKIS-3*):

    python path/to/TREKIS3_python_plots/run.py OUTPUT_Al2O3

From Python:

    from TREKIS3_python_plots import plot_file, save_figure
    fig = plot_file("Radial_electron_energy[eV_A^-3].txt", times=[1, 10, 100], xlim=(1, 300))
    save_figure(fig, "Radial_electron_energy[eV_A^-3].txt", out_dir="figures")

`example_usage.py` is a ready-made script with a CONFIG block.

## Expected folder structure

    OUTPUT_Al2O3/                        electron / hole / photon mean free paths, DOS
    |-- diff_CS/                         (ignored)
    `-- OUTPUT_Ar_in_Al2O3/              ion mean free path, energy loss, charge state, range
        `-- Ar_E_1.00_MeV_150.00_fs/     dynamics of one run (radial profiles, spectra, ...)

When a folder is given, it is searched recursively (use `--no-recursive` for that folder only).
Folders called `diff_CS` and `figures` are skipped. The figures of each folder are written to a
`figures` sub-folder next to its data (or to the folder given with `-o`, which then mirrors the
structure). A file without any non-zero data (e.g. holes of a shell that was never ionised) is
skipped with a warning, the other files are still plotted.

## Supported files

The kind of file is detected from its name; names with brackets (`...[eV_A^-3].txt`) and with
underscores (`..._eV_A_-3_.txt`, as after some downloads) both work.

| File | Plot | Default axes |
|---|---|---|
| `Total_energies.txt`, `Total_numbers.txt` | quantities vs. time; energies and particle numbers in separate panels; all-zero columns hidden; for `Total_energies.txt` the total energy from a `Total_numbers.txt` in the same folder is added | log x, linear y |
| `Radial_*.txt` (any particle or shell; density, energy, temperature) | radial profile, one curve per time | log-log |
| `*_distribution_vs_E[1_eV].txt` | energy spectrum, one curve per time | log-log |
| `*_theta_distribution.txt` | angular distribution, one curve per time | linear x (0-180 deg), log y |
| `OUTPUT_Electron_IMFPs*.dat`, `OUTPUT_Hole_IMFPs*.dat`, `OUTPUT_Photon_IMFPs*.dat` | inelastic mean free path (photons: attenuation length) vs. energy, per shell and total | log-log |
| `OUTPUT_Electron*EMFPs*.dat`, `OUTPUT_Hole*EMFPs*.dat` | elastic mean free path vs. energy | log-log |
| `OUTPUT_<material>_DOS_analysis.dat` | density of states, with the effective mass on a second axis | linear |
| `OUTPUT_<ion>*_IMFP.dat` | ion inelastic mean free path vs. ion energy (MeV), per shell and total | log-log |
| `OUTPUT_<ion>*_dEdx.dat` | electronic energy loss vs. ion energy (MeV), per shell and total | log x, linear y |
| `OUTPUT_<ion>*_effective_charges.dat` | effective ion charge vs. ion energy | log x, linear y |
| `OUTPUT_<ion>*_Range.dat` | ion range (micrometres) vs. ion energy (MeV) | log-log |

**Reading the files.** Dynamics files and the DOS file are read by the column names in their `#` headers,
never by column number. Numbers written by Fortran without the exponent letter (`1.0-100`) are
understood; mean free paths of `1e30` ("infinite") and non-positive values are not drawn.

**Shell names.** The mean-free-path and ion files have no header (columns: energy, one column per shell,
total). The shell names in their legends are taken from the header of a `Total_energies.txt` found near
the file (run folders contain one). If their number does not match the number of columns, generic names
("Shell 1", ...) are used and a note is printed. Pass `shell_names=[...]` to set them explicitly.

## Titles

Every figure gets a bold title above the axes and, if the folder names allow it, a grey line below it
with the context of the run:

    Radial profile: Electron
    Ar 1 MeV in Al2O3

The context comes from the folder names: `Ar_E_1.00_MeV_150.00_fs` (ion and energy),
`OUTPUT_Ar_in_Al2O3` (ion and material) and `OUTPUT_Al2O3` (material). Switch the titles off with
`--no-title` (Python: `title=False`) or replace the first line with `--title "text"` (`title="text"`).

## Command line options

`python -m TREKIS3_python_plots PATH [PATH ...] [options]`, where a path is a data file or a folder.

| Option | Meaning |
|---|---|
| `-o`, `--out DIR` | output folder (default: next to the data, in a `figures` sub-folder for folders) |
| `-f`, `--formats ...` | output formats, e.g. `png pdf svg` (default `png pdf`) |
| `--dpi N` | resolution of raster formats (default 300) |
| `-r`, `--no-recursive` | for folders: do not descend into sub-folders |
| `--only WORD ...` | for folders: only files whose name contains one of these words (e.g. `Radial IMFP`) |
| `-t`, `--times T ...` | only the time steps closest to these times in fs (radial, spectra, angles) |
| `--width single\|double` | figure size: one column (3.4 in) or full page width (7 in) |
| `--xlin`, `--ylin`, `--ylog` | linear x axis, linear y axis, logarithmic y axis |
| `--ydecades N` | dynamic range shown on a logarithmic y axis (default 8; 5 for mean free paths) |
| `--show-empty` | time series: also draw columns that are zero everywhere |
| `--no-title`, `--title TEXT` | no title / own title |
| `--show` | open the figures in a window |

## Python functions

All plot functions take the path of the data file and return a matplotlib `Figure`; `plot_file(path, **kw)`
chooses the right one. Keyword arguments (where they make sense): `times`, `xlog`, `ylog`, `xlim`, `ylim`,
`ydecades`, `figsize`, `title`, `shell_names`; unknown keywords are ignored when plotting folders.

| Function | Files |
|---|---|
| `plot_time_series(path, hide_empty=True, columns=None)` | `Total_energies.txt`, `Total_numbers.txt` |
| `plot_radial(path, times=None)` | `Radial_*.txt` |
| `plot_spectrum(path, times=None)` | `*_distribution_vs_E*.txt` |
| `plot_theta(path, times=None)` | `*_theta_distribution.txt` |
| `plot_mfp(path)`, `plot_emfp(path)` | electron / hole / photon inelastic and elastic mean free paths |
| `plot_dos(path)` | density of states |
| `plot_ion_imfp(path)`, `plot_ion_dedx(path)`, `plot_ion_zeff(path)`, `plot_ion_range(path)` | ion data |
| `save_figure(fig, path, out_dir=None, formats=("png","pdf"))` | save next to `path` or in `out_dir` |
| `plot_folder(folder, out_dir=None, formats=..., only=None, recursive=False, **kw)` | everything recognised in a folder |
| `detect_kind(path)` | kind of a file, or `None` |

## Look and customisation

The look is a clean, modern one: only left and bottom axis lines, outward ticks, a light grid, rounded
semi-transparent legend boxes, a bold left-aligned title, and a high-contrast, colour-blind-safe palette
(Paul Tol "vibrant" colours; curves ordered by time use a dark-purple-to-orange sequence). Everything is set
in `style.py`:

- `PALETTE`, `TOTAL_COLOUR`, `TIME_CMAP` - colours;
- `SINGLE`, `DOUBLE` - figure sizes in inches;
- `apply_style()` - fonts, line widths, grid, legend (`legend.framealpha`, ...) and output settings.

Long y labels are shrunk automatically only if they would not fit on the figure.

## Calling it from TREKIS-3

`run.py` is a launcher that works from any folder. `FORTRAN_INTEGRATION.md` describes the small changes
to the Fortran sources (keywords `python`, `python_cmd`, `python_script` in `INPUT_PARAMETERS.txt` and a
`Python_plots` subroutine) that start it after the material, ion and dynamics files have been written.

## Files of the package

| File | Content |
|---|---|
| `__main__.py` | command line (`python -m TREKIS3_python_plots`) |
| `run.py` | launcher for use from any folder / from TREKIS-3 |
| `__init__.py` | package interface (the functions above) |
| `plots.py` | time series, radial, spectrum and angular plots; dispatcher; saving; folder walking |
| `material_plots.py` | mean free paths, DOS, ion plots |
| `readers.py` | file readers and file-kind detection, shell-name lookup |
| `labels.py` | axis labels, legend names and titles from file and folder names |
| `style.py` | colours, sizes and the matplotlib style |
| `common.py` | helpers: empty-data check, titles, y-label fitting |
| `example_usage.py` | example script with a CONFIG block |
| `FORTRAN_INTEGRATION.md` | how to call the package from the Fortran code |

## Troubleshooting

- *"skipped ...: no non-zero data"* - the file contains only zeros (or only the `1e30` placeholder); nothing to plot.
- *Generic shell names ("Shell 1", ...)* - no `Total_energies.txt` with a matching number of shells was found
  near the file (e.g. in the first run of a new material); plot again after a run, or pass `shell_names`.
- *Figures too wide or titles cut* - use `--width double`, or a shorter `--title`.
