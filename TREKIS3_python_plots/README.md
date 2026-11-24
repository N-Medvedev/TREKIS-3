# TREKIS3_python_plots

Publication-quality matplotlib figures for TREKIS-3 output. Needs only `numpy` and `matplotlib`.

## Supported files (detected from the file name)

| File | Plot |
|---|---|
| `Total_energies.txt`, `Total_numbers.txt` | quantities vs. time (log time axis; energies and particle numbers in separate panels; all-zero columns hidden) |
| `Radial_*.txt` (any particle/shell, density/energy/temperature) | radial profile, one curve per time |
| `*_distribution_vs_E[1_eV].txt` | energy spectrum, one curve per time |
| `*_theta_distribution.txt` | angular distribution, one curve per time |
| `OUTPUT_Electron_IMFPs*.dat`, `OUTPUT_Hole_IMFPs*.dat`, `OUTPUT_Photon_IMFPs*.dat` | inelastic MFP vs. energy, per shell + total (log-log) |
| `OUTPUT_Electron*EMFPs*.dat`, `OUTPUT_Hole*EMFPs*.dat` | elastic MFP vs. energy |
| `OUTPUT_<material>_DOS_analysis.dat` | density of states and its integral |
| `OUTPUT_<ion>..._IMFP.dat`, `_dEdx.dat`, `_effective_charges.dat`, `_Range.dat` | ion MFP, energy loss, effective charge, range vs. ion energy (MeV) |

Dynamics files: column and time information is read from the file headers, never by column number.
The MFP and ion files have no header (columns: energy, one column per shell, total). Shell names for
their legends are taken from the header of a `Total_energies.txt` found near the file (the run folders
contain it); if their number does not match, generic names ("Shell 1", ...) are used. Pass
`shell_names=[...]` to set them explicitly.
File names with brackets (`...[eV_A^-3].txt`) or underscores (`..._eV_A_-3_.txt`) both work.

## Command line

Run from the folder that contains `TREKIS3_python_plots/`. A folder is searched recursively
(folders named `diff_CS` and `figures` are skipped) and the figures of each folder are written to
a `figures` sub-folder next to its data:

    python -m TREKIS3_python_plots OUTPUT_Al2O3                       # whole tree
    python -m TREKIS3_python_plots OUTPUT_Al2O3 --no-recursive        # only that folder
    python -m TREKIS3_python_plots path\to\run_folder
    python -m TREKIS3_python_plots Radial_electron_energy[eV_A^-3].txt -t 0.1 1 10 100 -o figures
    python -m TREKIS3_python_plots run_folder --only radial spectrum --width double -f pdf

Options: `-o` output dir, `-f` formats (default png pdf), `-t` times [fs], `--width single|double`,
`--xlin`, `--ylin`, `--ylog`, `--ydecades N`, `--show-empty`, `--only`, `--no-recursive`, `--show`.

## From Python

    from TREKIS3_python_plots import plot_file, save_figure
    fig = plot_file("Radial_electron_energy[eV_A^-3].txt", times=[1, 10, 100], xlim=(1, 300))
    save_figure(fig, "Radial_electron_energy[eV_A^-3].txt", out_dir="figures")

See `example_usage.py` for a script with a CONFIG block.

## Layout

- `readers.py` file parsing  - `style.py` look (fonts, sizes, colours)  - `labels.py` axis labels
- `plots.py` dynamics plots, dispatcher, folder walking  - `material_plots.py` MFP, DOS and ion plots  - `__main__.py` command line
