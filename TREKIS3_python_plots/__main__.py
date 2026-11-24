"""
! This file is part of TREKIS-3
Command line of TREKIS3_python_plots:

    python -m TREKIS3_python_plots <file(s) or folder(s)> [options]

(run from the folder that contains the TREKIS3_python_plots folder; run.py does the same from
any directory). The kind of plot is detected from the file name. Folders are searched
recursively unless --no-recursive is given; 'diff_CS' and 'figures' folders are skipped, and
files without non-zero data are skipped with a warning.

Options:
    -o, --out DIR          output directory (default: next to the data / in a 'figures' folder)
    -f, --formats ...      output formats, e.g. png pdf svg (default: png pdf)
    --dpi N                resolution of raster formats (default 300)
    -r, --no-recursive     for folders: do not descend into sub-folders
    --only WORD ...        for folders: only files whose name contains one of these words
    -t, --times T ...      only the time steps closest to these times [fs]
    --width single|double  figure size: one column (3.4 in) or full page width (7 in)
    --xlin, --ylin, --ylog linear x axis / linear y axis / logarithmic y axis
    --ydecades N           dynamic range of a logarithmic y axis
    --show-empty           time series: also draw columns that are zero everywhere
    --no-title             no title above the figures
    --title TEXT           use this title instead of the automatic one
    --show                 open the figures in a window
"""
import argparse
from pathlib import Path

import matplotlib

from .plots import accepted_kwargs, plot_file, plot_folder, save_figure
from .readers import detect_kind
from .style import DOUBLE, SINGLE


def main():
    """Parse the command line and make the plots."""
    ap = argparse.ArgumentParser(
        prog="python -m TREKIS3_python_plots",
        description="Plot TREKIS-3 output files. The plot type is detected from the file name.")
    ap.add_argument("paths", nargs="+", help="output file(s) and/or run folder(s)")
    ap.add_argument("-o", "--out", default=None,
                    help="output directory (default: next to the data file)")
    ap.add_argument("-f", "--formats", nargs="+", default=["png", "pdf"],
                    help="output formats (default: png pdf)")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("-r", "--no-recursive", action="store_true",
                    help="for folders: do not descend into sub-folders (default: descend; "
                         "'diff_CS' and 'figures' folders are skipped)")
    ap.add_argument("--only", nargs="+", default=None,
                    help="for folders: only files whose name contains one of these words "
                         "(e.g. Radial IMFP dEdx)")
    ap.add_argument("-t", "--times", type=float, nargs="+", default=None,
                    help="plot only the time steps closest to these times [fs] "
                         "(radial, spectra, angular distributions)")
    ap.add_argument("--width", choices=["single", "double"], default="single",
                    help="figure size: single (3.4 in) or double (7 in) column width")
    ap.add_argument("--xlin", action="store_true", help="linear x axis (default: log where sensible)")
    ap.add_argument("--ylin", action="store_true", help="linear y axis")
    ap.add_argument("--ylog", action="store_true", help="logarithmic y axis")
    ap.add_argument("--ydecades", type=float, default=None,
                    help="decades of dynamic range shown on a log y axis (default 8)")
    ap.add_argument("--show-empty", action="store_true",
                    help="time series: also plot columns that are zero everywhere")
    ap.add_argument("--no-title", action="store_true", help="no title above the figures")
    ap.add_argument("--title", default=None, help="title text instead of the automatic one")
    ap.add_argument("--show", action="store_true", help="open the figures in a window")
    args = ap.parse_args()

    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    kw = {}
    if args.times:
        kw["times"] = args.times
    if args.xlin:
        kw["xlog"] = False
    if args.ylin:
        kw["ylog"] = False
    if args.ylog:
        kw["ylog"] = True
    if args.ydecades is not None:
        kw["ydecades"] = args.ydecades
    if args.show_empty:
        kw["hide_empty"] = False
    if args.no_title:
        kw["title"] = False
    elif args.title:
        kw["title"] = args.title
    kw["figsize"] = DOUBLE if args.width == "double" else SINGLE

    saved = []
    for p in map(Path, args.paths):
        if p.is_dir():
            saved += plot_folder(p, args.out, args.formats, args.only,
                                 recursive=not args.no_recursive, **kw)
        elif detect_kind(p) is None:
            print(f"WARNING: {p.name}: not a recognised TREKIS-3 output file, skipped")
        else:
            try:
                fig = plot_file(p, **accepted_kwargs(detect_kind(p), kw))
                saved += save_figure(fig, p, args.out, args.formats, args.dpi)
            except Exception as exc:
                print(f"WARNING: skipped {p}: {exc}")
                continue
            if not args.show:
                plt.close(fig)
    for s in saved:
        print(f"Saved {s}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
