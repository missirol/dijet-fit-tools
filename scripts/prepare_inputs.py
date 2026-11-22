#!/usr/bin/env python3
"""Step 1: collect the m_jj histograms of all analysis categories into one ROOT file.

Input : one ROOT file per signal mass point (Z'->qq, mass parsed from '_M<mass>_' in the
        file name) and one data file.
Output: a ROOT file with, for every category CAT (rapidity-separation slice),
            CAT/data_obs          data, copied unchanged
            CAT/ZprimeQQ_<MASS>   signal, scaled by  lumi / N_processed
        where N_processed is the integral of the histogram 'eventsProcessed' of the same
        signal file and lumi is in nb^-1. The signal histograms are therefore the expected
        yields for sigma(pp->Z'->qq) = 1 nb.

Categories are selected with regular expressions (--categories, optionally --exclude): a
histogram of the data file is used if its whole name matches one of the --categories
patterns and none of the --exclude patterns. The category name CAT is the first capture group
of the matching pattern, or the full histogram name if the pattern has no capture group.
Every signal file must contain the same histograms.

Usage:
    prepare_inputs.py --data DATA.root --signal SIGNAL_M200_*.root ... \
        --categories 'REGEX' ['REGEX' ...] [--exclude 'REGEX' ...] -o inputs.root
"""
import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dijetfit.common import (DATA_NAME, mass_from_filename, select_categories,  # noqa: E402
                             setup_root, signal_hist_name)


def bin_edges(h):
    ax = h.GetXaxis()
    return [ax.GetBinLowEdge(i) for i in range(1, ax.GetNbins() + 2)]


def same_binning(a, b):
    ea, eb = bin_edges(a), bin_edges(b)
    return len(ea) == len(eb) and all(abs(x - y) <= 1e-9 * max(1.0, abs(x)) for x, y in zip(ea, eb))


def histogram_names(tfile):
    """Names of the 1D histograms at the top level of `tfile`."""
    return [k.GetName() for k in tfile.GetListOfKeys() if k.GetClassName().startswith("TH1")]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="data ROOT file")
    ap.add_argument("--signal", nargs="+", required=True, help="signal ROOT files (one per mass)")
    ap.add_argument("--lumi", type=float, default=680000.0,
                    help="integrated luminosity in nb^-1 (default: 680000)")
    ap.add_argument("--nevents-hist", default="eventsProcessed",
                    help="histogram with the number of processed signal events "
                         "(default: eventsProcessed)")
    ap.add_argument("--categories", nargs="+", required=True, metavar="REGEX",
                    help="regular expressions selecting the category histograms (whole-name "
                         "match); the first capture group, if any, is the category name")
    ap.add_argument("--exclude", nargs="*", default=[], metavar="REGEX",
                    help="regular expressions of histogram names to leave out")
    ap.add_argument("-o", "--output", default="inputs.root", help="output ROOT file")
    args = ap.parse_args()

    ROOT = setup_root()

    fdata = ROOT.TFile.Open(args.data)
    if not fdata or fdata.IsZombie():
        sys.exit(f"ERROR: cannot open {args.data}")
    try:
        cats, used = select_categories(histogram_names(fdata), args.categories, args.exclude)
    except (ValueError, re.error) as exc:  # bad regex or duplicate category names
        sys.exit(f"ERROR: category selection: {exc}")
    for pat, n in used.items():
        if n == 0:
            print(f"WARNING: --categories pattern '{pat}' selects no histogram in {args.data}",
                  file=sys.stderr)
    if not cats:
        sys.exit(f"ERROR: no histogram of {args.data} matches the --categories patterns")
    print(f"{len(cats)} categories: " + ", ".join(f"{c} <- {cats[c]}" for c in sorted(cats)))

    signals = {}
    for path in args.signal:
        mass = mass_from_filename(path)
        if mass is None:
            sys.exit(f"ERROR: cannot parse the mass ('_M<mass>_') from the file name {path}")
        if mass in signals:
            sys.exit(f"ERROR: two signal files for mass {mass}")
        signals[mass] = path

    out = ROOT.TFile(args.output, "RECREATE")
    yields = {}
    for cat, raw in sorted(cats.items()):
        d = out.mkdir(cat)
        h = fdata.Get(raw).Clone(DATA_NAME)
        h.SetDirectory(0)
        h.SetTitle(f"data, {cat}")
        d.WriteObject(h, DATA_NAME)
        yields[(cat, "data")] = h.Integral()

    for mass, path in sorted(signals.items()):
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            sys.exit(f"ERROR: cannot open {path}")
        hn = f.Get(args.nevents_hist)
        if not hn or not hn.InheritsFrom("TH1"):
            sys.exit(f"ERROR: no histogram '{args.nevents_hist}' in {path}")
        nproc = hn.Integral(0, hn.GetNbinsX() + 1)
        if nproc <= 0:
            sys.exit(f"ERROR: '{args.nevents_hist}' in {path} has integral {nproc}")
        scale = args.lumi / nproc
        missing = [cats[c] for c in cats if not f.Get(cats[c])]
        if missing:
            sys.exit(f"ERROR: {path} lacks the histograms {sorted(missing)}")
        print(f"M = {mass} GeV: {Path(path).name}, N_processed = {nproc:.8g}, "
              f"scale = lumi/N = {scale:.6g} (1 nb)")
        for cat in sorted(cats):
            src = f.Get(cats[cat])
            if not same_binning(src, fdata.Get(cats[cat])):
                sys.exit(f"ERROR: {path}: histogram {cats[cat]} has a different binning from the "
                         f"data histogram ({src.GetNbinsX()} bins in [{src.GetXaxis().GetXmin():g}, "
                         f"{src.GetXaxis().GetXmax():g}] vs {fdata.Get(cats[cat]).GetNbinsX()} bins in "
                         f"[{fdata.Get(cats[cat]).GetXaxis().GetXmin():g}, "
                         f"{fdata.Get(cats[cat]).GetXaxis().GetXmax():g}])")
            h = src.Clone(signal_hist_name(mass))
            h.SetDirectory(0)
            if h.GetSumw2N() == 0:
                h.Sumw2()
            h.Scale(scale)
            h.SetTitle(f"Z'->qq M={mass} GeV, sigma=1 nb, L={args.lumi:g}/nb, {cat}")
            out.Get(cat).WriteObject(h, signal_hist_name(mass))
            yields[(cat, mass)] = h.Integral()
        f.Close()

    # provenance
    meta = ROOT.TNamed("provenance",
                       f"lumi_invnb={args.lumi:g}; signal_xsec=1nb; nevents_hist={args.nevents_hist}; "
                       f"data={Path(args.data).name}; "
                       f"categories={' '.join(args.categories)}; exclude={' '.join(args.exclude)}; "
                       + "; ".join(f"M{m}={Path(p).name}" for m, p in sorted(signals.items())))
    out.cd()
    meta.Write()
    out.Close()

    masses = sorted(signals)
    print(f"\nYields in the full histogram range (signal: sigma = 1 nb, L = {args.lumi:g} nb^-1)")
    print(f"  {'category':<42s} {'data':>12s}" + "".join(f" {'M' + str(m):>10s}" for m in masses))
    for cat in sorted(cats):
        print(f"  {cat:<42s} {yields[(cat, 'data')]:>12.6g}"
              + "".join(f" {yields[(cat, m)]:>10.4g}" for m in masses))
    print(f"\nWritten {args.output}")


if __name__ == "__main__":
    main()
