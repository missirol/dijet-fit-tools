#!/usr/bin/env python3
"""Step 3: fit the signal templates with a double-sided Crystal Ball.

For every category and signal mass, the template CAT/ZprimeQQ_<MASS> (yield for 1 nb) is
fitted with an extended RooCrystalBall in the same fit range and binning as the background
of that category (read from background_fits.json). The fitted 'norm' is therefore the signal
yield inside the fit range. The templates are weighted histograms; by default the fit uses
SumW2-corrected uncertainties (--weights sumw2). --weights counts fits the underlying MC counts
and scales the yield back instead.

A fit is *failed* if it has status != 0, a covariance matrix that is not full and accurate
(covQual != 3), or a floating parameter at its limit (after the fixing rules of
dijetfit/signal.py). Every signal fit must pass in every category: the outputs are written (a
failed fit has the flag "usable": false) and the script then stops with an error, i.e. the
analysis case fails. The goodness of fit is reported as a quality figure and does
not reject fits: the p-value of the Poisson likelihood-ratio chi2 on the equivalent MC counts
for uniformly weighted templates (the normal case), of a Neyman chi2 with the bin errors
otherwise.

Mass points outside the fit range or too close to its edges are dropped: a mass M is used only
if XMIN + D <= M <= XMAX - D, with D = --min-distance (default 10 GeV); e.g. with the default
distance, M = 200 GeV is dropped when XMIN > 190 GeV.

Outputs (in --outdir):
  signal_fits.json    fit results per category and mass
  signal_summary.md   human-readable table
  plots/sig_<CAT>_M<MASS>.{png,pdf}
"""
import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dijetfit import signal as S  # noqa: E402
from dijetfit.common import (SIGNAL_PROCESS, category_labels, read_json, setup_root,  # noqa: E402
                             write_json)
from dijetfit.histtools import window_hist  # noqa: E402
from dijetfit.plotting import fit_plot  # noqa: E402


def signal_masses(tdir):
    out = []
    for k in tdir.GetListOfKeys():
        m = re.fullmatch(rf"{SIGNAL_PROCESS}_(\d+)", k.GetName())
        if m:
            out.append(int(m.group(1)))
    return sorted(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="ROOT file from prepare_inputs.py")
    ap.add_argument("--background", required=True,
                    help="background_fits.json (fit range and binning per category)")
    ap.add_argument("-o", "--outdir", default="output/signal")
    ap.add_argument("--categories", nargs="*")
    ap.add_argument("--min-distance", type=float, default=10.0,
                    help="drop masses M closer than this to either edge of the fit range, or "
                         "outside it [GeV] (default: 10)")
    ap.add_argument("--weights", choices=S.WEIGHT_MODES, default="sumw2",
                    help="treatment of the weighted templates: fit the weighted template with "
                         "SumW2Error (sumw2, default), fit the MC counts and scale the yield "
                         "back (counts), or fit the weighted template with AsymptoticError "
                         "(asymptotic) or without error correction (none, for comparisons)")
    ap.add_argument("--formats", nargs="+", default=["png", "pdf"])
    args = ap.parse_args()

    ROOT = setup_root()
    bkg = read_json(args.background)
    xmin, xmax = bkg["fit_range"]
    cats = args.categories or sorted(bkg["categories"])
    unknown = [c for c in cats if c not in bkg["categories"]]
    if unknown:
        sys.exit(f"ERROR: categories not in {args.background}: {', '.join(unknown)}")
    tfile = ROOT.TFile.Open(args.input)
    if not tfile or tfile.IsZombie():
        sys.exit(f"ERROR: cannot open {args.input}")

    outdir = Path(args.outdir)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)
    payload = dict(fit_range=[xmin, xmax], min_distance=args.min_distance, weights=args.weights,
                   categories={})
    md = ["# Signal fits (double-sided Crystal Ball)\n",
          f"Fit range [{xmin:g}, {xmax:g}] GeV; yields for sigma = 1 nb inside the fit range\n"]
    label = category_labels(cats)
    nbad = 0
    for c in cats:
        tdir = tfile.Get(c)
        if not tdir:
            sys.exit(f"ERROR: no category {c} in {args.input}")
        rebin = bkg["categories"][c]["rebin"]
        masses = signal_masses(tdir)
        used = [m for m in masses if xmin + args.min_distance <= m <= xmax - args.min_distance]
        dropped = [m for m in masses if m not in used]
        if dropped:
            print(f"{c}: dropping masses {dropped} (outside [{xmin:g} + {args.min_distance:g}, "
                  f"{xmax:g} - {args.min_distance:g}] GeV)")
        if len(used) < 2:
            sys.exit(f"ERROR: {c}: fewer than two signal masses inside the fit range")
        payload["categories"][c] = dict(rebin=rebin, masses=used, dropped=dropped, fits={})
        md += [f"## {c}\n", f"rebin {rebin}; dropped masses: {dropped or 'none'}\n",
               "| M | status | covQual | chi2/ndf | p | warn/err | at limit | mean | sigma | "
               "alphaL | nL | alphaR | nR | yield |", "|" + "---|" * 14]
        for m in used:
            raw = tdir.Get(f"{SIGNAL_PROCESS}_{m}")
            h = window_hist(ROOT, raw, xmin, xmax, rebin, name=f"{c}_M{m}")
            fit = S.fit_signal(ROOT, h, m, f"{c}_M{m}", weights=args.weights)
            ok = S.fit_ok(fit)
            nbad += not ok
            P = fit["params"]
            labels = [f"{label[c]}, M = {m} GeV",
                      f"status {fit['status']}, covQual {fit['cov_qual']}",
                      f"#chi^{{2}}/ndf = {fit['chi2']:.1f}/{fit['ndf']}, p = {fit['pvalue']:.2g}",
                      f"RooFit warnings/errors: {fit['n_warnings']}/{fit['n_errors']}"]
            labels += [f"{k} = {P[k]['value']:.4g} #pm {P[k]['error']:.2g}" for k in S.ALL_PARAMS]
            fit_plot(ROOT, str(outdir / "plots" / f"sig_{c}_M{m}"), args.formats, h,
                     fit["expected"], fit["curve"], labels, data_label=f"Z' M = {m} GeV (1 nb)",
                     fit_label="Crystal Ball fit")
            payload["categories"][c]["fits"][str(m)] = {
                k: v for k, v in fit.items() if not k.startswith("_") and k not in ("expected", "curve")}
            payload["categories"][c]["fits"][str(m)]["usable"] = ok
            md.append(f"| {m} | {fit['status']} | {fit['cov_qual']} | {fit['chi2']:.1f}/{fit['ndf']} | "
                      f"{fit['pvalue']:.3g} | {fit['n_warnings']}/{fit['n_errors']} | "
                      f"{','.join(fit['params_at_limit']) or '-'} | "
                      + " | ".join(f"{P[k]['value']:.4g}" for k in S.ALL_PARAMS) + " |")
            flag = "" if ok else "   <-- FAILED"
            print(f"{c} M{m}: status {fit['status']} covQual {fit['cov_qual']} "
                  f"chi2/ndf {fit['chi2']:.1f}/{fit['ndf']} warn/err {fit['n_warnings']}/"
                  f"{fit['n_errors']} at-limit {fit['params_at_limit'] or '-'}{flag}", flush=True)
        md.append("")
    write_json(outdir / "signal_fits.json", payload)
    (outdir / "signal_summary.md").write_text("\n".join(md) + "\n")
    print(f"\n{nbad} failed fit(s) (status != 0, covQual != 3 or a parameter at its limit)")
    print(f"Written {outdir / 'signal_fits.json'}, {outdir / 'signal_summary.md'}")
    if nbad:
        sys.exit(f"ERROR: {nbad} signal fit(s) failed: the analysis case fails (every signal fit must "
                 "pass in every category; see signal_summary.md)")


if __name__ == "__main__":
    main()
