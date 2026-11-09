#!/usr/bin/env python3
"""Step 2: background fits of data_obs in every category.

For each category, every candidate function (see dijetfit/background.py) is fitted to data_obs
in the fit range [XMIN, XMAX] (extended binned maximum likelihood). A function *passes* if
  - the fit has status 0 and a full, accurate covariance matrix (covQual 3),
  - no parameter is at its limit,
  - the goodness-of-fit p-value (Poisson likelihood-ratio chi2) exceeds --pvalue-threshold,
  - within its family, its order is not above the order selected by the Fisher F-test.
RooFit/Minuit warnings and errors are counted and reported, but do not fail a fit.
The passing functions of a category later enter its RooMultiPdf.

The fit range [--xmin, --xmax] is common to all categories and must be given explicitly (the
script has no built-in range). With --scan, the lower edge is chosen automatically, starting
from --xmin, following this priority:
  1. original binning: XMIN, XMIN + --scan-step, ... up to --scan-max-unbinned, until every
     category has at least one passing function;
  2. if some categories still fail at that XMIN: rebin only the failing categories, trying the
     factors --rebin-factors [default 2 5 10] in turn;
  3. if some categories still fail: raise XMIN further (by --scan-step) up to --scan-max,
     each category keeping the smallest rebin factor that works.
--scan-max-unbinned and --scan-max are required with --scan.
Rebinning merges consecutive bins starting from XMIN (the last bin may be narrower).
Every decision is logged in the output JSON and in the summary.

Outputs (in --outdir):
  background_fits.json   fit range, binning per category, all fit results, selection
  background_summary.md  human-readable summary (incl. the scan log)
  plots/bkg_<CAT>_<FUNC>.{png,pdf}   one plot per fit (data, fit, pulls)
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dijetfit import background as B  # noqa: E402
from dijetfit.common import DATA_NAME, category_labels, setup_root, write_json  # noqa: E402
from dijetfit.histtools import bin_edges_in_window, binning_label, window_hist  # noqa: E402
from dijetfit.plotting import fit_plot  # noqa: E402


def list_categories(tfile):
    return sorted(k.GetName() for k in tfile.GetListOfKeys() if k.GetClassName().startswith("TDirectory"))


def fit_category(ROOT, raw, cat, xmin, xmax, rebin, functions, pthr, fthr):
    h = window_hist(ROOT, raw, xmin, xmax, rebin, name=f"{cat}_data")
    fits = [B.fit_function(ROOT, f, h, cat) for f in functions]
    verdict = B.select_functions(fits, h.GetNbinsX(), pthr, fthr)
    for f in fits:
        f["verdict"] = verdict[f["name"]]
    passing = [f["name"] for f in fits if f["verdict"] == "pass"]
    return h, fits, passing


def public(fit):
    return {k: v for k, v in fit.items() if not k.startswith("_") and k not in ("curve", "expected")}


def parse_rebin(items):
    out = {}
    for it in items or []:
        try:
            cat, n = it.split("=")
            out[cat] = int(n)
        except ValueError:
            sys.exit(f"ERROR: bad --rebin '{it}' (expected CATEGORY=N)")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="ROOT file from prepare_inputs.py")
    ap.add_argument("-o", "--outdir", default="output/background")
    ap.add_argument("--xmin", type=float, required=True,
                    help="lower edge of the fit range [GeV] (the starting point with --scan)")
    ap.add_argument("--xmax", type=float, required=True, help="upper edge of the fit range [GeV]")
    ap.add_argument("--rebin", nargs="*", metavar="CAT=N",
                    help="merge N bins in category CAT (default: no rebinning)")
    ap.add_argument("--categories", nargs="*", help="only these categories (default: all)")
    ap.add_argument("--functions", nargs="+", default=B.function_names(),
                    help=f"candidate functions (default: {' '.join(B.function_names())})")
    ap.add_argument("--pvalue-threshold", type=float, default=0.05)
    ap.add_argument("--ftest-threshold", type=float, default=0.05)
    ap.add_argument("--scan", action="store_true",
                    help="choose XMIN (and rebinning) automatically, starting from --xmin")
    ap.add_argument("--scan-step", type=float, default=5.0, help="step of the scan [GeV] (default: 5)")
    ap.add_argument("--scan-max-unbinned", type=float,
                    help="largest XMIN tried before rebinning (required with --scan) [GeV]")
    ap.add_argument("--scan-max", type=float, help="largest XMIN tried at all (required with --scan) [GeV]")
    ap.add_argument("--rebin-factors", type=int, nargs="+", default=[2, 5, 10])
    ap.add_argument("--formats", nargs="+", default=["png", "pdf"])
    args = ap.parse_args()

    for f in args.functions:
        B.parse_function(f)
    if args.xmax <= args.xmin:
        sys.exit(f"ERROR: empty fit range [{args.xmin:g}, {args.xmax:g}]")
    if args.scan and (args.scan_max_unbinned is None or args.scan_max is None):
        sys.exit("ERROR: --scan needs --scan-max-unbinned and --scan-max")
    if args.scan and not (args.xmin <= args.scan_max_unbinned <= args.scan_max < args.xmax):
        sys.exit("ERROR: need XMIN <= --scan-max-unbinned <= --scan-max < XMAX")

    ROOT = setup_root()
    tfile = ROOT.TFile.Open(args.input)
    if not tfile or tfile.IsZombie():
        sys.exit(f"ERROR: cannot open {args.input}")
    cats = args.categories or list_categories(tfile)
    raw = {}
    for c in cats:
        h = tfile.Get(f"{c}/{DATA_NAME}")
        if not h:
            sys.exit(f"ERROR: no {c}/{DATA_NAME} in {args.input}")
        h.SetDirectory(0)
        raw[c] = h
    rebin = {c: 1 for c in cats}
    rebin.update(parse_rebin(args.rebin))

    outdir = Path(args.outdir)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)
    pthr, fthr = args.pvalue_threshold, args.ftest_threshold

    def run(xmin, which, factors):
        res = {}
        for c in which:
            res[c] = fit_category(ROOT, raw[c], c, xmin, args.xmax, factors[c], args.functions,
                                  pthr, fthr)
        return res

    scan_log = []

    def log(msg):
        print(msg, flush=True)
        scan_log.append(msg)

    def summary_line(xmin, res, factors):
        failing = [c for c in res if not res[c][2]]
        log(f"  xmin = {xmin:g} GeV: {len(res) - len(failing)}/{len(res)} categories with a "
            f"passing function" + (f"; failing: {', '.join(failing)}" if failing else "")
            + ("" if all(factors[c] == 1 for c in res) else
               "; rebin: " + ", ".join(f"{c}={factors[c]}" for c in res if factors[c] != 1)))
        return failing

    xmin = args.xmin
    if not args.scan:
        log(f"Fixed fit range [{xmin:g}, {args.xmax:g}] GeV")
        results = run(xmin, cats, rebin)
        summary_line(xmin, results, rebin)
    else:
        log(f"Scan of the lower edge of the fit range (upper edge {args.xmax:g} GeV)")
        log("Stage 1: original binning (no rebinning)")
        while True:
            results = run(xmin, cats, rebin)
            failing = summary_line(xmin, results, rebin)
            if not failing or xmin + args.scan_step > args.scan_max_unbinned + 1e-9:
                break
            xmin += args.scan_step
        if failing:
            log(f"Stage 2: rebin the failing categories at xmin = {xmin:g} GeV")
            for c in failing:
                for fac in args.rebin_factors:
                    rebin[c] = fac
                    results[c] = run(xmin, [c], rebin)[c]
                    ok = bool(results[c][2])
                    log(f"  {c}: rebin {fac} -> {'pass' if ok else 'fail'}")
                    if ok:
                        break
            failing = summary_line(xmin, results, rebin)
        if failing:
            log("Stage 3: raise the lower edge further")
            while failing and xmin + args.scan_step <= args.scan_max + 1e-9:
                xmin += args.scan_step
                results = run(xmin, cats, rebin)
                for c in [c for c in cats if not results[c][2]]:
                    for fac in [f for f in args.rebin_factors if f > rebin[c]]:
                        rebin[c] = fac
                        results[c] = run(xmin, [c], rebin)[c]
                        if results[c][2]:
                            break
                failing = summary_line(xmin, results, rebin)
        if failing:
            log(f"WARNING: no acceptable configuration found up to xmin = {xmin:g} GeV; "
                f"failing categories: {', '.join(failing)}")
        log(f"Selected fit range: [{xmin:g}, {args.xmax:g}] GeV")

    # ---- plots and outputs for the final configuration
    payload = dict(fit_range=[xmin, args.xmax], pvalue_threshold=pthr, ftest_threshold=fthr,
                   functions=args.functions,
                   function_definitions={f: B.describe(f) for f in args.functions},
                   integrator_eps_rel=B.INTEGRATOR_EPS_REL, scan_log=scan_log, categories={})
    md = ["# Background fits\n", f"Fit range: [{xmin:g}, {args.xmax:g}] GeV  ",
          f"Pass criteria: status 0, covQual 3, no parameter at limit, p > {pthr:g}, "
          f"F-test (threshold {fthr:g}) within each family\n", "## Scan / configuration log\n",
          "```", *scan_log, "```\n"]
    label = category_labels(cats)
    for c in cats:
        h, fits, passing = results[c]
        edges = bin_edges_in_window(raw[c], xmin, args.xmax, rebin[c])
        if rebin[c] != 1:
            print(f"NOTE: category {c} is rebinned by {rebin[c]}: {binning_label(edges)}")
        for f in fits:
            labels = [label[c],
                      f"status {f['status']}, covQual {f['cov_qual']}",
                      f"#chi^{{2}}/ndf = {f['chi2']:.1f}/{f['ndf']}, p = {f['pvalue']:.3g}",
                      f"RooFit warnings/errors: {f['n_warnings']}/{f['n_errors']}",
                      f"{f['verdict']}"]
            fit_plot(ROOT, str(outdir / "plots" / f"bkg_{c}_{f['name']}"), args.formats, h,
                     f["expected"], f["curve"], labels, logy=True, data_label="data_obs",
                     fit_label=f"{f['name']} fit")
        best = max((f for f in fits if f["name"] in passing), key=lambda f: f["pvalue"],
                   default=None)
        payload["categories"][c] = dict(
            rebin=rebin[c], bin_edges=edges, binning=binning_label(edges),
            n_data=h.Integral(), passing=passing, best=best["name"] if best else None,
            fits={f["name"]: public(f) for f in fits})
        md += [f"## {c}\n", f"Binning: {binning_label(edges)}"
               + (f" (**rebinned by {rebin[c]}**)" if rebin[c] != 1 else "") + "  ",
               f"Passing: {', '.join(passing) or 'NONE'}; best p-value: {best['name'] if best else '-'}\n",
               "| function | status | covQual | chi2/ndf | p-value | F-test p (next order) | "
               "warn/err | verdict |", "|---|---|---|---|---|---|---|---|"]
        for f in fits:
            fp = f.get("ftest_next")
            md.append(f"| {f['name']} | {f['status']} | {f['cov_qual']} | "
                      f"{f['chi2']:.1f}/{f['ndf']} | {f['pvalue']:.3g} | "
                      f"{'-' if fp is None else f'{fp:.3g}'} | {f['n_warnings']}/{f['n_errors']} | "
                      f"{f['verdict']} |")
        md.append("")
    write_json(outdir / "background_fits.json", payload)
    (outdir / "background_summary.md").write_text("\n".join(md) + "\n")

    print(f"\n{'category':<42s} {'rebin':>5s}  passing functions")
    for c in cats:
        print(f"{c:<42s} {rebin[c]:>5d}  {', '.join(results[c][2]) or 'NONE'}")
    nw = sum(f["n_warnings"] for c in cats for f in results[c][1])
    ne = sum(f["n_errors"] for c in cats for f in results[c][1])
    print(f"\nRooFit/Minuit messages in the final fits: {nw} warnings, {ne} errors")
    print(f"Written {outdir / 'background_fits.json'}, {outdir / 'background_summary.md'}")
    if any(not results[c][2] for c in cats):
        sys.exit("ERROR: some categories have no passing background function")


if __name__ == "__main__":
    main()
