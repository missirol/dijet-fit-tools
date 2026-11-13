#!/usr/bin/env python3
"""Step 2: background fits of data_obs in every category.

For each category, every candidate function (see dijetfit/background.py) is fitted to data_obs
in the fit range [XMIN, XMAX] (extended binned maximum likelihood). A fit is acceptable if
  - the fit has status 0 and a full, accurate covariance matrix (covQual 3),
  - no parameter is at its limit,
  - the goodness-of-fit p-value (Poisson likelihood-ratio chi2) exceeds --pvalue-threshold.
The functions of the category's RooMultiPdf are chosen with standard Fisher F-tests (see
dijetfit/background.py, select_functions):
  - in every family (dijet, moddijet, ua2), starting from the lowest acceptable order, move to
    the next order while it is significantly better (F-test p-value < --ftest-threshold) and
    its fit is acceptable: this gives the family's highest order; the family contributes its
    acceptable orders up to that one;
  - default: among the highest orders of the families, the one with the largest p-value;
  - alternates: the other contributed functions ranked by p-value, at most --max-alternates;
    with --alternates family-max, only the highest order of each family is used (at most one
    function per family).
RooFit/Minuit warnings and errors are counted and reported, but do not fail a fit.

The fit range [--xmin, --xmax] is common to all categories and must be given explicitly (the
script has no built-in range). With --scan, the lower edge is chosen automatically, starting
from --xmin and raised in steps of --scan-step up to --scan-max, until every category has at
least one passing function (see dijetfit/scan.py).

Rebinning is off by default, because build_workspace.py refuses categories with bins wider than
the original 1 GeV (Combine evaluates the pdfs at the bin centres; see docs/OUTPUT_FORMAT.md).
It can be switched on for studies with --rebin-factors (e.g. 2 5 10): the lower edge is then
raised without rebinning only up to --scan-max-unbinned; there, the failing categories are
rebinned with the smallest factor that works; beyond, at every new lower edge each category
starts again from no rebinning. --rebin fixes the rebinning of given categories without a scan.

If a category has no passing function (with --scan: still at --scan-max), the script writes
its outputs for inspection and stops with an error: every category of an analysis case needs at
least one passing function; categories are never dropped (a different categorisation is needed).
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
from dijetfit.scan import scan_fit_range  # noqa: E402


def list_categories(tfile):
    return sorted(k.GetName() for k in tfile.GetListOfKeys() if k.GetClassName().startswith("TDirectory"))


def fit_category(ROOT, raw, cat, xmin, xmax, rebin, functions, sel):
    """Fit all functions; returns (histogram, fits, passing) with passing = [default] +
    alternates, in this order (empty if no function is acceptable)."""
    h = window_hist(ROOT, raw, xmin, xmax, rebin, name=f"{cat}_data")
    fits = []
    for f in functions:
        try:
            fits.append(B.fit_function(ROOT, f, h, cat))
        except Exception as exc:  # e.g. a numerical failure of one function: that fit fails
            print(f"WARNING: {cat}: fit of {f} raised {type(exc).__name__}: {exc}", file=sys.stderr)
            fits.append(B.failed_fit(f, h.GetNbinsX(), exc))
    verdict, default = B.select_functions(fits, h.GetNbinsX(), **sel)
    for f in fits:
        f["verdict"] = verdict[f["name"]]
    alternates = sorted((f for f in fits if f["verdict"].startswith("pass (alternate")),
                        key=lambda f: int(f["verdict"].split()[-1].rstrip(")")))
    passing = ([default] + [f["name"] for f in alternates]) if default else []
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
    ap.add_argument("--pvalue-threshold", type=float, default=0.05,
                    help="minimum goodness-of-fit p-value of an acceptable fit (default: 0.05)")
    ap.add_argument("--ftest-threshold", type=float, default=0.05,
                    help="F-test p-value below which the next order is significantly better "
                         "(default: 0.05)")
    ap.add_argument("--max-alternates", type=int, default=4,
                    help="maximum number of alternate functions in the RooMultiPdf (default: 4)")
    ap.add_argument("--alternates", choices=["orders", "family-max"], default="orders",
                    help="'orders' (default): every acceptable order up to the family maximum; "
                         "'family-max': only the highest order of each family")
    ap.add_argument("--scan", action="store_true",
                    help="choose XMIN (and rebinning) automatically, starting from --xmin")
    ap.add_argument("--scan-step", type=float, default=5.0, help="step of the scan [GeV] (default: 5)")
    ap.add_argument("--scan-max-unbinned", type=float,
                    help="with --rebin-factors: largest XMIN tried before rebinning [GeV] "
                         "(default: --scan-max)")
    ap.add_argument("--scan-max", type=float, help="largest XMIN tried at all (required with --scan) [GeV]")
    ap.add_argument("--rebin-factors", type=int, nargs="*", default=[],
                    help="with --scan: rebin factors to try (default: none, i.e. no rebinning)")
    ap.add_argument("--formats", nargs="+", default=["png", "pdf"])
    args = ap.parse_args()

    for f in args.functions:
        B.parse_function(f)
    if args.xmax <= args.xmin:
        sys.exit(f"ERROR: empty fit range [{args.xmin:g}, {args.xmax:g}]")
    if args.scan and args.scan_max is None:
        sys.exit("ERROR: --scan needs --scan-max")
    if args.scan and args.scan_max_unbinned is None:
        args.scan_max_unbinned = args.scan_max
    if args.scan and args.rebin:
        sys.exit("ERROR: --rebin cannot be combined with --scan (the scan chooses the rebinning; "
                 "use --rebin-factors)")
    if args.rebin_factors and not args.scan:
        sys.exit("ERROR: --rebin-factors only applies with --scan")
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
    pthr = args.pvalue_threshold
    sel = dict(pvalue_threshold=pthr, ftest_threshold=args.ftest_threshold,
               max_alternates=args.max_alternates, family_max_only=args.alternates == "family-max")

    def run(xmin, which, factors):
        res = {}
        for c in which:
            res[c] = fit_category(ROOT, raw[c], c, xmin, args.xmax, factors[c], args.functions,
                                  sel)
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
        cache = {}

        def evaluate(x, c, fac):
            if (x, c, fac) not in cache:
                cache[(x, c, fac)] = fit_category(ROOT, raw[c], c, x, args.xmax, fac, args.functions, sel)
            res = cache[(x, c, fac)]
            return res, bool(res[2])

        xmin, rebin, results, failing = scan_fit_range(
            evaluate, cats, args.xmin, args.scan_step, args.scan_max_unbinned, args.scan_max,
            args.rebin_factors, log)
        if failing:
            log(f"WARNING: no acceptable configuration found up to xmin = {xmin:g} GeV; "
                f"failing categories: {', '.join(failing)}")
        log(f"Selected fit range: [{xmin:g}, {args.xmax:g}] GeV")

    # ---- plots and outputs for the final configuration
    payload = dict(fit_range=[xmin, args.xmax], pvalue_threshold=pthr,
                   ftest_threshold=args.ftest_threshold, max_alternates=args.max_alternates,
                   alternates=args.alternates,
                   selection="per family, F-tests from the lowest acceptable order give the highest "
                             "order; default: the family maximum with the largest p-value; "
                             + ("alternates: the other acceptable orders up to the family maxima, by p-value"
                                if args.alternates == "orders" else
                                "alternates: the highest orders of the other families, by p-value"),
                   functions=args.functions,
                   function_definitions={f: B.describe(f) for f in args.functions},
                   pdf_normalisation="sum over the bin centres (1 GeV bins); integral, adaptive Gauss-Kronrod "
                                     "(rebinned categories, pdf integrated over the bins)",
                   integrator_eps_rel=B.INTEGRATOR_EPS_REL, scan_log=scan_log, categories={})
    md = ["# Background fits\n", f"Fit range: [{xmin:g}, {args.xmax:g}] GeV  ",
          f"Acceptable fit: status 0, covQual 3, no parameter at limit, p > {pthr:g}. Per family, "
          f"F-tests (threshold {args.ftest_threshold:g}) from the lowest acceptable order give the "
          "highest order; default: the family maximum with the largest p-value; up to "
          f"{args.max_alternates} alternates ranked by p-value"
          + (" (only the highest order of each family)" if args.alternates == "family-max" else "")
          + ".\n", "## Scan / configuration log\n",
          "```", *scan_log, "```\n"]
    label = category_labels(cats)
    for c in cats:
        h, fits, passing = results[c]
        edges = bin_edges_in_window(raw[c], xmin, args.xmax, rebin[c])
        if rebin[c] != 1:
            print(f"WARNING: category {c} is rebinned by {rebin[c]} ({binning_label(edges)}); "
                  "build_workspace.py will refuse it (Combine evaluates the pdfs at the bin centres)",
                  file=sys.stderr)
        for f in fits:
            if f.get("exception"):
                continue  # nothing to plot
            labels = [label[c],
                      f"status {f['status']}, covQual {f['cov_qual']}",
                      f"#chi^{{2}}/ndf = {f['chi2']:.1f}/{f['ndf']}, p = {f['pvalue']:.3g}",
                      f"RooFit warnings/errors: {f['n_warnings']}/{f['n_errors']}",
                      f"{f['verdict']}"]
            fit_plot(ROOT, str(outdir / "plots" / f"bkg_{c}_{f['name']}"), args.formats, h,
                     f["expected"], f["curve"], labels, logy=True, data_label="data_obs",
                     fit_label=f"{f['name']} fit")
        default = next((f for f in fits if passing and f["name"] == passing[0]), None)
        entry = dict(
            rebin=rebin[c], bin_edges=edges, binning=binning_label(edges),
            n_data=h.Integral(), passing=passing, default=default["name"] if default else None,
            fits={f["name"]: public(f) for f in fits})
        payload["categories"][c] = entry
        md += [f"## {c}\n", f"Binning: {binning_label(edges)}"
               + (f" (**rebinned by {rebin[c]}**)" if rebin[c] != 1 else "") + "  ",
               f"RooMultiPdf: {', '.join(passing) or 'NONE'} (default: {default['name'] if default else '-'})\n",
               "| function | status | covQual | chi2/ndf | p-value | F-test p vs next order | "
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
    failing = [c for c in cats if not results[c][2]]
    if failing:
        sys.exit(f"ERROR: {len(failing)} categor{'y has' if len(failing) == 1 else 'ies have'} no passing "
                 f"background function ({', '.join(failing)}): the analysis case fails (every category "
                 "needs one; categories are not dropped, a different categorisation is needed)")


if __name__ == "__main__":
    main()
