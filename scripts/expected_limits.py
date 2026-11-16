#!/usr/bin/env python3
"""Step 7: expected upper limits on the signal strength with CMS Combine (all categories combined).

For every mass in a grid over the signal mass range of the workspace, Combine computes the
expected 95% CL upper limit on r (x sigma = 1 nb) with the asymptotic CLs method
(`combine -M AsymptoticLimits --run blind`): the Asimov data set of the background-only
hypothesis, built from the background model stored in the workspace (the fit to data), without
looking at the data themselves. The background function is profiled (RooMultiPdf, discrete
profiling with the recommended options); r is common to all categories, the background
parameters, yields and indices are per category.

The mass grid is every --step GeV from the lowest to the highest signal mass of the workspace
(the highest one is always included), or the masses given with --masses; masses between the
signal mass points use the interpolated signal model.

A mass point fails if Combine exits with an error or does not write all five quantiles; the
script then still writes its outputs and exits with an error.

Outputs (in --outdir):
  limits.json            per mass: the quantiles 0.025, 0.16, 0.5, 0.84, 0.975 of the expected
                         limit (`expected`), Combine's exit code (`rc`) and its log file
  limits_summary.md      table of the expected limits
  expected_limits.<fmt>  median expected limit vs mass with the 68% and 95% bands
  work/                  Combine outputs and logs (M<mass>.log)
"""
import argparse
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dijetfit import combine as C  # noqa: E402
from dijetfit.common import read_json, setup_root, write_json  # noqa: E402

QUANTILES = ("0.025", "0.16", "0.5", "0.84", "0.975")


def mass_grid(lo, hi, step):
    """Masses from lo to hi every `step` GeV; hi is always included."""
    if step <= 0:
        raise ValueError("the step must be positive")
    lo, hi = int(round(lo)), int(round(hi))
    masses = list(range(lo, hi + 1, step))
    if masses[-1] != hi:
        masses.append(hi)
    return masses


def combine_command(ws, mass):
    return ["combine", "-M", "AsymptoticLimits", str(ws), "-m", str(mass), "--run", "blind",
            *C.DISCRETE_PROFILING_OPTIONS, "-n", f".M{mass}"]


def read_limits(ROOT, path):
    """{quantile string: limit} from a Combine AsymptoticLimits output file ({} if unreadable)."""
    if not Path(path).exists():
        return {}
    f = ROOT.TFile.Open(str(path))
    t = f.Get("limit") if f and not f.IsZombie() else None
    out = {}
    for i in range(t.GetEntries() if t else 0):
        t.GetEntry(i)
        q = f"{t.quantileExpected:.3f}".rstrip("0").rstrip(".")
        if q in QUANTILES:
            out[q] = float(t.limit)
    if f:
        f.Close()
    return out


def is_complete(expected):
    """All five quantiles present, positive and ordered."""
    vals = [expected.get(q) for q in QUANTILES]
    return all(v is not None and v > 0 for v in vals) and all(a <= b for a, b in zip(vals, vals[1:]))


def summary_plot(ROOT, path_stem, masses, limits, formats, title):
    from array import array
    ROOT.gROOT.SetBatch(True)
    ok = [m for m in masses if is_complete(limits[m]["expected"])]
    if not ok:
        return
    q = {k: [limits[m]["expected"][k] for m in ok] for k in QUANTILES}
    x = array("d", ok)
    zero = array("d", [0.0] * len(ok))
    med = array("d", q["0.5"])

    def band(lo, hi):
        return ROOT.TGraphAsymmErrors(len(ok), x, med, zero, zero,
                                      array("d", [a - b for a, b in zip(q["0.5"], q[lo])]),
                                      array("d", [b - a for a, b in zip(q["0.5"], q[hi])]))

    c = ROOT.TCanvas("c_limits", "", 800, 600)
    c.SetLogy()
    c.SetLeftMargin(0.13)
    g95, g68 = band("0.025", "0.975"), band("0.16", "0.84")
    g95.SetFillColor(ROOT.kOrange)
    g68.SetFillColor(ROOT.kGreen + 1)
    gmed = ROOT.TGraph(len(ok), x, med)
    gmed.SetLineStyle(2)
    gmed.SetLineWidth(2)
    g95.SetTitle(f"{title};m_{{Z'}} [GeV];95% CL upper limit on r (#sigma = 1 nb)")
    g95.SetMinimum(0.5 * min(q["0.025"]))
    g95.SetMaximum(2.0 * max(q["0.975"]))
    g95.Draw("A3")
    g68.Draw("3 same")
    gmed.Draw("L same")
    leg = ROOT.TLegend(0.55, 0.72, 0.88, 0.88)
    leg.AddEntry(gmed, "Expected (median)", "l")
    leg.AddEntry(g68, "68% expected", "f")
    leg.AddEntry(g95, "95% expected", "f")
    leg.Draw()
    for fmt in formats:
        c.SaveAs(f"{path_stem}.{fmt}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("datacards", help="datacards.json from make_datacards.py")
    ap.add_argument("-o", "--outdir", default="output/limits")
    ap.add_argument("--step", type=int, default=25,
                    help="spacing of the mass grid [GeV] (default: 25)")
    ap.add_argument("--masses", nargs="+", type=int, help="masses (default: the grid of --step)")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1,
                    help="Combine processes in parallel (default: number of CPUs)")
    ap.add_argument("--formats", nargs="+", default=["png", "pdf"])
    args = ap.parse_args()

    C.require_tool("combine")
    meta_path = Path(args.datacards)
    meta = read_json(meta_path)
    if not meta.get("combined_workspace"):
        sys.exit("ERROR: no combined workspace (make_datacards.py was run with --no-text2workspace)")
    ws = (meta_path.parent / meta["combined_workspace"]).resolve()
    lo_m, hi_m = meta["mass_range"]
    try:
        masses = sorted(set(args.masses)) if args.masses else mass_grid(lo_m, hi_m, args.step)
    except ValueError as e:
        sys.exit(f"ERROR: {e}")
    bad = [m for m in masses if not lo_m <= m <= hi_m]
    if bad:
        sys.exit(f"ERROR: masses {bad} outside the mass range of the workspace [{lo_m}, {hi_m}]")
    outdir = Path(args.outdir)
    work = outdir / "work"
    work.mkdir(parents=True, exist_ok=True)
    ncat = len(meta["categories"])
    print(f"Expected limits (AsymptoticLimits, blind, envelope), {ncat} categories combined, "
          f"masses {masses}, {args.jobs} in parallel")

    def run(m):
        log = work / f"M{m}.log"
        log.unlink(missing_ok=True)
        try:
            C.run(combine_command(ws, m), cwd=work, log=log)
            return m, 0
        except RuntimeError:
            return m, 1

    with ThreadPoolExecutor(max(1, args.jobs)) as ex:
        rcs = dict(ex.map(run, masses))

    ROOT = setup_root(load_combine=True)  # ROOT I/O in the main thread only
    limits, failed = {}, []
    for m in masses:
        exp = read_limits(ROOT, work / f"higgsCombine.M{m}.AsymptoticLimits.mH{m}.root")
        limits[m] = dict(rc=rcs[m], expected=exp, log=str((work / f"M{m}.log").relative_to(outdir)))
        if rcs[m] or not is_complete(exp):
            failed.append(m)

    write_json(outdir / "limits.json", dict(
        datacards=str(meta_path.resolve()), workspace=str(ws), categories=list(meta["categories"]),
        fit_range=meta.get("fit_range"), method="AsymptoticLimits --run blind",
        options=C.DISCRETE_PROFILING_OPTIONS, masses=masses, quantiles=list(QUANTILES),
        limits={str(m): v for m, v in limits.items()}, failed_masses=failed))

    lines = ["# Expected 95% CL upper limits on the signal strength", "",
             f"{ncat} categories combined; r x sigma = 1 nb; asymptotic CLs, Asimov data of the "
             "background-only hypothesis (blind), discrete profiling of the background function.", "",
             "| M [GeV] | -2 sigma | -1 sigma | median | +1 sigma | +2 sigma |",
             "|---|---|---|---|---|---|"]
    for m in masses:
        e = limits[m]["expected"]
        cells = [f"{e[q]:.4g}" if q in e else "-" for q in QUANTILES]
        lines.append(f"| {m} | " + " | ".join(cells) + (" | **failed** |" if m in failed else " |"))
    if failed:
        lines += ["", f"Failed mass points: {failed} (see work/M<mass>.log)."]
    (outdir / "limits_summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[4:]))
    summary_plot(ROOT, outdir / "expected_limits", masses, limits, args.formats,
                 f"{ncat} categories combined")
    print(f"Written {outdir / 'limits.json'}, {outdir / 'limits_summary.md'}")
    if failed:
        sys.exit(f"ERROR: the expected limit failed at {failed} GeV")


if __name__ == "__main__":
    main()
