#!/usr/bin/env python3
"""Step 7: bias test of the signal strength with CMS Combine (all categories combined).

For every signal mass M and injected signal strength mu, Combine generates pseudo-data
(binned toys, `combine -M GenerateOnly --expectSignal mu`) from the combined workspace of
make_datacards.py: background = the functions selected by the index (by default the default
function of every category, with the parameters of the fit to data stored in the workspace),
signal = mu x the signal model for sigma = 1 nb at M. The background parameters and yields are
free and uncorrelated across the categories; r is common to all categories.

Every toy is fitted with signal + background, r floating:
  envelope  the background function is profiled (RooMultiPdf, discrete profiling, Combine's
            penalty of 0.5 per parameter on -log L), as in the analysis;
  fixed     the background is fixed to the generating function (its parameters float): a
            reference, which separates any bias of the envelope from that of the fit itself.
Best fits only (MultiDimFit --algo none); the functions chosen by the envelope are recorded.

The bias is measured by the mean of (r_fit - mu) in units of the standard deviation of r_fit
over the toys (the actual statistical uncertainty, including the effect of the envelope). Its
uncertainty is 1/sqrt(N toys). A common criterion is |bias| < 0.14 sigma, i.e. a bias that
would increase the total uncertainty by less than 1% if added in quadrature.

With --pull-toys K, the first K toys of each point are also fitted with uncertainties
(FitDiagnostics --robustFit 1, envelope; ~1 minute per toy): pull = (r_fit - mu) / sigma, with
the uncertainty on the side of the injected value. The default Hesse/Minos errors are not usable
with the envelope (the function index changes during the error scan).

The expected uncertainty sigma_r at each mass (Asimov data, mu = 1, envelope, --robustFit) is
measured first; the range of r in the fits is mu +- R_RANGE_SIGMAS x sigma_r.

Outputs (in --outdir):
  bias_results.json     per point: settings, per-toy results, summary statistics
  bias_summary.md       tables
  plots/rfit_<gen>_<fit>_M<M>_mu<mu>.png   distribution of r_fit with a Gaussian fit
  plots/pull_<gen>_M<M>_mu<mu>.png         pull distributions (with --pull-toys)
  plots/bias_summary_<gen>.png             bias vs mass for all mu
  work/                 Combine outputs and logs
"""
import argparse
import itertools
import math
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dijetfit import combine as C  # noqa: E402
from dijetfit.common import read_json, setup_root, write_json  # noqa: E402
from dijetfit.plotting import fit_plot  # noqa: E402

R_RANGE_SIGMAS = 20.0
BIAS_CRITERION = 0.14


def fmt_mu(mu):
    return f"{mu:g}".replace(".", "p").replace("-", "m")


def gen_indices(meta, gen):
    """Index of the generating function per category: 'default', or 'altK' (the K-th alternate
    of every category, or its default if it has fewer alternates)."""
    out = {}
    for c, v in meta["categories"].items():
        n = len(v["index_to_function"])
        if gen == "default":
            out[c] = v["default_index"]
        else:
            k = int(re.fullmatch(r"alt(\d+)", gen).group(1))
            out[c] = k if k < n else v["default_index"]
    return out


def read_tree(ROOT, path, tree, columns):
    """Columns of a TTree as numpy arrays (empty if the file or tree is missing)."""
    if not Path(path).exists():
        return {c: np.array([]) for c in columns}
    f = ROOT.TFile.Open(str(path))
    t = f.Get(tree) if f and not f.IsZombie() else None
    if not t or t.GetEntries() == 0:
        return {c: np.array([]) for c in columns}
    have = {b.GetName() for b in t.GetListOfBranches()}
    cols = [c for c in columns if c in have]
    df = ROOT.RDataFrame(t)
    arr = df.AsNumpy(cols)
    return {c: np.asarray(arr[c]) if c in arr else np.array([]) for c in columns}


class Runner:
    def __init__(self, args, meta, card_dir, work):
        self.args, self.meta, self.work = args, meta, work
        self.ws = (card_dir / meta["combined_workspace"]).resolve()
        self.cats = list(meta["categories"])
        self.idx_names = [meta["categories"][c]["index"] for c in self.cats]

    def combine(self, args, cwd, log):
        C.run(["combine", *map(str, args)], cwd, log)

    def asimov_sigma(self, mass):
        d = self.work / f"asimov_M{mass}"
        d.mkdir(parents=True, exist_ok=True)
        log = d / "combine.log"
        log.write_text("")
        self.combine(["-M", "FitDiagnostics", self.ws, "-m", mass, "-t", "-1", "--expectSignal", 1,
                      "--toysNoSystematics", "--rMin", -20, "--rMax", 20, "--robustFit", 1,
                      "--skipBOnlyFit", *C.DISCRETE_PROFILING_OPTIONS, "-n", ".asimov"], d, log)
        m = re.findall(r"Best fit r: ([-0-9.eE+]+)\s+-([0-9.eE+]+)/\+([0-9.eE+]+)", log.read_text())
        if not m:
            raise RuntimeError(f"no Asimov fit result in {log}")
        r, lo, hi = map(float, m[-1])
        return dict(r=r, err_lo=lo, err_hi=hi, sigma=0.5 * (lo + hi))

    def job(self, gen, mass, mu, j, ntoys, seed, rrange, pull_toys):
        """Generate ntoys toys and fit them; returns the per-toy results of this job."""
        d = self.work / f"{gen}_M{mass}_mu{fmt_mu(mu)}" / f"job{j}"
        d.mkdir(parents=True, exist_ok=True)
        log = d / "combine.log"
        log.write_text("")
        gidx = gen_indices(self.meta, gen)
        setp = ",".join(f"{self.meta['categories'][c]['index']}={gidx[c]}" for c in self.cats)
        common = [self.ws, "-m", mass, "--toysNoSystematics", "-s", seed]
        self.combine(["-M", "GenerateOnly", *common, "-t", ntoys, "--expectSignal", mu,
                      "--saveToys", "--setParameters", setp, "-n", ".gen"], d, log)
        toys = d / f"higgsCombine.gen.GenerateOnly.mH{mass}.{seed}.root"
        fit = [*common, "--toysFile", toys.name, "--rMin", rrange[0], "--rMax", rrange[1],
               *C.DISCRETE_PROFILING_OPTIONS]
        out = dict(seed=seed, ntoys=ntoys, fits={})
        for mode in self.args.fits:
            extra = ["--saveSpecifiedIndex", ",".join(self.idx_names)] if mode == "envelope" else \
                ["--setParameters", setp, "--freezeParameters", ",".join(self.idx_names)]
            self.combine(["-M", "MultiDimFit", *fit, "-t", ntoys, "--algo", "none", *extra,
                          "-n", f".{mode}"], d, log)
            out["fits"][mode] = d / f"higgsCombine.{mode}.MultiDimFit.mH{mass}.{seed}.root"
        if pull_toys:
            self.combine(["-M", "FitDiagnostics", *fit, "-t", pull_toys, "--robustFit", 1,
                          "--skipBOnlyFit", "-n", ".pull"], d, log)
            out["pull_file"] = d / "fitDiagnostics.pull.root"
        return out


def collect(ROOT, job, idx_names):
    """Per-toy results of one job from the Combine output files."""
    res = {}
    for mode, path in job["fits"].items():
        cols = ["iToy", "r", "quantileExpected"] + (idx_names if mode == "envelope" else [])
        a = read_tree(ROOT, path, "limit", cols)
        keep = a["quantileExpected"] == -1 if len(a["quantileExpected"]) else []
        toys = {}
        for i in np.flatnonzero(keep):
            t = int(a["iToy"][i])
            toys[t] = dict(r=float(a["r"][i]))
            if mode == "envelope" and len(a[idx_names[0]]):
                toys[t]["index"] = [int(a[n][i]) for n in idx_names]
        res[mode] = toys
    if job.get("pull_file"):
        # one row per toy, in order (the tree has no toy number)
        a = read_tree(ROOT, job["pull_file"], "tree_fit_sb", ["r", "rLoErr", "rHiErr", "fit_status"])
        res["pull"] = {t: dict(r=float(r), err_lo=float(lo), err_hi=float(hi), status=int(s))
                       for t, (r, lo, hi, s) in enumerate(zip(a["r"], a["rLoErr"], a["rHiErr"],
                                                              a["fit_status"]), start=1)}
    return res


def summarise(values, mu, rrange):
    v = [x for x in values if rrange[0] < x < rrange[1]]  # a fit at the range boundary failed
    st = C.mean_and_error(v)
    out = dict(n=st["n"], n_at_boundary=len(values) - len(v), mean=st["mean"],
               mean_err=st["mean_err"], std=st["std"], median=C.median(v))
    out["bias"] = st["mean"] - mu
    out["bias_over_sigma"] = (st["mean"] - mu) / st["std"] if st["std"] > 0 else float("nan")
    out["bias_over_sigma_err"] = 1 / math.sqrt(st["n"]) if st["n"] > 0 else float("nan")
    out["median_bias_over_sigma"] = (out["median"] - mu) / st["std"] if st["std"] > 0 else float("nan")
    return out


def gauss_plot(ROOT, path, values, title, xtitle, labels, formats, center=None, width=None):
    """Histogram of `values` with a binned maximum-likelihood Gaussian fit (pull sub-pad)."""
    v = np.asarray(values, dtype=float)
    if len(v) < 5:
        return None
    m, s = float(v.mean()), float(v.std(ddof=1)) or 1.0
    lo, hi = (m - 4 * s, m + 4 * s) if center is None else (center - 4 * width, center + 4 * width)
    nb = 32
    h = ROOT.TH1D(f"h_{Path(path).name}", f";{xtitle};toys", nb, lo, hi)
    h.SetDirectory(0)
    for x in v:
        h.Fill(x)
    for i in range(1, nb + 1):
        h.SetBinError(i, math.sqrt(h.GetBinContent(i)))
    f = ROOT.TF1(f"g_{Path(path).name}", "gausn", lo, hi)
    f.SetParameters(len(v) * (hi - lo) / nb, m, s)
    h.Fit(f, "QLN0")  # Poisson likelihood
    norm, mean, sigma = f.GetParameter(0), f.GetParameter(1), f.GetParameter(2)
    w = (hi - lo) / nb
    expected = [f.Integral(h.GetXaxis().GetBinLowEdge(i), h.GetXaxis().GetBinUpEdge(i)) / w * w
                for i in range(1, nb + 1)]
    xs = np.linspace(lo, hi, 400)
    curve = (xs.tolist(), [f.Eval(x) / w for x in xs])
    labels = labels + [f"Gaussian fit: mean {mean:.3g} #pm {f.GetParError(1):.2g}",
                       f"width {abs(sigma):.3g} #pm {f.GetParError(2):.2g}"]
    fit_plot(ROOT, str(path), formats, h, expected, curve, [title] + labels,
             data_label="toys", fit_label="Gaussian fit")
    return dict(mean=mean, mean_err=f.GetParError(1), width=abs(sigma), width_err=f.GetParError(2))


def summary_plot(ROOT, path, points, masses, mus, fits, formats):
    """bias / sigma vs mass, one graph per (mu, fit), with the +-0.14 band."""
    c = ROOT.TCanvas("c_sum", "", 900, 600)
    c.SetLeftMargin(0.12)
    ext = [abs(s["bias_over_sigma"]) + s["bias_over_sigma_err"] for p in points.values()
           for f, s in p["summary"].items() if f in fits]
    ymax = max([0.6] + [1.15 * e for e in ext])
    frame = c.DrawFrame(min(masses) - 40, -ymax, max(masses) + 40, 1.4 * ymax,
                        ";signal mass [GeV];(<r_{fit}> - #mu) / #sigma(r_{fit})")
    frame.GetYaxis().SetTitleOffset(1.3)
    box = ROOT.TBox(min(masses) - 40, -BIAS_CRITERION, max(masses) + 40, BIAS_CRITERION)
    box.SetFillColorAlpha(ROOT.kGreen - 9, 0.5)
    box.Draw()
    line = ROOT.TLine(min(masses) - 40, 0, max(masses) + 40, 0)
    line.SetLineStyle(2)
    line.Draw()
    leg = ROOT.TLegend(0.14, 0.70, 0.60, 0.89)
    leg.SetNColumns(2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    colors = [ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1,
              ROOT.kOrange + 7]
    keep = [box, line, leg]
    for (i, mu), (k, fit) in itertools.product(enumerate(mus), enumerate(fits)):
        xs, ys, es = [], [], []
        for m in masses:
            p = points.get((m, mu))
            if p and fit in p["summary"]:
                s = p["summary"][fit]
                xs.append(m + (i - (len(mus) - 1) / 2) * 6 + k * 2)
                ys.append(s["bias_over_sigma"])
                es.append(s["bias_over_sigma_err"])
        if not xs:
            continue
        g = ROOT.TGraphErrors(len(xs), np.array(xs, "d"), np.array(ys, "d"),
                              np.zeros(len(xs)), np.array(es, "d"))
        g.SetMarkerStyle(20 if fit == "envelope" else 24)
        g.SetMarkerColor(colors[i % len(colors)])
        g.SetLineColor(colors[i % len(colors)])
        g.Draw("P same")
        leg.AddEntry(g, f"#mu = {mu:g}, {fit}", "pl")
        keep.append(g)
    leg.Draw()
    for fmtx in formats:
        c.SaveAs(f"{path}.{fmtx}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("datacards", help="datacards.json from make_datacards.py")
    ap.add_argument("-o", "--outdir", default="output/bias")
    ap.add_argument("--masses", nargs="+", type=int,
                    help="signal masses (default: the signal mass points of the workspace)")
    ap.add_argument("--mu", nargs="+", type=float, default=[0.0, 0.5, 1.0, 5.0],
                    help="injected signal strengths (x sigma = 1 nb; default: 0 0.5 1 5)")
    ap.add_argument("--ntoys", type=int, default=300, help="toys per point (default: 300)")
    ap.add_argument("--toys-per-job", type=int, default=25)
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1,
                    help="Combine processes in parallel (default: number of CPUs)")
    ap.add_argument("--generate", nargs="+", default=["default"],
                    help="generating background functions: 'default' and/or 'altK' (the K-th "
                         "alternate of every category, its default if it has fewer)")
    ap.add_argument("--fits", nargs="+", choices=["envelope", "fixed"],
                    default=["envelope", "fixed"])
    ap.add_argument("--pull-toys", type=int, default=0,
                    help="toys per point also fitted with uncertainties (--robustFit; slow)")
    ap.add_argument("--pull-mu", nargs="*", type=float,
                    help="signal strengths for which the pulls are computed (default: all)")
    ap.add_argument("--seed", type=int, default=1000)
    ap.add_argument("--formats", nargs="+", default=["png", "pdf"])
    args = ap.parse_args()

    for g in args.generate:
        if g != "default" and not re.fullmatch(r"alt\d+", g):
            sys.exit(f"ERROR: unknown generator '{g}' (use 'default' or 'altK')")
    C.require_tool("combine")
    meta_path = Path(args.datacards)
    meta = read_json(meta_path)
    if not meta.get("combined_workspace"):
        sys.exit("ERROR: no combined workspace (make_datacards.py was run with --no-text2workspace)")
    masses = args.masses or sorted(set.intersection(
        *(set(v["masses"]) for v in meta["categories"].values())))
    lo_m, hi_m = meta["mass_range"]
    bad = [m for m in masses if not lo_m <= m <= hi_m]
    if bad:
        sys.exit(f"ERROR: masses {bad} outside the mass range of the workspace [{lo_m}, {hi_m}]")
    outdir = Path(args.outdir)
    work = outdir / "work"
    (outdir / "plots").mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)
    ROOT = setup_root(load_combine=True)
    runner = Runner(args, meta, meta_path.parent, work)
    idx_names = runner.idx_names
    pull_mu = args.pull_mu if args.pull_mu is not None else args.mu

    print(f"Expected uncertainty on r (Asimov, mu = 1, envelope) at {masses} ...", flush=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        asimov = dict(zip(masses, ex.map(runner.asimov_sigma, masses)))
    for m in masses:
        a = asimov[m]
        print(f"  M{m}: r = {a['r']:.4f} -{a['err_lo']:.4f}/+{a['err_hi']:.4f}")

    # ---- toys
    tasks = []
    for gi, gen in enumerate(args.generate):
        for mi, m in enumerate(masses):
            for ui, mu in enumerate(args.mu):
                half = R_RANGE_SIGMAS * asimov[m]["sigma"]
                rrange = (round(mu - half, 4), round(mu + half, 4))
                njobs = math.ceil(args.ntoys / args.toys_per_job)
                for j in range(njobs):
                    n = min(args.toys_per_job, args.ntoys - j * args.toys_per_job)
                    seed = args.seed + 100000 * gi + 10000 * mi + 1000 * ui + j
                    pt = args.pull_toys if (j == 0 and any(abs(mu - p) < 1e-9 for p in pull_mu)) else 0
                    tasks.append((gen, m, mu, j, n, seed, rrange, min(pt, n)))
    print(f"\n{len(tasks)} Combine jobs ({args.ntoys} toys per point, {len(args.generate)} "
          f"generator(s) x {len(masses)} masses x {len(args.mu)} signal strengths), "
          f"{args.jobs} in parallel", flush=True)
    results = {}
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(runner.job, *t): t for t in tasks}
        for k, fut in enumerate(as_completed(futs), 1):
            t = futs[fut]
            job = fut.result()
            results.setdefault(t[:3], []).append((t, collect(ROOT, job, idx_names)))
            if k % max(1, len(tasks) // 20) == 0 or k == len(tasks):
                print(f"  {k}/{len(tasks)} jobs done", flush=True)

    # ---- analysis
    payload = dict(datacards=str(meta_path), masses=masses, mu=args.mu, ntoys=args.ntoys,
                   generate=args.generate, fits=args.fits, bias_criterion=BIAS_CRITERION,
                   r_range_sigmas=R_RANGE_SIGMAS, discrete_profiling_options=C.DISCRETE_PROFILING_OPTIONS,
                   asimov=asimov, points=[])
    md = ["# Bias test of the signal strength (11 categories combined)\n",
          f"Toys generated by Combine from the combined workspace ({args.ntoys} per point); "
          "bias = (mean r_fit - mu) / std(r_fit) over the toys, uncertainty 1/sqrt(N). "
          f"Criterion |bias| < {BIAS_CRITERION} sigma.\n",
          "Expected uncertainty sigma_r (Asimov, mu = 1, envelope, robustFit):\n",
          "| M [GeV] | r (Asimov) | sigma_r |", "|---|---|---|"]
    md += [f"| {m} | {asimov[m]['r']:.4f} | -{asimov[m]['err_lo']:.4f} / +{asimov[m]['err_hi']:.4f} |"
           for m in masses]
    for gen in args.generate:
        gidx = gen_indices(meta, gen)
        pts = {}
        md += ["", f"## Generating functions: {gen}\n",
               ", ".join(f"{c}: {meta['categories'][c]['index_to_function'][str(gidx[c])]}"
                         for c in runner.cats) + "\n",
               "| M | mu | fit | toys ok | failed | at r boundary | mean r_fit | std r_fit | "
               "bias / sigma | median bias / sigma | envelope: same function as generated |",
               "|" + "---|" * 11]
        for m, mu in itertools.product(masses, args.mu):
            jobs = sorted(results.get((gen, m, mu), []), key=lambda x: x[0][3])
            half = R_RANGE_SIGMAS * asimov[m]["sigma"]
            rrange = (round(mu - half, 4), round(mu + half, 4))
            point = dict(generate=gen, mass=m, mu=mu, r_range=list(rrange), toys={}, summary={})
            for mode in args.fits:
                vals, same, chosen = [], [], {}
                for t, res in jobs:
                    for toy, rec in res.get(mode, {}).items():
                        vals.append(rec["r"])
                        if "index" in rec:
                            same.append(np.mean([i == gidx[c] for i, c in zip(rec["index"], runner.cats)]))
                            for i, c in zip(rec["index"], runner.cats):
                                f = meta["categories"][c]["index_to_function"][str(i)]
                                chosen.setdefault(c, {}).setdefault(f, 0)
                                chosen[c][f] += 1
                point["toys"][mode] = vals
                s = summarise(vals, mu, rrange)
                s["n_failed"] = args.ntoys - len(vals)
                if same:
                    s["fraction_same_function"] = float(np.mean(same))
                    s["chosen_functions"] = chosen
                point["summary"][mode] = s
                tag = f"{gen}_{mode}_M{m}_mu{fmt_mu(mu)}"
                gauss_plot(ROOT, outdir / "plots" / f"rfit_{tag}",
                           [v for v in vals if rrange[0] < v < rrange[1]],
                           f"{gen} bkg, {mode} fit, M = {m} GeV, #mu = {mu:g}", "r_{fit}",
                           [f"{s['n']} toys, bias = ({s['bias_over_sigma']:+.3f} #pm "
                            f"{s['bias_over_sigma_err']:.3f}) #sigma"], args.formats)
                md.append(f"| {m} | {mu:g} | {mode} | {s['n']} | {s['n_failed']} | "
                          f"{s['n_at_boundary']} | {s['mean']:.4f} | {s['std']:.4f} | "
                          f"{s['bias_over_sigma']:+.3f} ± {s['bias_over_sigma_err']:.3f} | "
                          f"{s['median_bias_over_sigma']:+.3f} | "
                          + (f"{100 * s['fraction_same_function']:.0f}%" if same else "-") + " |")
            pulls = []
            for t, res in jobs:
                for toy, rec in res.get("pull", {}).items():
                    p = C.pull(rec["r"], mu, rec["err_lo"], rec["err_hi"]) if rec["status"] == 0 else None
                    # an error that reaches the range boundary means the crossing was not found
                    if p is not None and rec["r"] - rec["err_lo"] > rrange[0] + 1e-6 and \
                            rec["r"] + rec["err_hi"] < rrange[1] - 1e-6:
                        pulls.append(p)
                    point["toys"].setdefault("pull_fits", []).append(rec)
            if pulls:
                st = C.mean_and_error(pulls)
                g = gauss_plot(ROOT, outdir / "plots" / f"pull_{gen}_M{m}_mu{fmt_mu(mu)}", pulls,
                               f"{gen} bkg, envelope fit, M = {m} GeV, #mu = {mu:g}",
                               "(r_{fit} - #mu) / #sigma_{r}",
                               [f"{len(pulls)} toys: mean {st['mean']:+.3f} #pm {st['mean_err']:.3f}, "
                                f"std {st['std']:.3f}"], args.formats, center=0.0, width=1.0)
                point["summary"]["pull"] = dict(st, n_requested=sum(len(r.get("pull", {})) for _, r in jobs),
                                                gauss=g)
            pts[(m, mu)] = point
            payload["points"].append(point)
        if any("pull" in p["summary"] for p in pts.values()):
            md += ["", "Pulls (envelope fits with --robustFit uncertainties):\n",
                   "| M | mu | toys | mean pull | std pull |", "|---|---|---|---|---|"]
            for (m, mu), p in pts.items():
                if "pull" in p["summary"]:
                    s = p["summary"]["pull"]
                    md.append(f"| {m} | {mu:g} | {s['n']} / {s['n_requested']} | "
                              f"{s['mean']:+.3f} ± {s['mean_err']:.3f} | {s['std']:.3f} ± {s['std_err']:.3f} |")
        summary_plot(ROOT, outdir / "plots" / f"bias_summary_{gen}", pts, masses, args.mu,
                     args.fits, args.formats)
    write_json(outdir / "bias_results.json", payload)
    (outdir / "bias_summary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    print(f"\nWritten {outdir / 'bias_results.json'}, {outdir / 'bias_summary.md'}")


if __name__ == "__main__":
    main()
