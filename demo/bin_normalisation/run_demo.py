#!/usr/bin/env python3
"""Standalone demonstration: normalisation of pdfs evaluated at the bin centres.

Artificial data are generated with a known truth, then fitted (signal + background, signal
strength r floating, background parameters floating) with four models that differ only in how
the expected counts per bin are computed:

  OLD       pdf at the bin centre, normalised by its integral, for background and signal
            (RooFit's and Combine's default for binned data)
  BKG_ONLY  background normalised by its sum over the bin centres, signal by its integral
  NEW       background and signal normalised by their sums over the bin centres (the workflow)
  EXACT     exact integrals over each bin (Simpson rule for the background, exact bin integrals
            of the fixed signal shape), i.e. the model of Combine's RooParametricShapeBinPdf

The truth is always the exact one (bin integrals of the generating functions), so a model is
unbiased if it returns the injected r. See demo/bin_normalisation/README.md.

Steps (run all with `all`, or one at a time):
  generate  write OUTDIR/data.root: per mass and mu, Asimov histograms (and --toys Poisson
            toys) for every category, plus the truth in OUTDIR/truth.json
  fit       fit every data set with every model; OUTDIR/results.json, OUTDIR/results.md
  plot      OUTDIR/plots: bias vs mass, differences of the fitted predictions, fit plots

Usage:  python3 demo/bin_normalisation/run_demo.py all -o demo_output
"""
import argparse
import json
import math
import sys
from array import array
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from dijetfit import background as B  # noqa: E402
from dijetfit.common import SQRTS  # noqa: E402

XMIN, XMAX, NBINS = 150.0, 600.0, 450  # 1 GeV bins, as in the analysis
MODELS = ("OLD", "BKG_ONLY", "NEW", "EXACT")
# Background truth: the dijet function (1-x)^p1 / x^(p2 + p3 ln x), x = m / 13.6 TeV, with
# parameters that reproduce the background fitted to one real |dy| category to within 4%
# (falling by a factor ~1200 between 150 and 600 GeV); ~2e7 events per category, as in the data.
BKG_P = (188.32, -11.379, -1.5612)
BKG_EVENTS = 2.06e7
# Signal truth: double-sided Crystal Ball, parameters and yield per category for mu = 1
# (sigma = 1 nb) from the fits of the real signal templates of one category.
SIGNALS = {
    200: dict(mean=203.2, sigma=28.41, alphaL=0.384, nL=100.0, alphaR=1.792, nR=1.627, events=8250.0),
    300: dict(mean=298.8, sigma=36.95, alphaL=0.467, nL=100.0, alphaR=2.538, nR=0.5, events=11149.0),
    400: dict(mean=398.3, sigma=39.03, alphaL=0.410, nL=100.0, alphaR=2.563, nR=0.5, events=12046.0),
    500: dict(mean=492.9, sigma=41.00, alphaL=0.346, nL=100.0, alphaR=7.679, nR=10.0, events=11838.0),
}
CB_KEYS = ("mean", "sigma", "alphaL", "nL", "alphaR", "nR")


# ---- truth (numpy) ----------------------------------------------------------------------------
def edges():
    return np.linspace(XMIN, XMAX, NBINS + 1)


def centres():
    e = edges()
    return 0.5 * (e[:-1] + e[1:])


def dijet(m, p1, p2, p3):
    x = np.asarray(m, float) / SQRTS
    return np.exp(p1 * np.log1p(-x) - (p2 + p3 * np.log(x)) * np.log(x))


def crystal_ball(m, mean, sigma, alphaL, nL, alphaR, nR):
    """Same definition as RooCrystalBall (unnormalised)."""
    t = (np.asarray(m, float) - mean) / sigma
    out = np.exp(-0.5 * t * t)
    left, right = t < -alphaL, t > alphaR
    out[left] = np.exp(nL * np.log(nL / alphaL) - 0.5 * alphaL ** 2 - nL * np.log(nL / alphaL - alphaL - t[left]))
    out[right] = np.exp(nR * np.log(nR / alphaR) - 0.5 * alphaR ** 2 - nR * np.log(nR / alphaR - alphaR + t[right]))
    return out


def bin_integrals(func, nsub=40):
    """Exact (to ~1e-12) integral of func over each bin: composite Simpson rule."""
    e = edges()
    lo, hi = e[:-1], e[1:]
    k = np.arange(nsub + 1)
    wgt = np.where((k == 0) | (k == nsub), 1.0, np.where(k % 2 == 1, 4.0, 2.0))
    pts = lo[:, None] + (hi - lo)[:, None] * k[None, :] / nsub
    return (func(pts) * wgt).sum(axis=1) * (hi - lo) / (3 * nsub)


def truth_bkg():
    """Expected background counts per bin (exact bin integrals), adding up to BKG_EVENTS."""
    s = bin_integrals(lambda m: dijet(m, *BKG_P))
    return BKG_EVENTS * s / s.sum()


def truth_sig(mass):
    """Expected signal counts per bin for mu = 1 (exact bin integrals; the yield is the number of
    signal events inside the fit range)."""
    sp = SIGNALS[mass]
    s = bin_integrals(lambda m: crystal_ball(m, **{k: sp[k] for k in CB_KEYS}))
    return sp["events"] * s / s.sum()


def midpoint_mismatch(func):
    """sum over the bin centres x bin width / integral - 1: what the OLD normalisation gets wrong."""
    s = func(centres()).sum() * (XMAX - XMIN) / NBINS
    return s / bin_integrals(func).sum() - 1.0


# ---- ROOT models ------------------------------------------------------------------------------
def setup_root():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    return ROOT


def to_hist(ROOT, name, counts):
    h = ROOT.TH1D(name, "", NBINS, XMIN, XMAX)
    h.SetDirectory(0)
    for i, v in enumerate(counts, start=1):
        h.SetBinContent(i, v)
        h.SetBinError(i, math.sqrt(max(v, 0.0)))
    return h


def simpson_formula(formula):
    """Formula whose value at a bin centre is the average of `formula` over the 1 GeV bin
    (Simpson rule; relative error ~ f''''/f / 2880 ~ 1e-10 here)."""
    shift = lambda d: formula.replace("@0", f"(@0{d:+.1f})")  # noqa: E731  (@0 = observable)
    return f"(({shift(-0.5)})+4*({formula})+({shift(+0.5)}))/6"


class Category:
    """One category: observable, data set and the four fit models (signal shape fixed)."""

    def __init__(self, ROOT, tag, hist, mass):
        self.ROOT, self.tag, self.mass = ROOT, tag, mass
        self.x = ROOT.RooRealVar(f"mjj_{tag}", "m_{jj} [GeV]", XMIN, XMAX)
        self.x.setBins(NBINS)
        self.data = ROOT.RooDataHist(f"data_{tag}", "", ROOT.RooArgList(self.x), hist)
        self.n = np.array([hist.GetBinContent(i) for i in range(1, NBINS + 1)])
        self.keep = [hist, self.data]
        # background: dijet function in the decorrelated basis of the workflow (start values
        # from the data; the same family as the truth, so any bias comes from the normalisation)
        m0 = B.reference_mass(XMIN, XMAX)
        pvals, mu = B.start_values("dijet4", hist, m0)
        M, tbar = B.decorrelated_basis("dijet4", hist, mu, m0)
        self.q0 = np.linalg.solve(M, pvals)
        self.formula = B.formula_of("dijet4", m0, M, tbar)
        sp = SIGNALS[mass]
        self.cb_consts = [ROOT.RooRealVar(f"{k}_{tag}", k, sp[k]) for k in CB_KEYS]
        self.sig_events = sp["events"]
        self.sig_template = to_hist(ROOT, f"sigtemplate_{tag}", truth_sig(mass))
        self.keep += self.cb_consts + [self.sig_template]
        self.models = {m: self.build(m) for m in MODELS}

    def build(self, model):
        R, tag, x = self.ROOT, f"{self.tag}_{model}", self.x
        params = [R.RooRealVar(f"q{k}_{tag}", "", v, v - 50, v + 50) for k, v in enumerate(self.q0, 1)]
        args = R.RooArgList(x)
        for p in params:
            args.add(p)
        if model == "OLD":
            bkg = B.make_pdf(R, f"bkg_{tag}", "dijet4", args, self.formula)  # integral normalisation
        elif model == "EXACT":
            bkg = B.make_pdf(R, f"bkg_{tag}", "dijet4", args, simpson_formula(self.formula), NBINS)
        else:
            bkg = B.make_pdf(R, f"bkg_{tag}", "dijet4", args, self.formula, NBINS)  # bin sum
        if model == "EXACT":  # exact bin integrals of the (fixed) signal shape
            dh = R.RooDataHist(f"sigdh_{tag}", "", R.RooArgList(x), self.sig_template)
            sig = R.RooHistPdf(f"sig_{tag}", "", R.RooArgSet(x), dh, 0)
            self.keep.append(dh)
        else:
            sig = R.RooCrystalBall(f"sig_{tag}", "", x, *self.cb_consts)
            if model == "NEW":
                B.set_bin_sum_normalisation(sig, NBINS)
        nbkg = R.RooRealVar(f"nbkg_{tag}", "", self.n.sum(), 0.5 * self.n.sum(), 1.5 * self.n.sum())
        self.keep += params + [args, bkg, sig, nbkg]
        return dict(params=params, bkg=bkg, sig=sig, nbkg=nbkg)

    def nll(self, model, r):
        R, mm = self.ROOT, self.models[model]
        nsig = R.RooFormulaVar(f"nsig_{self.tag}_{model}", "@0*@1",
                               R.RooArgList(r, R.RooFit.RooConst(self.sig_events)))
        pdf = R.RooAddPdf(f"sb_{self.tag}_{model}", "", R.RooArgList(mm["sig"], mm["bkg"]),
                          R.RooArgList(nsig, mm["nbkg"]))
        nll = pdf.createNLL(self.data, R.RooFit.Extended(True))
        self.keep += [nsig, pdf, nll]
        return nll

    def reset(self, model):
        mm = self.models[model]
        for p, v in zip(mm["params"], self.q0):
            p.setVal(v)
        mm["nbkg"].setVal(self.n.sum())

    def predictions(self, model, r):
        """Expected counts per bin as RooFit/Combine compute them: N * pdf(centre) * width."""
        mm = self.models[model]
        ns = self.ROOT.RooArgSet(self.x)
        sig, bkg = np.empty(NBINS), np.empty(NBINS)
        for i, c in enumerate(centres()):
            self.x.setVal(c)
            bkg[i] = mm["nbkg"].getVal() * mm["bkg"].getVal(ns)
            sig[i] = r * self.sig_events * mm["sig"].getVal(ns)
        return sig, bkg


def fit_dataset(ROOT, hists, mass, mu, tag):
    """Fit one data set (all categories) with every model. Returns per-model results."""
    cats = [Category(ROOT, f"{tag}_c{i}", h, mass) for i, h in enumerate(hists)]
    out = {}
    for model in MODELS:
        r = ROOT.RooRealVar(f"r_{tag}_{model}", "r", mu + 0.5, -100.0, 100.0)
        for c in cats:
            c.reset(model)
        nlls = ROOT.RooArgList()
        for c in cats:
            nlls.add(c.nll(model, r))
        total = ROOT.RooAddition(f"nll_{tag}_{model}", "", nlls)
        mn = ROOT.RooMinimizer(total)
        mn.setPrintLevel(-1)
        mn.setStrategy(2)
        mn.setOffsetting(True)
        mn.setEps(1e-6)
        mn.minimize("Minuit2", "migrad")
        mn.minimize("Minuit2", "migrad")
        mn.hesse()
        fr = mn.save()
        sig0, bkg0 = cats[0].predictions(model, r.getVal())
        out[model] = dict(
            r=r.getVal(), err=r.getError(), status=fr.status(), cov_qual=fr.covQual(),
            # sum of the per-bin predictions minus the yield used in the extended term
            mismatch_bkg=[float(c.predictions(model, r.getVal())[1].sum() - c.models[model]["nbkg"].getVal())
                          for c in cats],
            mismatch_sig=float(sig0.sum() - r.getVal() * cats[0].sig_events),
            cat0=dict(sig=sig0.tolist(), bkg=bkg0.tolist(),
                      q=[p.getVal() for p in cats[0].models[model]["params"]]))
        ROOT.SetOwnership(total, True)
    # OLD model at the parameters of the EXACT fit (= the truth for Asimov data), no refit: the
    # direct effect of the normalisation, before the fit distorts the shape to compensate for it
    c0, old, exact = cats[0], cats[0].models["OLD"], cats[0].models["EXACT"]
    for p_old, p_ex in zip(old["params"], exact["params"]):
        p_old.setVal(p_ex.getVal())
    old["nbkg"].setVal(exact["nbkg"].getVal())
    sig0, bkg0 = c0.predictions("OLD", out["EXACT"]["r"])
    out["OLD_no_refit"] = dict(cat0=dict(sig=sig0.tolist(), bkg=bkg0.tolist()))
    return out


# ---- steps ------------------------------------------------------------------------------------
def step_generate(args):
    ROOT = setup_root()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    f = ROOT.TFile(str(out / "data.root"), "RECREATE")
    b = truth_bkg()
    truth = dict(bkg_params=BKG_P, bkg_events=BKG_EVENTS, signals=SIGNALS, ncat=args.ncat,
                 masses=args.masses, mu=args.mu, toys=args.toys,
                 mismatch_bkg=midpoint_mismatch(lambda m: dijet(m, *BKG_P)),
                 mismatch_sig={m: midpoint_mismatch(lambda x, sp=SIGNALS[m]: crystal_ball(x, **{k: sp[k] for k in CB_KEYS}))
                               for m in args.masses})
    for mass in args.masses:
        s = truth_sig(mass)
        for mu in args.mu:
            expected = b + mu * s
            d = f.mkdir(f"M{mass}_mu{mu:g}")
            d.cd()
            for c in range(args.ncat):
                to_hist(ROOT, f"asimov_c{c}", expected).Write()
                for t in range(args.toys):
                    to_hist(ROOT, f"toy{t}_c{c}", rng.poisson(expected).astype(float)).Write()
    f.Close()
    (out / "truth.json").write_text(json.dumps(truth, indent=1))
    print(f"Written {out / 'data.root'}: {len(args.masses)} masses x {len(args.mu)} mu, {args.ncat} "
          f"categories, Asimov + {args.toys} toys each")
    print(f"Midpoint mismatch (sum over bin centres x width / integral - 1): background "
          f"{truth['mismatch_bkg']:+.2e}; signal " +
          ", ".join(f"M{m} {v:+.1e}" for m, v in truth["mismatch_sig"].items()))


def step_fit(args):
    ROOT = setup_root()
    out = Path(args.outdir)
    truth = json.loads((out / "truth.json").read_text())
    f = ROOT.TFile.Open(str(out / "data.root"))
    res = dict(truth=truth, fits={})
    md = ["# Bin-normalisation demo: fit results\n",
          f"{truth['ncat']} categories, ~{BKG_EVENTS:.2g} background events each; truth = exact bin "
          "integrals. Bias = (r_fit - mu) / sigma_r (sigma_r from the fit).\n",
          "| M | mu | data | " + " | ".join(MODELS) + " |", "|---|---|---|" + "---|" * len(MODELS)]
    for mass in truth["masses"]:
        for mu in truth["mu"]:
            d = f.Get(f"M{mass}_mu{mu:g}")
            sets = ["asimov"] + [f"toy{t}" for t in range(truth["toys"])]
            point = {}
            for s in sets:
                hists = [d.Get(f"{s}_c{c}") for c in range(truth["ncat"])]
                for h in hists:
                    h.SetDirectory(0)
                point[s] = fit_dataset(ROOT, hists, mass, mu, f"M{mass}_mu{mu:g}_{s}")
                row = point[s]
                print(f"M{mass} mu={mu:g} {s:7s} " + "  ".join(
                    f"{m}: r={row[m]['r']:+.4f} ({(row[m]['r'] - mu) / row[m]['err']:+.3f}s)" for m in MODELS),
                    flush=True)
            res["fits"][f"{mass}_{mu:g}"] = point
            a = point["asimov"]
            md.append(f"| {mass} | {mu:g} | Asimov | " + " | ".join(
                f"{(a[m]['r'] - mu) / a[m]['err']:+.3f}" for m in MODELS) + " |")
            if truth["toys"]:
                cells = []
                for m in MODELS:
                    b = np.array([(point[s][m]["r"] - mu) / point[s][m]["err"] for s in sets[1:]])
                    cells.append(f"{b.mean():+.3f} ± {b.std(ddof=1) / math.sqrt(len(b)):.3f}" if len(b) > 1
                                 else f"{b.mean():+.3f}")
                md.append(f"| {mass} | {mu:g} | {truth['toys']} toys (mean) | " + " | ".join(cells) + " |")
    md += ["", "Per-bin predictions minus the yield in the extended term (Asimov, mu = first value, "
           "category 0):", "", "| M | " + " | ".join(f"{m} bkg / sig" for m in MODELS) + " |",
           "|---|" + "---|" * len(MODELS)]
    for mass in truth["masses"]:
        a = res["fits"][f"{mass}_{truth['mu'][-1]:g}"]["asimov"]
        md.append(f"| {mass} | " + " | ".join(f"{a[m]['mismatch_bkg'][0]:+.2f} / {a[m]['mismatch_sig']:+.4f}"
                                             for m in MODELS) + " |")
    (out / "results.json").write_text(json.dumps(res))
    (out / "results.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


def step_plot(args):
    ROOT = setup_root()
    ROOT.gStyle.SetOptStat(0)
    from dijetfit.plotting import fit_plot
    out = Path(args.outdir)
    res = json.loads((out / "results.json").read_text())
    truth = res["truth"]
    pdir = out / "plots"
    pdir.mkdir(exist_ok=True)
    col = dict(OLD=ROOT.kRed + 1, BKG_ONLY=ROOT.kOrange + 1, NEW=ROOT.kBlue + 1, EXACT=ROOT.kGreen + 2)
    mk = dict(OLD=20, BKG_ONLY=22, NEW=21, EXACT=24)
    keep = []
    xs = centres()
    exact_n = truth_bkg()

    # 1) bias vs mass (Asimov, and toys if any)
    c = ROOT.TCanvas("cb", "", 1000, 650)
    c.SetGridy()
    c.SetLeftMargin(0.1)
    masses = truth["masses"]
    fr = c.DrawFrame(min(masses) - 30, -3, max(masses) + 30, 3, ";signal mass [GeV];(#hat{r} - #mu) / #sigma_{r}")
    leg = ROOT.TLegend(0.12, 0.14, 0.6, 0.30)
    leg.SetNColumns(2)
    leg.SetBorderSize(0)
    for k, m in enumerate(MODELS):
        for j, mu in enumerate(truth["mu"]):
            pts = []
            for mass in masses:
                a = res["fits"][f"{mass}_{mu:g}"]["asimov"][m]
                pts.append(((a["r"] - mu) / a["err"]))
            g = ROOT.TGraph(len(masses), array("d", [ms + (k - 1.5) * 6 + j * 1.5 for ms in masses]), array("d", pts))
            g.SetMarkerStyle(mk[m])
            g.SetMarkerColor(col[m])
            g.SetMarkerSize(1.2)
            g.Draw("P same")
            keep.append(g)
            if j == 0:
                leg.AddEntry(g, m, "p")
    leg.Draw()
    t = ROOT.TLatex(0.12, 0.93, f"Asimov data (truth: exact bin integrals), #mu = {', '.join(f'{u:g}' for u in truth['mu'])}, "
                                f"{truth['ncat']} categories")
    t.SetNDC()
    t.SetTextSize(0.035)
    t.Draw()
    keep += [fr, leg, t]
    c.SaveAs(str(pdir / "bias_vs_mass.png"))

    # 2) differences of the fitted background predictions (category 0, Asimov, first mass, mu=0)
    mass, mu = masses[0], truth["mu"][0]
    a = res["fits"][f"{mass}_{mu:g}"]["asimov"]
    c2 = ROOT.TCanvas("cd", "", 1000, 800)
    c2.Divide(1, 2)
    panels = (("(fitted prediction - EXACT) / #sigma_{stat}", 1.0, ("OLD", "OLD_no_refit", "BKG_ONLY", "NEW")),
              ("(fitted prediction - EXACT) / #sigma_{stat} #times 10^{3}", 1e3, ("BKG_ONLY", "NEW")))
    labels = dict(OLD="OLD (fitted)", OLD_no_refit="OLD at the true parameters (no refit)",
                  BKG_ONLY="BKG_ONLY (fitted)", NEW="NEW (fitted)")
    nref = np.array(a["EXACT"]["cat0"]["bkg"]) + np.array(a["EXACT"]["cat0"]["sig"])
    for ip, (ttl, scale, models) in enumerate(panels):
        p = c2.cd(ip + 1)
        p.SetGridy()
        p.SetLeftMargin(0.1)
        p.SetRightMargin(0.03)
        curves = {m: (np.array(a[m]["cat0"]["bkg"]) + np.array(a[m]["cat0"]["sig"]) - nref) / np.sqrt(exact_n) * scale
                  for m in models}
        hi = 1.15 * max(np.abs(v).max() for v in curves.values())
        frm = p.DrawFrame(XMIN, -hi, XMAX, 1.6 * hi, f";m_{{jj}} [GeV];{ttl}")
        frm.GetYaxis().SetTitleSize(0.045)
        frm.GetYaxis().SetTitleOffset(0.9)
        frm.GetXaxis().SetTitleSize(0.045)
        keep.append(frm)
        lg = ROOT.TLegend(0.12, 0.72, 0.97, 0.89)
        lg.SetNColumns(2)
        lg.SetBorderSize(0)
        lg.SetFillStyle(0)
        lg.SetTextSize(0.045)
        for m in models:
            g = ROOT.TGraph(NBINS, array("d", xs), array("d", curves[m]))
            g.SetLineColor(col.get(m, ROOT.kRed + 1))
            g.SetLineStyle(2 if m == "OLD_no_refit" else 1)
            g.SetLineWidth(2)
            g.Draw("L same")
            lg.AddEntry(g, labels[m], "l")
            keep.append(g)
        lg.Draw()
        keep.append(lg)
    c2.cd(1)
    tt = ROOT.TLatex(0.12, 0.93, f"Asimov, M = {mass} GeV, #mu = {mu:g}, category 0; #sigma_{{stat}} = #sqrt{{N}} per bin")
    tt.SetNDC()
    tt.SetTextSize(0.045)
    tt.Draw()
    keep.append(tt)
    c2.SaveAs(str(pdir / f"prediction_differences_M{mass}_mu{mu:g}.png"))

    # 3) fit plots with pull sub-pads (category 0, Asimov) for every model and data set
    f = ROOT.TFile.Open(str(out / "data.root"))
    for key, point in res["fits"].items():
        mass_s, mu_s = key.split("_")
        h = f.Get(f"M{mass_s}_mu{mu_s}/asimov_c0")
        for m in MODELS:
            fit = point["asimov"][m]
            nu = np.array(fit["cat0"]["bkg"]) + np.array(fit["cat0"]["sig"])
            fit_plot(ROOT, str(pdir / f"fit_M{mass_s}_mu{mu_s}_{m}"), ["png"], h, nu.tolist(),
                     (xs.tolist(), nu.tolist()),
                     [f"{m}: M = {mass_s} GeV, #mu = {mu_s} (Asimov, category 0)",
                      f"#hat{{r}} = {fit['r']:.4f} #pm {fit['err']:.4f}",
                      f"#Sigma#nu_{{bkg}} - N_{{bkg}} = {fit['mismatch_bkg'][0]:+.1f} events"],
                     logy=True, data_label="Asimov data (exact)", fit_label=f"S+B fit, {m}")
    print(f"Plots in {pdir}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["generate", "fit", "plot", "all"])
    ap.add_argument("-o", "--outdir", default="demo_output")
    ap.add_argument("--masses", nargs="+", type=int, default=[200, 300, 400, 500], choices=sorted(SIGNALS))
    ap.add_argument("--mu", nargs="+", type=float, default=[0.0, 1.0])
    ap.add_argument("--ncat", type=int, default=11,
                    help="identical, independent categories fitted together (r common); default 11 as "
                         "in the analysis")
    ap.add_argument("--toys", type=int, default=0, help="Poisson toys per point, besides Asimov")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    if args.step in ("generate", "all"):
        step_generate(args)
    if args.step in ("fit", "all"):
        step_fit(args)
    if args.step in ("plot", "all"):
        step_plot(args)


if __name__ == "__main__":
    main()
