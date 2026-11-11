#!/usr/bin/env python3
"""Step 4: parametrise the signal shape and yield as a function of the resonance mass.

Only usable signal fits enter (failed fits, see fit_signal.py, are left out and listed). For
every category, each Crystal Ball parameter (mean, sigma, alphaL, nL, alphaR, nR) and
the yield inside the fit range (norm, for 1 nb) is interpolated in the mass with a natural
cubic spline (default) or a polynomial of order N ('polN'), chosen with --default-model
and --model PARAM=MODEL. Each parameter is interpolated only through the mass points where it was
floating in the fit: a value fixed at a limit (or because the tail lies outside the fit range)
is set by the fitting procedure, not measured; with fewer than two floating points, all points
are used. Beyond the outermost points used, the parameter is held constant (no extrapolation). In the workspace the interpolated parameters are clamped (CLAMP below) so that
interpolation can never produce an invalid shape: alphaL, alphaR and nL, nR to the ranges
allowed in the signal fits, sigma to at least 1 GeV and the yield to non-negative values; the
mean is not clamped. Polynomial models are fits: their plots have a pull sub-pad.

Outputs (in --outdir):
  signal_parametrisation.json         interpolation models per category and parameter
  plots/param_<CAT>_<PARAM>.{png,pdf}  fitted values vs mass with the interpolation
  plots/model_<CAT>.{png,pdf}         the parametric shape at and between the mass points,
                                      on top of the normalised templates
"""
import argparse
import sys
from array import array
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dijetfit import interp  # noqa: E402
from dijetfit.common import (SIGNAL_PROCESS, category_labels, read_json, setup_root,  # noqa: E402
                             write_json)
from dijetfit.histtools import window_hist  # noqa: E402
from dijetfit.signal import ALL_PARAMS, RANGES  # noqa: E402

# interpolated values are clamped to these ranges in the workspace (None = no limit)
CLAMP = dict(mean=(None, None), sigma=(1.0, None), alphaL=RANGES["alpha"], alphaR=RANGES["alpha"],
             nL=RANGES["n"], nR=RANGES["n"], norm=(0.0, None))
LABELS = dict(mean="mean [GeV]", sigma="#sigma [GeV]", alphaL="#alpha_{L}", nL="n_{L}",
              alphaR="#alpha_{R}", nR="n_{R}", norm="yield in fit range (1 nb)")


def clamp(v, name):
    lo, hi = CLAMP[name]
    if lo is not None:
        v = max(lo, v)
    if hi is not None:
        v = min(hi, v)
    return v


def param_plot(ROOT, path, formats, cat, label, name, masses, vals, errs, fixed, model):
    g = ROOT.TGraphErrors(len(masses), array("d", masses), array("d", vals),
                          array("d", [0.0] * len(masses)), array("d", errs))
    xs = np.linspace(min(masses), max(masses), 200)
    ys = [clamp(interp.evaluate(model, v), name) for v in xs]
    curve = ROOT.TGraph(len(xs), array("d", xs), array("d", ys))
    is_fit = model["kind"] == "poly"  # a spline interpolates: its pulls would all be zero
    c = ROOT.TCanvas("c", "", 700, 700 if is_fit else 550)
    if is_fit:
        top = ROOT.TPad("top", "", 0, 0.3, 1, 1)
        bot = ROOT.TPad("bot", "", 0, 0, 1, 0.3)
        top.SetBottomMargin(0.03)
        bot.SetTopMargin(0.03)
        bot.SetBottomMargin(0.32)
        for pd in (top, bot):
            pd.SetLeftMargin(0.15)
            pd.SetGrid()
            pd.Draw()
        top.cd()
    else:
        c.SetLeftMargin(0.15)
        c.SetGrid()
    lo = min(min(v - e for v, e in zip(vals, errs)), min(ys))
    hi = max(max(v + e for v, e in zip(vals, errs)), max(ys))
    pad = 0.15 * (hi - lo) if hi > lo else 1.0
    frame = c.DrawFrame(min(masses) - 25, lo - pad, max(masses) + 25, hi + 2 * pad)
    frame.GetXaxis().SetTitle("m_{Z'} [GeV]")
    frame.GetYaxis().SetTitle(LABELS[name])
    frame.GetYaxis().SetTitleOffset(1.4)
    curve.SetLineColor(ROOT.kRed + 1)
    curve.SetLineWidth(2)
    curve.Draw("L")
    g.SetMarkerStyle(20)
    g.Draw("P")
    gf = None
    if any(fixed):
        fx = [m for m, f in zip(masses, fixed) if f]
        fy = [v for v, f in zip(vals, fixed) if f]
        gf = ROOT.TGraph(len(fx), array("d", fx), array("d", fy))
        gf.SetMarkerStyle(24)
        gf.SetMarkerSize(1.8)
        gf.SetMarkerColor(ROOT.kBlue + 1)
        gf.Draw("P")
    t = ROOT.TPaveText(0.18, 0.80, 0.90, 0.89, "NDC")
    t.SetFillStyle(0)
    t.SetBorderSize(0)
    t.SetTextAlign(12)
    t.SetTextSize(0.033)
    desc = "natural cubic spline" if model["kind"] == "spline" else f"pol{model['order']} fit"
    t.AddText(f"{label}: {name}, {desc}"
              + ("  (open circles: fixed in the fit, not used if enough points float)" if gf else ""))
    if is_fit:
        t.AddText(f"#chi^{{2}}/ndf = {model['chi2']:.2f}/{model['ndf']}"
                  + ("" if model.get("weighted", True) else " (equal weights: some points have no error)"))
    t.Draw()
    keep = []
    if is_fit:
        # pull = (fitted value - polynomial) / uncertainty, for points with an uncertainty
        bot.cd()
        pm = [m for m, e in zip(masses, errs) if e > 0]
        pv = [(v - interp.evaluate(model, m)) / e for m, v, e in zip(masses, vals, errs) if e > 0]
        fr = bot.DrawFrame(min(masses) - 25, -4.9, max(masses) + 25, 4.9)
        fr.GetXaxis().SetTitle("m_{Z'} [GeV]")
        fr.GetXaxis().SetTitleSize(0.12)
        fr.GetXaxis().SetLabelSize(0.10)
        fr.GetYaxis().SetTitle("#frac{value - fit}{#sigma}")
        fr.GetYaxis().SetTitleSize(0.10)
        fr.GetYaxis().SetTitleOffset(0.6)
        fr.GetYaxis().SetLabelSize(0.10)
        fr.GetYaxis().SetNdivisions(505)
        if pm:
            gp = ROOT.TGraph(len(pm), array("d", pm), array("d", pv))
            gp.SetMarkerStyle(20)
            gp.Draw("P")
            keep.append(gp)
    for fmt in formats:
        c.SaveAs(f"{path}.{fmt}")
    c.Close()


def model_plot(ROOT, path, formats, cat, models, templates, xmin, xmax):
    """Parametric shape (unit area) at the template masses (solid) and in between (dashed)."""
    x = ROOT.RooRealVar("x", "m_{jj} [GeV]", xmin, xmax)
    pars = {k: ROOT.RooRealVar(k, k, 1.0, 1e-3, 1e6) for k in ALL_PARAMS}
    pdf = ROOT.RooCrystalBall("cb", "", x, pars["mean"], pars["sigma"], pars["alphaL"], pars["nL"],
                              pars["alphaR"], pars["nR"])
    nset = ROOT.RooArgSet(x)
    masses = sorted(templates)
    mids = [0.5 * (a + b) for a, b in zip(masses[:-1], masses[1:])]
    c = ROOT.TCanvas("c", "", 900, 650)
    c.SetLeftMargin(0.12)
    frame = c.DrawFrame(xmin, 0, xmax, 1)
    frame.GetXaxis().SetTitle("m_{jj} [GeV]")
    frame.GetYaxis().SetTitle("probability density [1/GeV]")
    frame.GetYaxis().SetTitleOffset(1.4)
    colors = [ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1,
              ROOT.kOrange + 2, ROOT.kCyan + 2]
    keep, ymax = [], 0.0
    leg = ROOT.TLegend(0.70, 0.55, 0.92, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for i, m in enumerate(sorted(masses + mids)):
        for k in ALL_PARAMS:
            pars[k].setVal(clamp(interp.evaluate(models[k], m), k))
        xs = np.linspace(xmin, xmax, 400)
        ys = []
        for v in xs:
            x.setVal(v)
            ys.append(pdf.getVal(nset))
        g = ROOT.TGraph(len(xs), array("d", xs), array("d", ys))
        g.SetLineColor(colors[i % len(colors)])
        g.SetLineWidth(2)
        g.SetLineStyle(1 if m in masses else 2)
        keep.append(g)
        ymax = max(ymax, max(ys))
        leg.AddEntry(g, f"M = {m:g} GeV" + ("" if m in masses else " (interp.)"), "l")
        if m in templates:
            h = templates[m].Clone(f"t{m}")
            h.Scale(1.0 / h.Integral("width"))
            h.SetMarkerStyle(20)
            h.SetMarkerSize(0.4)
            h.SetMarkerColor(colors[i % len(colors)])
            h.SetLineColor(colors[i % len(colors)])
            keep.append(h)
    frame.SetMaximum(1.3 * ymax)
    for o in keep:
        o.Draw("L" if o.InheritsFrom("TGraph") else "P same")
    leg.Draw()
    t = ROOT.TLatex(0.14, 0.92, f"{cat}: parametric signal shape (points: templates)")
    t.SetNDC()
    t.SetTextSize(0.03)
    t.Draw()
    for fmt in formats:
        c.SaveAs(f"{path}.{fmt}")
    c.Close()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("signal_fits", help="signal_fits.json from fit_signal.py")
    ap.add_argument("--input", required=True, help="ROOT file from prepare_inputs.py (for plots)")
    ap.add_argument("-o", "--outdir", default="output/signal_param")
    ap.add_argument("--default-model", default="spline", help="'spline' (default) or 'polN'")
    ap.add_argument("--model", nargs="*", metavar="PARAM=MODEL",
                    help=f"per-parameter override, e.g. norm=pol2 (parameters: {', '.join(ALL_PARAMS)})")
    ap.add_argument("--formats", nargs="+", default=["png", "pdf"])
    args = ap.parse_args()

    choice = {k: args.default_model for k in ALL_PARAMS}
    for item in args.model or []:
        try:
            k, v = item.split("=")
        except ValueError:
            sys.exit(f"ERROR: bad --model '{item}'")
        if k not in choice:
            sys.exit(f"ERROR: unknown parameter '{k}'")
        choice[k] = v
    for v in choice.values():
        interp.parse_model(v)

    ROOT = setup_root()
    fits = read_json(args.signal_fits)
    xmin, xmax = fits["fit_range"]
    tfile = ROOT.TFile.Open(args.input)
    outdir = Path(args.outdir)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)

    payload = dict(fit_range=[xmin, xmax], models_requested=choice,
                   clamp={k: list(v) for k, v in CLAMP.items()}, categories={})
    label = category_labels(fits["categories"])
    for cat, cd in sorted(fits["categories"].items()):
        fitted = sorted(int(m) for m in cd["fits"])
        masses = [m for m in fitted if cd["fits"][str(m)].get("usable", True)]
        excluded = [m for m in fitted if m not in masses]
        if excluded:
            print(f"{cat}: failed signal fits at {excluded} GeV are not used", file=sys.stderr)
        if len(masses) < 2:
            sys.exit(f"ERROR: {cat}: fewer than 2 usable signal fits")
        models = {}
        for k in ALL_PARAMS:
            vals = [cd["fits"][str(m)]["params"][k]["value"] for m in masses]
            errs = [cd["fits"][str(m)]["params"][k]["error"] for m in masses]
            fixed = [cd["fits"][str(m)]["params"][k].get("constant", False) for m in masses]
            # interpolate through the points where the parameter was determined by the fit;
            # values fixed at a limit (or because the tail is outside the fit range) are set by
            # the fitting procedure, not measured (fewer than 2 such points: use all points)
            use = [i for i, f in enumerate(fixed) if not f]
            if len(use) < 2:
                use = list(range(len(masses)))
            models[k] = interp.build(choice[k], [masses[i] for i in use], [vals[i] for i in use],
                                     [errs[i] for i in use])
            used = [masses[i] for i in use]
            if used[0] > masses[0] or used[-1] < masses[-1]:
                # determined only at some masses: constant beyond them, no extrapolation
                models[k]["hold_outside"] = [used[0], used[-1]]
            models[k].update(points=dict(masses=masses, values=vals, errors=errs, fixed=fixed,
                                         used_for_interpolation=[masses[i] for i in use]),
                             clamp=list(CLAMP[k]))
            param_plot(ROOT, str(outdir / "plots" / f"param_{cat}_{k}"), args.formats, cat, label[cat], k,
                       masses, vals, errs, fixed, models[k])
        templates = {}
        for m in masses:
            raw = tfile.Get(f"{cat}/{SIGNAL_PROCESS}_{m}")
            templates[m] = window_hist(ROOT, raw, xmin, xmax, cd["rebin"], name=f"{cat}_{m}_t")
        model_plot(ROOT, str(outdir / "plots" / f"model_{cat}"), args.formats, cat, models,
                   templates, xmin, xmax)
        payload["categories"][cat] = dict(masses=masses, excluded_masses=excluded, rebin=cd["rebin"],
                                          models=models)
        print(f"{cat}: masses {masses}, models: "
              + ", ".join(f"{k}={choice[k]}" for k in ALL_PARAMS))
    write_json(outdir / "signal_parametrisation.json", payload)
    print(f"Written {outdir / 'signal_parametrisation.json'}")


if __name__ == "__main__":
    main()
