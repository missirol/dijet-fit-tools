#!/usr/bin/env python3
"""Step 4: parametrise the signal shape and yield as a function of the resonance mass.

For every category, each Crystal Ball parameter (mean, sigma, alphaL, nL, alphaR, nR) and
the yield inside the fit range (norm, for 1 nb) is interpolated in the mass with a natural
cubic spline (default) or a polynomial of order N ('polN'), chosen with --default-model
and --model PARAM=MODEL. Parameters fixed in the fits (see signal_fits.json) enter with their
fixed values. In the workspace the interpolated parameters are clamped to the ranges allowed
in the fits (e.g. n >= 0.5), so that interpolation can never produce an invalid shape.

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
    c = ROOT.TCanvas("c", "", 700, 550)
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
              + ("  (open circles: fixed in the fit)" if gf else ""))
    t.Draw()
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
        masses = sorted(int(m) for m in cd["fits"])
        if len(masses) < 2:
            sys.exit(f"ERROR: {cat}: fewer than 2 mass points")
        models = {}
        for k in ALL_PARAMS:
            vals = [cd["fits"][str(m)]["params"][k]["value"] for m in masses]
            errs = [cd["fits"][str(m)]["params"][k]["error"] for m in masses]
            fixed = [cd["fits"][str(m)]["params"][k].get("constant", False) for m in masses]
            models[k] = interp.build(choice[k], masses, vals, errs)
            models[k].update(points=dict(masses=masses, values=vals, errors=errs, fixed=fixed),
                             clamp=list(CLAMP[k]))
            param_plot(ROOT, str(outdir / "plots" / f"param_{cat}_{k}"), args.formats, cat, label[cat], k,
                       masses, vals, errs, fixed, models[k])
        templates = {}
        for m in masses:
            raw = tfile.Get(f"{cat}/{SIGNAL_PROCESS}_{m}")
            templates[m] = window_hist(ROOT, raw, xmin, xmax, cd["rebin"], name=f"{cat}_{m}_t")
        model_plot(ROOT, str(outdir / "plots" / f"model_{cat}"), args.formats, cat, models,
                   templates, xmin, xmax)
        payload["categories"][cat] = dict(masses=masses, rebin=cd["rebin"], models=models)
        print(f"{cat}: masses {masses}, models: "
              + ", ".join(f"{k}={choice[k]}" for k in ALL_PARAMS))
    write_json(outdir / "signal_parametrisation.json", payload)
    print(f"Written {outdir / 'signal_parametrisation.json'}")


if __name__ == "__main__":
    main()
