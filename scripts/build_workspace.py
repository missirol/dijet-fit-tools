#!/usr/bin/env python3
"""Step 5: write the RooWorkspace for the statistical analysis (one file, all categories).

Reads the prepared histograms, the background fit results and the signal parametrisation, and
writes a single ROOT file with one RooWorkspace ('w') holding, for every category CAT, the
objects described in docs/OUTPUT_FORMAT.md (CMS Combine conventions):
  mjj_CAT                      observable (fit range and binning of the category)
  data_obs_CAT                 RooDataHist
  ZprimeQQ_CAT, ZprimeQQ_CAT_norm   parametric signal pdf and yield (functions of MH)
  bkg_CAT, bkg_CAT_norm        RooMultiPdf of the passing background functions, and yield
  pdfindex_CAT                 RooCategory selecting the background function
  bkg_CAT_<func>               the individual background pdfs
After writing, the file is re-opened and checked: the signal parameters and yields in the
workspace must reproduce the parametrisation at the mass points, and the background pdfs the
fitted chi2 values. The deviation of the parametrised yield from the fitted yields (zero for
splines, non-zero for polynomial fits) is reported.

RooMultiPdf is a class of CMS Combine; the Combine library must be loadable. With
--no-multipdf the RooMultiPdf is not written (all other objects are).

The background and signal pdfs are normalised by their sum over the bin centres (RooBinIntegrator, stored
with each pdf; see dijetfit/background.make_pdf), as in the fits: Combine predicts N x pdf(bin
centre) x bin width per bin, and these predictions then add up exactly to the yield N.

Categories with bins wider than the original 1 GeV (rebinned by fit_background.py) are refused:
for such bins the fits integrate the pdfs over each bin, while Combine evaluates them at the bin
centres, so the fitted parameters would not describe the model that Combine uses.

The mass variable MH covers the mass range common to all categories (no extrapolation of any
category's signal model).
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dijetfit import background as B  # noqa: E402
from dijetfit import interp  # noqa: E402
from dijetfit.common import (DATA_NAME, SIGNAL_PROCESS, capture_output, count_messages,  # noqa: E402
                             message_summary, read_json, setup_root, write_json)
from dijetfit.histtools import window_hist  # noqa: E402
from dijetfit.signal import ALL_PARAMS, CB_PARAMS  # noqa: E402

MASS_VAR = "MH"


def wimport(w, obj, ROOT):
    if obj.InheritsFrom("RooAbsData"):
        getattr(w, "import")(obj, ROOT.RooFit.Silence())
    else:
        getattr(w, "import")(obj, ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="ROOT file from prepare_inputs.py")
    ap.add_argument("--background", required=True, help="background_fits.json")
    ap.add_argument("--signal", required=True, help="signal_parametrisation.json")
    ap.add_argument("-o", "--output", default="output/workspace/dijet_workspace.root")
    ap.add_argument("--workspace-name", default="w")
    ap.add_argument("--no-multipdf", action="store_true",
                    help="do not write the RooMultiPdf (no CMS Combine library needed)")
    args = ap.parse_args()

    use_multipdf = not args.no_multipdf
    ROOT = setup_root(load_combine=use_multipdf)
    bkg = read_json(args.background)
    sig = read_json(args.signal)
    if bkg["fit_range"] != sig["fit_range"]:
        sys.exit(f"ERROR: fit ranges differ: background {bkg['fit_range']}, signal {sig['fit_range']}")
    xmin, xmax = bkg["fit_range"]
    cats = sorted(bkg["categories"])
    if sorted(sig["categories"]) != cats:
        sys.exit("ERROR: the background and signal files have different categories")
    tfile = ROOT.TFile.Open(args.input)

    wide = {}
    for c in cats:
        e = bkg["categories"][c]["bin_edges"]
        lo, hi = np.array(e[:-1]), np.array(e[1:])
        if B.needs_bin_integration(lo, hi):
            wide[c] = bkg["categories"][c]["binning"]
    if wide:
        sys.exit("ERROR: categories with bins wider than 1 GeV: "
                 + "; ".join(f"{c} ({b})" for c, b in wide.items())
                 + ". The background and signal fits integrated the pdfs over these bins, but the "
                 "RooWorkspace written by this script is used with the bin-centre evaluation that "
                 "Combine applies to binned data, i.e. with a different model than the one fitted. "
                 "Only 1 GeV bins are supported for the workspace.")

    # MH: the mass range covered by every category (no extrapolation of any signal model)
    lo_m = max(min(sig["categories"][c]["masses"]) for c in cats)
    hi_m = min(max(sig["categories"][c]["masses"]) for c in cats)
    if hi_m <= lo_m:
        sys.exit(f"ERROR: the categories have no common signal mass range ({lo_m} - {hi_m})")
    if any(sig["categories"][c]["masses"] != sig["categories"][cats[0]]["masses"] for c in cats):
        print(f"WARNING: the categories have different signal mass points; MH is limited to the "
              f"common range [{lo_m}, {hi_m}] GeV", file=sys.stderr)
    all_masses = [lo_m, hi_m]
    mh = ROOT.RooRealVar(MASS_VAR, "m_{Z'} [GeV]", 0.5 * (lo_m + hi_m), lo_m, hi_m)
    # the mass hypothesis is a fixed parameter: Combine sets it with -m; if it were floating,
    # Combine would fit it as a free parameter
    mh.setConstant(True)
    w = ROOT.RooWorkspace(args.workspace_name, "dijet resonance search: signal and background models")
    keep = [mh]
    summary = dict(file=str(args.output), workspace=args.workspace_name, mass_variable=MASS_VAR,
                   mass_range=[all_masses[0], all_masses[-1]], fit_range=[xmin, xmax],
                   multipdf=use_multipdf, categories={})

    with capture_output() as log:
        for c in cats:
            cb, cs = bkg["categories"][c], sig["categories"][c]
            if cb["rebin"] != cs["rebin"]:
                sys.exit(f"ERROR: {c}: background and signal binning differ")
            if not cb["passing"]:
                sys.exit(f"ERROR: {c}: no passing background function")
            edges = cb["bin_edges"]
            nbins = B.bin_sum_nbins(np.array(edges[:-1]), np.array(edges[1:]))
            if not nbins:
                sys.exit(f"ERROR: {c}: the bins are not uniform")
            x = ROOT.RooRealVar(f"mjj_{c}", "m_{jj} [GeV]", xmin, xmax)
            x.setBinning(ROOT.RooBinning(len(edges) - 1, np.array(edges, dtype="d")))
            h = window_hist(ROOT, tfile.Get(f"{c}/{DATA_NAME}"), xmin, xmax, cb["rebin"],
                            name=f"{c}_data_ws")
            data = ROOT.RooDataHist(f"{DATA_NAME}_{c}", f"data, {c}", ROOT.RooArgList(x), h)
            wimport(w, data, ROOT)
            keep += [x, h, data]

            # ---- signal
            sp = f"{SIGNAL_PROCESS}_{c}"
            fv = {}
            for k in ALL_PARAMS:
                model = cs["models"][k]
                lo, hi = model["clamp"]
                expr = interp.clamped_formula(model, lo, hi, "@0")
                name = f"{sp}_norm" if k == "norm" else f"{sp}_{k}"
                fv[k] = ROOT.RooFormulaVar(name, f"{k}(MH), {c}", expr, ROOT.RooArgList(mh))
                keep.append(fv[k])
            spdf = ROOT.RooCrystalBall(sp, f"Z'->qq signal, {c}", x, *[fv[k] for k in CB_PARAMS])
            B.set_bin_sum_normalisation(spdf, nbins)  # as the background and the signal fits
            wimport(w, spdf, ROOT)
            wimport(w, fv["norm"], ROOT)
            keep.append(spdf)

            # ---- background
            pdfs = ROOT.RooArgList()
            cat_index = ROOT.RooCategory(
                f"pdfindex_{c}", "background function index: "
                + ", ".join(f"{i}={f}" for i, f in enumerate(cb["passing"])))
            bnames = []
            for i, func in enumerate(cb["passing"]):
                fr = cb["fits"][func]
                arglist = ROOT.RooArgList(x)
                for k in range(1, fr["npar"]):
                    pr = fr["params"][f"p{k}"]
                    v = ROOT.RooRealVar(f"{func}_p{k}_{c}", f"{func} q{k}", pr["value"],
                                        pr["value"] - B.QRANGE, pr["value"] + B.QRANGE)
                    v.setError(pr["error"])
                    arglist.add(v)
                    keep.append(v)
                # normalised by the sum over the bin centres, as in the fits (see B.make_pdf)
                pdf = B.make_pdf(ROOT, f"bkg_{c}_{func}", func, arglist, fr["formula"], nbins)
                pdfs.add(pdf)
                if not use_multipdf:  # RooMultiPdf defines the category states itself
                    cat_index.defineType(func, i)
                bnames.append(func)
                keep += [pdf, arglist]
            default_fn = cb["default"]
            nb = cb["fits"][default_fn]["norm"]
            bnorm = ROOT.RooRealVar(f"bkg_{c}_norm", f"background yield, {c}", nb["value"],
                                    0.0, 2.0 * nb["value"])
            bnorm.setError(nb["error"])
            keep += [cat_index, bnorm]
            if use_multipdf:
                mp = ROOT.RooMultiPdf(f"bkg_{c}", f"background, {c}", cat_index, pdfs)
                keep.append(mp)
                cat_index.setIndex(bnames.index(default_fn))
                wimport(w, mp, ROOT)
            else:
                cat_index.setIndex(bnames.index(default_fn))
                for i in range(pdfs.getSize()):
                    wimport(w, pdfs.at(i), ROOT)
                wimport(w, cat_index, ROOT)
            wimport(w, bnorm, ROOT)
            summary["categories"][c] = dict(
                binning=cb["binning"], rebin=cb["rebin"], n_data=float(h.Integral()),
                signal=dict(pdf=sp, norm=f"{sp}_norm", masses=cs["masses"]),
                background=dict(multipdf=f"bkg_{c}" if use_multipdf else None,
                                norm=f"bkg_{c}_norm", index=f"pdfindex_{c}",
                                functions={f: f"bkg_{c}_{f}" for f in bnames},
                                index_to_function={i: f for i, f in enumerate(bnames)},
                                default=default_fn, default_index=bnames.index(default_fn)))
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        w.writeToFile(str(args.output), True)
    nwarn, nerr = count_messages(log["text"])
    for line in message_summary(log["text"]):
        print("  ", line)
    print(f"Written {args.output} (RooFit warnings/errors while building: {nwarn}/{nerr})")

    # ---- re-open and check
    print("\nChecks on the written workspace:")
    f = ROOT.TFile.Open(str(args.output))
    w2 = f.Get(args.workspace_name)
    mh2 = w2.var(MASS_VAR)
    worst_param = worst_chi2 = worst_yield_vs_fit = 0.0
    for c in cats:
        x2 = w2.var(f"mjj_{c}")
        for name in (f"{DATA_NAME}_{c}", f"{SIGNAL_PROCESS}_{c}", f"{SIGNAL_PROCESS}_{c}_norm",
                     f"bkg_{c}_norm", f"pdfindex_{c}") + ((f"bkg_{c}",) if use_multipdf else ()):
            if not w2.obj(name):
                sys.exit(f"ERROR: {name} missing from the workspace")
        # signal parameters and yield at the fitted masses
        models = sig["categories"][c]["models"]
        for k in ALL_PARAMS:
            name = f"{SIGNAL_PROCESS}_{c}_norm" if k == "norm" else f"{SIGNAL_PROCESS}_{c}_{k}"
            lo_, hi_ = models[k]["clamp"]
            for m, fitted in zip(models[k]["points"]["masses"], models[k]["points"]["values"]):
                if not mh2.getMin() <= m <= mh2.getMax():
                    continue  # outside the common mass range
                mh2.setVal(m)
                got = w2.function(name).getVal()
                want = interp.evaluate(models[k], m)
                want = max(lo_, want) if lo_ is not None else want
                want = min(hi_, want) if hi_ is not None else want
                worst_param = max(worst_param, abs(got - want) / max(abs(want), 1e-12))
                if k == "norm":
                    worst_yield_vs_fit = max(worst_yield_vs_fit, abs(got / fitted - 1))
        # background: chi2 of each stored pdf against data_obs
        d = w2.data(f"{DATA_NAME}_{c}")
        n = np.array([d.weight(i) if d.get(i) else 0 for i in range(d.numEntries())])
        edges = bkg["categories"][c]["bin_edges"]
        lo, hi = np.array(edges[:-1]), np.array(edges[1:])
        for func in bkg["categories"][c]["passing"]:
            pdf = w2.pdf(f"bkg_{c}_{func}")
            fr = bkg["categories"][c]["fits"][func]
            # bin-centre evaluation with the bin-sum normalisation, as in the fits and in Combine
            nu = B.expected_counts(ROOT, x2, pdf, fr["norm"]["value"], lo, hi, False)
            chi2 = B.baker_cousins(n, nu)
            worst_chi2 = max(worst_chi2, abs(chi2 - fr["chi2"]))
    print(f"  max relative difference, workspace vs parametrisation (signal): {worst_param:.2e}")
    print(f"  max relative difference, parametrised vs fitted signal yield: {worst_yield_vs_fit:.2e}")
    print(f"  max |chi2(workspace) - chi2(fit)| of the background functions: {worst_chi2:.3g}")
    summary["checks"] = dict(max_rel_diff_signal_workspace_vs_parametrisation=worst_param,
                             max_rel_diff_signal_yield_parametrised_vs_fitted=worst_yield_vs_fit,
                             max_abs_diff_bkg_chi2=worst_chi2)
    write_json(Path(args.output).with_suffix(".json"), summary)
    print(f"Written {Path(args.output).with_suffix('.json')} (content summary)")
    if worst_param > 1e-6 or worst_chi2 > 0.1:
        sys.exit("ERROR: the workspace does not reproduce the fit results")


if __name__ == "__main__":
    main()
