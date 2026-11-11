#!/usr/bin/env python3
"""Generate synthetic input files with a known truth, for the closure tests of the workflow.

The files have the same format as the real inputs, so the workflow (run_all.sh) can run on
them unchanged:
  data.root                         one m_jj histogram per category (Poisson-fluctuated background)
  ZprimeQQ_M<mass>_synthetic.root   per mass: the same histograms filled with signal MC counts
                                    (multinomially sampled), plus 'eventsProcessed'
  truth.json                        the generating functions and parameters, and the numbers
                                    of generated events, to compare the fit results with
Some settings deliberately differ from the real analysis (category names, histogram range,
luminosity), so that the tests also check that nothing depends on the real setup:
  * 3 categories TestJet_dY<A>to<B>_mjj, an inclusive histogram (to be excluded) and an
    unrelated histogram (to be ignored);
  * histogram range 100-800 GeV in 1 GeV bins (wider than any fit range used in the tests);
  * background: a 4-parameter dijet function per category; optionally multiplied by a
    low-mass turn-on, efficiency 1 / (1 + exp(-(m - mid) / width)), to test the fit-range scan;
  * signal: double-sided Crystal Ball shapes whose parameters depend linearly on the mass.
    Masses VALIDATION_MASSES are written to a separate directory: they are not given to the
    workflow, and serve to check the interpolation of the signal model.
The fit range is not a property of the inputs: it is passed to the workflow (workflow_config),
and the tests take it from their command line (pytest --fit-xmin/--fit-xmax, default 150-600).

Usage:  python tests/synthetic.py -o OUTDIR [--seed N] [--xmin 150 --xmax 600]
                                  [--turnon-mid M --turnon-width W]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

SQRTS = 13600.0
HIST_RANGE = (100.0, 800.0)
BIN_WIDTH = 1.0
LUMI = 1000.0            # nb^-1 given to prepare_inputs.py
N_GENERATED = 1_000_000  # signal events processed per mass point
FIT_MASSES = (200, 300, 400, 500)
VALIDATION_MASSES = (350,)
CATEGORY_PREFIX = "TestJet_dY"
CATEGORY_SUFFIX = "_mjj"
INCLUSIVE = "0p0to1p5"

# category -> background truth: dijet4 parameters (p1, p2, p3), expected events in the histogram
# range, and the fraction of the selected signal events in this category
CATEGORIES = {
    "0p0to0p5": dict(p=(10.0, 13.0, 1.00), n_bkg=5.0e6, sig_frac=0.45),
    "0p5to1p0": dict(p=(10.0, 12.8, 0.97), n_bkg=4.0e6, sig_frac=0.35),
    "1p0to1p5": dict(p=(12.0, 12.5, 0.95), n_bkg=3.0e6, sig_frac=0.20),
}


def signal_params(mass):
    """Crystal Ball parameters (RooCrystalBall convention, common sigma) at a given mass."""
    return dict(mean=0.99 * mass + 1.0, sigma=15.0 + 0.04 * mass, alphaL=0.5 + 0.0002 * mass,
                nL=6.0 + 0.004 * mass, alphaR=1.6, nR=4.0)


def signal_efficiency(mass):
    """Fraction of the processed signal events that end up in the histograms (all categories)."""
    return 0.15 + 0.0001 * mass


# ---------------------------------------------------------------------------------------------
def dijet4(m, p1, p2, p3):
    x = np.asarray(m, float) / SQRTS
    return np.exp(p1 * np.log1p(-x) - (p2 + p3 * np.log(x)) * np.log(x))


def turnon_efficiency(m, turnon):
    """turnon = (mid, width), or None for no turn-on."""
    if turnon is None:
        return np.ones_like(np.asarray(m, float))
    mid, width = turnon
    return 1.0 / (1.0 + np.exp(-(np.asarray(m, float) - mid) / width))


def crystal_ball(x, mean, sigma, alphaL, nL, alphaR, nR):
    """Unnormalised double-sided Crystal Ball, same definition as RooCrystalBall."""
    t = (np.asarray(x, float) - mean) / sigma
    out = np.exp(-0.5 * t * t)
    left, right = t < -alphaL, t > alphaR
    bl = nL / alphaL - alphaL
    out[left] = np.exp(nL * np.log(nL / alphaL) - 0.5 * alphaL ** 2 - nL * np.log(bl - t[left]))
    br = nR / alphaR - alphaR
    out[right] = np.exp(nR * np.log(nR / alphaR) - 0.5 * alphaR ** 2 - nR * np.log(br + t[right]))
    return out


def bin_integrals(func, edges, nsub=20):
    """Integral of func over each bin (composite Simpson rule, nsub even)."""
    edges = np.asarray(edges, float)
    lo, hi = edges[:-1], edges[1:]
    k = np.arange(nsub + 1)
    w = np.where((k == 0) | (k == nsub), 1.0, np.where(k % 2 == 1, 4.0, 2.0))
    pts = lo[:, None] + (hi - lo)[:, None] * k[None, :] / nsub
    return (func(pts) * w).sum(axis=1) * (hi - lo) / (3 * nsub)


def edges():
    n = int(round((HIST_RANGE[1] - HIST_RANGE[0]) / BIN_WIDTH))
    return np.linspace(HIST_RANGE[0], HIST_RANGE[1], n + 1)


def background_expected(cat, turnon=None):
    """Expected background counts per bin (histogram range) for a category."""
    p = CATEGORIES[cat]["p"]
    shape = bin_integrals(lambda m: dijet4(m, *p) * turnon_efficiency(m, turnon), edges())
    return CATEGORIES[cat]["n_bkg"] * shape / shape.sum()


def signal_probabilities(mass):
    """Probability per bin of a selected signal event (histogram range)."""
    sp = signal_params(mass)
    prob = bin_integrals(lambda m: crystal_ball(m, **sp), edges())
    return prob / prob.sum()


def category_names():
    return {c: f"{CATEGORY_PREFIX}{c}{CATEGORY_SUFFIX}" for c in CATEGORIES}


def make_hist(ROOT, name, counts):
    from array import array
    e = edges()
    h = ROOT.TH1D(name, name, len(e) - 1, array("d", e))
    h.SetDirectory(0)
    h.Sumw2()
    for i, v in enumerate(counts, start=1):
        h.SetBinContent(i, float(v))
        h.SetBinError(i, float(np.sqrt(v)))
    return h


def generate(outdir, seed=12345, turnon=None):
    """Write the synthetic inputs to outdir and return the truth (also written to truth.json).
    turnon: (mid, width) of a low-mass turn-on applied to the background, or None."""
    import ROOT
    ROOT.gROOT.SetBatch(True)
    rng = np.random.default_rng(seed)
    outdir = Path(outdir)
    (outdir / "validation").mkdir(parents=True, exist_ok=True)
    names = category_names()
    truth = dict(seed=seed, sqrts=SQRTS, hist_range=list(HIST_RANGE), bin_width=BIN_WIDTH,
                 lumi=LUMI, n_generated=N_GENERATED, turnon=list(turnon) if turnon else None,
                 categories={}, signal={}, files=dict(signal=[], validation=[]))

    # ---- data
    f = ROOT.TFile(str(outdir / "data.root"), "RECREATE")
    inclusive = np.zeros(len(edges()) - 1)
    for cat, name in names.items():
        n = rng.poisson(background_expected(cat, turnon))
        inclusive += n
        make_hist(ROOT, name, n).Write()
        truth["categories"][cat] = dict(histogram=name, dijet4=list(CATEGORIES[cat]["p"]),
                                        n_bkg_expected=CATEGORIES[cat]["n_bkg"], n_data=int(n.sum()),
                                        sig_frac=CATEGORIES[cat]["sig_frac"])
    make_hist(ROOT, f"{CATEGORY_PREFIX}{INCLUSIVE}{CATEGORY_SUFFIX}", inclusive).Write()  # excluded
    make_hist(ROOT, f"{CATEGORY_PREFIX}0p0to0p5_pt", inclusive).Write()                  # unrelated
    f.Close()
    truth["files"]["data"] = str(outdir / "data.root")
    truth["inclusive_histogram"] = f"{CATEGORY_PREFIX}{INCLUSIVE}{CATEGORY_SUFFIX}"

    # ---- signal
    for mass in FIT_MASSES + VALIDATION_MASSES:
        sub = outdir / ("validation" if mass in VALIDATION_MASSES else ".")
        path = sub / f"ZprimeQQ_M{mass}_synthetic.root"
        prob = signal_probabilities(mass)
        n_sel_total = rng.binomial(N_GENERATED, signal_efficiency(mass))
        n_cat = rng.multinomial(n_sel_total, [CATEGORIES[c]["sig_frac"] for c in CATEGORIES])
        f = ROOT.TFile(str(path), "RECREATE")
        sig = dict(params=signal_params(mass), efficiency=signal_efficiency(mass), categories={})
        for (cat, name), ncat in zip(names.items(), n_cat):
            make_hist(ROOT, name, rng.multinomial(int(ncat), prob)).Write()
            sig["categories"][cat] = dict(n_selected=int(ncat))
        make_hist(ROOT, f"{CATEGORY_PREFIX}{INCLUSIVE}{CATEGORY_SUFFIX}", np.zeros(len(prob))).Write()
        ev = ROOT.TH1D("eventsProcessed", "", 1, 0, 1)
        ev.SetBinContent(1, float(N_GENERATED))
        ev.Write()
        f.Close()
        truth["signal"][str(mass)] = sig
        truth["files"]["validation" if mass in VALIDATION_MASSES else "signal"].append(str(path))
    (outdir / "truth.json").write_text(json.dumps(truth, indent=1))
    return truth


def workflow_config(xmin, xmax, scan_max=None):
    """Settings for run_all.sh matching the synthetic inputs (environment variables). The scan
    may raise the lower edge up to scan_max (default: xmin + 60 GeV)."""
    return dict(
        LUMI=str(LUMI), XMIN=f"{xmin:g}", XMAX=f"{xmax:g}",
        SCAN_MAX=f"{scan_max if scan_max is not None else xmin + 60.0:g}",
        CATEGORIES=f"({CATEGORY_PREFIX}\\d+p\\d+to\\d+p\\d+){CATEGORY_SUFFIX}",
        EXCLUDE=f"{CATEGORY_PREFIX}{INCLUSIVE}{CATEGORY_SUFFIX}", FORMATS="png")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--xmin", type=float, default=150.0, help="fit range for the printed command")
    ap.add_argument("--xmax", type=float, default=600.0, help="fit range for the printed command")
    ap.add_argument("--turnon-mid", type=float, help="add a background turn-on centred here [GeV]")
    ap.add_argument("--turnon-width", type=float, default=4.0)
    args = ap.parse_args()
    turnon = (args.turnon_mid, args.turnon_width) if args.turnon_mid is not None else None
    truth = generate(args.outdir, args.seed, turnon)
    print(f"Written synthetic inputs to {args.outdir} (truth.json); run the workflow with:")
    env = " ".join(f"{k}='{v}'" for k, v in workflow_config(args.xmin, args.xmax).items())
    print(f"  {env} ./run_all.sh {truth['files']['data']} {' '.join(truth['files']['signal'])}")


if __name__ == "__main__":
    sys.exit(main())
