"""Signal shape fits: double-sided Crystal Ball (RooCrystalBall with a common sigma)."""
import numpy as np

from .background import (BIN_INTEGRATION_PRECISION, baker_cousins, bin_sum_nbins,
                         expected_counts, hist_arrays, needs_bin_integration,
                         set_bin_sum_normalisation)
from .common import capture_output, count_messages, message_summary

CB_PARAMS = ["mean", "sigma", "alphaL", "nL", "alphaR", "nR"]
ALL_PARAMS = CB_PARAMS + ["norm"]

# allowed ranges (the starting values are derived from the histogram)
RANGES = dict(alpha=(0.1, 8.0), n=(0.5, 100.0))
START = dict(alphaL=0.5, alphaR=2.0, n=10.0)  # starting values of the tail parameters
N_CURVE = 600  # points of the drawn curve


def start_values(hist):
    """Peak position and a Gaussian width estimate (from the full width at half maximum)."""
    lo, hi, n = hist_arrays(hist)
    mc = 0.5 * (lo + hi)
    dens = n / (hi - lo)
    # smooth over 5 bins to be robust against fluctuations
    k = np.ones(5) / 5
    sm = np.convolve(dens, k, mode="same")
    i = int(np.argmax(sm))
    half = 0.5 * sm[i]
    left = i
    while left > 0 and sm[left] > half:
        left -= 1
    right = i
    while right < len(sm) - 1 and sm[right] > half:
        right += 1
    fwhm = max(mc[right] - mc[left], 2 * (hi[i] - lo[i]))
    return float(mc[i]), float(fwhm / 2.355)


def common_weight(hist):
    """The per-event weight w if all entries of `hist` have the same weight (errors^2 = w * content
    in every bin), else None."""
    ws = [hist.GetBinError(i) ** 2 / hist.GetBinContent(i)
          for i in range(1, hist.GetNbinsX() + 1) if hist.GetBinContent(i) > 0]
    if not ws:
        return None
    w = sum(ws) / len(ws)
    return w if all(abs(x - w) < 1e-6 * w for x in ws) else None


WEIGHT_MODES = ("sumw2", "counts", "asymptotic", "none")


def fit_signal(ROOT, hist, mass, tag, strategy=2, weights="sumw2"):
    """Extended binned fit of the signal template `hist` with a double-sided Crystal Ball.
    The fit range is the histogram range (= the analysis fit range); the fitted norm is the
    signal yield inside that range.

    The templates are MC counts scaled by a common weight w = lumi * (1 nb) / N_processed.
    weights selects how this is handled:
      sumw2      : fit the weighted template, errors corrected with SumW2Error (default);
                   the weights are first rescaled to a mean of 1 (the yield is scaled back):
                   neither the SumW2 errors nor the parameters depend on this overall scale,
                   but Minuit's convergence criterion does (a likelihood scaled by a small
                   weight converges much less precisely);
      counts     : fit the raw counts (hist / w) with Poisson errors and scale the yield back
                   by w (exact MC uncertainties, no covariance correction needed); falls back
                   to 'sumw2' if the weights are not uniform;
      asymptotic : fit the weighted template, errors from AsymptoticError (weights rescaled to a
                   mean of 1 as for sumw2);
      none       : fit the weighted template without error correction (errors not valid for
                   weighted data; for comparisons only)."""
    if weights not in WEIGHT_MODES:
        raise ValueError(f"unknown weights mode '{weights}'")
    w_common = common_weight(hist) if weights == "counts" else None
    if weights == "counts" and w_common is None:
        weights = "sumw2"
    scale_back = 1.0
    if w_common is not None:
        hist = hist.Clone(hist.GetName() + "_counts")
        hist.SetDirectory(0)
        hist.Scale(1.0 / w_common)
        for i in range(1, hist.GetNbinsX() + 1):
            hist.SetBinError(i, max(hist.GetBinContent(i), 0.0) ** 0.5)
        scale_back = w_common
    elif weights in ("sumw2", "asymptotic"):
        sumw = sum(hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1))
        sumw2 = sum(hist.GetBinError(i) ** 2 for i in range(1, hist.GetNbinsX() + 1))
        if sumw > 0 and sumw2 > 0:
            mean_w = sumw2 / sumw  # mean event weight (sum w^2 / sum w)
            hist = hist.Clone(hist.GetName() + "_unitw")
            hist.SetDirectory(0)
            hist.Scale(1.0 / mean_w)  # scales contents and errors: weights w / mean_w
            scale_back = mean_w
    lo, hi, n = hist_arrays(hist)
    xmin, xmax = float(lo[0]), float(hi[-1])
    edges = list(lo) + [xmax]
    x = ROOT.RooRealVar(f"mjj_{tag}", "m_{jj} [GeV]", xmin, xmax)
    x.setBinning(ROOT.RooBinning(len(edges) - 1, np.array(edges, dtype="d")))

    peak, sig0 = start_values(hist)
    sig0 *= 0.75  # the FWHM-based width overestimates the core width (tails)
    p = dict(
        mean=ROOT.RooRealVar(f"mean_{tag}", "mean", peak, peak - 3 * sig0, peak + 3 * sig0),
        sigma=ROOT.RooRealVar(f"sigma_{tag}", "sigma", sig0, 0.2 * sig0, 3.0 * sig0),
        alphaL=ROOT.RooRealVar(f"alphaL_{tag}", "alphaL", START["alphaL"], *RANGES["alpha"]),
        nL=ROOT.RooRealVar(f"nL_{tag}", "nL", START["n"], *RANGES["n"]),
        alphaR=ROOT.RooRealVar(f"alphaR_{tag}", "alphaR", START["alphaR"], *RANGES["alpha"]),
        nR=ROOT.RooRealVar(f"nR_{tag}", "nR", START["n"], *RANGES["n"]),
    )
    pdf = ROOT.RooCrystalBall(f"cb_{tag}", "", x, p["mean"], p["sigma"], p["alphaL"], p["nL"],
                              p["alphaR"], p["nR"])
    nbins = bin_sum_nbins(lo, hi)
    if nbins:  # bins evaluated at their centres: normalise by the bin sum, as the background
        set_bin_sum_normalisation(pdf, nbins)
    ntot = float(n.sum())
    p["norm"] = ROOT.RooRealVar(f"norm_{tag}", "N_{sig}", ntot, 0.5 * ntot, 1.5 * ntot)
    model = ROOT.RooExtendPdf(f"model_{tag}", "", pdf, p["norm"])
    data = ROOT.RooDataHist(f"data_{tag}", "", ROOT.RooArgList(x), hist)

    integrate = needs_bin_integration(lo, hi)  # rebinned category: integrate over the bins

    texts = []  # RooFit/Minuit output of every stage (messages are counted over all of them)

    def fit():
        with capture_output() as log:
            r = model.fitTo(data, ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1),
                            ROOT.RooFit.Strategy(strategy),
                            ROOT.RooFit.SumW2Error(weights == "sumw2"),
                            ROOT.RooFit.AsymptoticError(weights == "asymptotic"),
                            ROOT.RooFit.Extended(True),
                            ROOT.RooFit.IntegrateBins(BIN_INTEGRATION_PRECISION if integrate else -1))
        texts.append(log["text"])
        return r, log["text"]

    # The parameters are released in stages, from the core to the tails:
    #   1: mean, sigma, norm   2: + alphaL, alphaR   3: + nL, nR
    # Then tail parameters that the template does not determine are fixed and the fit is
    # repeated (up to 3 times):
    #   - a tail whose onset (mean -/+ alpha*sigma) lies outside the fit range has no effect
    #     inside it: its alpha and n are fixed;
    #   - an alpha or n that runs to its limit is fixed there (for n -> infinity the power-law
    #     tail becomes exponential, exp(-alpha|t|), and n is not determined).
    stages = []
    tails = ("alphaL", "nL", "alphaR", "nR")
    for t in tails:
        p[t].setConstant(True)
    res, text = fit()
    stages.append(dict(stage=1, status=int(res.status()), counts=count_messages(text)))
    for t in ("alphaL", "alphaR"):
        p[t].setConstant(False)
    res, text = fit()
    stages.append(dict(stage=2, status=int(res.status()), counts=count_messages(text)))
    for t in ("nL", "nR"):
        p[t].setConstant(False)
    res, text = fit()
    stages.append(dict(stage=3, status=int(res.status()), counts=count_messages(text)))

    fixed = {}
    for it in range(3):
        new = {}
        mu_, sg_ = p["mean"].getVal(), p["sigma"].getVal()
        for side, sign, edge in (("L", -1, xmin), ("R", +1, xmax)):
            a, nn = p["alpha" + side], p["n" + side]
            onset = mu_ + sign * a.getVal() * sg_
            if (sign < 0 and onset <= edge) or (sign > 0 and onset >= edge):
                for v in (a, nn):
                    if not v.isConstant():
                        new[v.GetName().split("_")[0]] = ("tail outside the fit range", v)
        for t in tails:
            v = p[t]
            if v.isConstant() or t in new:
                continue
            span = v.getMax() - v.getMin()
            if v.getMax() - v.getVal() < 0.01 * span:
                v.setVal(v.getMax())
                new[t] = ("at upper limit", v)
            elif v.getVal() - v.getMin() < 0.01 * span:
                v.setVal(v.getMin())
                new[t] = ("at lower limit", v)
        if not new:
            break
        for t, (why, v) in new.items():
            v.setConstant(True)
            fixed[t] = dict(value=v.getVal(), reason=why)
        res, text = fit()
        stages.append(dict(stage=4 + it, status=int(res.status()), counts=count_messages(text),
                           fixed=sorted(new)))
    all_text = "\n".join(texts)
    nwarn, nerr = count_messages(all_text)

    nu = expected_counts(ROOT, x, pdf, p["norm"].getVal(), lo, hi, integrate)
    if w_common is not None:  # raw MC counts: Poisson likelihood-ratio chi2
        chi2 = baker_cousins(n, nu)
        chi2_kind = "Poisson likelihood ratio (MC counts)"
    elif common_weight(hist) is not None:  # weighted fit of a uniformly weighted template:
        wu = common_weight(hist)            # same statistic, on the equivalent MC counts
        chi2 = baker_cousins(n / wu, nu / wu)
        chi2_kind = "Poisson likelihood ratio (equivalent MC counts)"
    else:
        err = np.array([hist.GetBinError(i) for i in range(1, hist.GetNbinsX() + 1)])
        ok = err > 0
        chi2 = float(np.sum(((n[ok] - nu[ok]) / err[ok]) ** 2))
        chi2_kind = "Neyman (bin errors)"
    nfloat = res.floatParsFinal().getSize()
    ndf = int(len(n) - nfloat)
    pvalue = float(ROOT.TMath.Prob(chi2, ndf)) if ndf > 0 else 0.0

    at_limit = []
    for name, v in p.items():
        if v.isConstant():
            continue
        span = v.getMax() - v.getMin()
        if v.getVal() - v.getMin() < 1e-3 * span or v.getMax() - v.getVal() < 1e-3 * span:
            at_limit.append(name)

    wscale = scale_back
    nu = nu * wscale

    nset = ROOT.RooArgSet(x)
    xs = np.linspace(xmin, xmax, N_CURVE)
    curve = []
    for m in xs:
        x.setVal(m)
        curve.append(wscale * p["norm"].getVal() * pdf.getVal(nset))

    return dict(
        mass=mass, status=int(res.status()), cov_qual=int(res.covQual()), edm=float(res.edm()),
        chi2=chi2, chi2_kind=chi2_kind, ndf=ndf, pvalue=pvalue, n_warnings=nwarn, n_errors=nerr,
        messages=message_summary(all_text), params_at_limit=at_limit,
        # a parameter fixed by the procedure has no uncertainty (RooFit would keep the stale
        # error of the fit before it was fixed)
        params={k: dict(value=v.getVal() * (wscale if k == "norm" else 1.0),
                        error=0.0 if v.isConstant() else v.getError() * (wscale if k == "norm" else 1.0),
                        constant=bool(v.isConstant()))
                for k, v in p.items()},
        weights_mode=weights,
        mc_counts_in_range=float(n.sum()), event_weight=wscale, bin_integration=integrate,
        start=dict(mean=peak, sigma=sig0), fixed=fixed, stages=stages,
        expected=nu.tolist(), curve=(xs.tolist(), curve),
        _keep=(x, p, pdf, model, data, res),
    )


def fit_ok(fit):
    return fit["status"] == 0 and fit["cov_qual"] == 3 and not fit["params_at_limit"]

