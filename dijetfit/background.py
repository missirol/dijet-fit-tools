"""Background functions for the dijet m_jj spectrum and their fits.

All functions are of the form  f(m) = exp( sum_k p_k * T_k(m) ),  i.e. log-linear in their
parameters, with the normalisation given by a separate yield parameter. The families are the
standard ones of dijet resonance searches, written with x = m / sqrt(s):

  dijet    (CMS)          f = (1-x)^p1 / x^(p2 + p3 ln x + p4 ln^2 x)
  moddijet (modified)     f = (1-x^(1/3))^p1 / x^(p2 + p3 ln x + p4 ln^2 x)
  ua2      (UA2/CDF)      f = x^(-p2) exp(-p1 x - p3 x^2 - p4 x^3)

<family><N> denotes N parameters *including* the normalisation (the usual convention in CMS
papers), e.g. dijet4 = the 4-parameter dijet function (p1, p2, p3 + yield). Orders 3-6 are
available (the 6-parameter forms add p5 ln^3 x, resp. p5 x^4).

For numerical stability the polynomial in ln x (resp. in x) is written in the centred
variable L = ln(m / M0) (resp. z = m / M0 - 1), with the reference mass M0 = the geometric
mean of the fit range (so nothing depends on a particular mass range). A polynomial
of the same degree in a shifted variable spans exactly the same set of functions (the
constant term is absorbed in the normalisation), so the families are identical to the
published ones; only the meaning of the individual parameters differs.
"""
import math

import numpy as np

from .common import SQRTS, capture_output, count_messages, message_summary

INTEGRATOR_EPS_REL = 1e-10
BIN_INTEGRATION_ABOVE = 1.0 + 1e-6  # GeV: bins wider than this are integrated over
BIN_INTEGRATION_PRECISION = 1e-10
QRANGE = 50.0  # half-range of the decorrelated parameters (their errors are ~1)
MAX_RECENTRE = 3  # re-centring attempts for a parameter that ends at its limit
NORM_RANGE = (0.5, 1.5)  # allowed range of the yield, in units of the events in the fit range

_X = f"(@0/{SQRTS!r})"


def term_defs(m0):
    """term name -> (numpy function of m, RooFit expression of @0), for reference mass m0."""
    r = float(m0)
    return {
        "ln1mx": (lambda m: np.log1p(-m / SQRTS), f"log(1-{_X})"),
        "ln1mx13": (lambda m: np.log(1 - np.cbrt(m / SQRTS)), f"log(1-pow({_X},1./3.))"),
        "L1": (lambda m: np.log(m / r), f"log(@0/{r!r})"),
        "L2": (lambda m: np.log(m / r) ** 2, f"pow(log(@0/{r!r}),2)"),
        "L3": (lambda m: np.log(m / r) ** 3, f"pow(log(@0/{r!r}),3)"),
        "L4": (lambda m: np.log(m / r) ** 4, f"pow(log(@0/{r!r}),4)"),
        "z1": (lambda m: m / r - 1, f"(@0/{r!r}-1)"),
        "z2": (lambda m: (m / r - 1) ** 2, f"pow(@0/{r!r}-1,2)"),
        "z3": (lambda m: (m / r - 1) ** 3, f"pow(@0/{r!r}-1,3)"),
        "z4": (lambda m: (m / r - 1) ** 4, f"pow(@0/{r!r}-1,4)"),
    }


def reference_mass(xmin, xmax):
    """M0: geometric mean of the fit range."""
    return float(np.sqrt(xmin * xmax))


FAMILIES = {
    "dijet": ["ln1mx", "L1", "L2", "L3", "L4"],
    "moddijet": ["ln1mx13", "L1", "L2", "L3", "L4"],
    "ua2": ["L1", "z1", "z2", "z3", "z4"],
}
FAMILY_FORMULAS = {
    "dijet": "(1-x)^p1 / x^(p2 + p3 ln x + p4 ln^2 x [+ p5 ln^3 x])",
    "moddijet": "(1-x^(1/3))^p1 / x^(p2 + p3 ln x + p4 ln^2 x [+ p5 ln^3 x])",
    "ua2": "x^(-p2) exp(-p1 x - p3 x^2 - p4 x^3 [- p5 x^4])",
}
ORDERS = (3, 4, 5, 6)  # number of parameters including the normalisation


def function_names(families=None, orders=ORDERS):
    return [f"{fam}{n}" for fam in (families or FAMILIES) for n in orders]


def parse_function(name):
    for fam in FAMILIES:
        if name.startswith(fam) and name[len(fam):].isdigit():
            n = int(name[len(fam):])
            if n in ORDERS:
                return fam, n
    raise ValueError(f"unknown background function '{name}'")


def terms_of(name):
    fam, n = parse_function(name)
    return FAMILIES[fam][:n - 1]


def formula_of(name, m0, M=None, tbar=None):
    """RooGenericPdf expression with @0 = m and @1.. = the shape parameters.

    Without a basis: exp(sum_k p_k T_k). With a basis (M, tbar), see decorrelated_basis:
    exp(sum_j q_j U_j) with U_j = sum_k M_kj (T_k - tbar_k)."""
    terms = terms_of(name)
    T = term_defs(m0)
    if M is None:
        return "exp(" + " + ".join(f"@{k + 1}*{T[t][1]}" for k, t in enumerate(terms)) + ")"
    centred = [f"({T[t][1]}-({float(tbar[k])!r}))" for k, t in enumerate(terms)]
    parts = []
    for j in range(len(terms)):
        u = "+".join(f"({float(M[k, j])!r})*{centred[k]}" for k in range(len(terms)))
        parts.append(f"@{j + 1}*({u})")
    return "exp(" + " + ".join(parts) + ")"


def describe(name):
    fam, n = parse_function(name)
    terms = terms_of(name)
    return (f"{FAMILY_FORMULAS[fam]} with {n - 1} shape parameters; implemented as "
            f"exp(sum_k a_k T_k) with T = {terms}, L = ln(m/M0), z = m/M0 - 1, "
            f"x = m/{SQRTS:g}, M0 = geometric mean of the fit range")


def hist_arrays(hist):
    nb = hist.GetNbinsX()
    ax = hist.GetXaxis()
    lo = np.array([ax.GetBinLowEdge(i) for i in range(1, nb + 1)])
    hi = np.array([ax.GetBinUpEdge(i) for i in range(1, nb + 1)])
    n = np.array([hist.GetBinContent(i) for i in range(1, nb + 1)])
    return lo, hi, n


def start_values(name, hist, m0, max_iter=100):
    """Maximum-likelihood estimate of the shape parameters, used as starting values.

    The binned model  nu_i = w_i * exp(c0 + sum_k p_k T_k(m_i))  is a Poisson generalised
    linear model with log link, so its maximum-likelihood solution is found reliably by
    iteratively reweighted least squares. It is very close to the minimum of the RooFit extended
    binned likelihood (where c0 is traded for the yield and the pdf integral; the two differ
    only by the bin-centre vs integral normalisation), so Minuit starts next to the minimum.
    Returns (p_1..p_k, expected counts per bin at the solution).
    """
    lo, hi, n = hist_arrays(hist)
    mc, w = 0.5 * (lo + hi), hi - lo
    T = term_defs(m0)
    X = np.column_stack([np.ones(len(n))] + [T[t][0](mc) for t in terms_of(name)])
    off = np.log(w)
    ok = n > 0
    # initial guess: least squares on ln(n/w)
    sw = np.sqrt(n[ok])
    c, *_ = np.linalg.lstsq(X[ok] * sw[:, None], np.log(n[ok] / w[ok]) * sw, rcond=None)

    def deviance(c):
        nu = np.exp(X @ c + off)
        return baker_cousins(n, nu)

    dev = deviance(c)
    for _ in range(max_iter):
        mu = np.exp(X @ c + off)
        z = X @ c + (n - mu) / mu
        A = X * np.sqrt(mu)[:, None]
        c_new, *_ = np.linalg.lstsq(A, z * np.sqrt(mu), rcond=None)
        step = c_new - c
        for _half in range(30):  # step halving if the deviance does not decrease
            dev_new = deviance(c + step)
            if dev_new <= dev + 1e-9 * abs(dev):
                break
            step *= 0.5
        c = c + step
        converged = abs(dev - dev_new) < 1e-8 * max(1.0, dev)
        dev = dev_new
        if converged:
            break
    mu = np.exp(X @ c + off)
    return c[1:], mu


def decorrelated_basis(name, hist, mu, m0):
    """Linear combinations U_j = sum_k M_kj (T_k(m) - tbar_k) of the terms of `name` that are
    uncorrelated, with unit variance, at the expected bin contents `mu` (the starting point).
    In this basis Minuit sees a nearly diagonal, well-scaled problem; the family of functions
    is unchanged (the constants tbar_k are absorbed in the normalisation).
    Returns (M, tbar)."""
    lo, hi, _ = hist_arrays(hist)
    mc = 0.5 * (lo + hi)
    defs = term_defs(m0)
    T = np.column_stack([defs[t][0](mc) for t in terms_of(name)])
    tbar = (mu[:, None] * T).sum(axis=0) / mu.sum()
    Tc = T - tbar
    F = (Tc * mu[:, None]).T @ Tc
    L = np.linalg.cholesky(F)
    M = np.linalg.inv(L).T  # U = Tc @ M  ->  U^T diag(mu) U = 1
    return M, tbar


def baker_cousins(n, nu):
    """Poisson likelihood-ratio chi2, summed over bins."""
    n = np.asarray(n, float)
    nu = np.asarray(nu, float)
    if np.any(nu <= 0):
        return float("inf")
    term = nu - n
    pos = n > 0
    term[pos] += n[pos] * np.log(n[pos] / nu[pos])
    return float(2.0 * term.sum())


def make_pdf(ROOT, pdf_name, func, arglist, formula, nbins=None):
    """RooGenericPdf for background function `func` (arglist = observable, p1, p2, ...).

    nbins given (uniform bins evaluated at their centres, i.e. not rebinned): the pdf is
    normalised by the sum of its values at the nbins bin centres times the bin width
    (RooBinIntegrator), not by its integral. The expected counts per bin, N * pdf(centre) *
    width, then add up exactly to the yield N, i.e. the extended likelihood of RooFit and Combine
    is the binned Poisson likelihood of these expectations. With the integral as normalisation
    the per-bin expectations add up to N * (1 + O(1e-5)) (midpoint-rule error), and the
    parameter dependence of that mismatch is enough, with O(10^7) events, to pull the shape
    parameters by O(1) standard deviations and the signal strength in Combine by up to ~1
    standard deviation (also between Combine's own toys and its fits).

    nbins None (bins integrated over, see fit_function): the integral is used, with an adaptive
    Gauss-Kronrod integrator of relative precision 1e-10 (with the default 1e-7 the likelihood
    is noisy at the O(1) level for these sample sizes).
    The integrator setting is stored with the pdf in the workspace."""
    pdf = ROOT.RooGenericPdf(pdf_name, func, formula, arglist)
    if nbins:
        set_bin_sum_normalisation(pdf, nbins)
    else:
        cfg = pdf.specialIntegratorConfig(True)
        cfg.setEpsRel(INTEGRATOR_EPS_REL)
        cfg.setEpsAbs(0.0)
        cfg.method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D")
    return pdf


def set_bin_sum_normalisation(pdf, nbins):
    """Normalise `pdf` by the sum of its values at the centres of nbins uniform bins over the
    range of its observable, times the bin width (RooBinIntegrator; the setting is stored with
    the pdf). Analytical integrals (e.g. of RooCrystalBall) are switched off for this."""
    pdf.forceNumInt(True)
    cfg = pdf.specialIntegratorConfig(True)
    cfg.method1D().setLabel("RooBinIntegrator")
    cfg.getConfigSection("RooBinIntegrator").setRealValue("numBins", int(nbins))
    return pdf


def bin_sum_nbins(lo, hi):
    """Number of bins for the bin-sum normalisation of make_pdf, or None if the bins are
    integrated over (wider than 1 GeV) or not uniform."""
    w = hi - lo
    if needs_bin_integration(lo, hi) or np.max(w) - np.min(w) > 1e-6 * np.max(w):
        return None
    return len(lo)


def at_limit_of(var, frac=1e-3):
    """True if a RooRealVar is within frac x (its range) of either end of its range."""
    span = var.getMax() - var.getMin()
    return var.getVal() - var.getMin() < frac * span or var.getMax() - var.getVal() < frac * span


def needs_bin_integration(lo, hi):
    return bool(np.max(hi - lo) > BIN_INTEGRATION_ABOVE)


def expected_counts(ROOT, x, pdf, norm, lo, hi, integrate, npts=16):
    """Expected events per bin: norm * pdf (normalised on the range of x) at the bin centre times
    the bin width, or, if `integrate`, integrated over the bin (composite Simpson rule)."""
    nset = ROOT.RooArgSet(x)

    def dens(m):
        x.setVal(float(m))
        return pdf.getVal(nset)

    out = []
    for a, b in zip(lo, hi):
        if not integrate:
            out.append(dens(0.5 * (a + b)) * (b - a))
            continue
        ts = np.linspace(a, b, npts + 1)
        fs = np.array([dens(t) for t in ts])
        h = (b - a) / npts
        out.append(h / 3 * (fs[0] + fs[-1] + 4 * fs[1:-1:2].sum() + 2 * fs[2:-1:2].sum()))
    return norm * np.array(out)


def fit_function(ROOT, name, hist, tag, strategy=1):
    """Extended binned maximum-likelihood fit of function `name` to `hist` (whose axis is the
    fit range). Returns a dict with the results (python values) and the RooFit objects.

    With 1 GeV bins the pdf is evaluated at the bin centres (RooFit's default for binned data)
    and normalised by its sum over the bin centres (see make_pdf), so that the expected counts
    add up to the fitted yield. With wider bins (rebinned categories) the pdf is integrated over
    each bin: RooFit's default binned likelihood (bin-centre density, but total yield from the
    exact integral) is then biased at the level of many standard deviations for these sample
    sizes."""
    lo, hi, n = hist_arrays(hist)
    xmin, xmax = float(lo[0]), float(hi[-1])
    edges = list(lo) + [xmax]

    x = ROOT.RooRealVar(f"mjj_{tag}", "m_{jj} [GeV]", xmin, xmax)
    x.setBinning(ROOT.RooBinning(len(edges) - 1, np.array(edges, dtype="d")))
    m0 = reference_mass(xmin, xmax)
    pvals, mu = start_values(name, hist, m0)
    M, tbar = decorrelated_basis(name, hist, mu, m0)
    qvals = np.linalg.solve(M, pvals)  # p = M q
    formula = formula_of(name, m0, M, tbar)
    params = []
    arglist = ROOT.RooArgList(x)
    for k, v in enumerate(qvals, start=1):
        # in the decorrelated basis the parameter uncertainties are ~1
        p = ROOT.RooRealVar(f"{name}_p{k}_{tag}", f"q{k}", v, v - QRANGE, v + QRANGE)
        params.append(p)
        arglist.add(p)
    pdf = make_pdf(ROOT, f"{name}_{tag}", name, arglist, formula, bin_sum_nbins(lo, hi))
    ntot = float(n.sum())
    norm = ROOT.RooRealVar(f"{name}_{tag}_norm", "N_{bkg}", ntot, NORM_RANGE[0] * ntot, NORM_RANGE[1] * ntot)
    model = ROOT.RooExtendPdf(f"{name}_{tag}_ext", "", pdf, norm)
    data = ROOT.RooDataHist(f"data_{tag}", "", ROOT.RooArgList(x), hist)

    integrate = needs_bin_integration(lo, hi)
    opts = [ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.Strategy(strategy),
            ROOT.RooFit.Offset(True), ROOT.RooFit.Extended(True),
            ROOT.RooFit.IntegrateBins(BIN_INTEGRATION_PRECISION if integrate else -1)]
    # A shape parameter that ends at its limit: re-centre its range and fit again (at most
    # MAX_RECENTRE times). The at-limit check below is made on the result of the last fit, so a
    # parameter that still ends at a limit after the last attempt is reported as such.
    texts = []  # RooFit/Minuit output of every attempt (messages are counted over all of them)
    with capture_output() as log:
        res = model.fitTo(data, *opts)
    texts.append(log["text"])
    for _ in range(MAX_RECENTRE):
        hit = [p for p in params if at_limit_of(p)]
        if not hit:
            break
        for p in hit:
            p.setRange(p.getVal() - QRANGE, p.getVal() + QRANGE)
        with capture_output() as log:
            res = model.fitTo(data, *opts)
        texts.append(log["text"])
    all_text = "\n".join(texts)
    nwarn, nerr = count_messages(all_text)

    # expected events per bin and the goodness of fit
    nu = expected_counts(ROOT, x, pdf, norm.getVal(), lo, hi, integrate)
    chi2 = baker_cousins(n, nu)
    nfloat = res.floatParsFinal().getSize()
    ndf = len(n) - nfloat
    pvalue = ROOT.TMath.Prob(chi2, ndf) if (ndf > 0 and math.isfinite(chi2)) else 0.0

    at_limit = [p.GetName() for p in params + [norm] if at_limit_of(p)]

    nset = ROOT.RooArgSet(x)
    xs = np.linspace(xmin, xmax, 400)
    curve = []
    for m in xs:
        x.setVal(m)
        curve.append(norm.getVal() * pdf.getVal(nset))

    return dict(
        name=name, family=parse_function(name)[0], npar=parse_function(name)[1],
        status=int(res.status()), cov_qual=int(res.covQual()), edm=float(res.edm()),
        min_nll=float(res.minNll()), chi2=chi2, ndf=int(ndf), pvalue=float(pvalue),
        n_warnings=nwarn, n_errors=nerr, messages=message_summary(all_text), n_attempts=len(texts),
        params_at_limit=at_limit,
        params={f"p{k}": dict(value=p.getVal(), error=p.getError(), min=p.getMin(),
                              max=p.getMax()) for k, p in enumerate(params, start=1)},
        norm=dict(value=norm.getVal(), error=norm.getError()),
        expected=nu.tolist(), curve=(xs.tolist(), curve), bin_integration=integrate,
        formula=formula,
        basis=dict(terms=terms_of(name), m0=m0, M=M.tolist(), tbar=tbar.tolist()),
        params_original=original_params(M, params, res),
        _keep=(x, params, pdf, norm, model, data, res),
    )


def failed_fit(name, nbins, exc):
    """Result record for a fit that could not be performed (an exception was raised)."""
    fam, npar = parse_function(name)
    return dict(name=name, family=fam, npar=npar, status=-1, cov_qual=-1, edm=float("nan"),
                min_nll=float("nan"), chi2=float("inf"), ndf=int(nbins - npar), pvalue=0.0,
                n_warnings=0, n_errors=1, messages=[f"{type(exc).__name__}: {exc}"],
                params_at_limit=[], params={}, norm=dict(value=float("nan"), error=float("nan")),
                expected=[], curve=([], []), bin_integration=False, formula=None, basis=None,
                params_original={}, exception=f"{type(exc).__name__}: {exc}")


def original_params(M, params, res):
    """The fitted coefficients p_k of the terms T_k (p = M q), with errors from the fit
    covariance, for reporting."""
    q = np.array([p.getVal() for p in params])
    pv = M @ q
    try:
        names = [p.GetName() for p in params]
        fl = res.floatParsFinal()
        idx = [fl.index(fl.find(nm)) for nm in names]
        C = np.array([[res.covarianceMatrix()[i][j] for j in idx] for i in idx])
        pe = np.sqrt(np.clip(np.diag(M @ C @ M.T), 0, None))
    except Exception:
        pe = np.full(len(pv), float("nan"))
    return {f"p{k}": dict(value=float(v), error=float(e))
            for k, (v, e) in enumerate(zip(pv, pe), start=1)}


def fit_quality_ok(fit):
    return fit["status"] == 0 and fit["cov_qual"] == 3 and not fit["params_at_limit"]


def ftest(fit_lo, fit_hi, nbins):
    """p-value of the Fisher F-test: is fit_hi (more parameters) significantly better?"""
    import ROOT
    dp = fit_hi["npar"] - fit_lo["npar"]
    ndf_hi = nbins - fit_hi["npar"]
    if dp <= 0 or ndf_hi <= 0 or not math.isfinite(fit_lo["chi2"]) or fit_hi["chi2"] <= 0:
        return None
    F = ((fit_lo["chi2"] - fit_hi["chi2"]) / dp) / (fit_hi["chi2"] / ndf_hi)
    if F <= 0:
        return 1.0
    return float(1.0 - ROOT.Math.fdistribution_cdf(F, dp, ndf_hi))


def acceptable(fit, pvalue_threshold):
    """A valid background fit: status 0, full covariance matrix (covQual 3), no parameter at its
    limit, and goodness-of-fit p-value above the threshold."""
    return fit_quality_ok(fit) and fit["pvalue"] > pvalue_threshold


def select_functions(fits, nbins, pvalue_threshold=0.05, ftest_threshold=0.05, max_alternates=4):
    """Choose the functions of the RooMultiPdf: a default and up to `max_alternates` alternates.

    1. Only acceptable fits are used (see `acceptable`; p-value > pvalue_threshold).
    2. In every family, a sequence of standard Fisher F-tests fixes the highest order: start
       from the lowest acceptable order; move to the next order if it is significantly better
       (F-test p-value < ftest_threshold) and its own fit is acceptable; otherwise stop. The
       family contributes its acceptable orders from the lowest one up to this highest order.
    3. Default: among the highest orders of the families, the one with the largest p-value.
    4. Alternates: the other contributed functions, ranked by p-value, at most max_alternates.

    Every fit gets fit['ftest_next'] (F-test p-value against the next order of its family,
    when both fits have a good quality) and fit['ftest_max'] (the highest order of its family).
    Returns (verdict {name: 'pass (default)', 'pass (alternate k)' or the reason it is not used},
    default name or None)."""
    verdict, candidates, family_max = {}, [], {}
    for fam in sorted({f["family"] for f in fits}):
        ff = sorted((f for f in fits if f["family"] == fam), key=lambda f: f["npar"])
        by_npar = {f["npar"]: f for f in ff}
        for f in ff:
            nxt = by_npar.get(f["npar"] + 1)
            f["ftest_next"] = (ftest(f, nxt, nbins) if nxt and fit_quality_ok(f) and fit_quality_ok(nxt)
                               else None)
        start = next((f for f in ff if acceptable(f, pvalue_threshold)), None)
        top = start
        while top is not None:
            nxt = by_npar.get(top["npar"] + 1)
            if (nxt is None or top["ftest_next"] is None or top["ftest_next"] >= ftest_threshold
                    or not acceptable(nxt, pvalue_threshold)):
                break
            top = nxt
        for f in ff:
            f["ftest_max"] = top["name"] if top else None
        if top is not None:
            family_max[fam] = top
        for f in ff:
            if f.get("exception"):
                verdict[f["name"]] = "fail: the fit raised an exception"
            elif not fit_quality_ok(f):
                why = []
                if f["status"] != 0:
                    why.append(f"status {f['status']}")
                if f["cov_qual"] != 3:
                    why.append(f"covQual {f['cov_qual']}")
                if f["params_at_limit"]:
                    why.append("parameter at limit")
                verdict[f["name"]] = "fail: " + ", ".join(why)
            elif f["pvalue"] <= pvalue_threshold:
                verdict[f["name"]] = f"fail: p-value {f['pvalue']:.3g}"
            elif f["npar"] > top["npar"]:
                if top["npar"] + 1 not in by_npar:  # e.g. --functions dijet3 dijet5
                    verdict[f["name"]] = (f"not used: the {fam}{top['npar'] + 1} fit is missing, so no "
                                          f"F-test beyond {top['name']}")
                else:
                    verdict[f["name"]] = f"not used: above the F-test maximum ({top['name']})"
            else:
                candidates.append(f)
    if not family_max:
        return verdict, None
    default = max(family_max.values(), key=lambda f: f["pvalue"])
    verdict[default["name"]] = "pass (default)"
    others = sorted((f for f in candidates if f is not default), key=lambda f: -f["pvalue"])
    for k, f in enumerate(others, start=1):
        verdict[f["name"]] = (f"pass (alternate {k})" if k <= max_alternates
                              else f"not used: beyond the {max_alternates} alternates")
    return verdict, default["name"]
