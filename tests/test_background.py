"""Unit tests for dijetfit.background: function definitions, starting values, fits on Asimov
histograms (contents equal to the expectation, no fluctuations), goodness-of-fit and selection."""
import numpy as np
import pytest

import synthetic as syn
from conftest import make_th1
from dijetfit import background as B

TRUTH = (10.0, 13.0, 1.0)  # dijet4 (p1, p2, p3) in the published form
N_EVENTS = 2.0e7


def published_to_terms(p1, p2, p3, m0):
    """Coefficients of the terms (ln(1-x), L, L^2) of dijet4 for the published parameters, with
    L = ln(m/m0): ln f = p1 ln(1-x) - (p2 + p3 ln x) ln x, ln x = L + c, c = ln(m0/sqrt(s))."""
    c = np.log(m0 / syn.SQRTS)
    return np.array([p1, -(p2 + 2 * p3 * c), -p3])


def asimov(ROOT, fit_range, width=1.0, integrate=False, name="asimov"):
    """Histogram over the fit range whose contents are the dijet4 expectation: bin-centre density
    x width (the model RooFit uses without bin integration) or the integral over each bin.
    The number of bins is rounded so that they fill the fit range exactly."""
    XMIN, XMAX = fit_range
    nb = max(1, int(round((XMAX - XMIN) / width)))
    edges = np.linspace(XMIN, XMAX, nb + 1)
    lo, hi = edges[:-1], edges[1:]
    f = lambda m: syn.dijet4(m, *TRUTH)  # noqa: E731
    norm = syn.bin_integrals(f, [XMIN, XMAX], nsub=2000)[0]
    nu = syn.bin_integrals(f, edges) if integrate else f(0.5 * (lo + hi)) * (hi - lo)
    return make_th1(ROOT, name, edges, N_EVENTS * nu / norm)


def test_function_names_and_parsing():
    assert B.function_names(["dijet"]) == ["dijet3", "dijet4", "dijet5", "dijet6"]
    assert B.function_names(["ua2"], orders=(3, 4)) == ["ua23", "ua24"]
    assert B.parse_function("ua25") == ("ua2", 5)
    assert B.terms_of("dijet4") == ["ln1mx", "L1", "L2"]
    for bad in ("dijet2", "dijet7", "expo3", "modexp6", "polypow6", "polyext4", "ua22"):
        with pytest.raises(ValueError):
            B.parse_function(bad)
    # the families of the analysis note, with their orders
    assert B.function_names(["modexp", "polypow", "polyext"]) == [
        "modexp3", "modexp4", "modexp5", "polypow3", "polypow4", "polypow5", "polyext5", "polyext6"]
    assert len(B.function_names()) == 20 and B.parse_function("polyext6") == ("polyext", 6)
    assert B.is_loglinear("ua23") and not B.is_loglinear("modexp4") and B.terms_of("polypow4") == []


def test_reference_mass_is_geometric_mean():  # pure arithmetic, independent of the configured range
    assert B.reference_mass(150.0, 600.0) == pytest.approx(300.0)
    assert B.reference_mass(120.0, 700.0) == pytest.approx(np.sqrt(120 * 700))


def test_formula_matches_numpy_terms(ROOT, fit_range):
    """The RooFit expression and the numpy terms describe the same function."""
    XMIN, XMAX = fit_range
    m0 = 280.0
    coef = np.array([7.0, -4.0, -0.8, 0.3])
    x = ROOT.RooRealVar("x_f", "", XMIN, XMAX)
    pars = [ROOT.RooRealVar(f"c{k}", "", v) for k, v in enumerate(coef)]
    al = ROOT.RooArgList(x)
    for p in pars:
        al.add(p)
    fv = ROOT.RooFormulaVar("fv", B.formula_of("dijet5", m0), al)
    T = B.term_defs(m0)
    for m in np.linspace(XMIN, XMAX, 7)[1:-1]:
        x.setVal(m)
        expected = np.exp(sum(c * T[t][0](np.array([m]))[0] for c, t in zip(coef, B.terms_of("dijet5"))))
        assert fv.getVal() == pytest.approx(expected, rel=1e-10)


def test_start_values_recover_truth_on_asimov(ROOT, fit_range):
    """IRLS on the bin-centre Asimov histogram gives back the true coefficients exactly."""
    h = asimov(ROOT, fit_range)
    m0 = B.reference_mass(*fit_range)
    p, mu = B.start_values("dijet4", h, m0)
    assert p == pytest.approx(published_to_terms(*TRUTH, m0), rel=1e-7)
    assert mu.sum() == pytest.approx(h.Integral(), rel=1e-9)


def test_decorrelated_basis_is_orthonormal(ROOT, fit_range):
    h = asimov(ROOT, fit_range)
    m0 = B.reference_mass(*fit_range)
    _, mu = B.start_values("dijet5", h, m0)
    M, tbar = B.decorrelated_basis("dijet5", h, mu, m0)
    lo, hi, _ = B.hist_arrays(h)
    T = np.column_stack([B.term_defs(m0)[t][0](0.5 * (lo + hi)) for t in B.terms_of("dijet5")])
    U = (T - tbar) @ M
    assert (U * mu[:, None]).T @ U == pytest.approx(np.eye(4), abs=1e-6)


def _pulls(r, truth):
    return np.array([(r["params_original"][f"p{k}"]["value"] - t) / r["params_original"][f"p{k}"]["error"]
                     for k, t in enumerate(truth, start=1)])


def test_fit_recovers_truth_on_asimov_with_bin_integration(ROOT, fit_range, monkeypatch):
    """With the pdf integrated over the bins, the fit of the integrated Asimov histogram returns
    the true coefficients and yield exactly."""
    monkeypatch.setattr(B, "BIN_INTEGRATION_ABOVE", 0.5)  # integrate also 1 GeV bins
    h = asimov(ROOT, fit_range, integrate=True, name="asimov_int")
    r = B.fit_function(ROOT, "dijet4", h, "t_asimov_int")
    assert r["bin_integration"]
    assert r["status"] == 0 and r["cov_qual"] == 3 and not r["params_at_limit"]
    assert r["n_errors"] == 0
    # the only message allowed: RooFit flags the non-integer Asimov contents as weighted data
    assert all("weighted data" in m for m in r["messages"])
    assert r["chi2"] < 1e-3
    assert np.abs(_pulls(r, published_to_terms(*TRUTH, B.reference_mass(*fit_range)))) == pytest.approx(0, abs=0.01)
    assert r["norm"]["value"] == pytest.approx(h.Integral(), rel=1e-6)
    assert np.array(r["expected"]) == pytest.approx(B.hist_arrays(h)[2], rel=1e-5)


def test_bin_centre_evaluation_at_1gev(ROOT, fit_range):
    """Default for 1 GeV bins: the pdf is evaluated at the bin centres. This documents the size of
    the resulting approximation for 2x10^7 events: the predicted curve is essentially unchanged
    (chi2 against the true expectation < 2 for the ~450 bins of the default 150-600 GeV range),
    while the coefficients move along a
    nearly flat direction by less than 1 standard deviation."""
    h = asimov(ROOT, fit_range, integrate=True, name="asimov_cen")
    r = B.fit_function(ROOT, "dijet4", h, "t_asimov_cen")
    assert not r["bin_integration"]
    assert r["status"] == 0 and r["cov_qual"] == 3 and r["n_errors"] == 0
    assert r["chi2"] < 2.0
    assert np.all(np.abs(_pulls(r, published_to_terms(*TRUTH, B.reference_mass(*fit_range)))) < 1.0)


def test_wide_bins_are_integrated(ROOT, fit_range, monkeypatch):
    """Regression test, 10 GeV bins and ~10^7 events (truth: the function integrated over the
    bins). The fit integrates over the bins and reproduces the truth. Evaluating the pdf at the
    bin centres also works if it is normalised by its sum over the bin centres (the expected
    counts then add up to the yield), but is biased by many standard deviations with RooFit's
    default normalisation by the integral (the expected counts do not add up to the yield)."""
    h = asimov(ROOT, fit_range, width=10.0, integrate=True, name="asimov10")
    r = B.fit_function(ROOT, "dijet4", h, "t_int")
    assert r["bin_integration"] and r["status"] == 0 and r["n_errors"] == 0
    assert r["chi2"] < 0.05
    monkeypatch.setattr(B, "BIN_INTEGRATION_ABOVE", 1e9)  # force bin-centre evaluation
    r_c = B.fit_function(ROOT, "dijet4", h, "t_cen")
    assert not r_c["bin_integration"] and r_c["chi2"] < 0.05
    monkeypatch.setattr(B, "bin_sum_nbins", lambda lo, hi: None)  # normalise by the integral
    r_i = B.fit_function(ROOT, "dijet4", h, "t_cen_int")
    assert not r_i["bin_integration"]
    assert r_i["chi2"] > 100 * max(r["chi2"], 1e-3)


def test_pdf_carries_precise_integrator(ROOT, fit_range):
    x = ROOT.RooRealVar("x_i", "", *fit_range)
    p = ROOT.RooRealVar("p_i", "", 1.0)
    pdf = B.make_pdf(ROOT, "pdf_i", "dijet3", ROOT.RooArgList(x, p, p), "exp(@1*log(@0))")
    cfg = pdf.specialIntegratorConfig()
    assert cfg.epsRel() == pytest.approx(B.INTEGRATOR_EPS_REL)
    assert cfg.method1D().getCurrentLabel() == "RooAdaptiveGaussKronrodIntegrator1D"


def test_bin_sum_normalisation(ROOT, fit_range):
    """With uniform bins, the pdf is normalised by its sum over the bin centres: the expected
    counts N * pdf(centre) * width add up exactly to N (with the integral as normalisation they
    differ by the midpoint-rule error, ~1e-5 here, which biases fits with ~1e7 events)."""
    xmin, xmax = fit_range
    nb = int(round(xmax - xmin))
    lo = np.arange(nb) + xmin
    hi = lo + 1.0
    assert B.bin_sum_nbins(lo, hi) == nb
    assert B.bin_sum_nbins(np.array([0.0, 2.0]), np.array([2.0, 4.0])) is None  # integrated
    assert B.bin_sum_nbins(np.array([0.0, 1.0]), np.array([1.0, 1.5])) is None  # not uniform
    x = ROOT.RooRealVar("x_bs", "", xmin, xmax)
    p = ROOT.RooRealVar("p_bs", "", -0.03)
    formula = "exp(@1*@0)"
    pdf = B.make_pdf(ROOT, "pdf_bs", "test", ROOT.RooArgList(x, p), formula, nb)
    assert pdf.specialIntegratorConfig().method1D().getCurrentLabel() == "RooBinIntegrator"
    nu = B.expected_counts(ROOT, x, pdf, 1e7, lo, hi, False)
    assert nu.sum() == pytest.approx(1e7, rel=1e-12)
    ref = B.make_pdf(ROOT, "pdf_bs_int", "test", ROOT.RooArgList(x, p), formula)
    nu_int = B.expected_counts(ROOT, x, ref, 1e7, lo, hi, False)
    assert abs(nu_int.sum() / 1e7 - 1) > 1e-6  # the mismatch that the bin sum removes


def test_baker_cousins():
    n = np.array([0.0, 3.0, 10.0])
    assert B.baker_cousins(n, n + (n == 0) * 1e-300) == pytest.approx(0.0, abs=1e-12)
    assert B.baker_cousins([5.0], [4.0]) == pytest.approx(2 * (4 - 5 + 5 * np.log(5 / 4)))
    assert B.baker_cousins([0.0], [2.0]) == pytest.approx(4.0)
    assert B.baker_cousins([1.0], [0.0]) == float("inf")


def _fake(name, chi2, ok=True, p=0.5):
    fam, n = B.parse_function(name)
    return dict(name=name, family=fam, npar=n, chi2=chi2, pvalue=p, status=0 if ok else 3,
                cov_qual=3, params_at_limit=[])


def test_ftest():
    nb = 450
    # adding one parameter improves chi2 by 0.5: not significant; by 70: very significant
    assert B.ftest(_fake("dijet4", 450.5), _fake("dijet5", 450.0), nb) > 0.4
    assert B.ftest(_fake("dijet4", 520.0), _fake("dijet5", 450.0), nb) < 1e-10
    assert B.ftest(_fake("dijet4", 450.0), _fake("dijet5", 451.0), nb) == 1.0  # worse: no gain


def test_selection_follows_standard_ftest_sequence():
    nb = 450
    fits = [
        # dijet: 3 fails p; 4 acceptable; 5 significantly better than 4 and acceptable; 6 not
        # significantly better than 5 -> family maximum dijet5
        _fake("dijet3", 900, p=1e-20), _fake("dijet4", 470, p=0.06), _fake("dijet5", 440, p=0.6),
        _fake("dijet6", 439.8, p=0.6),
        # moddijet: 3 acceptable, 4 not significantly better -> maximum moddijet3 (4 not used even
        # though its p-value is larger)
        _fake("moddijet3", 452, p=0.3), _fake("moddijet4", 451.5, p=0.31),
        # ua2: 3 fails quality; 4 acceptable; 5 significantly better but its p-value fails -> stop
        _fake("ua23", 470, p=0.2, ok=False), _fake("ua24", 470, p=0.08), _fake("ua25", 430, p=0.04),
    ]
    v, default = B.select_functions(fits, nb, 0.05, 0.05, 4)
    assert v["dijet3"].startswith("fail: p-value")
    assert v["dijet6"] == "not used: above the F-test maximum (dijet5)"
    assert v["moddijet4"] == "not used: above the F-test maximum (moddijet3)"
    assert v["ua23"].startswith("fail: status") and v["ua25"].startswith("fail: p-value")
    f = {x["name"]: x for x in fits}
    assert f["dijet4"]["ftest_next"] < 0.05 and f["dijet5"]["ftest_next"] > 0.05
    assert f["dijet3"]["ftest_max"] == "dijet5" and f["ua24"]["ftest_max"] == "ua24"
    # default: the family maximum with the largest p-value (dijet5 0.6 > moddijet3 0.3 > ua24 0.08)
    assert default == "dijet5" and v["dijet5"] == "pass (default)"
    # alternates: the other contributed functions (orders up to the family maxima) by p-value
    assert v["moddijet3"] == "pass (alternate 1)"
    assert v["ua24"] == "pass (alternate 2)"
    assert v["dijet4"] == "pass (alternate 3)"


def test_selection_with_only_the_family_maxima():
    """--alternates family-max: one function per family, its F-test maximum."""
    nb = 450
    fits = [_fake("dijet4", 470, p=0.06), _fake("dijet5", 440, p=0.6), _fake("dijet6", 439.8, p=0.6),
            _fake("moddijet3", 452, p=0.3), _fake("moddijet4", 451.5, p=0.31),
            _fake("ua24", 470, p=0.08), _fake("ua25", 440, p=0.5)]
    v, default = B.select_functions(fits, nb, 0.05, 0.05, 4, family_max_only=True)
    assert default == "dijet5"
    assert v["ua25"] == "pass (alternate 1)" and v["moddijet3"] == "pass (alternate 2)"
    assert v["dijet4"] == "not used: below the family maximum (dijet5)"
    assert v["ua24"] == "not used: below the family maximum (ua25)"
    assert v["dijet6"] == "not used: above the F-test maximum (dijet5)"
    # the default selection keeps the lower orders
    v, _ = B.select_functions(fits, nb, 0.05, 0.05, 4)
    assert v["dijet4"].startswith("pass (alternate") and v["ua24"].startswith("pass (alternate")


def test_selection_caps_alternates_and_handles_no_acceptable_fit():
    nb = 450
    fits = [_fake(f"{fam}3", 450, p=p) for fam, p in (("dijet", 0.9), ("moddijet", 0.5), ("ua2", 0.4))]
    fits += [_fake("dijet4", 449.9, p=0.88), _fake("moddijet4", 449.9, p=0.45)]
    v, default = B.select_functions(fits, nb, 0.05, 0.05, max_alternates=1)
    assert default == "dijet3"
    assert v["moddijet3"] == "pass (alternate 1)" and v["ua23"] == "not used: beyond the 1 alternates"
    v, default = B.select_functions([_fake("dijet3", 900, p=1e-9)], nb)
    assert default is None and v["dijet3"].startswith("fail")


def test_parameter_at_limit_is_reported_after_the_last_attempt(ROOT, fit_range, monkeypatch):
    """Regression test: a parameter that still ends at its limit after the last re-centring must
    be reported, and the fit must fail the quality requirement, even though Minuit reports a
    good fit (status 0, covQual 3). The starting values are moved by 3 standard deviations
    (in the decorrelated basis) and a narrow parameter range (QRANGE, in units of the parameter
    uncertainties) makes the parameters hit their limits."""
    h = asimov(ROOT, fit_range, integrate=True, name="asimov_lim")
    orig = B.start_values

    def shifted(name, hist, m0):
        p, mu = orig(name, hist, m0)
        M, _ = B.decorrelated_basis(name, hist, mu, m0)
        return p + M @ np.full(len(p), 3.0), mu

    monkeypatch.setattr(B, "start_values", shifted)
    monkeypatch.setattr(B, "QRANGE", 0.05)
    monkeypatch.setattr(B, "MAX_RECENTRE", 1)
    r = B.fit_function(ROOT, "dijet4", h, "t_lim")
    assert r["params_at_limit"], r["params"]
    assert r["status"] == 0 and r["cov_qual"] == 3  # Minuit alone would not notice
    assert not B.fit_quality_ok(r)
    # with enough re-centring attempts the same fit ends inside the (moved) ranges
    monkeypatch.setattr(B, "QRANGE", 1.0)
    monkeypatch.setattr(B, "MAX_RECENTRE", 3)
    r = B.fit_function(ROOT, "dijet4", h, "t_lim2")
    assert not r["params_at_limit"] and B.fit_quality_ok(r)


def test_yield_at_limit_is_reported(ROOT, fit_range, monkeypatch):
    """The yield parameter is also checked: a yield range that excludes the number of events in
    the fit range makes the fitted yield end at its limit, which must be flagged."""
    h = asimov(ROOT, fit_range, integrate=True, name="asimov_norm")
    r = B.fit_function(ROOT, "dijet4", h, "t_norm")
    assert not r["params_at_limit"]
    monkeypatch.setattr(B, "NORM_RANGE", (1.01, 1.5))  # fitted yield ~1.0 x events: below the range
    r = B.fit_function(ROOT, "dijet4", h, "t_norm_lim")
    assert any(n.endswith("_norm") for n in r["params_at_limit"]), r["params_at_limit"]
    assert not B.fit_quality_ok(r)
    x = ROOT.RooRealVar("x_n", "", 0.0, 0.0, 10.0)
    assert B.at_limit_of(x) and not B.at_limit_of(ROOT.RooRealVar("y_n", "", 5.0, 0.0, 10.0))


def test_selection_with_a_missing_order_and_a_failed_fit():
    """An order above a gap in the fitted orders gets a distinct verdict (no F-test was possible);
    a fit that raised an exception is marked as such."""
    fits = [_fake("dijet3", 452, p=0.3), _fake("dijet5", 440, p=0.6),
            B.failed_fit("ua23", 450, RuntimeError("boom"))]
    v, default = B.select_functions(fits, 450)
    assert default == "dijet3"
    assert v["dijet5"].startswith("not used: the dijet4 fit is missing")
    assert v["ua23"] == "fail: the fit raised an exception"


# ---- all families: shared formula, the families that are not log-linear --------------------------
# typical parameters of the non-linear functions: fits to a category of the real data (150-600 GeV)
NL_TRUTH = {
    "modexp3": (-54.82, 0.1767), "modexp4": (26.81, -0.09235, 2.417),
    "modexp5": (-551.18, 1.55835, 9.79123, 93.05),
    "polypow3": (259.0, 5.816), "polypow4": (1831.0, 2.995e4, 3.967), "polypow5": (578.9, -5784.0, 9.629e4, 5.687),
    "polyext5": (-61.11, 9.102, 0.4708, -19.07), "polyext6": (-152.7, 13.78, 0.9088, -38.62, 380.6),
}


def _nl_truth(fit_range, name):
    """NL_TRUTH, or skip if it is not a valid falling spectrum over the configured fit range."""
    p = np.array(NL_TRUTH[name])
    m = np.linspace(*fit_range, 50)
    lf = B.nonlinear_log_f(name, m / B.SQRTS, p)
    if not np.all(np.isfinite(lf)) or not np.all(np.diff(lf) < 0):
        pytest.skip(f"{name}: the typical parameters are not valid over the fit range {fit_range}")
    return p


def test_formula_with_constants_matches_numpy_for_every_function(ROOT, fit_range):
    """The RooFit formula of the decorrelated parametrisation (per-fit constants as arguments)
    describes the same function as numpy, for every function of every family."""
    XMIN, XMAX = fit_range
    rng = np.random.default_rng(3)
    m0 = B.reference_mass(XMIN, XMAX)
    x = ROOT.RooRealVar("x_all", "", XMIN, XMAX)
    for name in B.function_names():
        k = B.parse_function(name)[1] - 1
        M = np.eye(k) + 0.01 * rng.normal(size=(k, k))
        q = 0.1 * rng.normal(size=k)
        if B.is_loglinear(name):
            centre = rng.normal(size=k)
            T = B.term_defs(m0)
            lnf = lambda m: sum(q[j] * sum(M[i, j] * (T[t][0](m) - centre[i])  # noqa: E731
                                           for i, t in enumerate(B.terms_of(name))) for j in range(k))
        else:
            centre = np.array(NL_TRUTH[name])
            lnf = lambda m: B.nonlinear_log_f(name, m / B.SQRTS, centre + M @ q)  # noqa: E731
        formula, consts = B.formula_and_constants(name, k, m0, M, centre)
        al = ROOT.RooArgList(x)
        keep = [ROOT.RooRealVar(f"q{j}_{name}", "", v) for j, v in enumerate(q)]
        for v in keep:
            al.add(v)
        keep += B.add_constants(ROOT, al, consts, "c", name)
        fv = ROOT.RooFormulaVar(f"fv_{name}", formula, al)
        ms = np.linspace(XMIN, XMAX, 6)[1:-1]
        ref = np.array([lnf(np.array([m]))[0] for m in ms])
        got = []
        for m in ms:
            x.setVal(m)
            got.append(np.log(fv.getVal()))
        # ratios between masses (the normalisation is free)
        assert np.array(got) - got[0] == pytest.approx(ref - ref[0], abs=1e-8), name


def test_fits_of_a_function_share_one_compiled_formula(ROOT, fit_range):
    """The formula string does not depend on the data or the fit range (the per-fit constants
    are arguments), so ROOT compiles it once per function."""
    XMIN, XMAX = fit_range
    h1 = asimov(ROOT, fit_range, name="share1")
    h2 = asimov(ROOT, (XMIN + 10.0, XMAX), name="share2")
    for name in ("dijet5", "polyext5"):
        r1 = B.fit_function(ROOT, name, h1, f"share1_{name}")
        r2 = B.fit_function(ROOT, name, h2, f"share2_{name}")
        assert r1["formula"] == r2["formula"]
        assert r1["constants"] != r2["constants"]


# At these typical parameters the likelihood of modexp4 and polypow4 has an almost flat, curved
# direction (e.g. (1 + p1 x + p2 x^2)^-p3 for large p1, p2); Minuit's covariance matrix may then be
# forced positive-definite (covQual 2), and such a fit fails the covQual 3 requirement of the
# selection, as for any function.
NEARLY_FLAT = {"modexp4", "polypow4"}


@pytest.mark.parametrize("name", B.function_names(["modexp", "polypow", "polyext"]))
def test_nonlinear_fit_recovers_its_own_asimov(ROOT, fit_range, name):
    """Asimov histogram (bin-centre expectation) of a non-linear function: the start search and
    the RooFit fit reproduce it exactly (covQual 3, except for NEARLY_FLAT: at least 2)."""
    XMIN, XMAX = fit_range
    p = _nl_truth(fit_range, name)
    edges = np.arange(XMIN, XMAX + 0.5, 1.0)
    c = 0.5 * (edges[:-1] + edges[1:])
    f = np.exp(B.nonlinear_log_f(name, c / B.SQRTS, p))
    h = make_th1(ROOT, f"asimov_nl_{name}", edges, N_EVENTS * f / f.sum())
    p_start, mu = B.start_values_nonlinear(name, h)
    assert B.baker_cousins(B.hist_arrays(h)[2], mu) < 1e-3
    r = B.fit_function(ROOT, name, h, f"t_nl_{name}")
    assert r["status"] == 0 and not r["params_at_limit"], r["messages"]
    assert r["cov_qual"] == 3 or (name in NEARLY_FLAT and r["cov_qual"] == 2), r["cov_qual"]
    assert r["chi2"] < 1e-2 and r["pvalue"] > 0.99
    assert np.array(r["expected"]) == pytest.approx(B.hist_arrays(h)[2], rel=1e-4)
