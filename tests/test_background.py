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
    assert B.function_names(["dijet"]) == ["dijet3", "dijet4", "dijet5"]
    assert B.parse_function("ua25") == ("ua2", 5)
    assert B.terms_of("dijet4") == ["ln1mx", "L1", "L2"]
    for bad in ("dijet2", "dijet7", "expo3"):
        with pytest.raises(ValueError):
            B.parse_function(bad)


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
    (chi2 against the true expectation < 2 for 450 bins), while the coefficients move along a
    nearly flat direction by less than 1 standard deviation."""
    h = asimov(ROOT, fit_range, integrate=True, name="asimov_cen")
    r = B.fit_function(ROOT, "dijet4", h, "t_asimov_cen")
    assert not r["bin_integration"]
    assert r["status"] == 0 and r["cov_qual"] == 3 and r["n_errors"] == 0
    assert r["chi2"] < 2.0
    assert np.all(np.abs(_pulls(r, published_to_terms(*TRUTH, B.reference_mass(*fit_range)))) < 1.0)


def test_wide_bins_are_integrated(ROOT, fit_range, monkeypatch):
    """Regression test: with 10 GeV bins and ~10^7 events, evaluating the pdf at the bin centres
    is biased by many standard deviations; the fit must integrate over the bins."""
    h = asimov(ROOT, fit_range, width=10.0, integrate=True, name="asimov10")
    r = B.fit_function(ROOT, "dijet4", h, "t_int")
    assert r["bin_integration"] and r["status"] == 0 and r["n_errors"] == 0
    assert r["chi2"] < 0.05
    monkeypatch.setattr(B, "BIN_INTEGRATION_ABOVE", 1e9)  # force bin-centre evaluation
    r_c = B.fit_function(ROOT, "dijet4", h, "t_cen")
    assert not r_c["bin_integration"]
    assert r_c["chi2"] > 100 * max(r["chi2"], 1e-3)


def test_pdf_carries_precise_integrator(ROOT, fit_range):
    x = ROOT.RooRealVar("x_i", "", *fit_range)
    p = ROOT.RooRealVar("p_i", "", 1.0)
    pdf = B.make_pdf(ROOT, "pdf_i", "dijet3", ROOT.RooArgList(x, p, p), "exp(@1*log(@0))")
    cfg = pdf.specialIntegratorConfig()
    assert cfg.epsRel() == pytest.approx(B.INTEGRATOR_EPS_REL)
    assert cfg.method1D().getCurrentLabel() == "RooAdaptiveGaussKronrodIntegrator1D"


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


def test_ftest_and_selection():
    nb = 450
    # adding one parameter improves chi2 by 0.5: not significant (F-test reported for information)
    assert B.ftest(_fake("dijet4", 450.5), _fake("dijet5", 450.0), nb) > 0.4
    assert B.ftest(_fake("dijet4", 520.0), _fake("dijet5", 450.0), nb) < 1e-10
    fits = [_fake("dijet3", 900, p=1e-20), _fake("dijet4", 455, p=0.4), _fake("dijet5", 454.8, p=0.4),
            _fake("ua23", 470, p=0.2, ok=False), _fake("ua24", 460, p=0.3), _fake("ua25", 459, p=0.3),
            _fake("moddijet3", 480, p=0.06), _fake("moddijet4", 440, p=0.8)]
    v = B.select_functions(fits, nb, 0.05)
    # per family: the acceptable function with the fewest parameters
    assert v["dijet3"].startswith("fail: p-value")
    assert v["dijet4"] == "pass"
    assert v["dijet5"] == "not used: dijet4 has fewer parameters"
    assert v["ua23"].startswith("fail: status")
    assert v["ua24"] == "pass" and v["ua25"].startswith("not used: ua24")
    # moddijet3 passes p > 0.05 although the next order is much better (by design)
    assert v["moddijet3"] == "pass" and v["moddijet4"].startswith("not used: moddijet3")
    assert sorted(n for n, x in v.items() if x == "pass") == ["dijet4", "moddijet3", "ua24"]
    assert fits[6]["ftest_next"] < 1e-6  # F-test vs moddijet4 stored for information


def test_parameter_at_limit_is_reported_after_the_last_attempt(ROOT, fit_range, monkeypatch):
    """Regression test: a parameter that still ends at its limit after the last re-centring must
    be reported, and the fit must fail the quality requirement, even though Minuit reports a
    good fit (status 0, covQual 3). A narrow parameter range (QRANGE, in units of the parameter
    uncertainties) makes the parameters hit their limits: the minimum of the bin-centre fit is
    ~0.7 sigma away from the starting values for this histogram."""
    h = asimov(ROOT, fit_range, integrate=True, name="asimov_lim")
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


def test_yield_at_limit_is_reported(ROOT, fit_range):
    """The yield parameter is also checked: an expectation outside its range is flagged."""
    h = asimov(ROOT, fit_range, integrate=True, name="asimov_norm")
    r = B.fit_function(ROOT, "dijet4", h, "t_norm")
    assert not r["params_at_limit"]
    x = ROOT.RooRealVar("x_n", "", 0.0, 0.0, 10.0)
    assert B.at_limit_of(x) and not B.at_limit_of(ROOT.RooRealVar("y_n", "", 5.0, 0.0, 10.0))
