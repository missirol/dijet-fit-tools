"""Unit tests for dijetfit.signal (Crystal Ball fits) and the plotting helper."""
import numpy as np
import pytest

import synthetic as syn
from conftest import make_th1
from dijetfit import signal as S
from dijetfit.plotting import fit_plot

TRUTH = dict(mean=300.0, sigma=25.0, alphaL=0.8, nL=5.0, alphaR=1.5, nR=4.0)  # mean: see truth_for
WEIGHT = 0.01
N_MC = 2.0e5


def truth_for(fit_range):
    """Signal truth for the Asimov tests: TRUTH, centred in the configured fit range."""
    return dict(TRUTH, mean=0.5 * (fit_range[0] + fit_range[1]))


def test_numpy_crystal_ball_matches_roofit(ROOT):  # pure function test, no fit range involved
    x = ROOT.RooRealVar("x_cb", "", 100.0, 700.0)
    pars = {k: ROOT.RooRealVar(f"{k}_cb", "", v) for k, v in TRUTH.items()}
    cb = ROOT.RooCrystalBall("cb_t", "", x, *[pars[k] for k in S.CB_PARAMS[:6]])
    ms = np.array([155.0, 230.0, 279.0, 300.0, 333.0, 400.0, 590.0])
    ref = []
    for m in ms:
        x.setVal(m)
        ref.append(cb.getVal())  # unnormalised value
    assert syn.crystal_ball(ms, **TRUTH) == pytest.approx(np.array(ref), rel=1e-10)


def fit_edges(fit_range):
    xmin, xmax = fit_range
    return np.linspace(xmin, xmax, int(round(xmax - xmin)) + 1)  # 1 GeV bins


def asimov_template(ROOT, fit_range, name):
    """Weighted Asimov template over the fit range: MC counts = N_MC x probability integrated over
    each bin, each event with weight WEIGHT (contents = WEIGHT x counts, errors =
    sqrt(WEIGHT^2 x counts)). Fitted with bin integration (fixture `integrate_bins`), the truth
    is recovered exactly; with bin-centre evaluation it is recovered up to the small O(w^2)
    approximation (see test_background.test_bin_centre_evaluation_at_1gev)."""
    truth = truth_for(fit_range)
    edges = fit_edges(fit_range)
    f = lambda m: syn.crystal_ball(m, **truth)  # noqa: E731
    counts = N_MC * syn.bin_integrals(f, edges) / syn.bin_integrals(f, list(fit_range), nsub=4000)[0]
    return make_th1(ROOT, name, edges, WEIGHT * counts, errors=WEIGHT * np.sqrt(counts))


@pytest.fixture
def integrate_bins(monkeypatch):
    """Integrate the pdf over the bins also for 1 GeV bins."""
    from dijetfit import background as B
    monkeypatch.setattr(B, "BIN_INTEGRATION_ABOVE", 0.5)


def test_common_weight(ROOT, fit_range):
    h = asimov_template(ROOT, fit_range, "cw")
    assert S.common_weight(h) == pytest.approx(WEIGHT)
    h.SetBinError(100, 2 * h.GetBinError(100))
    assert S.common_weight(h) is None


@pytest.mark.parametrize("mode", ["sumw2", "counts"])
def test_fit_recovers_truth_on_asimov(ROOT, fit_range, integrate_bins, mode):
    h = asimov_template(ROOT, fit_range, f"asimov_{mode}")
    truth = truth_for(fit_range)
    r = S.fit_signal(ROOT, h, truth["mean"], f"t_{mode}", weights=mode)
    assert r["status"] == 0 and r["cov_qual"] == 3 and not r["params_at_limit"]
    assert r["n_errors"] == 0 and not r["fixed"]
    assert r["weights_mode"] == mode and r["bin_integration"]
    # Minuit stops at EDM < 1e-3, i.e. within ~0.05 sigma of the exact minimum
    for k, v in truth.items():
        assert abs(r["params"][k]["value"] - v) < 0.1 * r["params"][k]["error"], k
    assert abs(r["params"]["norm"]["value"] - h.Integral()) < 0.1 * r["params"]["norm"]["error"]
    assert r["chi2"] < 0.05  # zero up to Minuit's convergence tolerance


def test_sumw2_and_counts_give_the_same_uncertainties(ROOT, fit_range, integrate_bins):
    """For uniform weights, SumW2-corrected errors equal the MC-statistics errors."""
    h = asimov_template(ROOT, fit_range, "asimov_err")
    m = truth_for(fit_range)["mean"]
    a = S.fit_signal(ROOT, h, m, "t_e1", weights="sumw2")
    b = S.fit_signal(ROOT, h, m, "t_e2", weights="counts")
    c = S.fit_signal(ROOT, h, m, "t_e3", weights="none")
    # equal up to the numerical precision of Hesse's finite-difference second derivatives
    for k in S.ALL_PARAMS:
        assert a["params"][k]["error"] == pytest.approx(b["params"][k]["error"], rel=0.05), k
    # without correction the errors are wrong by ~1/sqrt(w): exactly for the (parabolic) yield and
    # mean; for the other parameters the flat weighted likelihood also distorts Hesse's estimate
    for k in ("norm", "mean"):
        assert c["params"][k]["error"] == pytest.approx(b["params"][k]["error"] / np.sqrt(WEIGHT), rel=0.02), k
    for k in S.ALL_PARAMS:
        assert c["params"][k]["error"] > 0.8 * b["params"][k]["error"] / np.sqrt(WEIGHT), k
    # norm uncertainty = MC statistics: weight x sqrt(N_MC in the range)
    assert b["params"]["norm"]["error"] == pytest.approx(WEIGHT * np.sqrt(h.Integral() / WEIGHT), rel=0.02)


def test_undetermined_tail_is_fixed(ROOT, fit_range):
    """A right tail starting beyond the fit range is fixed, not left floating."""
    edges = fit_edges(fit_range)
    # peak 40 GeV below the upper edge; right tail onset at mean + 3 sigma = upper edge + 35 GeV
    truth = dict(TRUTH, mean=fit_range[1] - 40.0, alphaR=3.0)
    c = 0.5 * (edges[:-1] + edges[1:])
    counts = 1e5 * syn.crystal_ball(c, **truth) / syn.crystal_ball(c, **truth).sum()
    h = make_th1(ROOT, "tail", edges, counts)
    r = S.fit_signal(ROOT, h, truth["mean"], "t_tail")
    assert r["status"] == 0 and r["cov_qual"] == 3
    assert {"alphaR", "nR"} <= set(r["fixed"])
    assert r["fixed"]["alphaR"]["reason"] == "tail outside the fit range"


def test_fit_plot_writes_files(ROOT, tmp_path):
    edges = [0.0, 1.0, 2.0, 3.0]
    h = make_th1(ROOT, "plot", edges, [10.0, 20.0, 15.0])
    fit_plot(ROOT, str(tmp_path / "p"), ["png", "pdf"], h, [11.0, 19.0, 15.0],
             ([0.0, 3.0], [12.0, 14.0]), ["label"])
    assert (tmp_path / "p.png").stat().st_size > 1000
    assert (tmp_path / "p.pdf").stat().st_size > 1000
