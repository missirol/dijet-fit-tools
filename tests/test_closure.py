"""End-to-end closure tests on synthetic inputs with a known truth (tests/synthetic.py).

* Main test: the complete workflow (run_all.sh) runs with the configured fit range (pytest
  --fit-xmin/--fit-xmax, default 150-600 GeV) on a background without turn-on; the range must
  stay as configured, and every step must give back the truth.
* Scan test: a background with a turn-on just above the configured lower edge; the background
  step with --scan must raise the lower edge past the turn-on.

The workflow runs once per test module (a few minutes); select or skip these tests with
`pytest -m closure` / `pytest -m "not closure"`. Pass `--basetemp DIR` to keep the outputs.
"""
import json
import math
import os
import re
import subprocess
import sys

import numpy as np
import pytest

import synthetic as syn
from conftest import REPO
from dijetfit import background as B
from dijetfit import combine as C
from dijetfit.histtools import window_hist

pytestmark = pytest.mark.closure
MAX_PULL = 4.0


@pytest.fixture(scope="module")
def run(tmp_path_factory, fit_range):
    base = tmp_path_factory.mktemp("closure")
    truth = syn.generate(base / "inputs", seed=12345)
    out = base / "output"
    env = dict(os.environ, OUT=str(out), **syn.workflow_config(*fit_range))
    proc = subprocess.run([str(REPO / "run_all.sh"), truth["files"]["data"], *truth["files"]["signal"]],
                          env=env, capture_output=True, text=True, cwd=base)
    (base / "run_all.log").write_text(proc.stdout + proc.stderr)
    assert proc.returncode == 0, f"run_all.sh failed, see {base / 'run_all.log'}"
    load = lambda p: json.loads((out / p).read_text())  # noqa: E731
    return dict(base=base, out=out, truth=truth, fit_range=fit_range, bkg=load("background/background_fits.json"),
                sig=load("signal/signal_fits.json"), param=load("signal_param/signal_parametrisation.json"),
                ws=load("workspace/dijet_workspace.json"))


def cat_name(cat):
    return syn.CATEGORY_PREFIX + cat


def expected_masses(fit_range, min_distance=10.0):
    """Signal masses the workflow should use: those at least min_distance inside both edges of
    the fit range (fit_signal.py --min-distance default)."""
    xmin, xmax = fit_range
    return [m for m in syn.FIT_MASSES if xmin + min_distance <= m <= xmax - min_distance]


# ---- step 1 -------------------------------------------------------------------------------------
def test_inputs_categories_and_normalisation(ROOT, run):
    """Only the selected categories are written (inclusive excluded, unrelated histogram ignored),
    the data are copied unchanged and the signal is scaled by lumi / N_processed."""
    f = ROOT.TFile.Open(str(run["out"] / "inputs.root"))
    dirs = sorted(k.GetName() for k in f.GetListOfKeys() if k.GetClassName().startswith("TDirectory"))
    assert dirs == sorted(cat_name(c) for c in syn.CATEGORIES)
    scale = syn.LUMI / syn.N_GENERATED
    for cat, t in run["truth"]["categories"].items():
        assert f.Get(f"{cat_name(cat)}/data_obs").Integral() == t["n_data"]
        for m in syn.FIT_MASSES:  # step 1 copies all masses; masses are dropped in step 3
            h = f.Get(f"{cat_name(cat)}/ZprimeQQ_{m}")
            n_sel = run["truth"]["signal"][str(m)]["categories"][cat]["n_selected"]
            assert h.Integral() == pytest.approx(n_sel * scale, rel=1e-12)
            b = h.GetMaximumBin()  # a well-populated bin: weight^2 x counts = scale x content
            assert h.GetBinError(b) ** 2 == pytest.approx(scale * h.GetBinContent(b), rel=1e-9)


# ---- step 2 -------------------------------------------------------------------------------------
def test_fit_range_is_the_configured_one(run):
    """Without a turn-on, the fits pass at the configured lower edge: although the scan may
    raise it (up to XMIN + 60 GeV), the fit range must stay exactly the configured one."""
    assert float(syn.workflow_config(*run["fit_range"])["SCAN_MAX"]) > run["fit_range"][0]
    assert run["bkg"]["fit_range"] == list(run["fit_range"])
    assert all(c["rebin"] == 1 for c in run["bkg"]["categories"].values())


def expected_counts_of_fit(ROOT, fit, edges, tag):
    x = ROOT.RooRealVar(f"x_{tag}", "", edges[0], edges[-1])
    al = ROOT.RooArgList(x)
    keep = []
    for k in range(1, fit["npar"]):
        v = ROOT.RooRealVar(f"q{k}_{tag}", "", fit["params"][f"p{k}"]["value"])
        al.add(v)
        keep.append(v)
    keep += B.add_constants(ROOT, al, fit.get("constants"), "c", tag)  # the per-fit constants
    lo, hi = np.array(edges[:-1]), np.array(edges[1:])
    pdf = B.make_pdf(ROOT, f"pdf_{tag}", fit["name"], al, fit["formula"], B.bin_sum_nbins(lo, hi))
    return B.expected_counts(ROOT, x, pdf, fit["norm"]["value"], lo, hi, B.needs_bin_integration(lo, hi))


def test_background_closure(ROOT, run):
    """The fit of the generating function (dijet4) reproduces the true expected spectrum in the
    fit range within the statistical precision, in every category; and the selection never
    rejects it for its fit quality or goodness of fit. The RooMultiPdf follows the F-test
    selection: a default (the family maximum with the largest p-value) and up to 4 alternates,
    none above its family's F-test maximum.

    Other functions are only required to describe the data, which is what the selection tests:
    a passing function may still differ from the truth by more than the statistics. That is the
    modelling uncertainty the RooMultiPdf / discrete profiling is meant to cover. The distances
    from the truth of all passing functions are printed for information."""
    xmin, xmax = run["bkg"]["fit_range"]
    e = syn.edges()
    sel = (e[:-1] >= xmin - 1e-9) & (e[1:] <= xmax + 1e-9)
    for cat in syn.CATEGORIES:
        cb = run["bkg"]["categories"][cat_name(cat)]
        # the log-linear families (dijet4 generates the data) always converge; a non-linear function
        # that does not (it may sit in an almost flat valley on these data) must be failed by the
        # selection and not be used
        bad = {n: (f["status"], f["n_errors"], f["messages"][:2]) for n, f in cb["fits"].items()
               if B.is_loglinear(n) and not (f["n_errors"] == 0 and f["status"] == 0)}
        assert not bad, (cat, bad)
        for n, f in cb["fits"].items():
            if not B.fit_quality_ok(f):
                assert f["verdict"].startswith("fail") and n not in cb["passing"], (cat, n, f["verdict"])
        assert cb["passing"], cat
        gen = cb["fits"]["dijet4"]
        assert not gen["verdict"].startswith("fail"), (cat, gen["verdict"])
        # RooMultiPdf: default + at most 4 alternates; default = the family maximum with the largest
        # p-value; alternates ordered by p-value; nothing above a family's F-test maximum
        assert 1 <= len(cb["passing"]) <= 5 and cb["default"] == cb["passing"][0]
        maxima = {cb["fits"][f]["ftest_max"] for f in cb["fits"]} - {None}
        assert cb["default"] in maxima
        assert cb["fits"][cb["default"]]["pvalue"] == max(cb["fits"][f]["pvalue"] for f in maxima)
        alt_p = [cb["fits"][f]["pvalue"] for f in cb["passing"][1:]]
        assert alt_p == sorted(alt_p, reverse=True)
        for f in cb["passing"]:
            top = cb["fits"][f]["ftest_max"]
            assert cb["fits"][f]["npar"] <= cb["fits"][top]["npar"], (f, top)
        nu_true = syn.background_expected(cat)[sel]
        for name in sorted(set(cb["passing"]) | {"dijet4"}):
            fit = cb["fits"][name]
            nu_fit = expected_counts_of_fit(ROOT, fit, cb["bin_edges"], f"{cat}_{name}")
            # distance between the fitted and the true expectation in units of the statistical
            # precision: for an unbiased fit, chi2-distributed with ~npar degrees of freedom
            dist = float(np.sum((nu_fit - nu_true) ** 2 / nu_true))
            pull = (fit["norm"]["value"] - nu_true.sum()) / fit["norm"]["error"]
            print(f"{cat} {name} ({fit['verdict']}): distance from truth {dist:.1f} (npar {fit['npar']}), "
                  f"yield pull {pull:+.2f}")
            if name == "dijet4":
                assert dist < 3 * fit["npar"], (cat, name, dist)
                assert abs(pull) < MAX_PULL, (cat, name, pull)


# ---- step 3 -------------------------------------------------------------------------------------
def true_yield_in_range(run, cat, mass, xmin, xmax):
    e = syn.edges()
    sel = (e[:-1] >= xmin - 1e-9) & (e[1:] <= xmax + 1e-9)
    frac = syn.signal_probabilities(mass)[sel].sum()
    n_sel = run["truth"]["signal"][str(mass)]["categories"][cat]["n_selected"]
    return n_sel * frac * syn.LUMI / syn.N_GENERATED


def fitted_signal_counts(fit, edges, scale):
    """Expected MC counts per bin from a signal fit (bin-centre density x width, normalised to
    the fitted yield in the fit range), in units of generated events."""
    p = {k: fit["params"][k]["value"] for k in ("mean", "sigma", "alphaL", "nL", "alphaR", "nR")}
    lo, hi = np.array(edges[:-1]), np.array(edges[1:])
    dens = syn.crystal_ball(0.5 * (lo + hi), **p) * (hi - lo)
    return fit["params"]["norm"]["value"] / scale * dens / dens.sum()


def test_signal_closure(ROOT, run):
    """Each signal fit converges and closes on the generating model, judged with the likelihood
    ratio between the fitted and the true expectation for the MC template:
        dchi2 = chi2(truth) - chi2(fit)   (Poisson likelihood-ratio chi2 on the MC counts)
    The fit must be at least as good as the truth (dchi2 >= 0 up to numerical precision: else
    the minimisation failed), and the truth must be compatible with the fit (dchi2 below the
    chi2 quantile at p = 1e-4 for the number of floating parameters). This accounts for the
    strong correlations between the tail parameters (e.g. alphaR and nR), whose individual
    pulls are not standard-normal. The peak position and the yield are also checked directly."""
    xmin, xmax = run["sig"]["fit_range"]
    assert run["sig"]["weights"] == "sumw2"
    scale = syn.LUMI / syn.N_GENERATED
    e = syn.edges()
    sel = (e[:-1] >= xmin - 1e-9) & (e[1:] <= xmax + 1e-9)
    f = ROOT.TFile.Open(str(run["out"] / "inputs.root"))
    for cat in syn.CATEGORIES:
        cs = run["sig"]["categories"][cat_name(cat)]
        assert cs["masses"] == expected_masses((xmin, xmax))
        edges = run["bkg"]["categories"][cat_name(cat)]["bin_edges"]
        for m in cs["masses"]:
            fit = cs["fits"][str(m)]
            assert fit["status"] == 0 and fit["cov_qual"] == 3 and fit["n_errors"] == 0
            h = window_hist(ROOT, f.Get(f"{cat_name(cat)}/ZprimeQQ_{m}"), xmin, xmax,
                            run["bkg"]["categories"][cat_name(cat)]["rebin"], name=f"sc_{cat}_{m}")
            n = np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX() + 1)]) / scale
            n_sel = run["truth"]["signal"][str(m)]["categories"][cat]["n_selected"]
            nu_true = n_sel * syn.signal_probabilities(m)[sel]
            nu_fit = fitted_signal_counts(fit, edges, scale)
            nfloat = sum(not v["constant"] for v in fit["params"].values())
            dchi2 = B.baker_cousins(n, nu_true) - B.baker_cousins(n, nu_fit)
            limit = ROOT.Math.chisquared_quantile_c(1e-4, nfloat)
            pm = (fit["params"]["mean"]["value"] - syn.signal_params(m)["mean"]) / fit["params"]["mean"]["error"]
            ty = true_yield_in_range(run, cat, m, xmin, xmax)
            nn = fit["params"]["norm"]
            print(f"{cat} M{m}: chi2(truth) - chi2(fit) = {dchi2:.1f} (limit {limit:.1f}, {nfloat} floating "
                  f"parameters), mean pull {pm:+.2f}, yield {nn['value']:.4g} vs {ty:.4g}")
            assert dchi2 > -0.5, (cat, m, dchi2)
            assert dchi2 < limit, (cat, m, dchi2, limit)
            assert abs(pm) < MAX_PULL, (cat, m, pm)
            # the template's own fluctuation is part of the truth here (fixed n_selected), so the
            # fitted yield must reproduce the expected yield in range closely
            assert abs(nn["value"] - ty) < MAX_PULL * nn["error"], (cat, m, nn, ty)


# ---- steps 4-5 ----------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def workspace(ROOT, run):
    assert ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit") >= 0
    f = ROOT.TFile.Open(str(run["out"] / "workspace" / "dijet_workspace.root"))
    yield f.Get("w")
    f.Close()


def test_workspace_content(run, workspace):
    w = workspace
    mh = w.var("MH")
    assert mh and mh.isConstant()
    masses = expected_masses(run["fit_range"])
    assert (mh.getMin(), mh.getMax()) == (min(masses), max(masses))
    for cat in syn.CATEGORIES:
        c = cat_name(cat)
        for name in (f"data_obs_{c}", f"ZprimeQQ_{c}", f"ZprimeQQ_{c}_norm", f"bkg_{c}", f"bkg_{c}_norm",
                     f"pdfindex_{c}", f"mjj_{c}"):
            assert w.obj(name), name
        assert w.pdf(f"bkg_{c}").getNumPdfs() == len(run["bkg"]["categories"][c]["passing"])
    # the functions that are not log-linear reach the workspace, with their constant arguments
    nl = [(c, f) for c, v in run["bkg"]["categories"].items() for f in v["passing"] if not B.is_loglinear(f)]
    assert nl, "no non-linear background function in any RooMultiPdf of the synthetic inputs"
    for c, f in nl:
        assert w.pdf(f"bkg_{c}_{f}") and w.obj(f"{f}_c0_{c}"), (c, f)
    # each parameter is interpolated through the masses where it was floating (if >= 2)
    for c, cd in run["param"]["categories"].items():
        for k, md in cd["models"].items():
            pts = md["points"]
            floating = [m for m, f in zip(pts["masses"], pts["fixed"]) if not f]
            assert pts["used_for_interpolation"] == (floating if len(floating) >= 2 else pts["masses"]), (c, k)
    checks = run["ws"]["checks"]
    assert checks["max_rel_diff_signal_workspace_vs_parametrisation"] < 1e-9
    assert checks["max_abs_diff_bkg_chi2"] < 1e-6


def test_workspace_pdfs_sum_to_one_over_the_bin_centres(ROOT, run, workspace):
    """Every pdf in the workspace (signal at several masses, every background function) is
    normalised by its sum over the bin centres: Combine's per-bin predictions N x pdf(centre) x
    width then add up exactly to the yield. A mismatch of only ~1e-6 (the signal pdf with its
    analytical integral) biased r by 0.2-0.3 sigma at the lowest mass with ~1e7 background events."""
    w = workspace
    mh = w.var("MH")
    for c in (cat_name(cat) for cat in syn.CATEGORIES):
        x = w.var(f"mjj_{c}")
        ns = ROOT.RooArgSet(x)
        edges = run["bkg"]["categories"][c]["bin_edges"]
        centres = [0.5 * (a + b) for a, b in zip(edges[:-1], edges[1:])]
        width = edges[1] - edges[0]

        def bin_sum(pdf):
            tot = 0.0
            for xc in centres:
                x.setVal(xc)
                tot += pdf.getVal(ns) * width
            return tot

        for m in (mh.getMin(), 0.5 * (mh.getMin() + mh.getMax()), mh.getMax()):
            mh.setVal(m)
            assert bin_sum(w.pdf(f"ZprimeQQ_{c}")) == pytest.approx(1.0, abs=1e-12), (c, m)
        for func in run["bkg"]["categories"][c]["passing"]:
            assert bin_sum(w.pdf(f"bkg_{c}_{func}")) == pytest.approx(1.0, abs=1e-12), (c, func)


def test_signal_model_at_a_mass_not_used_in_the_fits(ROOT, run, workspace):
    """The interpolated signal model at M = 350 GeV (not given to the workflow) describes the
    independently generated M = 350 GeV template: shape (goodness of fit) and yield."""
    w = workspace
    xmin, xmax = run["bkg"]["fit_range"]
    fval = ROOT.TFile.Open(run["truth"]["files"]["validation"][0])
    scale = syn.LUMI / syn.N_GENERATED
    masses = expected_masses((xmin, xmax))
    validation = [m for m in syn.VALIDATION_MASSES if masses and masses[0] < m < masses[-1]]
    if not validation:
        pytest.skip("no validation mass inside the interpolation range for this fit range")
    for m in validation:
        w.var("MH").setVal(m)
        for cat in syn.CATEGORIES:
            c = cat_name(cat)
            h = fval.Get(syn.category_names()[cat])
            x = w.var(f"mjj_{c}")
            nset = ROOT.RooArgSet(x)
            edges = run["bkg"]["categories"][c]["bin_edges"]
            lo, hi = np.array(edges[:-1]), np.array(edges[1:])
            hw = window_hist(ROOT, h, xmin, xmax, run["bkg"]["categories"][c]["rebin"], name=f"val_{cat}_{m}")
            counts = np.array([hw.GetBinContent(i) for i in range(1, hw.GetNbinsX() + 1)])
            norm = w.function(f"ZprimeQQ_{c}_norm").getVal() / scale  # in MC counts
            dens = []
            for a, b in zip(lo, hi):
                x.setVal(0.5 * (a + b))
                dens.append(w.pdf(f"ZprimeQQ_{c}").getVal(nset))
            nu = norm * np.array(dens) * (hi - lo)
            chi2 = B.baker_cousins(counts, nu)
            pval = float(ROOT.TMath.Prob(chi2, len(counts)))
            assert pval > 1e-3, (cat, m, chi2, len(counts))
            ty = true_yield_in_range(run, cat, m, xmin, xmax) / scale
            assert abs(norm - ty) < MAX_PULL * np.sqrt(ty) + 0.02 * ty, (cat, m, norm, ty)


def test_every_fit_has_a_plot(run):
    out = run["out"]
    nbkg = sum(len(c["fits"]) for c in run["bkg"]["categories"].values())
    nsig = sum(len(c["fits"]) for c in run["sig"]["categories"].values())
    assert len(list((out / "background" / "plots").glob("bkg_*.png"))) == nbkg
    assert len(list((out / "signal" / "plots").glob("sig_*.png"))) == nsig
    ncat = len(run["param"]["categories"])
    assert len(list((out / "signal_param" / "plots").glob("param_*.png"))) == 7 * ncat
    assert len(list((out / "signal_param" / "plots").glob("model_*.png"))) == ncat


# ---- step 6 and the bias test (CMS Combine) -------------------------------------------------------
def combine_cmd(run, *args):
    return subprocess.run(["combine", *map(str, args)], cwd=run["out"] / "datacards",
                          capture_output=True, text=True)


def best_fit_r(run, mass, mu, rrange, tag, *extra):
    """r of a Combine fit (MultiDimFit, best fit only) to the Asimov data set for mu."""
    p = combine_cmd(run, "-M", "MultiDimFit", "datacard_combined.root", "-m", mass, "-t", -1,
                    "--expectSignal", mu, "--toysNoSystematics", "--rMin", rrange[0],
                    "--rMax", rrange[1],
                    "--algo", "none", *C.DISCRETE_PROFILING_OPTIONS, "-n", tag, *extra)
    assert p.returncode == 0 and "MultiDimFit failed" not in p.stdout, p.stdout[-3000:]
    found = re.search(r"best fit parameter values:\s*\n\s*r :\s*([-+0-9.eE]+)", p.stdout)
    assert found, p.stdout[-2000:]
    return float(found.group(1))


def test_datacards(run):
    """One datacard per category, taking everything from the workspace; the combined card and
    workspace exist (make_datacards.py also checks the expected yields in the latter)."""
    d = run["out"] / "datacards"
    meta = json.loads((d / "datacards.json").read_text())
    assert sorted(meta["categories"]) == sorted(run["ws"]["categories"])
    for c in meta["categories"]:
        card = (d / meta["cards"][c]).read_text()
        assert f"w:bkg_{c}" in card and f"w:ZprimeQQ_{c}" in card and f"w:data_obs_{c}" in card
        assert f"pdfindex_{c}" in card
    combined = (d / "datacard_combined.txt").read_text()
    assert re.search(rf"imax {len(meta['categories'])}\b", combined)
    assert (d / "datacard_combined.root").exists()


def sigma_r_estimate(ROOT, run, mass):
    """Lower bound on the uncertainty of r (background fixed): 1 / sqrt(sum s_i^2 / b_i), with
    s_i the signal (sigma = 1 nb) and b_i the data in each bin of every category."""
    f = ROOT.TFile.Open(str(run["out"] / "workspace" / "dijet_workspace.root"))
    w = f.Get("w")
    w.var("MH").setVal(mass)
    info = 0.0
    for c in run["ws"]["categories"]:
        x, d = w.var(f"mjj_{c}"), w.data(f"data_obs_{c}")
        pdf, norm = w.pdf(f"ZprimeQQ_{c}"), w.function(f"ZprimeQQ_{c}_norm").getVal()
        for i in range(d.numEntries()):
            x.setVal(d.get(i).getRealValue(x.GetName()))
            s = norm * pdf.getVal(ROOT.RooArgSet(x)) * 1.0
            info += s * s / max(d.weight(), 1.0)
    return 1.0 / math.sqrt(info)


def test_combine_asimov_fit_returns_the_injected_signal(ROOT, run):
    """Regression test for the normalisation of the background pdfs: a Combine fit to its own
    Asimov data set, with the generating background function, must return the injected signal
    strength (within 2% of its uncertainty). With the pdfs normalised by their integral instead
    of the sum over the bin centres, the expected counts per bin do not add up to the yield, and
    with ~10^7 events per category the fitted r was off by up to one standard deviation."""
    masses = expected_masses(run["fit_range"])
    meta = json.loads((run["out"] / "datacards" / "datacards.json").read_text())
    idx = [v["index"] for v in meta["categories"].values()]
    fix = ["--setParameters", ",".join(f"{v['index']}={v['default_index']}"
                                       for v in meta["categories"].values()),
           "--freezeParameters", ",".join(idx)]
    for m in (masses[0], masses[-1]):
        sigma = sigma_r_estimate(ROOT, run, m)
        for mu in (0.0, 3.0 * sigma):
            r = best_fit_r(run, m, mu, (mu - 20 * sigma, mu + 20 * sigma), f".asimov{m}_{mu:g}", *fix)
            assert abs(r - mu) < 0.02 * sigma, (m, mu, r, sigma)


def test_combine_mass_is_set_by_m(run):
    """Regression test: -m sets the signal mass (MH must be a constant parameter, otherwise
    Combine fits it and -m has no effect)."""
    masses = expected_masses(run["fit_range"])
    if len(masses) < 2:
        pytest.skip("fewer than two signal masses for this fit range")
    limits = {}
    for m in (masses[0], masses[-1]):
        p = combine_cmd(run, "-M", "AsymptoticLimits", "datacard_combined.root", "-m", m,
                        "--run", "blind", "--cminDefaultMinimizerStrategy", "0", "-n", f".m{m}")
        found = re.search(r"Expected 50\.0%: r < ([0-9.eE+-]+)", p.stdout)
        assert found, p.stdout[-2000:]
        limits[m] = float(found.group(1))
    assert abs(limits[masses[0]] / limits[masses[-1]] - 1) > 0.1, limits


def test_bias_test_script(run):
    """bias_test.py end to end on a small configuration (10 toys in 3 jobs, two generators):
    every toy is fitted and counted in the Combine log, the seeds are unique, the pull toys are
    taken from several jobs, the per-toy function indices are saved, and the fixed-function fits
    close (mean r compatible with the injected mu)."""
    masses = expected_masses(run["fit_range"])
    meta = json.loads((run["out"] / "datacards" / "datacards.json").read_text())
    # a family that is an alternate (not the default) in at least one category
    fams = {re.sub(r"\d+$", "", f) for v in meta["categories"].values()
            for i, f in v["index_to_function"].items() if int(i) != v["default_index"]}
    gen2 = f"family-{sorted(fams)[0]}"
    out = run["base"] / "bias"
    mu = 1.0
    p = subprocess.run([sys.executable, str(REPO / "scripts" / "bias_test.py"),
                        run["out"] / "datacards" / "datacards.json", "-o", out, "--masses",
                        str(masses[-1]), "--mu", f"{mu:g}", "--ntoys", "10", "--toys-per-job", "4",
                        "--generate", "default", gen2, "--pull-toys", "6", "--jobs", "3",
                        "--formats", "png"],
                       capture_output=True, text=True)
    (run["base"] / "bias_test.log").write_text(p.stdout + p.stderr)
    assert p.returncode == 0, p.stderr[-3000:]
    assert "WARNING" not in p.stdout, p.stdout[-3000:]
    res = json.loads((out / "bias_results.json").read_text())
    assert res["asimov"][str(masses[-1])]["sigma"] > 0 and res["failed_jobs"] == []
    assert res["categories"] == list(meta["categories"])
    assert {pt["generate"] for pt in res["points"]} == {"default", gen2}
    for point in res["points"]:
        for mode in ("envelope", "fixed"):
            s = point["summary"][mode]
            assert s["n"] == 10 and s["n_failed"] == 0 and s["n_reported_ok"] == 10
            assert s["n_failed_fit"] == 0 and s["n_failed_unexplained"] == 0
            assert (out / "plots" / f"rfit_{point['generate']}_{mode}_M{masses[-1]}_mu1.png").exists()
        idx = point["toys"]["envelope_index"]
        assert len(idx) == 10 and all(len(i) == len(meta["categories"]) for i in idx)
        assert len(point["toys"]["pull_fits"]) == 6  # 4 from the first job, 2 from the second
        fixed = point["summary"]["fixed"]
        assert abs(fixed["mean"] - mu) < 4 * fixed["std"] / math.sqrt(fixed["n"])
    seeds = [int(re.search(r"\.(\d+)\.root$", f.name).group(1))
             for f in (out / "work").glob("*/job*/higgsCombine.gen.GenerateOnly.*.root")]
    assert len(seeds) == 6 and len(set(seeds)) == 6
    assert (out / "plots" / "bias_summary_default.png").exists()
    # verdict of every point against the criterion, and the list of failing points
    crit = res["bias_criterion"]
    failing = {(f["generate"], f["fit"]) for f in res["failing_points"]}
    for point in res["points"]:
        for mode in ("envelope", "fixed"):
            s = point["summary"][mode]
            assert s["within_criterion"] == (abs(s["bias_over_sigma"]) < crit)
            assert ((point["generate"], mode) in failing) == (not s["within_criterion"])
    assert f"## Verdict (|bias| < {crit} sigma)" in (out / "bias_summary.md").read_text()


def test_bias_test_survives_a_failing_job_and_resumes(run, tmp_path, monkeypatch):
    """A Combine job that fails is recorded (failed_jobs, its toys counted as failed); the other
    jobs are analysed and the outputs are written. --resume then reruns only the failed job."""
    sys.path.insert(0, str(REPO / "scripts"))
    import bias_test as BT
    masses = expected_masses(run["fit_range"])
    original = BT.Runner.job

    def job(self, gen, mass, mu, j, *a):
        if j == 1:
            raise RuntimeError("simulated failure")
        return original(self, gen, mass, mu, j, *a)

    monkeypatch.setattr(BT.Runner, "job", job)
    monkeypatch.setattr(sys, "argv", ["bias_test.py", str(run["out"] / "datacards" / "datacards.json"),
                                      "-o", str(tmp_path), "--masses", str(masses[-1]), "--mu", "0",
                                      "--ntoys", "6", "--toys-per-job", "3", "--fits", "fixed",
                                      "--jobs", "2", "--formats", "png"])
    BT.main()
    res = json.loads((tmp_path / "bias_results.json").read_text())
    assert len(res["failed_jobs"]) == 1 and "simulated failure" in res["failed_jobs"][0]["error"]
    s = res["points"][0]["summary"]["fixed"]
    assert s["n"] == 3 and s["n_failed"] == 3 and s["n_lost_in_failed_jobs"] == 3
    assert s["n_failed_unexplained"] == 0
    assert "1 Combine job(s) failed" in (tmp_path / "bias_summary.md").read_text()

    # --resume with the same arguments: only the failed job runs again, the finished job and the
    # Asimov fit are reused, and the result is complete
    first = res["points"][0]["toys"]["fixed"]  # r of the toys of job 0
    monkeypatch.setattr(BT.Runner, "job", original)
    calls = []
    monkeypatch.setattr(BT.C, "run", lambda cmd, *a, _run=BT.C.run: (calls.append(cmd[2]), _run(cmd, *a)))
    monkeypatch.setattr(sys, "argv", sys.argv + ["--resume"])
    BT.main()
    res2 = json.loads((tmp_path / "bias_results.json").read_text())
    assert res2["failed_jobs"] == [] and res2["resumed"] == {"asimov": 1, "toys": 1}
    assert calls == ["GenerateOnly", "MultiDimFit"]  # one job, one fit mode; no Asimov fit
    s = res2["points"][0]["summary"]["fixed"]
    assert s["n"] == 6 and s["n_failed"] == 0 and s["n_failed_unexplained"] == 0
    assert res2["points"][0]["toys"]["fixed"][:3] == first  # job 0 taken from the first run

    # different settings (another seed): nothing is reused
    calls.clear()
    monkeypatch.setattr(sys, "argv", sys.argv + ["--seed", "777"])
    BT.main()
    assert json.loads((tmp_path / "bias_results.json").read_text())["resumed"] == {"asimov": 1, "toys": 0}
    assert calls.count("GenerateOnly") == 2

    # same arguments, but other Combine options: nothing is reused (the commands are part of the
    # settings), not even the Asimov fit
    calls.clear()
    monkeypatch.setattr(BT.C, "DISCRETE_PROFILING_OPTIONS",
                        BT.C.DISCRETE_PROFILING_OPTIONS + ["--X-rtd", "MINIMIZER_freezeDisassociatedParams"])
    BT.main()
    assert json.loads((tmp_path / "bias_results.json").read_text())["resumed"] == {"asimov": 0, "toys": 0}
    assert calls.count("GenerateOnly") == 2 and calls.count("FitDiagnostics") == 1


def test_bias_test_refuses_duplicate_values(run, tmp_path, monkeypatch):
    """Duplicate values (e.g. --mu 1 1.0) would make two jobs share a directory."""
    sys.path.insert(0, str(REPO / "scripts"))
    import bias_test as BT
    monkeypatch.setattr(sys, "argv", ["bias_test.py", str(run["out"] / "datacards" / "datacards.json"),
                                      "-o", str(tmp_path), "--mu", "1", "1.0"])
    with pytest.raises(SystemExit, match="duplicate values in --mu"):
        BT.main()


def test_make_datacards_detects_a_shape_mismatch(ROOT, run):
    """The shape check of make_datacards.py passes on the combined workspace and fails when a
    background shape parameter or the signal differs from the source workspace."""
    sys.path.insert(0, str(REPO / "scripts"))
    import make_datacards as MD
    summary = run["ws"]
    src_f = ROOT.TFile.Open(str(run["out"] / "workspace" / "dijet_workspace.root"))
    src = src_f.Get("w")
    f = ROOT.TFile.Open(str(run["out"] / "datacards" / "datacard_combined.root"))
    w = f.Get("w")
    c = next(iter(summary["categories"]))
    masses = summary["categories"][c]["signal"]["masses"]
    mh, mh_src = w.var(summary["mass_variable"]), src.var(summary["mass_variable"])
    assert MD.check_shapes(ROOT, w, src, summary, c, mh, mh_src, masses) == []
    default = summary["categories"][c]["background"]["functions"][summary["categories"][c]["background"]["default"]]
    par = next(v for v in w.pdf(default).getVariables()
               if not v.isConstant() and "mjj" not in v.GetName() and v.GetName() != summary["mass_variable"])
    par.setVal(par.getVal() + 0.5)
    problems = MD.check_shapes(ROOT, w, src, summary, c, mh, mh_src, masses)
    assert any("background function" in x for x in problems), problems


# ---- fit-range scan ------------------------------------------------------------------------------
def test_fit_range_scan_cuts_a_turnon(ROOT, fit_range, tmp_path_factory):
    """A background with a turn-on centred just above the configured lower edge: the fits fail
    there, and fit_background.py --scan must raise the lower edge, in steps of 5 GeV and without
    rebinning, to a point where the dijet functions pass again."""
    xmin, xmax = fit_range
    base = tmp_path_factory.mktemp("scan")
    turnon = (xmin + 2.0, 4.0)  # efficiency 0.38 at the configured lower edge
    truth = syn.generate(base / "inputs", seed=54321, turnon=turnon)
    cfg = syn.workflow_config(xmin, xmax, scan_max=xmin + 40.0)
    step = lambda *a: subprocess.run([sys.executable, *map(str, a)], capture_output=True, text=True)  # noqa: E731
    p1 = step(REPO / "scripts" / "prepare_inputs.py", "--data", truth["files"]["data"],
              "--signal", *truth["files"]["signal"], "--categories", cfg["CATEGORIES"],
              "--exclude", cfg["EXCLUDE"], "-o", base / "inputs.root")
    assert p1.returncode == 0, p1.stderr[-2000:]
    p2 = step(REPO / "scripts" / "fit_background.py", base / "inputs.root", "-o", base / "background",
              "--xmin", xmin, "--xmax", xmax, "--scan", "--scan-max", cfg["SCAN_MAX"], "--formats", "png")
    (base / "fit_background.log").write_text(p2.stdout + p2.stderr)
    assert p2.returncode == 0, p2.stdout[-2000:] + p2.stderr[-2000:]
    bkg = json.loads((base / "background" / "background_fits.json").read_text())
    new_xmin = bkg["fit_range"][0]
    print(f"configured lower edge {xmin:g} GeV, turn-on centred at {turnon[0]:g} GeV: "
          f"scan selected {new_xmin:g} GeV")
    assert xmin < new_xmin <= xmin + 40.0
    steps = (new_xmin - xmin) / 5.0
    assert steps == pytest.approx(round(steps), abs=1e-9)
    assert all(c["rebin"] == 1 and c["passing"] for c in bkg["categories"].values())
    # the scan log shows that the fits failed at the configured lower edge
    first = next(line for line in bkg["scan_log"] if "xmin = " in line)
    assert f"xmin = {xmin:g} GeV" in first and "failing" in first, first


# ---- options and error paths, on the closure inputs ------------------------------------------
def script(name, *args):
    return subprocess.run([sys.executable, str(REPO / "scripts" / name), *map(str, args)],
                          capture_output=True, text=True)


def first_category():
    return cat_name(next(iter(syn.CATEGORIES)))


def test_weights_counts_option_agrees_with_sumw2(run, tmp_path):
    """--weights counts (fit the MC counts, scale the yield back) gives the same results as the
    default SumW2 fit of the weighted templates: uniform weights."""
    c = first_category()
    p = script("fit_signal.py", run["out"] / "inputs.root", "--background",
               run["out"] / "background" / "background_fits.json", "-o", tmp_path,
               "--categories", c, "--weights", "counts", "--formats", "png")
    assert p.returncode == 0, p.stderr[-2000:]
    counts = json.loads((tmp_path / "signal_fits.json").read_text())
    assert counts["weights"] == "counts"
    for m, fc in counts["categories"][c]["fits"].items():
        fs = run["sig"]["categories"][c]["fits"][m]
        assert fc["status"] == 0 and fc["cov_qual"] == 3
        for k in ("mean", "sigma", "norm"):
            a, b = fs["params"][k], fc["params"][k]
            assert abs(a["value"] - b["value"]) < 0.2 * a["error"], (m, k, a, b)
            assert b["error"] == pytest.approx(a["error"], rel=0.1), (m, k, a, b)


def test_polynomial_parametrisation_and_no_multipdf(ROOT, run, tmp_path):
    """Polynomial signal models (with pull plots) and a workspace without RooMultiPdf."""
    p = script("parametrize_signal.py", run["out"] / "signal" / "signal_fits.json", "--input",
               run["out"] / "inputs.root", "-o", tmp_path / "param", "--default-model", "pol1",
               "--model", "mean=pol2", "--formats", "png")
    assert p.returncode == 0, p.stderr[-2000:]
    par = json.loads((tmp_path / "param" / "signal_parametrisation.json").read_text())
    for cd in par["categories"].values():
        assert cd["models"]["mean"]["kind"] == "poly" and cd["models"]["mean"]["order"] == 2
        assert cd["models"]["sigma"]["kind"] == "poly" and cd["models"]["sigma"]["order"] == 1
    assert len(list((tmp_path / "param" / "plots").glob("param_*.png"))) == 7 * len(par["categories"])
    ws = tmp_path / "ws" / "w.root"
    p = script("build_workspace.py", "--input", run["out"] / "inputs.root", "--background",
               run["out"] / "background" / "background_fits.json", "--signal",
               tmp_path / "param" / "signal_parametrisation.json", "-o", ws, "--no-multipdf")
    assert p.returncode == 0, p.stdout[-2000:] + p.stderr[-2000:]
    f = ROOT.TFile.Open(str(ws))
    w = f.Get("w")
    classes = {o.ClassName() for o in w.components()}
    assert "RooMultiPdf" not in classes
    for c in par["categories"]:
        cat = w.cat(f"pdfindex_{c}")
        assert cat and cat.size() == len(run["bkg"]["categories"][c]["passing"])
        for fn in run["bkg"]["categories"][c]["passing"]:
            assert w.pdf(f"bkg_{c}_{fn}")
    f.Close()


def test_workspace_refuses_bins_wider_than_1gev(run, tmp_path):
    bkg = json.loads((run["out"] / "background" / "background_fits.json").read_text())
    c = first_category()
    e = bkg["categories"][c]["bin_edges"]
    bkg["categories"][c]["bin_edges"] = e[::2] + ([e[-1]] if (len(e) - 1) % 2 else [])
    bkg["categories"][c]["binning"] = "2 GeV bins (test)"
    (tmp_path / "bkg.json").write_text(json.dumps(bkg))
    p = script("build_workspace.py", "--input", run["out"] / "inputs.root", "--background",
               tmp_path / "bkg.json", "--signal", run["out"] / "signal_param" / "signal_parametrisation.json",
               "-o", tmp_path / "w.root")
    assert p.returncode != 0
    assert "bins wider than 1 GeV" in p.stderr and c in p.stderr and "bin-centre" in p.stderr
    assert not (tmp_path / "w.root").exists()


def test_a_category_without_a_passing_function_fails_the_case(run, tmp_path):
    """With a goodness-of-fit threshold that one category cannot meet, fit_background.py writes
    its outputs and stops with an error naming that category: an analysis case needs a passing
    function in every category (categories are never dropped)."""
    xmin, xmax = run["fit_range"]
    pvals = {c: v["fits"]["dijet4"]["pvalue"] for c, v in run["bkg"]["categories"].items()}
    worst = min(pvals, key=pvals.get)
    others = sorted(v for c, v in pvals.items() if c != worst)
    thr = 0.5 * (pvals[worst] + others[0])
    p = script("fit_background.py", run["out"] / "inputs.root", "--xmin", xmin, "--xmax", xmax,
               "--functions", "dijet4", "--pvalue-threshold", thr, "--formats", "png", "-o", tmp_path / "a")
    assert p.returncode != 0 and "the analysis case fails" in p.stderr and worst in p.stderr
    out = json.loads((tmp_path / "a" / "background_fits.json").read_text())
    assert sorted(out["categories"]) == sorted(pvals) and "dropped_categories" not in out
    assert out["categories"][worst]["passing"] == []


def test_signal_binning_must_match_data(ROOT, run, tmp_path):
    """prepare_inputs.py refuses signal histograms whose binning differs from the data."""
    from array import array
    path = tmp_path / "ZprimeQQ_M300_rebinned.root"
    f = ROOT.TFile(str(path), "RECREATE")
    lo, hi = syn.HIST_RANGE
    for name in syn.category_names().values():
        h = ROOT.TH1D(name, "", int((hi - lo) / 2), array("d", np.arange(lo, hi + 1, 2.0)))
        h.Fill(300.0)
        h.Write()
    ev = ROOT.TH1D("eventsProcessed", "", 1, 0, 1)
    ev.SetBinContent(1, 100.0)
    ev.Write()
    f.Close()
    cfg = syn.workflow_config(*run["fit_range"])
    p = script("prepare_inputs.py", "--data", run["truth"]["files"]["data"], "--signal", path,
               "--categories", cfg["CATEGORIES"], "--exclude", cfg["EXCLUDE"], "-o", tmp_path / "in.root")
    assert p.returncode != 0 and "different binning" in p.stderr, p.stderr[-2000:]


def test_masses_near_the_upper_edge_are_dropped(run, tmp_path):
    """A signal mass closer than --min-distance to the upper edge of the fit range is dropped."""
    bkg = json.loads((run["out"] / "background" / "background_fits.json").read_text())
    xmin = bkg["fit_range"][0]
    top = max(syn.FIT_MASSES)
    bkg["fit_range"] = [xmin, top + 5.0]  # 5 GeV above the highest mass < default distance 10
    if top + 5.0 > syn.HIST_RANGE[1] or xmin + 10.0 > min(syn.FIT_MASSES):
        pytest.skip("fit range not suitable for this check")
    (tmp_path / "bkg.json").write_text(json.dumps(bkg))
    c = first_category()
    p = script("fit_signal.py", run["out"] / "inputs.root", "--background", tmp_path / "bkg.json",
               "-o", tmp_path / "sig", "--categories", c, "--formats", "png")
    assert p.returncode == 0, p.stderr[-2000:]
    out = json.loads((tmp_path / "sig" / "signal_fits.json").read_text())
    assert top in out["categories"][c]["dropped"] and top not in out["categories"][c]["masses"]


def test_fit_signal_stops_on_a_failed_fit(run, tmp_path, monkeypatch):
    """A signal fit that fails (here: one mass made to fail the quality check) makes
    fit_signal.py write its outputs, with the fit flagged, and stop: the analysis case fails."""
    sys.path.insert(0, str(REPO / "scripts"))
    import fit_signal as FS
    c = first_category()
    masses = run["sig"]["categories"][c]["masses"] if "sig" in run else \
        json.loads((run["out"] / "signal" / "signal_fits.json").read_text())["categories"][c]["masses"]
    bad = masses[0]
    real = FS.S.fit_ok
    monkeypatch.setattr(FS.S, "fit_ok", lambda fit: real(fit) and fit["mass"] != bad)
    monkeypatch.setattr(sys, "argv", ["fit_signal.py", str(run["out"] / "inputs.root"), "--background",
                                      str(run["out"] / "background" / "background_fits.json"),
                                      "-o", str(tmp_path), "--categories", c, "--formats", "png"])
    with pytest.raises(SystemExit, match="the analysis case fails"):
        FS.main()
    out = json.loads((tmp_path / "signal_fits.json").read_text())
    fits = out["categories"][c]["fits"]
    assert fits[str(bad)]["usable"] is False
    assert all(f["usable"] for m, f in fits.items() if m != str(bad))


def test_a_failed_signal_fit_fails_the_case(run, tmp_path):
    """A signal fit flagged as failed stops the parametrisation: every signal fit must pass in
    every category of an analysis case."""
    sig = json.loads((run["out"] / "signal" / "signal_fits.json").read_text())
    c = first_category()
    masses = sig["categories"][c]["masses"]
    bad = masses[len(masses) // 2]
    assert all(f["usable"] for f in sig["categories"][c]["fits"].values())
    sig["categories"][c]["fits"][str(bad)]["usable"] = False
    (tmp_path / "sig.json").write_text(json.dumps(sig))
    p = script("parametrize_signal.py", tmp_path / "sig.json", "--input", run["out"] / "inputs.root",
               "-o", tmp_path / "param", "--formats", "png")
    assert p.returncode != 0
    assert "the analysis case fails" in p.stderr and c in p.stderr and str(bad) in p.stderr
