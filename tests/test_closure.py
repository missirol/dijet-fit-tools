"""End-to-end closure tests on synthetic inputs with a known truth (tests/synthetic.py).

* Main test: the complete workflow (run_all.sh) runs with the configured fit range (pytest
  --fit-xmin/--fit-xmax, default 150-600 GeV) on a background without turn-on; the range must
  stay as configured, and every step must give back the truth.
* Scan test: a background with a turn-on just above the configured lower edge; the background
  step with --scan must raise the lower edge past the turn-on.

The workflow runs once per test module (a few tens of seconds); select or skip these tests with
`pytest -m closure` / `pytest -m "not closure"`. Pass `--basetemp DIR` to keep the outputs.
"""
import json
import os
import re
import shutil
import subprocess
import sys

import numpy as np
import pytest

import synthetic as syn
from conftest import REPO
from dijetfit import background as B

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


def expected_masses(xmin, min_distance=10.0):
    """Signal masses the workflow should use: those at least min_distance above the lower edge
    (fit_signal.py --min-distance default)."""
    return [m for m in syn.FIT_MASSES if m - xmin >= min_distance]


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
            assert h.GetBinError(300) ** 2 == pytest.approx(scale * h.GetBinContent(300), rel=1e-9)


# ---- step 2 -------------------------------------------------------------------------------------
def test_fit_range_is_the_configured_one(run):
    """Without a turn-on, the fits pass at the configured lower edge: the fit range must be
    exactly the configured one, without rebinning."""
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
    pdf = B.make_pdf(ROOT, f"pdf_{tag}", fit["name"], al, fit["formula"])
    lo, hi = np.array(edges[:-1]), np.array(edges[1:])
    return B.expected_counts(ROOT, x, pdf, fit["norm"]["value"], lo, hi, B.needs_bin_integration(lo, hi))


def test_background_closure(ROOT, run):
    """The fit of the generating function (dijet4) reproduces the true expected spectrum in the
    fit range within the statistical precision, in every category; and the selection never
    rejects it for its fit quality or goodness of fit. (It may be excluded by the F-test in
    favour of a simpler member of its family that describes the data equally well, e.g. over a
    short fit range; at least one function must pass.)

    Other functions are only required to describe the data, which is what the selection tests:
    a passing function may still differ from the truth by more than the statistics. That is the
    modelling uncertainty the RooMultiPdf / discrete profiling is meant to cover. The distances
    from the truth of all passing functions are printed for information."""
    xmin, xmax = run["bkg"]["fit_range"]
    e = syn.edges()
    sel = (e[:-1] >= xmin - 1e-9) & (e[1:] <= xmax + 1e-9)
    for cat in syn.CATEGORIES:
        cb = run["bkg"]["categories"][cat_name(cat)]
        assert all(f["n_errors"] == 0 and f["status"] == 0 for f in cb["fits"].values())
        assert cb["passing"], cat
        gen = cb["fits"]["dijet4"]
        assert gen["verdict"] == "pass" or gen["verdict"].startswith("fail: F-test"), (cat, gen["verdict"])
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
        assert cs["masses"] == expected_masses(xmin)
        edges = run["bkg"]["categories"][cat_name(cat)]["bin_edges"]
        lo, hi = np.array(edges[:-1]), np.array(edges[1:])
        for m in cs["masses"]:
            fit = cs["fits"][str(m)]
            assert fit["status"] == 0 and fit["cov_qual"] == 3 and fit["n_errors"] == 0
            h = f.Get(f"{cat_name(cat)}/ZprimeQQ_{m}")
            n = np.array([h.GetBinContent(h.FindFixBin(0.5 * (a + b))) for a, b in zip(lo, hi)]) / scale
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
    masses = expected_masses(run["fit_range"][0])
    assert (mh.getMin(), mh.getMax()) == (min(masses), max(masses))
    for cat in syn.CATEGORIES:
        c = cat_name(cat)
        for name in (f"data_obs_{c}", f"ZprimeQQ_{c}", f"ZprimeQQ_{c}_norm", f"bkg_{c}", f"bkg_{c}_norm",
                     f"pdfindex_{c}", f"mjj_{c}"):
            assert w.obj(name), name
        assert w.pdf(f"bkg_{c}").getNumPdfs() == len(run["bkg"]["categories"][c]["passing"])
    checks = run["ws"]["checks"]
    assert checks["max_rel_diff_signal_workspace_vs_parametrisation"] < 1e-9
    assert checks["max_abs_diff_bkg_chi2"] < 1e-6


def test_signal_model_at_a_mass_not_used_in_the_fits(ROOT, run, workspace):
    """The interpolated signal model at M = 350 GeV (not given to the workflow) describes the
    independently generated M = 350 GeV template: shape (goodness of fit) and yield."""
    w = workspace
    xmin, xmax = run["bkg"]["fit_range"]
    fval = ROOT.TFile.Open(run["truth"]["files"]["validation"][0])
    scale = syn.LUMI / syn.N_GENERATED
    masses = expected_masses(xmin)
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
            counts = np.array([h.GetBinContent(h.FindFixBin(0.5 * (a + b))) for a, b in zip(lo, hi)])
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


@pytest.mark.skipif(shutil.which("text2workspace.py") is None, reason="CMS Combine not available")
def test_combine_reads_workspace_and_mass_is_set_by_m(run, tmp_path):
    """Regression test: Combine reads the workspace, and -m sets the signal mass (MH must be a
    constant parameter, otherwise Combine fits it and -m has no effect)."""
    c = cat_name(next(iter(syn.CATEGORIES)))
    ws = run["out"] / "workspace" / "dijet_workspace.root"
    (tmp_path / "card.txt").write_text(f"""imax 1
jmax 1
kmax *
---------------
shapes ZprimeQQ  {c}  {ws} w:ZprimeQQ_{c}
shapes bkg       {c}  {ws} w:bkg_{c}
shapes data_obs  {c}  {ws} w:data_obs_{c}
---------------
bin          {c}
observation  -1
---------------
bin          {c}   {c}
process      ZprimeQQ  bkg
process      0    1
rate         1    1
---------------
pdfindex_{c} discrete
""")
    sh = lambda *a: subprocess.run(a, cwd=tmp_path, capture_output=True, text=True)  # noqa: E731
    assert sh("text2workspace.py", "card.txt", "-o", "card.root").returncode == 0
    masses = expected_masses(run["fit_range"][0])
    if len(masses) < 2:
        pytest.skip("fewer than two signal masses for this fit range")
    limits = {}
    for m in (masses[0], masses[-1]):
        p = sh("combine", "-M", "AsymptoticLimits", "card.root", "-m", str(m), "--run", "blind",
               "--cminDefaultMinimizerStrategy", "0", "-n", f".m{m}")
        found = re.search(r"Expected 50\.0%: r < ([0-9.eE+-]+)", p.stdout)
        assert found, p.stdout[-2000:]
        limits[m] = float(found.group(1))
    assert abs(limits[masses[0]] / limits[masses[-1]] - 1) > 0.1, limits


# ---- fit-range scan ------------------------------------------------------------------------------
def test_fit_range_scan_cuts_a_turnon(ROOT, fit_range, tmp_path_factory):
    """A background with a turn-on centred just above the configured lower edge: the fits fail
    there, and fit_background.py --scan must raise the lower edge, in steps of 5 GeV and without
    rebinning, to a point where the dijet functions pass again."""
    xmin, xmax = fit_range
    base = tmp_path_factory.mktemp("scan")
    turnon = (xmin + 2.0, 4.0)  # efficiency 0.38 at the configured lower edge
    truth = syn.generate(base / "inputs", seed=54321, turnon=turnon)
    cfg = syn.workflow_config(xmin, xmax, scan_max_unbinned=xmin + 40.0, scan_max=xmin + 60.0)
    step = lambda *a: subprocess.run([sys.executable, *map(str, a)], capture_output=True, text=True)  # noqa: E731
    p1 = step(REPO / "scripts" / "prepare_inputs.py", "--data", truth["files"]["data"],
              "--signal", *truth["files"]["signal"], "--categories", cfg["CATEGORIES"],
              "--exclude", cfg["EXCLUDE"], "-o", base / "inputs.root")
    assert p1.returncode == 0, p1.stderr[-2000:]
    p2 = step(REPO / "scripts" / "fit_background.py", base / "inputs.root", "-o", base / "background",
              "--xmin", xmin, "--xmax", xmax, "--scan", "--scan-max-unbinned", cfg["SCAN_MAX_UNBINNED"],
              "--scan-max", cfg["SCAN_MAX"], "--formats", "png")
    (base / "fit_background.log").write_text(p2.stdout + p2.stderr)
    assert p2.returncode == 0, p2.stdout[-2000:] + p2.stderr[-2000:]
    bkg = json.loads((base / "background" / "background_fits.json").read_text())
    new_xmin = bkg["fit_range"][0]
    print(f"configured lower edge {xmin:g} GeV, turn-on centred at {turnon[0]:g} GeV: "
          f"scan selected {new_xmin:g} GeV")
    assert xmin < new_xmin <= xmin + 40.0
    assert (new_xmin - xmin) % 5.0 == pytest.approx(0.0)
    assert all(c["rebin"] == 1 and c["passing"] for c in bkg["categories"].values())
    assert any("xmin = " in line for line in bkg["scan_log"])
