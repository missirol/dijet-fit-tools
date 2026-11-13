"""The category study (studies/category_optimisation) on small synthetic histograms: it runs, the
schemes never exceed the upper bound, the ranked scheme improves with more groups, and the
mass-window correction is measured from the |dy| slices."""
import json
import subprocess
import sys

import numpy as np
import pytest

from conftest import REPO

pytestmark = pytest.mark.slow


def _write(ROOT, path, hist2d, slices, events):
    f = ROOT.TFile(str(path), "RECREATE")
    h = ROOT.TH2D("dydphi", "", 30, 0, 3.0, 20, 0, np.pi)
    for i in range(30):
        for j in range(20):
            h.SetBinContent(i + 1, j + 1, hist2d[i, j])
    h.Write()
    for name, (centre, content) in slices.items():
        hs = ROOT.TH1D(name, "", 100, 150, 650)
        for k in range(100):
            hs.SetBinContent(k + 1, content[k])
        hs.Write()
    ev = ROOT.TH1D("eventsProcessed", "", 1, 0, 1)
    ev.SetBinContent(1, events)
    ev.Write()
    f.Close()


def test_category_study_runs_on_synthetic_histograms(ROOT, tmp_path):
    dy = (np.arange(30) + 0.5) * 0.1
    dphi = (np.arange(20) + 0.5) * np.pi / 20
    DY, DPHI = np.meshgrid(dy, dphi, indexing="ij")
    bkg = 1e6 * np.exp(0.4 * DY) * (1 + np.exp(-(np.pi - DPHI)))      # rises with |dy|
    m = np.arange(150, 650, 5.0) + 2.5
    slice_names = {f"mjj_dY{a:.1f}to{a + 0.5:.1f}".replace(".", "p"): a + 0.25 for a in (0.0, 0.5, 1.0)}
    data_slices = {n: (c, 1e5 * np.exp(-m / (60 + 10 * c))) for n, c in slice_names.items()}
    _write(ROOT, tmp_path / "data.root", bkg, data_slices, 1.0)
    for mass in (200, 400):
        sig = 1e3 * np.exp(-DY / 1.0) * np.exp(-(np.pi - DPHI) / 0.3)    # small |dy|, back-to-back
        sl = {n: (c, 50 * np.exp(-0.5 * ((m - mass) / (0.1 * mass)) ** 2)) for n, c in slice_names.items()}
        _write(ROOT, tmp_path / f"sig_M{mass}_x.root", sig, sl, 1e5)
    out = tmp_path / "out"
    p = subprocess.run([sys.executable, str(REPO / "studies" / "category_optimisation" / "optimise_categories.py"),
                        "--data", tmp_path / "data.root", "--signal", tmp_path / "sig_M200_x.root",
                        tmp_path / "sig_M400_x.root", "--hist", "dydphi", "--ranked", "1", "3", "5",
                        "--slice-hists", r"mjj_dY(?P<lo>\d+p\d+)to(?P<hi>\d+p\d+)", "-o", out,
                        "--formats", "png"], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-3000:]
    res = json.loads((out / "results.json").read_text())
    assert res["masses"] == [200, 400] and set(res["window_correction"]["slopes"]) == {"200", "400"}
    for kind in ("raw", "corrected"):
        rows = {name: r[kind] for name, r in res["schemes"].items()}
        for name, r in rows.items():
            assert all(f <= 1 + 1e-9 for f in r["frac"]), (kind, name)  # never above the upper bound
        cur = next(n for n in rows if n.startswith("A "))
        assert rows[cur]["rel"] == pytest.approx([1.0, 1.0])
        ranked = [rows[f"R ranked s/b, {n} groups"]["frac"] for n in (1, 3, 5)]
        assert all(np.all(np.array(b) >= np.array(a) - 1e-9) for a, b in zip(ranked, ranked[1:]))
    assert (out / "summary.md").exists() and (out / "signal_regions.png").exists()
