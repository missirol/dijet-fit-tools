"""The bin-normalisation demo (demo/bin_normalisation) runs and shows what it should: on exact
Asimov data, the model normalised by the bin sum (NEW) returns the injected signal strength like
exact bin integration, while the integral normalisation (OLD; and BKG_ONLY for the signal at
200 GeV) is biased."""
import json
import subprocess
import sys

from conftest import REPO


def test_bin_normalisation_demo(tmp_path):
    p = subprocess.run([sys.executable, str(REPO / "demo" / "bin_normalisation" / "run_demo.py"), "all",
                        "-o", str(tmp_path), "--masses", "200", "--mu", "0"],
                       capture_output=True, text=True)
    (tmp_path / "demo.log").write_text(p.stdout + p.stderr)
    assert p.returncode == 0, p.stderr[-3000:]
    fit = json.loads((tmp_path / "results.json").read_text())["fits"]["200_0"]["asimov"]
    bias = {m: (fit[m]["r"] - 0.0) / fit[m]["err"] for m in ("OLD", "BKG_ONLY", "NEW", "EXACT")}
    assert abs(bias["EXACT"]) < 0.005 and abs(bias["NEW"]) < 0.02, bias
    assert abs(bias["OLD"]) > 1.0 and abs(bias["BKG_ONLY"]) > 0.1, bias
    assert abs(sum(fit["NEW"]["mismatch_bkg"])) < 1e-3 and abs(sum(fit["OLD"]["mismatch_bkg"])) > 100
    for name in ("bias_vs_mass.png", "prediction_differences_M200_mu0.png", "fit_M200_mu0_NEW.png"):
        assert (tmp_path / "plots" / name).exists(), name
