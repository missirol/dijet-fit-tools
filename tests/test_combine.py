"""Unit tests for dijetfit.combine: datacard text and bias-test statistics."""
import math

import pytest

from dijetfit import combine as C


def test_datacard_text():
    card = C.datacard_text("CAT_A", "../ws/w.root")
    lines = [ln.split() for ln in card.splitlines()]
    assert ["imax", "1"] == lines[1][:2] and ["jmax", "1"] == lines[2][:2]
    shapes = [ln for ln in lines if ln and ln[0] == "shapes"]
    assert shapes == [["shapes", "ZprimeQQ", "CAT_A", "../ws/w.root", "w:ZprimeQQ_CAT_A"],
                      ["shapes", "bkg", "CAT_A", "../ws/w.root", "w:bkg_CAT_A"],
                      ["shapes", "data_obs", "CAT_A", "../ws/w.root", "w:data_obs_CAT_A"]]
    rows = {ln[0]: ln[1:] for ln in lines if ln and ln[0] in ("observation", "rate")}
    assert rows == {"observation": ["-1"], "rate": ["1", "1"]}
    proc = [ln[1:] for ln in lines if ln and ln[0] == "process"]
    assert proc == [["ZprimeQQ", "bkg"], ["0", "1"]]  # signal = process 0 (scaled by r)
    assert ["pdfindex_CAT_A", "discrete"] in lines
    assert not any(ln and ln[-1] in ("lnN", "shape", "param") for ln in lines)  # no systematics


def test_pull_uses_the_error_on_the_side_of_the_truth():
    assert C.pull(1.2, 1.0, err_lo=0.1, err_hi=0.4) == pytest.approx(2.0)   # above: lower error
    assert C.pull(0.6, 1.0, err_lo=0.1, err_hi=0.4) == pytest.approx(-1.0)  # below: upper error
    assert C.pull(0.6, 1.0, err_lo=0.1, err_hi=0.0) is None
    assert C.pull(0.6, 1.0, err_lo=0.1, err_hi=float("nan")) is None


def test_mean_and_error_and_median():
    s = C.mean_and_error([1.0, 2.0, 3.0, 4.0])
    assert s["n"] == 4 and s["mean"] == pytest.approx(2.5)
    assert s["std"] == pytest.approx(math.sqrt(5 / 3))
    assert s["mean_err"] == pytest.approx(s["std"] / 2)
    assert C.median([3.0, 1.0, 2.0]) == 2.0 and C.median([4.0, 1.0, 2.0, 3.0]) == 2.5
    assert math.isnan(C.median([]))


def test_discrete_profiling_options_are_pairs():
    o = C.DISCRETE_PROFILING_OPTIONS
    assert o[0] == "--cminDefaultMinimizerStrategy"
    assert all(o[i] == "--X-rtd" for i in range(2, len(o), 2))


def test_bias_test_generators():
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    import bias_test as BT
    meta = dict(categories=dict(
        a=dict(default_index=0, index_to_function={"0": "dijet6", "1": "ua26", "2": "ua25"}),
        b=dict(default_index=0, index_to_function={"0": "moddijet5", "1": "dijet5", "2": "ua24"}),
        c=dict(default_index=0, index_to_function={"0": "dijet5", "1": "moddijet5"})))
    assert BT.gen_indices(meta, "default") == dict(a=0, b=0, c=0)
    assert BT.gen_indices(meta, "alt2") == dict(a=2, b=2, c=0)
    # highest order of the family in each envelope, the default if the family is missing
    assert BT.gen_indices(meta, "family-ua2") == dict(a=1, b=2, c=0)
    assert BT.gen_indices(meta, "family-dijet") == dict(a=0, b=1, c=0)
    assert BT.gen_indices(meta, "family-moddijet") == dict(a=0, b=0, c=1)
