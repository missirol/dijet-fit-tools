"""Unit tests for dijetfit.interp: splines, polynomials and their RooFit formulas."""
import numpy as np
import pytest

from dijetfit import interp

XS = [200.0, 300.0, 400.0, 500.0]
YS = [1.0, 3.0, 2.0, 5.0]


def test_spline_interpolates_and_is_natural():
    sp = interp.natural_cubic_spline(XS, YS)
    for x, y in zip(XS, YS):
        assert interp.evaluate(sp, x) == pytest.approx(y)
    # natural boundary: second derivative 0 at both ends (numerically)
    h = 1e-3
    for x0, side in ((XS[0], +1), (XS[-1], -1)):
        xa, xb, xc = x0, x0 + side * h, x0 + 2 * side * h
        d2 = (interp.evaluate(sp, xa) - 2 * interp.evaluate(sp, xb) + interp.evaluate(sp, xc)) / h ** 2
        assert abs(d2) < 1e-3
    # continuity of value and slope at the inner knots
    for k in XS[1:-1]:
        assert interp.evaluate(sp, k - 1e-7) == pytest.approx(interp.evaluate(sp, k + 1e-7), abs=1e-5)


def test_spline_of_linear_data_is_linear_also_outside():
    ys = [2.0 + 0.01 * x for x in XS]
    sp = interp.natural_cubic_spline(XS, ys)
    for x in (150.0, 250.0, 333.3, 480.0, 550.0):
        assert interp.evaluate(sp, x) == pytest.approx(2.0 + 0.01 * x)


@pytest.mark.parametrize("model", ["spline", "pol1", "pol2"])
def test_formula_matches_evaluate(ROOT, model):
    m = interp.build(model, XS, YS, [0.1] * 4)
    f = ROOT.TFormula(f"f_{model}", interp.formula(m, "x"))
    for x in (150.0, 200.0, 260.0, 300.0, 410.0, 500.0, 560.0):
        assert f.Eval(x) == pytest.approx(interp.evaluate(m, x), rel=1e-12, abs=1e-12)


def test_clamped_formula(ROOT):
    m = interp.build("spline", XS, YS, [0.1] * 4)
    f = ROOT.TFormula("fc", interp.clamped_formula(m, 1.5, 4.0, "x"))
    assert f.Eval(200.0) == pytest.approx(1.5)
    assert f.Eval(500.0) == pytest.approx(4.0)
    assert f.Eval(300.0) == pytest.approx(3.0)


def test_poly_fit_recovers_polynomial_and_handles_zero_errors():
    xs = np.array(XS)
    ys = 1.0 - 0.02 * xs + 3e-5 * xs ** 2
    m = interp.poly_fit(xs, ys, [0.1] * 4, 2)
    assert m["coefficients"] == pytest.approx([1.0, -0.02, 3e-5], rel=1e-6)
    assert m["weighted"] and m["chi2"] == pytest.approx(0.0, abs=1e-12)
    # points without uncertainty (parameters fixed in the fit, stored with error 0; they enter a
    # polynomial fit only when fewer than two points float) -> all points get equal weights
    m0 = interp.poly_fit(xs, [1.0, 1.0, 1.0, 5.0], [0.1, 0.0, 0.1, 0.1], 1)
    assert not m0["weighted"]
    assert m0["coefficients"] == pytest.approx(np.polyfit(xs, [1.0, 1.0, 1.0, 5.0], 1)[::-1])


def test_parse_model():
    assert interp.parse_model("spline") == ("spline", None)
    assert interp.parse_model("POL3") == ("poly", 3)
    with pytest.raises(ValueError):
        interp.parse_model("cubic")


def test_hold_outside_keeps_the_model_constant_beyond_its_points(ROOT):
    """A parameter determined only at some masses is held constant beyond them (no
    extrapolation), identically in python and in the RooFit formula."""
    m = interp.build("spline", XS[:3], YS[:3], [0.1] * 3)
    m["hold_outside"] = [XS[0], XS[2]]
    f = ROOT.TFormula("fh", interp.formula(m, "x"))
    assert interp.evaluate(m, 500.0) == pytest.approx(YS[2])
    assert interp.evaluate(m, 150.0) == pytest.approx(YS[0])
    for x in (150.0, 250.0, 380.0, 450.0, 500.0):
        assert f.Eval(x) == pytest.approx(interp.evaluate(m, x), rel=1e-12)
