"""Interpolation of the signal parameters as a function of the resonance mass.

Two models are available for each parameter:
  spline : natural cubic spline through the fitted values (default). It interpolates exactly;
           the point uncertainties are not used. Outside the knots it is continued linearly.
  polN   : polynomial of order N, fitted to the points (weighted by their uncertainties).
A model may carry 'hold_outside': [lo, hi]: outside this interval it is held constant at its
value at the nearest end (used when a parameter is determined only at some mass points).
The result can be evaluated in python and turned into a RooFit formula of the mass variable.
"""
import bisect

import numpy as np


def natural_cubic_spline(xs, ys):
    """Coefficients of the natural cubic spline (second derivative 0 at both ends).
    On [xs[i], xs[i+1]) the spline is a + b*dx + c*dx^2 + d*dx^3 with dx = x - xs[i]."""
    xs = [float(v) for v in xs]
    ys = [float(v) for v in ys]
    n = len(xs)
    if n < 2:
        raise ValueError("a spline needs at least 2 points")
    h = [xs[i + 1] - xs[i] for i in range(n - 1)]
    if any(v <= 0 for v in h):
        raise ValueError("spline knots must be strictly increasing")
    M = [0.0] * n  # second derivatives
    m = n - 2
    if m > 0:
        A = np.zeros((m, m))
        r = np.zeros(m)
        for k in range(m):
            A[k, k] = 2 * (h[k] + h[k + 1])
            if k > 0:
                A[k, k - 1] = h[k]
            if k < m - 1:
                A[k, k + 1] = h[k + 1]
            r[k] = 6 * ((ys[k + 2] - ys[k + 1]) / h[k + 1] - (ys[k + 1] - ys[k]) / h[k])
        M[1:n - 1] = list(np.linalg.solve(A, r))
    segs = []
    for i in range(n - 1):
        a = ys[i]
        b = (ys[i + 1] - ys[i]) / h[i] - h[i] * (2 * M[i] + M[i + 1]) / 6
        c = M[i] / 2
        d = (M[i + 1] - M[i]) / (6 * h[i])
        segs.append((a, b, c, d))
    a, b, c, d = segs[-1]
    slope_right = b + 2 * c * h[-1] + 3 * d * h[-1] ** 2
    return dict(kind="spline", knots=xs, values=ys, segments=segs,
                slope_left=segs[0][1], slope_right=slope_right)


def poly_fit(xs, ys, errs, order):
    """Least-squares polynomial fit, weighted by 1/error. If some points have no uncertainty
    (parameters fixed in the shape fits), all points get equal weights (the chi2 is then not
    a goodness of fit)."""
    xs, ys, errs = map(lambda v: np.asarray(v, float), (xs, ys, errs))
    order = min(order, len(xs) - 1)
    weighted = bool(np.all(errs > 0))
    w = 1.0 / errs if weighted else np.ones_like(ys)
    coef = np.polyfit(xs, ys, order, w=w)[::-1]  # c0 + c1 x + ...
    chi2 = float(np.sum(((ys - np.polyval(coef[::-1], xs)) * w) ** 2))
    ndf = len(xs) - (order + 1)
    return dict(kind="poly", order=int(order), coefficients=[float(c) for c in coef],
                weighted=weighted, chi2=chi2, ndf=int(ndf))


def evaluate(model, x):
    if model.get("hold_outside"):  # constant beyond the knots/points used (no extrapolation)
        lo, hi = model["hold_outside"]
        x = min(max(x, lo), hi)
    if model["kind"] == "poly":
        return float(sum(c * x ** k for k, c in enumerate(model["coefficients"])))
    xs = model["knots"]
    if x <= xs[0]:
        return model["values"][0] + model["slope_left"] * (x - xs[0])
    if x >= xs[-1]:
        return model["values"][-1] + model["slope_right"] * (x - xs[-1])
    i = min(bisect.bisect_right(xs, x) - 1, len(xs) - 2)
    a, b, c, d = model["segments"][i]
    dx = x - xs[i]
    return a + dx * (b + dx * (c + dx * d))


def formula(model, var="@0"):
    """Closed-form expression of the model in the variable `var` (for a RooFormulaVar)."""
    def num(v):
        return f"({float(v)!r})"

    if model.get("hold_outside"):
        lo, hi = model["hold_outside"]
        var = f"TMath::Min(TMath::Max({var},{num(lo)}),{num(hi)})"

    if model["kind"] == "poly":
        return "+".join(f"{num(c)}*pow({var},{k})" if k else num(c)
                        for k, c in enumerate(model["coefficients"]))
    xs, ys, segs = model["knots"], model["values"], model["segments"]

    def seg(i):
        a, b, c, d = segs[i]
        dx = f"({var}-{num(xs[i])})"
        return f"{num(a)}+{dx}*({num(b)}+{dx}*({num(c)}+{dx}*{num(d)}))"

    expr = f"({num(ys[-1])}+{num(model['slope_right'])}*({var}-{num(xs[-1])}))"
    for i in range(len(segs) - 1, -1, -1):
        expr = f"(({var}<{num(xs[i + 1])})?({seg(i)}):{expr})"
    left = f"{num(ys[0])}+{num(model['slope_left'])}*({var}-{num(xs[0])})"
    return f"(({var}<{num(xs[0])})?({left}):{expr})"


def clamped_formula(model, lo=None, hi=None, var="@0"):
    """formula() limited to [lo, hi] (either may be None), with TMath::Max/Min."""
    f = formula(model, var)
    if lo is not None:
        f = f"TMath::Max({float(lo)!r},{f})"
    if hi is not None:
        f = f"TMath::Min({float(hi)!r},{f})"
    return f


def parse_model(text):
    """'spline' or 'polN' -> ('spline', None) / ('poly', N)."""
    t = text.strip().lower()
    if t == "spline":
        return ("spline", None)
    if t.startswith("pol") and t[3:].isdigit():
        return ("poly", int(t[3:]))
    raise ValueError(f"unknown model '{text}' (use 'spline' or 'polN')")


def build(text, xs, ys, errs):
    kind, order = parse_model(text)
    if kind == "spline":
        return natural_cubic_spline(xs, ys)
    return poly_fit(xs, ys, errs, order)
