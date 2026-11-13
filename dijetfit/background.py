"""Background functions for the dijet m_jj spectrum and their fits.

The families are those of the analysis note's background table (x = m / sqrt(s), p0 = the
normalisation, a separate yield parameter in the fits):

  dijet    (3-6 par.)  p0 (1-x)^p1 / x^(p2 + p3 ln x + p4 ln^2 x + p5 ln^3 x)
  moddijet (3-6 par.)  p0 (1-x^(1/3))^p1 / x^(p2 + p3 ln x + p4 ln^2 x + p5 ln^3 x)
  ua2      (3-6 par.)  p0 x^(-p1) exp(-(p2 x + p3 x^2 + p4 x^3 + p5 x^4))
  modexp   (3-5 par.)  p0 exp(p1 x^p2),  p0 exp(p1 x^p2 + p1 (1-x)^p3),
                       p0 exp(p1 x^p2 + p3 (1-x)^p4)
  polypow  (3-5 par.)  p0 (1 + p1 x + p2 x^2 + p3 x^3)^(-p_last)
  polyext  (5-6 par.)  p0 (1-x)^p1 (1 + p4 x + p5 x^2) x^-(p2 + p3 ln x)

<family><N> denotes N parameters *including* the normalisation (the usual convention in CMS
papers), e.g. dijet4 = the 4-parameter dijet function (p1, p2, p3 + yield). Not in the note's
table: ua23 (x^-p1 exp(-p2 x)); the note's ua2 family starts at 4 parameters.

dijet, moddijet and ua2 are log-linear in their parameters, f = exp(sum_k p_k T_k(m)). For
numerical stability the polynomial in ln x (resp. in x) is written in the centred variable
L = ln(m / M0) (resp. z = m / M0 - 1), with the reference mass M0 = the geometric mean of the
fit range. A polynomial of the same degree in a shifted variable spans exactly the same set of
functions (the constant term is absorbed in the normalisation), so the families are identical to
the published ones; only the meaning of the individual parameters differs. Their starting values
come from a Poisson generalised linear model (exact maximum likelihood).

modexp, polypow and polyext are not log-linear; they are written exactly as in the note. Their
starting values come from a grid over the parameters that enter non-linearly, with the others
solved for exactly at each grid point, followed by a full Poisson maximum-likelihood fit (numpy).

All functions are fitted in a decorrelated parametrisation: p = p_start + M q (log-linear:
p = M q, terms centred), with M such that the q are uncorrelated with unit variance at the start.
The per-fit constants (M, the centres, M0) are passed to the RooFit formula as constant
arguments, so every fit of a function shares one compiled formula (ROOT compiles each distinct
formula string once; thousands of distinct ones exhaust the process memory maps).
"""
import math

import numpy as np

from .common import SQRTS, capture_output, count_messages, message_summary

INTEGRATOR_EPS_REL = 1e-10
BIN_INTEGRATION_ABOVE = 1.0 + 1e-6  # GeV: bins wider than this are integrated over
BIN_INTEGRATION_PRECISION = 1e-10
QRANGE = 50.0  # half-range of the decorrelated parameters (their errors are ~1)
MAX_RECENTRE = 3  # re-centring attempts for a parameter that ends at its limit
NORM_RANGE = (0.5, 1.5)  # allowed range of the yield, in units of the events in the fit range

_X = f"(@0/{SQRTS!r})"


_TERM_EXPR = {  # {M0}: the reference mass (a constant argument of the formula), {X}: x
    "ln1mx": "log(1-{X})",
    "ln1mx13": "log(1-pow({X},1./3.))",
    "L1": "log(@0/{M0})", "L2": "pow(log(@0/{M0}),2)", "L3": "pow(log(@0/{M0}),3)",
    "L4": "pow(log(@0/{M0}),4)",
    "z1": "(@0/{M0}-1)", "z2": "pow(@0/{M0}-1,2)", "z3": "pow(@0/{M0}-1,3)", "z4": "pow(@0/{M0}-1,4)",
}


def term_defs(m0):
    """term name -> (numpy function of m, RooFit expression of @0), for reference mass m0."""
    r = float(m0)
    fn = {
        "ln1mx": lambda m: np.log1p(-m / SQRTS),
        "ln1mx13": lambda m: np.log(1 - np.cbrt(m / SQRTS)),
        "L1": lambda m: np.log(m / r), "L2": lambda m: np.log(m / r) ** 2,
        "L3": lambda m: np.log(m / r) ** 3, "L4": lambda m: np.log(m / r) ** 4,
        "z1": lambda m: m / r - 1, "z2": lambda m: (m / r - 1) ** 2,
        "z3": lambda m: (m / r - 1) ** 3, "z4": lambda m: (m / r - 1) ** 4,
    }
    return {t: (fn[t], _TERM_EXPR[t].format(M0=repr(r), X=_X)) for t in fn}


def reference_mass(xmin, xmax):
    """M0: geometric mean of the fit range."""
    return float(np.sqrt(xmin * xmax))


LOGLINEAR = {  # family -> its terms, in order (N parameters: the first N-1)
    "dijet": ["ln1mx", "L1", "L2", "L3", "L4"],
    "moddijet": ["ln1mx13", "L1", "L2", "L3", "L4"],
    "ua2": ["L1", "z1", "z2", "z3", "z4"],
}
NONLINEAR = ("modexp", "polypow", "polyext")
FAMILIES = {**LOGLINEAR, **{f: None for f in NONLINEAR}}
FAMILY_ORDERS = {"dijet": (3, 4, 5, 6), "moddijet": (3, 4, 5, 6), "ua2": (3, 4, 5, 6),
                 "modexp": (3, 4, 5), "polypow": (3, 4, 5), "polyext": (5, 6)}
FAMILY_FORMULAS = {
    "dijet": "(1-x)^p1 / x^(p2 + p3 ln x + p4 ln^2 x [+ p5 ln^3 x])",
    "moddijet": "(1-x^(1/3))^p1 / x^(p2 + p3 ln x + p4 ln^2 x [+ p5 ln^3 x])",
    "ua2": "x^(-p1) exp(-(p2 x + p3 x^2 + p4 x^3 [+ p5 x^4]))",
    "modexp": "exp(p1 x^p2) | exp(p1 x^p2 + p1 (1-x)^p3) | exp(p1 x^p2 + p3 (1-x)^p4)",
    "polypow": "(1 + p1 x [+ p2 x^2 [+ p3 x^3]])^(-p_last)",
    "polyext": "(1-x)^p1 (1 + p4 x [+ p5 x^2]) x^-(p2 + p3 ln x)",
}
ORDERS = (3, 4, 5, 6)  # all orders used by any family (number of parameters incl. normalisation)


def function_names(families=None, orders=None):
    """All functions of the given families (default: all), each at its orders (optionally only
    those in `orders`)."""
    return [f"{fam}{n}" for fam in (families or FAMILIES) for n in FAMILY_ORDERS[fam]
            if orders is None or n in orders]


def parse_function(name):
    for fam in FAMILIES:
        if name.startswith(fam) and name[len(fam):].isdigit():
            n = int(name[len(fam):])
            if n in FAMILY_ORDERS[fam]:
                return fam, n
    raise ValueError(f"unknown background function '{name}'")


def is_loglinear(name):
    return parse_function(name)[0] in LOGLINEAR


def terms_of(name):
    """Terms of a log-linear function (empty for the other families)."""
    fam, n = parse_function(name)
    return LOGLINEAR[fam][:n - 1] if fam in LOGLINEAR else []


def formula_of(name, m0, M=None, tbar=None):
    """RooGenericPdf expression with @0 = m and @1.. = the shape parameters.

    Without a basis: exp(sum_k p_k T_k). With a basis (M, tbar), see decorrelated_basis:
    exp(sum_j q_j U_j) with U_j = sum_k M_kj (T_k - tbar_k)."""
    terms = terms_of(name)
    T = term_defs(m0)
    if M is None:
        return "exp(" + " + ".join(f"@{k + 1}*{T[t][1]}" for k, t in enumerate(terms)) + ")"
    centred = [f"({T[t][1]}-({float(tbar[k])!r}))" for k, t in enumerate(terms)]
    parts = []
    for j in range(len(terms)):
        u = "+".join(f"({float(M[k, j])!r})*{centred[k]}" for k in range(len(terms)))
        parts.append(f"@{j + 1}*({u})")
    return "exp(" + " + ".join(parts) + ")"


def describe(name):
    fam, n = parse_function(name)
    if fam in NONLINEAR:
        return (f"{FAMILY_FORMULAS[fam]} with {n - 1} shape parameters: {nonlinear_formula_text(name)}; "
                f"x = m/{SQRTS:g}")
    terms = terms_of(name)
    return (f"{FAMILY_FORMULAS[fam]} with {n - 1} shape parameters; implemented as "
            f"exp(sum_k a_k T_k) with T = {terms}, L = ln(m/M0), z = m/M0 - 1, "
            f"x = m/{SQRTS:g}, M0 = geometric mean of the fit range")


# ---- families that are not log-linear: ln f(x; p), x = m / sqrt(s), written as in the note
def nonlinear_log_f(name, x, p):
    """ln f (without the normalisation) of a non-linear function, numpy; p = (p1, ..., p_{N-1})."""
    fam, n = parse_function(name)
    x = np.asarray(x, float)
    with np.errstate(all="ignore"):
        if fam == "modexp":
            if n == 3:
                return p[0] * x ** p[1]
            if n == 4:
                return p[0] * (x ** p[1] + (1 - x) ** p[2])
            return p[0] * x ** p[1] + p[2] * (1 - x) ** p[3]
        if fam == "polypow":
            poly = 1 + sum(p[k] * x ** (k + 1) for k in range(n - 2))
            return -p[n - 2] * np.log(poly)
        if fam == "polyext":
            poly = 1 + p[3] * x + (p[4] * x * x if n == 6 else 0)
            return p[0] * np.log1p(-x) - (p[1] + p[2] * np.log(x)) * np.log(x) + np.log(poly)
    raise ValueError(name)


def nonlinear_expr(name, P, X):
    """RooFit expression of ln f for parameter expressions P[0..N-2] and the expression X of x."""
    fam, n = parse_function(name)
    if fam == "modexp":
        if n == 3:
            return f"{P[0]}*pow({X},{P[1]})"
        if n == 4:
            return f"{P[0]}*(pow({X},{P[1]})+pow(1-{X},{P[2]}))"
        return f"{P[0]}*pow({X},{P[1]})+{P[2]}*pow(1-{X},{P[3]})"
    if fam == "polypow":
        poly = "+".join(["1"] + [f"{P[k]}*pow({X},{k + 1})" for k in range(n - 2)])
        return f"-{P[n - 2]}*log({poly})"
    if fam == "polyext":
        poly = f"1+{P[3]}*{X}" + (f"+{P[4]}*pow({X},2)" if n == 6 else "")
        return f"{P[0]}*log(1-{X})-({P[1]}+{P[2]}*log({X}))*log({X})+log({poly})"
    raise ValueError(name)


def nonlinear_formula_text(name):
    n = parse_function(name)[1]
    return "exp(" + nonlinear_expr(name, [f"p{k}" for k in range(1, n)], "x") + ")"


# variable projection: for fixed values of the parameters that enter non-linearly (theta), ln f is
# linear in the others: ln f = offset(x; theta) + sum_l a_l B_l(x; theta).
_EXP = np.concatenate([-np.logspace(-2, 0.6, 10)[::-1], np.logspace(-2, 0.6, 10)])  # x^p exponents
_LIN = np.concatenate([-np.logspace(0, 4, 9)[::-1], np.logspace(0, 4, 9)])           # x, (1-x)^p
_SQ = np.concatenate([-np.logspace(1, 7, 7)[::-1], [0.0], np.logspace(1, 7, 7)])      # x^2
_CUB = np.concatenate([-np.logspace(2, 9, 8)[::-1], [0.0], np.logspace(2, 9, 8)])     # x^3


def _projection(name):
    """(theta grid, function theta -> (columns B_l, offset, p from (theta, a)))."""
    fam, n = parse_function(name)
    if fam == "modexp" and n == 3:
        return ([(e,) for e in _EXP],
                lambda x, t: ([x ** t[0]], 0.0, lambda a: [a[0], t[0]]))
    if fam == "modexp" and n == 4:
        return ([(e, b) for e in _EXP for b in _LIN],
                lambda x, t: ([x ** t[0] + (1 - x) ** t[1]], 0.0, lambda a: [a[0], t[0], t[1]]))
    if fam == "modexp":
        return ([(e, b) for e in _EXP for b in _LIN],
                lambda x, t: ([x ** t[0], (1 - x) ** t[1]], 0.0, lambda a: [a[0], t[0], a[1], t[1]]))
    if fam == "polypow":
        grids = {3: [(b,) for b in np.concatenate([_LIN, [-10.0, -5.0]])],
                 4: [(b, q) for b in _LIN for q in _SQ],
                 5: [(b, q, c) for b in _LIN[::2] for q in _SQ for c in _CUB]}[n]

        def proj(x, t):
            poly = 1 + sum(t[k] * x ** (k + 1) for k in range(len(t)))
            return [-np.log(poly)], 0.0, lambda a: [*t, a[0]]
        return grids, proj
    if fam == "polyext":
        grid = [(b,) for b in np.concatenate([_LIN, [0.0]])] if n == 5 else \
               [(b, q) for b in np.concatenate([_LIN, [0.0]]) for q in _SQ]

        def proj(x, t):
            poly = 1 + t[0] * x + (t[1] * x * x if n == 6 else 0)
            lx = np.log(x)
            return [np.log1p(-x), -lx, -lx * lx], np.log(poly), lambda a: [a[0], a[1], a[2], *t]
        return grid, proj
    raise ValueError(name)


def _irls(X, n, off, c=None, max_iter=100):
    """Poisson GLM (log link) maximum likelihood: nu = exp(X c + off). Returns (c, deviance)."""
    ok = n > 0
    if c is None:
        sw = np.sqrt(n[ok])
        c, *_ = np.linalg.lstsq(X[ok] * sw[:, None], (np.log(n[ok]) - off[ok]) * sw, rcond=None)
    dev = baker_cousins(n, np.exp(np.clip(X @ c + off, -700, 700)))
    for _ in range(max_iter):
        mu = np.exp(np.clip(X @ c + off, -700, 700))
        z = X @ c + (n - mu) / mu
        A = X * np.sqrt(mu)[:, None]
        c_new, *_ = np.linalg.lstsq(A, z * np.sqrt(mu), rcond=None)
        step = c_new - c
        dev_new = dev
        for _half in range(30):  # step halving if the deviance does not decrease
            dev_new = baker_cousins(n, np.exp(np.clip(X @ (c + step) + off, -700, 700)))
            if dev_new <= dev + 1e-9 * abs(dev):
                break
            step *= 0.5
        c = c + step
        converged = abs(dev - dev_new) < 1e-8 * max(1.0, dev)
        dev = dev_new
        if converged:
            break
    return c, dev


def _jacobian(name, x, p):
    """d ln f / d p (numerical, central differences), shape (len(x), len(p))."""
    J = np.empty((len(x), len(p)))
    for k in range(len(p)):
        h = 1e-6 * max(1.0, abs(p[k]))
        up, dn = np.array(p, float), np.array(p, float)
        up[k] += h
        dn[k] -= h
        J[:, k] = (nonlinear_log_f(name, x, up) - nonlinear_log_f(name, x, dn)) / (2 * h)
    return J


def _poisson_refine(name, x, n, off, p, c0, max_iter=300):
    """Full Poisson maximum likelihood over (c0, p) with Levenberg-Marquardt steps:
    nu = exp(c0 + off + ln f(x; p)). Returns (p, c0, deviance)."""
    def dev_of(p, c0):
        lf = nonlinear_log_f(name, x, p)
        if not np.all(np.isfinite(lf)):
            return np.inf
        return baker_cousins(n, np.exp(np.clip(c0 + off + lf, -700, 700)))

    p = np.array(p, float)
    dev, lam = dev_of(p, c0), 1e-3
    for _ in range(max_iter):
        nu = np.exp(np.clip(c0 + off + nonlinear_log_f(name, x, p), -700, 700))
        J = np.column_stack([np.ones(len(x)), _jacobian(name, x, p)])
        g = J.T @ (n - nu)
        H = (J * nu[:, None]).T @ J
        improved = False
        while lam < 1e12:
            try:
                d = np.linalg.solve(H + lam * np.diag(np.diag(H) + 1e-300), g)
            except np.linalg.LinAlgError:
                lam *= 10
                continue
            new = dev_of(p + d[1:], c0 + d[0])
            if new < dev:
                improved = abs(dev - new) > 1e-9 * max(1.0, dev)
                p, c0, dev = p + d[1:], c0 + d[0], new
                lam = max(lam / 3, 1e-9)
                break
            lam *= 4
        if not improved:
            break
    return p, c0, dev


PROFILE_STARTS, FINAL_STARTS = 30, 6  # start search of the non-linear functions
TIE_DEVIANCE = 0.1
_THETA_SCALE = {"exp": 0.01, "lin": 1.0, "sq": 10.0, "cub": 100.0}


def _theta_scales(name):
    """Scale of each non-linear parameter, for the asinh transformation of the search."""
    fam, n = parse_function(name)
    kinds = {("modexp", 3): ["exp"], ("modexp", 4): ["exp", "lin"], ("modexp", 5): ["exp", "lin"],
             ("polypow", 3): ["lin"], ("polypow", 4): ["lin", "sq"], ("polypow", 5): ["lin", "sq", "cub"],
             ("polyext", 5): ["lin"], ("polyext", 6): ["lin", "sq"]}[(fam, n)]
    return np.array([_THETA_SCALE[k] for k in kinds])


def _nelder_mead(f, u0, step=0.5, max_eval=400, tol=1e-6):
    """Minimise f from u0 (numpy Nelder-Mead). Returns (u, f(u))."""
    k = len(u0)
    pts = [np.array(u0, float)] + [np.array(u0, float) + step * np.eye(k)[i] for i in range(k)]
    vals = [f(u) for u in pts]
    neval = k + 1
    while neval < max_eval:
        order = np.argsort(vals)
        pts, vals = [pts[i] for i in order], [vals[i] for i in order]
        if np.isfinite(vals[0]) and abs(vals[-1] - vals[0]) < tol * max(1.0, abs(vals[0])):
            break
        c = np.mean(pts[:-1], axis=0)
        xr = c + (c - pts[-1]); fr = f(xr); neval += 1
        if fr < vals[0]:
            xe = c + 2 * (c - pts[-1]); fe = f(xe); neval += 1
            pts[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            xc = c + 0.5 * (pts[-1] - c); fc = f(xc); neval += 1
            if fc < vals[-1]:
                pts[-1], vals[-1] = xc, fc
            else:  # shrink
                pts = [pts[0]] + [pts[0] + 0.5 * (q - pts[0]) for q in pts[1:]]
                vals = [vals[0]] + [f(q) for q in pts[1:]]
                neval += k
    i = int(np.argmin(vals))
    return pts[i], vals[i]


def start_values_nonlinear(name, hist):
    """Maximum-likelihood estimate of p (numpy) for a non-linear function, used as the starting
    point. The parameters that enter non-linearly (theta, 1-3 of them) are searched with the
    others profiled: for fixed theta the model is a Poisson GLM, solved exactly (IRLS). A grid
    over theta gives the starting points of Nelder-Mead minimisations of the profiled deviance
    (in asinh(theta / scale)); the best results are refined by a full Poisson fit of all
    parameters. Returns (p, expected counts per bin at the solution)."""
    lo, hi, n = hist_arrays(hist)
    mc, w = 0.5 * (lo + hi), hi - lo
    x, off_w = mc / SQRTS, np.log(w)
    grid, proj = _projection(name)
    scale = _theta_scales(name)
    cache = {}

    def solve(t):
        key = tuple(np.round(t, 12))
        if key in cache:
            return cache[key]
        out = None
        with np.errstate(all="ignore"):
            cols, off, to_p = proj(x, t)
            off = off_w + off
            if np.all(np.isfinite(off)) and all(np.all(np.isfinite(b)) for b in cols):
                X = np.column_stack([np.ones(len(x))] + cols)
                try:
                    c, dev = _irls(X, n, off, max_iter=20)
                    if np.isfinite(dev):
                        out = (dev, to_p(c[1:]), c[0])
                except np.linalg.LinAlgError:
                    pass
        cache[key] = out
        return out

    def prof(u):
        r = solve(tuple(np.sinh(np.clip(u, -40, 40)) * scale))
        return r[0] if r else np.inf

    starts = []
    for t in grid:
        r = solve(t)
        if r:
            starts.append((r[0], t))
    if not starts:
        raise RuntimeError(f"{name}: no valid starting point")
    starts.sort(key=lambda z: z[0])
    # diverse starting points: the best grid point of each region (points at least 1 unit apart
    # in asinh(theta / scale) from all those already taken), so that separate valleys are explored
    chosen = []
    for d, t in starts:
        u = np.arcsinh(np.asarray(t) / scale)
        if all(np.max(np.abs(u - v)) > 0.7 for _, _, v in chosen):
            chosen.append((d, t, u))
        if len(chosen) == PROFILE_STARTS:
            break
    found = []
    for _dev, t, _u in chosen:
        u, d = _nelder_mead(prof, np.arcsinh(np.asarray(t) / scale))
        if np.isfinite(d):  # map back exactly as prof() did (clipped), so solve(t) is the cached point
            found.append((d, tuple(np.sinh(np.clip(u, -40, 40)) * scale)))
    found.sort(key=lambda z: z[0])
    results = []
    for d, t in found[:FINAL_STARTS]:
        _dev, p, c0 = solve(t)
        with np.errstate(all="ignore"):
            p, c0, dev = _poisson_refine(name, x, n, off_w, p, c0, max_iter=100)
        if np.isfinite(dev):
            mu = np.exp(c0 + off_w + nonlinear_log_f(name, x, p))
            results.append((dev, _condition(name, x, p, mu), np.asarray(p, float), mu))
    if not results:
        raise RuntimeError(f"{name}: no valid starting point")
    # equally good solutions (deviance within TIE_DEVIANCE of the best): take the best-conditioned
    # one. E.g. polyext5 on a dijet4-like spectrum: p4 -> infinity is as good as p4 ~ 0, since
    # (1 + p4 x) -> p4 x only shifts p2, but the fit cannot work there.
    dmin = min(r[0] for r in results)
    _dev, _cond, p, mu = min((r for r in results if r[0] <= dmin + TIE_DEVIANCE), key=lambda r: r[1])
    return p, mu


def _condition(name, x, p, mu):
    """Condition number of the Fisher matrix of the shape parameters (inf if singular)."""
    with np.errstate(all="ignore"):
        F = _fisher(name, x, p, mu)
    if not np.all(np.isfinite(F)):
        return np.inf
    ev = np.linalg.eigvalsh(F)
    return ev[-1] / ev[0] if ev[0] > 0 else np.inf


def _fisher(name, x, p, mu):
    J = _jacobian(name, x, p)
    Jc = J - (mu[:, None] * J).sum(axis=0) / mu.sum()
    return (Jc * mu[:, None]).T @ Jc


def decorrelated_nonlinear(name, hist, p, mu):
    """M such that p = p_start + M q gives uncorrelated q with unit variance at the start
    (linearised: Fisher information of the centred Jacobian d ln f / d p)."""
    lo, hi, _ = hist_arrays(hist)
    x = 0.5 * (lo + hi) / SQRTS
    F = _fisher(name, x, p, mu)
    try:
        return np.linalg.inv(np.linalg.cholesky(F)).T
    except np.linalg.LinAlgError:
        # a direction without information (degenerate start): eigen-decomposition with a floor,
        # so that such a parameter moves in steps of its own size; the fit will tell
        ev, V = np.linalg.eigh(F)
        floor = max(ev[-1], 1.0) * 1e-12
        return V / np.sqrt(np.maximum(ev, floor))


def formula_and_constants(name, k, m0, M, centre):
    """RooFit formula of a function in the decorrelated parametrisation, with @0 = m,
    @1..@k = q, and the per-fit constants as further arguments @(k+1).. (so that all fits of the
    function share one compiled formula). Returns (formula, constants).
      log-linear: ln f = sum_j q_j sum_i M_ij (T_i - tbar_i), centre = tbar, constants M0, tbar, M
      non-linear: ln f = ln f(x; p) - ln f(x0; p),  p_i = p_start_i + sum_j M_ij q_j,
                  x0 = M0/sqrt(s), centre = p_start, constants p_start, M, M0"""
    M = np.asarray(M, float)
    if is_loglinear(name):
        i_m0, i_c, i_M = k + 1, k + 2, 2 * k + 2
        T = [_TERM_EXPR[t].format(M0=f"@{i_m0}", X=_X) for t in terms_of(name)]
        parts = []
        for j in range(k):
            u = "+".join(f"@{i_M + i * k + j}*({T[i]}-@{i_c + i})" for i in range(k))
            parts.append(f"@{j + 1}*({u})")
        return "exp(" + "+".join(parts) + ")", [float(m0)] + [float(v) for v in centre] + M.ravel().tolist()
    i_c, i_M, i_m0 = k + 1, 2 * k + 1, 2 * k + 1 + k * k
    P = ["(@{}+{})".format(i_c + i, "+".join(f"@{i_M + i * k + j}*@{j + 1}" for j in range(k)))
         for i in range(k)]
    # ln f(x) - ln f(x0), x0 = M0/sqrt(s): the same shape, but without the huge constant that
    # exp(ln f) can carry (it would overflow; the normalisation is separate anyway)
    X0 = f"(@{i_m0}/{SQRTS!r})"
    expr = f"({nonlinear_expr(name, P, _X)})-({nonlinear_expr(name, P, X0)})"
    return "exp(" + expr + ")", [float(v) for v in centre] + M.ravel().tolist() + [float(m0)]


def hist_arrays(hist):
    nb = hist.GetNbinsX()
    ax = hist.GetXaxis()
    lo = np.array([ax.GetBinLowEdge(i) for i in range(1, nb + 1)])
    hi = np.array([ax.GetBinUpEdge(i) for i in range(1, nb + 1)])
    n = np.array([hist.GetBinContent(i) for i in range(1, nb + 1)])
    return lo, hi, n


def start_values(name, hist, m0, max_iter=100):
    """Maximum-likelihood estimate of the shape parameters, used as starting values.

    The binned model  nu_i = w_i * exp(c0 + sum_k p_k T_k(m_i))  is a Poisson generalised
    linear model with log link, so its maximum-likelihood solution is found reliably by
    iteratively reweighted least squares. It is very close to the minimum of the RooFit extended
    binned likelihood (where c0 is traded for the yield and the pdf integral; the two differ
    only by the bin-centre vs integral normalisation), so Minuit starts next to the minimum.
    Returns (p_1..p_k, expected counts per bin at the solution).
    """
    lo, hi, n = hist_arrays(hist)
    mc, w = 0.5 * (lo + hi), hi - lo
    T = term_defs(m0)
    X = np.column_stack([np.ones(len(n))] + [T[t][0](mc) for t in terms_of(name)])
    off = np.log(w)
    c, _ = _irls(X, n, off, max_iter=max_iter)
    mu = np.exp(X @ c + off)
    return c[1:], mu


def decorrelated_basis(name, hist, mu, m0):
    """Linear combinations U_j = sum_k M_kj (T_k(m) - tbar_k) of the terms of `name` that are
    uncorrelated, with unit variance, at the expected bin contents `mu` (the starting point).
    In this basis Minuit sees a nearly diagonal, well-scaled problem; the family of functions
    is unchanged (the constants tbar_k are absorbed in the normalisation).
    Returns (M, tbar)."""
    lo, hi, _ = hist_arrays(hist)
    mc = 0.5 * (lo + hi)
    defs = term_defs(m0)
    T = np.column_stack([defs[t][0](mc) for t in terms_of(name)])
    tbar = (mu[:, None] * T).sum(axis=0) / mu.sum()
    Tc = T - tbar
    F = (Tc * mu[:, None]).T @ Tc
    L = np.linalg.cholesky(F)
    M = np.linalg.inv(L).T  # U = Tc @ M  ->  U^T diag(mu) U = 1
    return M, tbar


def baker_cousins(n, nu):
    """Poisson likelihood-ratio chi2, summed over bins."""
    n = np.asarray(n, float)
    nu = np.asarray(nu, float)
    if np.any(nu <= 0):
        return float("inf")
    term = nu - n
    pos = n > 0
    term[pos] += n[pos] * np.log(n[pos] / nu[pos])
    return float(2.0 * term.sum())


def make_pdf(ROOT, pdf_name, func, arglist, formula, nbins=None):
    """RooGenericPdf for background function `func` (arglist = observable, p1, p2, ...).

    nbins given (uniform bins evaluated at their centres, i.e. not rebinned): the pdf is
    normalised by the sum of its values at the nbins bin centres times the bin width
    (RooBinIntegrator), not by its integral. The expected counts per bin, N * pdf(centre) *
    width, then add up exactly to the yield N, i.e. the extended likelihood of RooFit and Combine
    is the binned Poisson likelihood of these expectations. With the integral as normalisation
    the per-bin expectations add up to N * (1 + O(1e-5)) (midpoint-rule error), and the
    parameter dependence of that mismatch is enough, with O(10^7) events, to pull the shape
    parameters by O(1) standard deviations and the signal strength in Combine by up to ~1
    standard deviation (also between Combine's own toys and its fits).

    nbins None (bins integrated over, see fit_function): the integral is used, with an adaptive
    Gauss-Kronrod integrator of relative precision 1e-10 (with the default 1e-7 the likelihood
    is noisy at the O(1) level for these sample sizes).
    The integrator setting is stored with the pdf in the workspace."""
    pdf = ROOT.RooGenericPdf(pdf_name, func, formula, arglist)
    if nbins:
        set_bin_sum_normalisation(pdf, nbins)
    else:
        cfg = pdf.specialIntegratorConfig(True)
        cfg.setEpsRel(INTEGRATOR_EPS_REL)
        cfg.setEpsAbs(0.0)
        cfg.method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D")
    return pdf


def set_bin_sum_normalisation(pdf, nbins):
    """Normalise `pdf` by the sum of its values at the centres of nbins uniform bins over the
    range of its observable, times the bin width (RooBinIntegrator; the setting is stored with
    the pdf). Analytical integrals (e.g. of RooCrystalBall) are switched off for this."""
    pdf.forceNumInt(True)
    cfg = pdf.specialIntegratorConfig(True)
    cfg.method1D().setLabel("RooBinIntegrator")
    cfg.getConfigSection("RooBinIntegrator").setRealValue("numBins", int(nbins))
    return pdf


def bin_sum_nbins(lo, hi):
    """Number of bins for the bin-sum normalisation of make_pdf, or None if the bins are
    integrated over (wider than 1 GeV) or not uniform."""
    w = hi - lo
    if needs_bin_integration(lo, hi) or np.max(w) - np.min(w) > 1e-6 * np.max(w):
        return None
    return len(lo)


def at_limit_of(var, frac=1e-3):
    """True if a RooRealVar is within frac x (its range) of either end of its range."""
    span = var.getMax() - var.getMin()
    return var.getVal() - var.getMin() < frac * span or var.getMax() - var.getVal() < frac * span


def needs_bin_integration(lo, hi):
    return bool(np.max(hi - lo) > BIN_INTEGRATION_ABOVE)


def expected_counts(ROOT, x, pdf, norm, lo, hi, integrate, npts=16):
    """Expected events per bin: norm * pdf (normalised on the range of x) at the bin centre times
    the bin width, or, if `integrate`, integrated over the bin (composite Simpson rule)."""
    nset = ROOT.RooArgSet(x)

    def dens(m):
        x.setVal(float(m))
        return pdf.getVal(nset)

    out = []
    for a, b in zip(lo, hi):
        if not integrate:
            out.append(dens(0.5 * (a + b)) * (b - a))
            continue
        ts = np.linspace(a, b, npts + 1)
        fs = np.array([dens(t) for t in ts])
        h = (b - a) / npts
        out.append(h / 3 * (fs[0] + fs[-1] + 4 * fs[1:-1:2].sum() + 2 * fs[2:-1:2].sum()))
    return norm * np.array(out)


def fit_function(ROOT, name, hist, tag, strategy=1):
    """Extended binned maximum-likelihood fit of function `name` to `hist` (whose axis is the
    fit range). Returns a dict with the results (python values) and the RooFit objects.

    With 1 GeV bins the pdf is evaluated at the bin centres (RooFit's default for binned data)
    and normalised by its sum over the bin centres (see make_pdf), so that the expected counts
    add up to the fitted yield. With wider bins (rebinned categories) the pdf is integrated over
    each bin: RooFit's default binned likelihood (bin-centre density, but total yield from the
    exact integral) is then biased at the level of many standard deviations for these sample
    sizes."""
    lo, hi, n = hist_arrays(hist)
    xmin, xmax = float(lo[0]), float(hi[-1])
    edges = list(lo) + [xmax]

    x = ROOT.RooRealVar(f"mjj_{tag}", "m_{jj} [GeV]", xmin, xmax)
    x.setBinning(ROOT.RooBinning(len(edges) - 1, np.array(edges, dtype="d")))
    m0 = reference_mass(xmin, xmax)
    if is_loglinear(name):
        pvals, mu = start_values(name, hist, m0)
        M, centre = decorrelated_basis(name, hist, mu, m0)
        qvals = np.linalg.solve(M, pvals)  # p = M q
        offset = np.zeros(len(pvals))
    else:
        pvals, mu = start_values_nonlinear(name, hist)
        M = decorrelated_nonlinear(name, hist, pvals, mu)
        centre = offset = pvals  # p = p_start + M q
        qvals = np.zeros(len(pvals))
    formula, constants = formula_and_constants(name, len(qvals), m0, M, centre)
    params = []
    arglist = ROOT.RooArgList(x)
    for k, v in enumerate(qvals, start=1):
        # in the decorrelated basis the parameter uncertainties are ~1
        p = ROOT.RooRealVar(f"{name}_p{k}_{tag}", f"q{k}", v, v - QRANGE, v + QRANGE)
        params.append(p)
        arglist.add(p)
    consts = add_constants(ROOT, arglist, constants, f"{name}_c", tag)
    pdf = make_pdf(ROOT, f"{name}_{tag}", name, arglist, formula, bin_sum_nbins(lo, hi))
    ntot = float(n.sum())
    norm = ROOT.RooRealVar(f"{name}_{tag}_norm", "N_{bkg}", ntot, NORM_RANGE[0] * ntot, NORM_RANGE[1] * ntot)
    model = ROOT.RooExtendPdf(f"{name}_{tag}_ext", "", pdf, norm)
    data = ROOT.RooDataHist(f"data_{tag}", "", ROOT.RooArgList(x), hist)

    integrate = needs_bin_integration(lo, hi)
    opts = [ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.Strategy(strategy),
            ROOT.RooFit.Offset(True), ROOT.RooFit.Extended(True),
            ROOT.RooFit.IntegrateBins(BIN_INTEGRATION_PRECISION if integrate else -1)]
    # A shape parameter that ends at its limit: re-centre its range and fit again (at most
    # MAX_RECENTRE times). The at-limit check below is made on the result of the last fit, so a
    # parameter that still ends at a limit after the last attempt is reported as such.
    texts = []  # RooFit/Minuit output of every attempt (messages are counted over all of them)
    with capture_output() as log:
        res = model.fitTo(data, *opts)
    texts.append(log["text"])
    for _ in range(MAX_RECENTRE):
        hit = [p for p in params if at_limit_of(p)]
        if not hit:
            break
        for p in hit:
            p.setRange(p.getVal() - QRANGE, p.getVal() + QRANGE)
        with capture_output() as log:
            res = model.fitTo(data, *opts)
        texts.append(log["text"])
    all_text = "\n".join(texts)
    nwarn, nerr = count_messages(all_text)

    # expected events per bin and the goodness of fit
    nu = expected_counts(ROOT, x, pdf, norm.getVal(), lo, hi, integrate)
    chi2 = max(baker_cousins(n, nu), 0.0)  # rounding can give -1e-10 on exact Asimov data (p would be 0)
    nfloat = res.floatParsFinal().getSize()
    ndf = len(n) - nfloat
    pvalue = ROOT.TMath.Prob(chi2, ndf) if (ndf > 0 and math.isfinite(chi2)) else 0.0

    at_limit = [p.GetName() for p in params + [norm] if at_limit_of(p)]

    nset = ROOT.RooArgSet(x)
    xs = np.linspace(xmin, xmax, 400)
    curve = []
    for m in xs:
        x.setVal(m)
        curve.append(norm.getVal() * pdf.getVal(nset))

    return dict(
        name=name, family=parse_function(name)[0], npar=parse_function(name)[1],
        status=int(res.status()), cov_qual=int(res.covQual()), edm=float(res.edm()),
        min_nll=float(res.minNll()), chi2=chi2, ndf=int(ndf), pvalue=float(pvalue),
        n_warnings=nwarn, n_errors=nerr, messages=message_summary(all_text), n_attempts=len(texts),
        params_at_limit=at_limit,
        params={f"p{k}": dict(value=p.getVal(), error=p.getError(), min=p.getMin(),
                              max=p.getMax()) for k, p in enumerate(params, start=1)},
        norm=dict(value=norm.getVal(), error=norm.getError()),
        expected=nu.tolist(), curve=(xs.tolist(), curve), bin_integration=integrate,
        formula=formula, constants=constants,
        basis=dict(kind="loglinear" if is_loglinear(name) else "nonlinear", terms=terms_of(name),
                   m0=m0, M=M.tolist(), centre=np.asarray(centre).tolist()),
        params_original=original_params(M, params, res, offset),
        _keep=(x, params, consts, pdf, norm, model, data, res),
    )


def failed_fit(name, nbins, exc):
    """Result record for a fit that could not be performed (an exception was raised)."""
    fam, npar = parse_function(name)
    return dict(name=name, family=fam, npar=npar, status=-1, cov_qual=-1, edm=float("nan"),
                min_nll=float("nan"), chi2=float("inf"), ndf=int(nbins - npar), pvalue=0.0,
                n_warnings=0, n_errors=1, messages=[f"{type(exc).__name__}: {exc}"],
                params_at_limit=[], params={}, norm=dict(value=float("nan"), error=float("nan")),
                expected=[], curve=([], []), bin_integration=False, formula=None, constants=None,
                basis=None,
                params_original={}, exception=f"{type(exc).__name__}: {exc}")


def add_constants(ROOT, arglist, constants, prefix, tag):
    """Append the per-fit constants of a formula (see formula_and_constants) to arglist as
    RooConstVar; returns them (to keep them alive)."""
    out = []
    for i, v in enumerate(constants or []):
        c = ROOT.RooConstVar(f"{prefix}{i}_{tag}", "", float(v))
        arglist.add(c)
        out.append(c)
    return out


def original_params(M, params, res, offset=None):
    """The fitted parameters p_k in the published form (p = offset + M q; log-linear: the
    coefficients of the terms T_k, offset 0), with errors from the fit covariance, for reporting."""
    q = np.array([p.getVal() for p in params])
    pv = M @ q + (0 if offset is None else np.asarray(offset))
    try:
        names = [p.GetName() for p in params]
        fl = res.floatParsFinal()
        idx = [fl.index(fl.find(nm)) for nm in names]
        C = np.array([[res.covarianceMatrix()[i][j] for j in idx] for i in idx])
        pe = np.sqrt(np.clip(np.diag(M @ C @ M.T), 0, None))
    except Exception:
        pe = np.full(len(pv), float("nan"))
    return {f"p{k}": dict(value=float(v), error=float(e))
            for k, (v, e) in enumerate(zip(pv, pe), start=1)}


def fit_quality_ok(fit):
    return fit["status"] == 0 and fit["cov_qual"] == 3 and not fit["params_at_limit"]


def ftest(fit_lo, fit_hi, nbins):
    """p-value of the Fisher F-test: is fit_hi (more parameters) significantly better?"""
    import ROOT
    dp = fit_hi["npar"] - fit_lo["npar"]
    ndf_hi = nbins - fit_hi["npar"]
    if dp <= 0 or ndf_hi <= 0 or not math.isfinite(fit_lo["chi2"]) or fit_hi["chi2"] <= 0:
        return None
    F = ((fit_lo["chi2"] - fit_hi["chi2"]) / dp) / (fit_hi["chi2"] / ndf_hi)
    if F <= 0:
        return 1.0
    return float(1.0 - ROOT.Math.fdistribution_cdf(F, dp, ndf_hi))


def acceptable(fit, pvalue_threshold):
    """A valid background fit: status 0, full covariance matrix (covQual 3), no parameter at its
    limit, and goodness-of-fit p-value above the threshold."""
    return fit_quality_ok(fit) and fit["pvalue"] > pvalue_threshold


def select_functions(fits, nbins, pvalue_threshold=0.05, ftest_threshold=0.05, max_alternates=4,
                     family_max_only=False):
    """Choose the functions of the RooMultiPdf: a default and up to `max_alternates` alternates.

    1. Only acceptable fits are used (see `acceptable`; p-value > pvalue_threshold).
    2. In every family, a sequence of standard Fisher F-tests fixes the highest order: start
       from the lowest acceptable order; move to the next order if it is significantly better
       (F-test p-value < ftest_threshold) and its own fit is acceptable; otherwise stop. The
       family contributes its acceptable orders from the lowest one up to this highest order.
    3. Default: among the highest orders of the families, the one with the largest p-value.
    4. Alternates: the other contributed functions, ranked by p-value, at most max_alternates.
       With family_max_only, a family contributes only its highest order (at most one function
       per family in the RooMultiPdf).

    Every fit gets fit['ftest_next'] (F-test p-value against the next order of its family,
    when both fits have a good quality) and fit['ftest_max'] (the highest order of its family).
    Returns (verdict {name: 'pass (default)', 'pass (alternate k)' or the reason it is not used},
    default name or None)."""
    verdict, candidates, family_max = {}, [], {}
    for fam in sorted({f["family"] for f in fits}):
        ff = sorted((f for f in fits if f["family"] == fam), key=lambda f: f["npar"])
        by_npar = {f["npar"]: f for f in ff}
        for f in ff:
            nxt = by_npar.get(f["npar"] + 1)
            f["ftest_next"] = (ftest(f, nxt, nbins) if nxt and fit_quality_ok(f) and fit_quality_ok(nxt)
                               else None)
        start = next((f for f in ff if acceptable(f, pvalue_threshold)), None)
        top = start
        while top is not None:
            nxt = by_npar.get(top["npar"] + 1)
            if (nxt is None or top["ftest_next"] is None or top["ftest_next"] >= ftest_threshold
                    or not acceptable(nxt, pvalue_threshold)):
                break
            top = nxt
        for f in ff:
            f["ftest_max"] = top["name"] if top else None
        if top is not None:
            family_max[fam] = top
        for f in ff:
            if f.get("exception"):
                verdict[f["name"]] = "fail: the fit raised an exception"
            elif not fit_quality_ok(f):
                why = []
                if f["status"] != 0:
                    why.append(f"status {f['status']}")
                if f["cov_qual"] != 3:
                    why.append(f"covQual {f['cov_qual']}")
                if f["params_at_limit"]:
                    why.append("parameter at limit")
                verdict[f["name"]] = "fail: " + ", ".join(why)
            elif f["pvalue"] <= pvalue_threshold:
                verdict[f["name"]] = f"fail: p-value {f['pvalue']:.3g}"
            elif f["npar"] > top["npar"]:
                if top["npar"] + 1 not in by_npar:  # e.g. --functions dijet3 dijet5
                    verdict[f["name"]] = (f"not used: the {fam}{top['npar'] + 1} fit is missing, so no "
                                          f"F-test beyond {top['name']}")
                else:
                    verdict[f["name"]] = f"not used: above the F-test maximum ({top['name']})"
            elif family_max_only and f is not top:
                verdict[f["name"]] = f"not used: below the family maximum ({top['name']})"
            else:
                candidates.append(f)
    if not family_max:
        return verdict, None
    default = max(family_max.values(), key=lambda f: f["pvalue"])
    verdict[default["name"]] = "pass (default)"
    others = sorted((f for f in candidates if f is not default), key=lambda f: -f["pvalue"])
    for k, f in enumerate(others, start=1):
        verdict[f["name"]] = (f"pass (alternate {k})" if k <= max_alternates
                              else f"not used: beyond the {max_alternates} alternates")
    return verdict, default["name"]
