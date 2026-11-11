"""Choice of the common fit range (lower edge) and of the per-category rebinning.

The decision logic is kept separate from the fits so that it can be tested on its own: it only
needs a function evaluate(xmin, category, rebin) -> (result, passes).

Priority (agreed analysis convention):
  1. no rebinning: raise the lower edge from `start` in steps of `step`, up to `max_unbinned`,
     until every category passes;
  2. if some categories still fail there: try the rebin factors, in increasing order, for the
     failing categories only; each keeps the smallest factor that works (a category for which
     none works goes back to no rebinning);
  3. if some still fail: raise the lower edge further, up to `scan_max`. At every new lower edge
     each category starts again from no rebinning and keeps the smallest factor that works, so
     that no category is rebinned when it is not needed.
"""


def smallest_working_factor(evaluate, xmin, cat, factors):
    """Try `factors` in increasing order; return (factor, result, passes) for the first that
    passes, or for the smallest factor (the least rebinning) if none does."""
    tried = []
    for fac in sorted(factors):
        result, ok = evaluate(xmin, cat, fac)
        if ok:
            return fac, result, True
        tried.append((fac, result))
    return tried[0][0], tried[0][1], False


def scan_fit_range(evaluate, cats, start, step, max_unbinned, scan_max, rebin_factors, log=print):
    """Run the scan; returns (xmin, {cat: rebin}, {cat: result}, [failing categories])."""
    eps = 1e-9

    def report(xmin, rebin, failing):
        log(f"  xmin = {xmin:g} GeV: {len(cats) - len(failing)}/{len(cats)} categories with a passing "
            "function" + (f"; failing: {', '.join(failing)}" if failing else "")
            + ("" if all(rebin[c] == 1 for c in cats) else
               "; rebin: " + ", ".join(f"{c}={rebin[c]}" for c in cats if rebin[c] != 1)))

    log("Stage 1: original binning (no rebinning)")
    xmin = start
    while True:
        rebin = {c: 1 for c in cats}
        results, ok = {}, {}
        for c in cats:
            results[c], ok[c] = evaluate(xmin, c, 1)
        failing = [c for c in cats if not ok[c]]
        report(xmin, rebin, failing)
        if not failing or xmin + step > max_unbinned + eps:
            break
        xmin += step
    if not failing:
        return xmin, rebin, results, failing

    log(f"Stage 2: rebin the failing categories at xmin = {xmin:g} GeV")
    for c in failing:  # no rebinning was just tried in stage 1: try the factors only
        fac, res, good = smallest_working_factor(evaluate, xmin, c, list(rebin_factors))
        if good:
            rebin[c], results[c], ok[c] = fac, res, True
        log(f"  {c}: " + (f"passes with rebin {rebin[c]}" if ok[c] else
                          "fails with every rebin factor (kept without rebinning)"))
    failing = [c for c in cats if not ok[c]]
    report(xmin, rebin, failing)

    if failing:
        log("Stage 3: raise the lower edge further (each category with the smallest rebin factor that works)")
    while failing and xmin + step <= scan_max + eps:
        xmin += step
        for c in cats:
            rebin[c], results[c], ok[c] = smallest_working_factor(evaluate, xmin, c, [1] + list(rebin_factors))
        failing = [c for c in cats if not ok[c]]
        report(xmin, rebin, failing)
    return xmin, rebin, results, failing
