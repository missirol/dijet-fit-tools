"""Unit tests for the fit-range scan logic (dijetfit/scan.py), with a fake fit: each category
passes for the (xmin, rebin) combinations given by a rule, so that every stage can be forced."""
from dijetfit.scan import scan_fit_range, smallest_working_factor

START, STEP, MAX_UNBINNED, MAX = 150.0, 5.0, 170.0, 200.0
FACTORS = [2, 5, 10]


def run(rules, cats=("A", "B")):
    """rules: {cat: function(xmin, rebin) -> passes}; returns the scan result and the calls."""
    calls = []

    def evaluate(xmin, cat, fac):
        calls.append((xmin, cat, fac))
        ok = rules[cat](xmin, fac)
        return (xmin, cat, fac, ok), ok

    log = []
    out = scan_fit_range(evaluate, list(cats), START, STEP, MAX_UNBINNED, MAX, FACTORS, log.append)
    return out, calls, log


def always(x, f):
    return True


def test_everything_passes_at_the_start():
    (xmin, rebin, results, failing), calls, _ = run({"A": always, "B": always})
    assert xmin == START and rebin == {"A": 1, "B": 1} and failing == []
    assert calls == [(START, "A", 1), (START, "B", 1)]


def test_stage1_raises_the_lower_edge_without_rebinning():
    (xmin, rebin, _, failing), _, _ = run({"A": always, "B": lambda x, f: x >= 160})
    assert xmin == 160 and rebin == {"A": 1, "B": 1} and failing == []


def test_stage2_uses_the_smallest_factor_that_works_for_failing_categories_only():
    (xmin, rebin, results, failing), _, _ = run({"A": always, "B": lambda x, f: f >= 5})
    assert xmin == MAX_UNBINNED and failing == []
    assert rebin == {"A": 1, "B": 5}
    assert results["B"] == (MAX_UNBINNED, "B", 5, True)


def test_stage3_retries_without_rebinning_at_a_higher_edge():
    """Regression: a category that fails every factor in stage 2 must not stay at the largest
    factor; at the next edge it passes without rebinning."""
    (xmin, rebin, results, failing), _, _ = run({"A": always, "B": lambda x, f: x >= 175 and f == 1})
    assert xmin == 175 and failing == []
    assert rebin == {"A": 1, "B": 1}
    assert results["B"] == (175, "B", 1, True)


def test_stage3_removes_rebinning_that_is_no_longer_needed():
    """A category rebinned in stage 2 goes back to no rebinning when the edge is raised, if it
    passes there without it."""
    rules = {"A": lambda x, f: f >= 2 or x >= 175,   # rebinned at 170, unbinned OK from 175
             "B": lambda x, f: x >= 175}             # forces stage 3
    (xmin, rebin, _, failing), _, _ = run(rules)
    assert xmin == 175 and failing == [] and rebin == {"A": 1, "B": 1}


def test_nothing_works():
    (xmin, rebin, results, failing), _, log = run({"A": always, "B": lambda x, f: False})
    assert xmin == MAX and failing == ["B"]
    assert rebin["B"] == 1  # the least rebinning is reported for a failing category
    assert results["B"][3] is False
    assert any("Stage 3" in line for line in log)


def test_smallest_working_factor():
    seen = []

    def ev(x, c, f):
        seen.append(f)
        return f, f in (5, 10)

    assert smallest_working_factor(ev, 0, "c", [10, 1, 5, 2]) == (5, 5, True)
    assert seen == [1, 2, 5]
    assert smallest_working_factor(lambda x, c, f: (f, False), 0, "c", [1, 2]) == (1, 1, False)


def test_stage_boundaries():
    """Without rebinning the edge never goes beyond --scan-max-unbinned; rebinning is first tried
    exactly there; the edge never goes beyond --scan-max."""
    seq = []

    def ev(x, c, f):
        seq.append((x, f))
        return None, False

    xmin, _, _, failing = scan_fit_range(ev, ["A"], START, STEP, MAX_UNBINNED, MAX, FACTORS, lambda s: None)
    assert failing == ["A"] and xmin == MAX
    first_rebin = next(i for i, (x, f) in enumerate(seq) if f > 1)
    assert seq[first_rebin][0] == MAX_UNBINNED
    assert [x for x, f in seq[:first_rebin]] == [150.0, 155.0, 160.0, 165.0, 170.0]
    assert max(x for x, f in seq) == MAX
