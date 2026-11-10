"""Shared pytest fixtures and options. Run the tests inside the analysis environment:
    pytest [--fit-xmin 150] [--fit-xmax 600]
The fit range used by the tests is configurable (default 150-600 GeV); it has to lie inside
the synthetic histogram range (tests/synthetic.py, HIST_RANGE)."""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tests"))


def pytest_addoption(parser):
    parser.addoption("--fit-xmin", type=float, default=150.0,
                     help="lower edge of the fit range used by the tests [GeV] (default: 150)")
    parser.addoption("--fit-xmax", type=float, default=600.0,
                     help="upper edge of the fit range used by the tests [GeV] (default: 600)")


@pytest.fixture(scope="session")
def fit_range(request):
    import synthetic
    xmin, xmax = request.config.getoption("--fit-xmin"), request.config.getoption("--fit-xmax")
    lo, hi = synthetic.HIST_RANGE
    if not lo <= xmin < xmax <= hi:
        pytest.exit(f"fit range [{xmin:g}, {xmax:g}] must lie inside the synthetic histogram range "
                    f"[{lo:g}, {hi:g}]", returncode=4)
    return xmin, xmax


@pytest.fixture(scope="session")
def ROOT():
    from dijetfit.common import setup_root
    return setup_root()


def make_th1(ROOT, name, edges, contents, errors=None):
    """TH1D with the given bin edges, contents and (default: sqrt) errors."""
    from array import array
    import numpy as np
    h = ROOT.TH1D(name, name, len(edges) - 1, array("d", [float(e) for e in edges]))
    h.SetDirectory(0)
    h.Sumw2()
    errors = np.sqrt(np.abs(contents)) if errors is None else errors
    for i, (c, e) in enumerate(zip(contents, errors), start=1):
        h.SetBinContent(i, float(c))
        h.SetBinError(i, float(e))
    return h
