"""Unit tests for dijetfit.histtools: fit-range windows and rebinning."""
import numpy as np
import pytest

from conftest import make_th1
from dijetfit.histtools import bin_edges_in_window, binning_label, window_hist


@pytest.fixture
def hist(ROOT):
    edges = np.arange(100.0, 121.0, 1.0)  # 20 bins of 1 GeV
    contents = np.arange(1.0, 21.0)
    return make_th1(ROOT, "h_hist", edges, contents)


def test_edges_window_and_rebin(hist):
    assert bin_edges_in_window(hist, 102, 110) == list(np.arange(102.0, 111.0))
    assert bin_edges_in_window(hist, 102, 110, rebin=2) == [102.0, 104.0, 106.0, 108.0, 110.0]
    # 9 bins merged by 4: the last bin holds the single remaining bin
    assert bin_edges_in_window(hist, 101, 110, rebin=4) == [101.0, 105.0, 109.0, 110.0]


def test_window_edges_must_be_bin_edges(hist):
    with pytest.raises(ValueError):
        bin_edges_in_window(hist, 102.5, 110)
    with pytest.raises(ValueError):
        bin_edges_in_window(hist, 110, 102)


def test_window_hist_contents_and_errors(ROOT, hist):
    h = window_hist(ROOT, hist, 101, 110, rebin=4, name="w")
    # original bin i (1-based) holds i; window starts at bin 2 (101-102 GeV)
    assert [h.GetBinContent(i) for i in (1, 2, 3)] == [2 + 3 + 4 + 5, 6 + 7 + 8 + 9, 10]
    assert h.GetBinError(1) == pytest.approx(np.sqrt(2 + 3 + 4 + 5))
    assert h.GetXaxis().GetXmin() == 101 and h.GetXaxis().GetXmax() == 110
    assert h.Integral() == pytest.approx(sum(range(2, 11)))


def test_binning_label():
    assert binning_label([0.0, 2.0, 4.0]) == "2 bins of 2 GeV"
    assert "last bin narrower" in binning_label([0.0, 2.0, 3.0])
