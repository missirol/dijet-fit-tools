"""Unit tests for dijetfit.histtools: fit-range windows and rebinning."""
import numpy as np
import pytest

from conftest import make_th1
from dijetfit.histtools import bin_edges_in_window, binning_label, first_filled_edge, window_hist


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


def test_first_filled_edge(ROOT):
    """Lower edge of the first non-empty bin inside the window: xmin itself if that bin is filled,
    the edge of the first filled bin if lower bins are empty (a cut on m_jj above xmin), None if
    the window is empty."""
    edges = np.arange(100.0, 121.0, 1.0)
    contents = np.where(edges[:-1] >= 107.0, 5.0, 0.0)  # filled from 107 (a cut at 107)
    h = make_th1(ROOT, "h_cut", edges, contents)
    assert first_filled_edge(h, 108, 120) == 108.0
    assert first_filled_edge(h, 107, 120) == 107.0
    assert first_filled_edge(h, 100, 120) == 107.0
    assert first_filled_edge(h, 100, 105) is None

