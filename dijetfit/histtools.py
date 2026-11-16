"""Histogram windowing/rebinning helpers."""
from array import array


def bin_edges_in_window(hist, xmin, xmax, rebin=1):
    """Bin edges of `hist` restricted to [xmin, xmax], merging `rebin` consecutive bins
    starting from xmin. If the number of bins is not a multiple of `rebin`, the last
    (highest-mass) bin is narrower: it contains the remaining original bins.
    xmin and xmax must coincide with bin edges of `hist`."""
    ax = hist.GetXaxis()
    edges = [ax.GetBinLowEdge(i) for i in range(1, ax.GetNbins() + 2)]
    tol = 1e-6 * (edges[-1] - edges[0])

    def index(v):
        for i, e in enumerate(edges):
            if abs(e - v) < tol:
                return i
        raise ValueError(f"{v:g} is not a bin edge of '{hist.GetName()}' "
                         f"(range [{edges[0]:g}, {edges[-1]:g}])")

    i0, i1 = index(xmin), index(xmax)
    if i1 <= i0:
        raise ValueError(f"empty window [{xmin:g}, {xmax:g}]")
    out = edges[i0:i1 + 1:rebin]
    if abs(out[-1] - edges[i1]) > tol:
        out.append(edges[i1])
    return out


def first_filled_edge(hist, xmin, xmax):
    """Lower edge of the first non-empty bin of `hist` in [xmin, xmax] (None if all are empty).
    With xmin below a selection cut on m_jj (e.g. histograms filled from 170 GeV and xmin = 150),
    this is above xmin: the bins in between are empty because of the cut, not of the spectrum."""
    ax = hist.GetXaxis()
    for i in range(1, ax.GetNbins() + 1):
        lo, hi = ax.GetBinLowEdge(i), ax.GetBinUpEdge(i)
        if hi <= xmin + 1e-9 or lo < xmin - 1e-9:
            continue
        if lo >= xmax - 1e-9:
            break
        if hist.GetBinContent(i) != 0:
            return lo
    return None


def window_hist(ROOT, hist, xmin, xmax, rebin=1, name=None):
    """New TH1D with the content of `hist` in [xmin, xmax], bins merged by `rebin`
    (see bin_edges_in_window). Bin errors are combined in quadrature."""
    edges = bin_edges_in_window(hist, xmin, xmax, rebin)
    out = ROOT.TH1D(name or f"{hist.GetName()}_win", hist.GetTitle(), len(edges) - 1,
                    array("d", edges))
    out.SetDirectory(0)
    out.Sumw2()
    ax = hist.GetXaxis()
    for j in range(1, out.GetNbinsX() + 1):
        lo, hi = out.GetXaxis().GetBinLowEdge(j), out.GetXaxis().GetBinUpEdge(j)
        b0 = ax.FindFixBin(lo + 1e-6 * (hi - lo))
        b1 = ax.FindFixBin(hi - 1e-6 * (hi - lo))
        content = sum(hist.GetBinContent(b) for b in range(b0, b1 + 1))
        err2 = sum(hist.GetBinError(b) ** 2 for b in range(b0, b1 + 1))
        out.SetBinContent(j, content)
        out.SetBinError(j, err2 ** 0.5)
    return out


def binning_label(edges):
    widths = sorted({round(b - a, 6) for a, b in zip(edges[:-1], edges[1:])})
    if len(widths) == 1:
        return f"{len(edges) - 1} bins of {widths[0]:g} GeV"
    return (f"{len(edges) - 1} bins, widths {', '.join(f'{w:g}' for w in widths)} GeV "
            "(last bin narrower)")
