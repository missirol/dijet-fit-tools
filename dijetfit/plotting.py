"""Fit plot: histogram + fitted curve, with a pull sub-pad (data - fit) / data_error."""
from array import array


def fit_plot(ROOT, path_stem, formats, hist, expected, curve, labels, logy=False,
             data_label="data", fit_label="fit", curve_color=None):
    """Draw `hist` (with its bin errors) and the fit.

    expected : list of fitted expectations per bin (same binning as `hist`); the pulls are
               (content - expected) / error, for bins with error > 0.
    curve    : (xs, ys) of the fitted function as event density per GeV; it is drawn scaled
               to events per bin for the *typical* bin width when all bins have the same
               width, otherwise the histogram is drawn as events / GeV.
    labels   : lines of text for the information box.
    """
    nb = hist.GetNbinsX()
    widths = [hist.GetXaxis().GetBinWidth(i) for i in range(1, nb + 1)]
    uniform = max(widths) - min(widths) < 1e-6 * max(widths)
    hdraw = hist.Clone(hist.GetName() + "_draw")
    hdraw.SetDirectory(0)
    if uniform:
        scale = widths[0]
        ytitle = f"Events / {widths[0]:g} GeV"
    else:
        scale = 1.0
        ytitle = "Events / GeV"
        for i in range(1, nb + 1):
            hdraw.SetBinContent(i, hist.GetBinContent(i) / widths[i - 1])
            hdraw.SetBinError(i, hist.GetBinError(i) / widths[i - 1])
    xs, ys = curve
    g = ROOT.TGraph(len(xs), array("d", xs), array("d", [y * scale for y in ys]))

    pulls = ROOT.TH1D(hist.GetName() + "_pull", "", nb,
                      array("d", [hist.GetXaxis().GetBinLowEdge(i) for i in range(1, nb + 2)]))
    pulls.SetDirectory(0)
    for i in range(1, nb + 1):
        err = hist.GetBinError(i)
        if err > 0:
            pulls.SetBinContent(i, (hist.GetBinContent(i) - expected[i - 1]) / err)

    c = ROOT.TCanvas("c_" + hist.GetName(), "", 800, 800)
    top = ROOT.TPad("top", "", 0, 0.3, 1, 1)
    bot = ROOT.TPad("bot", "", 0, 0, 1, 0.3)
    for p in (top, bot):
        p.SetLeftMargin(0.14)
        p.SetRightMargin(0.04)
    top.SetBottomMargin(0.02)
    bot.SetTopMargin(0.03)
    bot.SetBottomMargin(0.33)
    bot.SetGridy(True)
    top.Draw()
    bot.Draw()

    top.cd()
    top.SetLogy(logy)
    hdraw.SetMarkerStyle(20)
    hdraw.SetMarkerSize(0.5 if nb > 100 else 0.8)
    hdraw.SetLineColor(ROOT.kBlack)
    hdraw.GetYaxis().SetTitle(ytitle)
    hdraw.GetYaxis().SetTitleSize(0.05)
    hdraw.GetYaxis().SetLabelSize(0.045)
    hdraw.GetYaxis().SetTitleOffset(1.35)
    hdraw.GetXaxis().SetLabelSize(0)
    if logy:
        positive = [hdraw.GetBinContent(i) for i in range(1, nb + 1) if hdraw.GetBinContent(i) > 0]
        hdraw.SetMinimum(0.3 * min(positive) if positive else 0.1)
        hdraw.SetMaximum(3 * hdraw.GetMaximum())
    else:
        hdraw.SetMinimum(0)
        hdraw.SetMaximum(1.8 * hdraw.GetMaximum())
    hdraw.Draw("E1 X0")
    g.SetLineColor(curve_color or ROOT.kRed + 1)
    g.SetLineWidth(2)
    g.Draw("L same")

    leg = ROOT.TLegend(0.55, 0.80, 0.94, 0.92)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    leg.AddEntry(hdraw, data_label, "lep")
    leg.AddEntry(g, fit_label, "l")
    leg.Draw()

    box = ROOT.TPaveText(0.55, 0.78 - 0.034 * len(labels), 0.94, 0.79, "NDC")
    box.SetBorderSize(0)
    box.SetFillStyle(0)
    box.SetTextAlign(12)
    box.SetTextSize(0.028)
    box.SetTextFont(42)
    for t in labels:
        box.AddText(t)
    box.Draw()

    bot.cd()
    pulls.SetFillColor(ROOT.kAzure - 9)
    pulls.SetLineColor(ROOT.kAzure + 2)
    pulls.GetYaxis().SetTitle("#frac{data - fit}{#sigma_{data}}")
    pulls.GetYaxis().SetRangeUser(-4.9, 4.9)
    pulls.GetYaxis().SetNdivisions(505)
    pulls.GetYaxis().SetTitleSize(0.10)
    pulls.GetYaxis().SetTitleOffset(0.55)
    pulls.GetYaxis().SetLabelSize(0.10)
    pulls.GetXaxis().SetTitle("m_{jj} [GeV]")
    pulls.GetXaxis().SetTitleSize(0.12)
    pulls.GetXaxis().SetLabelSize(0.10)
    pulls.GetXaxis().SetTitleOffset(1.1)
    pulls.Draw("HIST")

    c.Update()
    for fmt in formats:
        c.SaveAs(f"{path_stem}.{fmt}")
    c.Close()
