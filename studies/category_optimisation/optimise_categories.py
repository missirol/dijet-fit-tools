#!/usr/bin/env python3
"""Study: categories in |dy| and |dphi| of the two jets, compared by S/sqrt(B).

From a 2D histogram of |dy| (x) vs |dphi| (y) of the two leading jets, in the data file and in
the signal files (one per mass), this script compares ways of splitting the plane into
categories by the expected sensitivity

    Z(M) = sqrt( sum over categories of S(M)^2 / B ),

S = signal yield for sigma = 1 nb (x --lumi / eventsProcessed), B = data yield, both summed over
the 2D bins of the category (S << B). The upper bound Z_max(M) has every 2D bin as its own
category. The schemes compared (all with |dphi| > --dphi-min, as the current categories):

  A  the current categories: |dy| slices of --current-width up to --current-dy-max
  B  one category, |dy| < --current-dy-max
  C  |dy| slices with edges --wide-edges
  D  |dy| slices of 0.3 up to 2.1
  E  C, each split into back-to-back (pi - |dphi| < one bin of the histogram) and the rest
  F  D, split in the same way
  G  A, split in the same way
  R  the 2D bins ranked by the mass-averaged s/b (sum over masses of S(M)/S_tot(M), divided by
     B), grouped into N categories with the boundaries that maximise the mean of
     (Z/Z_max)^2 over the masses (dynamic programming on 300 groups of equal B); the same
     definition for every mass. The resulting regions follow the s/b contours (2D lookup).

The 2D histogram has no m_jj axis: B is the data yield for all m_jj, not the background under a
given signal peak. If --slice-hists is given (a regular expression matching the m_jj histograms
of |dy| slices in the data and signal files, with named groups `lo` and `hi`, e.g.
`...dY(?P<lo>\\d+p\\d+)to(?P<hi>\\d+p\\d+)_J1J2_mass`), the ratio (S/sqrt(B)) in the mass window
[--window] x M over (S/sqrt(B)) for all m_jj is measured per slice; its linear dependence on |dy|
(1 - k_M |dy|, relative to the first slice) is fitted and applied to S in every 2D bin, also
beyond the slices (an extrapolation). Results are given without and with this correction.

Outputs (in --outdir): summary.md (tables, ranked regions), results.json, and the figure
signal_regions.<fmt> (mass-averaged s/b map with the current and the C/E boundaries; gain of each
scheme over the current categories).

Example (the study of November 2026, on the second input set):
  optimise_categories.py --data L1ScoutingSelection_Run2026D.root \\
      --signal ZprimeToQQ_M{200,300,400,500}_TuneCP5_13p6TeV_pythia8.root \\
      --hist L1SCaloJet_Std_Mjj150_J1J2_deltaRapAbs_deltaPhiAbs \\
      --slice-hists 'L1SCaloJet_Std_Mjj150dPhi1p0dY(?P<lo>\\d+p\\d+)to(?P<hi>\\d+p\\d+)_J1J2_mass' \\
      --exclude-slices 0p0to1p1 -o category_study
"""
import argparse
import json
import re
import sys
from array import array
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from dijetfit.common import mass_from_filename, setup_root, write_json  # noqa: E402


def read_2d(ROOT, path, name):
    f = ROOT.TFile.Open(str(path))
    h = f.Get(name) if f else None
    if not h:
        sys.exit(f"ERROR: no histogram {name} in {path}")
    nx, ny = h.GetNbinsX(), h.GetNbinsY()
    a = np.array([[h.GetBinContent(i, j) for j in range(1, ny + 1)] for i in range(1, nx + 1)])
    xe = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, nx + 2)])
    ye = np.array([h.GetYaxis().GetBinLowEdge(j) for j in range(1, ny + 2)])
    ep = f.Get("eventsProcessed")
    return a, xe, ye, (ep.Integral() if ep else None)


def zvals(lab, B, S, masses):
    """Z(M) of the categories given by the labels (-1: not used)."""
    out = []
    for m in masses:
        z2 = 0.0
        for k in np.unique(lab[lab >= 0]):
            sel = lab == k
            b = B[sel].sum()
            if b > 0:
                z2 += S[m][sel].sum() ** 2 / b
        out.append(np.sqrt(z2))
    return np.array(out)


def ranked_groups(N, B, S, masses, zmax, ngroups=300):
    """Optimal split of the 2D bins, ranked by the mass-averaged s/b, into N categories."""
    ok = B > 0
    d = np.where(ok, sum(S[m] / S[m].sum() for m in masses) / np.where(ok, B, 1), 0.0)
    idx = np.flatnonzero(ok.ravel())
    order = idx[np.argsort(-d.ravel()[idx])]
    cb = np.cumsum(B.ravel()[order])
    grp = np.minimum((cb / cb[-1] * ngroups).astype(int), ngroups - 1)
    cB = np.concatenate([[0], np.cumsum(np.bincount(grp, B.ravel()[order], ngroups))])
    cS = {m: np.concatenate([[0], np.cumsum(np.bincount(grp, S[m].ravel()[order], ngroups))]) for m in masses}

    def seg(i, j):  # mean over masses of (Z/Z_max)^2 of the group [i, j)
        b = cB[j] - cB[i]
        return 0.0 if b <= 0 else float(np.mean([(cS[m][j] - cS[m][i]) ** 2 / b / zmax[m] ** 2 for m in masses]))

    G = ngroups
    best = np.full((N + 1, G + 1), -np.inf)
    arg = np.zeros((N + 1, G + 1), int)
    best[0, 0] = 0.0
    for k in range(1, N + 1):
        for j in range(1, G + 1):
            v = [best[k - 1, i] + seg(i, j) for i in range(j)]
            i = int(np.argmax(v))
            best[k, j], arg[k, j] = v[i], i
    j, cuts = int(np.argmax(best[N])), []
    for k in range(N, 0, -1):
        i = arg[k, j]
        cuts.append((i, j))
        j = i
    lab = -np.ones(B.size, int)
    for gi, (i, j) in enumerate(sorted(cuts)):
        lab[order[(grp >= i) & (grp < j)]] = gi
    return lab.reshape(B.shape)


def window_slopes(ROOT, args, data_path, sig_paths, masses):
    """k_M of (S/sqrtB)_window / (S/sqrtB)_all ~ 1 - k_M |dy| (relative to the first slice),
    measured on the |dy| slices matched by --slice-hists."""
    pat = re.compile(args.slice_hists)
    fd = ROOT.TFile.Open(str(data_path))
    names = [k.GetName() for k in fd.GetListOfKeys()]
    slices = []
    for n in names:
        mt = pat.fullmatch(n)
        if mt and not any(x in n for x in args.exclude_slices):
            lo, hi = (float(mt.group(g).replace("p", ".")) for g in ("lo", "hi"))
            slices.append((0.5 * (lo + hi), n))
    if len(slices) < 2:
        sys.exit(f"ERROR: --slice-hists matches fewer than 2 histograms in {data_path}")
    slices.sort()
    k, detail = {}, {}
    for m in masses:
        fs = ROOT.TFile.Open(str(sig_paths[m]))
        ratio = []
        for dy, n in slices:
            hb, hs = fd.Get(n), fs.Get(n)
            a, b = hb.FindBin(args.window[0] * m + 1e-6), hb.FindBin(args.window[1] * m - 1e-6)
            nb_all, ns_all = hb.Integral(0, hb.GetNbinsX() + 1), hs.Integral(0, hs.GetNbinsX() + 1)
            nb_w, ns_w = hb.Integral(a, b), hs.Integral(a, b)
            ratio.append((ns_w / np.sqrt(nb_w)) / (ns_all / np.sqrt(nb_all)))
        dys = np.array([s[0] for s in slices])
        rel = np.array(ratio) / ratio[0]
        # 1 - k (dy - dy0), least squares through (dy0, 1)
        x = dys - dys[0]
        k[m] = float(-np.sum(x * (rel - 1)) / np.sum(x * x))
        detail[m] = dict(dy=dys.tolist(), relative_ratio=rel.tolist())
    return k, detail


def plot(ROOT, path, formats, B, S, xe, ye, masses, table, current, wide):
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    s = sum(S[m] / S[m].sum() for m in masses) / len(masses)
    with np.errstate(all="ignore"):
        r = np.where(B > 0, s / (B / B.sum()), 0.0)
    nx = min(len(xe) - 1, int(np.searchsorted(xe, 3.0 - 1e-9)) + 1)
    ny_all = len(ye) - 1
    ny = min(ny_all, int(np.searchsorted(ye, 1.26 - 1e-9)) + 1)
    keep = []
    c = ROOT.TCanvas("c_sr", "", 1600, 700)
    c.Divide(2, 1)
    p = c.cd(1)
    p.SetRightMargin(0.15)
    p.SetLeftMargin(0.12)
    u_edges = np.pi - ye[::-1]  # pi - |dphi|, increasing
    h = ROOT.TH2D("h_sr", "mass-averaged signal / data density ratio;|#Deltay|;#pi - |#Delta#phi|",
                  nx, array("d", xe[:nx + 1]), ny, array("d", u_edges[:ny + 1]))
    for i in range(nx):
        for j in range(ny):
            h.SetBinContent(i + 1, j + 1, r[i, ny_all - 1 - j])
    h.Draw("COLZ")
    keep.append(h)

    def line(x1, y1, x2, y2, col, w=3, st=1):
        ln = ROOT.TLine(x1, y1, x2, y2)
        ln.SetLineColor(col)
        ln.SetLineWidth(w)
        ln.SetLineStyle(st)
        ln.Draw()
        keep.append(ln)
    top = u_edges[ny]
    for x in wide[1:]:
        if x <= xe[nx]:
            line(x, 0, x, top, ROOT.kRed + 1)
    line(0, u_edges[1], min(wide[-1], xe[nx]), u_edges[1], ROOT.kRed + 1)
    dymax, width, _ = current
    for x in np.arange(width, dymax + 1e-9, width):
        line(x, 0, x, top, ROOT.kGray + 2, 1, 2)
    t = ROOT.TLatex()
    t.SetTextSize(0.03)
    t.SetTextColor(ROOT.kRed + 1)
    t.DrawLatexNDC(0.14, 0.92, "red: scheme E (|#Deltay| slices #times back-to-back split)")
    t.SetTextColor(ROOT.kBlack)
    t.DrawLatexNDC(0.14, 0.88, "dashed: current |#Deltay| slices")
    p = c.cd(2)
    p.SetLeftMargin(0.12)
    p.SetGridy()
    rel = [row["corrected" if "corrected" in row else "raw"]["rel"] for row in table.values()]
    lo = min(min(v) for v in rel)
    hi = max(max(v) for v in rel)
    fr = p.DrawFrame(min(masses) - 50, min(0.9, lo - 0.02), max(masses) + 50, hi + 0.08,
                     "sensitivity relative to the current categories;signal mass [GeV];Z / Z_{current}")
    keep.append(fr)
    leg = ROOT.TLegend(0.14, 0.62, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.026)
    cols = [ROOT.kBlack, ROOT.kGray + 1, ROOT.kAzure + 1, ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kRed - 7,
            ROOT.kOrange + 1, ROOT.kGreen + 1, ROOT.kGreen + 2, ROOT.kGreen + 3, ROOT.kTeal + 2, ROOT.kCyan + 2]
    for k, (name, row) in enumerate(table.items()):
        v = row["corrected" if "corrected" in row else "raw"]["rel"]
        g = ROOT.TGraph(len(masses), array("d", masses), array("d", v))
        g.SetLineColor(cols[k % len(cols)])
        g.SetMarkerColor(cols[k % len(cols)])
        g.SetMarkerStyle(20)
        g.SetLineWidth(3 if name.startswith("E") else 2)
        g.Draw("LP same")
        leg.AddEntry(g, f"{name} ({row['N']} cat.)", "lp")
        keep.append(g)
    leg.Draw()
    keep.append(leg)
    for fmt in formats:
        c.SaveAs(f"{path}.{fmt}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="data ROOT file")
    ap.add_argument("--signal", nargs="+", required=True, help="signal ROOT files (mass from the file name)")
    ap.add_argument("--hist", required=True, help="2D histogram: x = |dy|, y = |dphi|")
    ap.add_argument("--lumi", type=float, default=680000.0, help="nb^-1 (default: 680000)")
    ap.add_argument("--dphi-min", type=float, default=1.0, help="|dphi| cut of all schemes (default: 1.0)")
    ap.add_argument("--current", nargs=2, type=float, default=[1.1, 0.1], metavar=("DY_MAX", "WIDTH"),
                    help="current categories: |dy| slices of WIDTH up to DY_MAX (default: 1.1 0.1)")
    ap.add_argument("--wide-edges", nargs="+", type=float, default=[0, 0.5, 1.0, 1.5, 2.0, 2.5],
                    help="|dy| edges of schemes C and E (default: 0 0.5 1 1.5 2 2.5)")
    ap.add_argument("--ranked", nargs="+", type=int, default=[3, 5, 8, 11],
                    help="numbers of categories of the ranked scheme R (default: 3 5 8 11)")
    ap.add_argument("--slice-hists", help="regex of the m_jj histograms of |dy| slices (groups lo, hi) "
                                          "for the mass-window correction")
    ap.add_argument("--exclude-slices", nargs="*", default=[], help="substrings of slice names to skip")
    ap.add_argument("--window", nargs=2, type=float, default=[0.85, 1.10], metavar=("LO", "HI"),
                    help="mass window in units of M for the correction (default: 0.85 1.10)")
    ap.add_argument("-o", "--outdir", default="category_study")
    ap.add_argument("--formats", nargs="+", default=["png", "pdf"])
    args = ap.parse_args()

    ROOT = setup_root()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    B, xe, ye, _ = read_2d(ROOT, args.data, args.hist)
    S0, sig_paths = {}, {}
    for f in args.signal:
        m = mass_from_filename(f)
        a, xe2, ye2, ep = read_2d(ROOT, f, args.hist)
        if not (np.allclose(xe, xe2) and np.allclose(ye, ye2)):
            sys.exit(f"ERROR: {f}: binning of {args.hist} differs from the data")
        if not ep:
            sys.exit(f"ERROR: {f}: no eventsProcessed histogram")
        S0[m], sig_paths[m] = a * args.lumi / ep, f
    masses = sorted(S0)
    xc, yc = 0.5 * (xe[:-1] + xe[1:]), 0.5 * (ye[:-1] + ye[1:])
    DY, DPHI = np.meshgrid(xc, yc, indexing="ij")
    U = np.pi - DPHI
    b2b = (U < (ye[-1] - ye[-2])).astype(int)  # the last |dphi| bin: back-to-back
    dymax, width = args.current
    use = DPHI > args.dphi_min

    def slices(edges):
        return np.where((DY < edges[-1]) & use, np.digitize(DY, edges) - 1, -1)

    wide = list(args.wide_edges)
    schemes = {
        f"A current: |dy| slices of {width:g} up to {dymax:g}": np.where((DY < dymax) & use, np.floor(DY / width + 1e-9), -1),
        f"B one category |dy| < {dymax:g}": np.where((DY < dymax) & use, 0, -1),
        f"C |dy| edges {wide}": slices(wide),
        "D |dy| slices of 0.3 up to 2.1": np.where((DY < 2.1) & use, np.floor(DY / 0.3 + 1e-9), -1),
    }
    schemes["E C x back-to-back split"] = np.where(schemes[f"C |dy| edges {wide}"] >= 0,
                                                   schemes[f"C |dy| edges {wide}"] * 2 + b2b, -1)
    schemes["F D x back-to-back split"] = np.where(schemes["D |dy| slices of 0.3 up to 2.1"] >= 0,
                                                   schemes["D |dy| slices of 0.3 up to 2.1"] * 2 + b2b, -1)
    cur = list(schemes)[0]
    schemes["G A x back-to-back split"] = np.where(schemes[cur] >= 0, schemes[cur] * 2 + b2b, -1)
    schemes = {k: v.astype(int) for k, v in schemes.items()}

    ok = B > 0
    zmax0 = {m: float(np.sqrt((S0[m][ok] ** 2 / B[ok]).sum())) for m in masses}
    ranked = {}
    for N in args.ranked:
        lab = ranked_groups(N, B, S0, masses, zmax0)
        schemes[f"R ranked s/b, {N} groups"] = lab
        ranked[N] = lab

    corrections = {"raw": {m: 0.0 for m in masses}}
    window_detail = None
    if args.slice_hists:
        k, window_detail = window_slopes(ROOT, args, args.data, sig_paths, masses)
        corrections["corrected"] = k
    table = {name: {"N": int(len(np.unique(lab[lab >= 0])))} for name, lab in schemes.items()}
    zmax_all = {}
    for kind, k in corrections.items():
        S = {m: S0[m] * np.clip(1 - k[m] * DY, 0, None) for m in masses}
        zmax = {m: float(np.sqrt((S[m][ok] ** 2 / B[ok]).sum())) for m in masses}
        zmax_all[kind] = zmax
        zcur = zvals(schemes[cur], B, S, masses)
        for name, lab in schemes.items():
            z = zvals(lab, B, S, masses)
            table[name][kind] = dict(z=z.tolist(), frac=(z / np.array([zmax[m] for m in masses])).tolist(),
                                     rel=(z / zcur).tolist())

    md = ["# Categories in |dy| and |dphi|: S/sqrt(B) study\n",
          f"Data: `{args.data}`; signal: {', '.join(f'M{m}' for m in masses)} (sigma = 1 nb); histogram "
          f"`{args.hist}`; B = data yield for all m_jj; |dphi| > {args.dphi_min:g} in every scheme.\n"]
    for kind in corrections:
        title = "without" if kind == "raw" else "with"
        md += [f"## Z / Z_max and Z / Z_current, {title} the mass-window correction along |dy|\n",
               "| scheme | N | " + " | ".join(f"Z/Zmax M{m}" for m in masses) + " | "
               + " | ".join(f"Z/Zcur M{m}" for m in masses) + " |",
               "|---|---|" + "---|" * (2 * len(masses))]
        for name, row in table.items():
            r = row[kind]
            md.append(f"| {name} | {row['N']} | " + " | ".join(f"{x:.3f}" for x in r["frac"]) + " | "
                      + " | ".join(f"{x:.3f}" for x in r["rel"]) + " |")
        md.append("")
    if window_detail:
        md += ["Mass-window correction: k_M (S/sqrt(B) in the window relative to all m_jj, falls as "
               "1 - k_M |dy|): " + ", ".join(f"M{m}: {corrections['corrected'][m]:.3f}" for m in masses) + "\n"]
    for N, lab in ranked.items():
        md += [f"## Ranked scheme, {N} groups: category of each 2D bin\n",
               "Rows: |dy| bins from 0 (up to 3); columns: pi - |dphi| from 0, one character per bin; "
               "'.': not used.\n", "```"]
        ny_all = len(ye) - 1
        for i in range(min(len(xc), int(np.searchsorted(xe, 3.0 - 1e-9)))):
            md.append(f"{xe[i]:4.1f} " + "".join(
                "." if lab[i, ny_all - 1 - j] < 0 else "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[lab[i, ny_all - 1 - j] % 36]
                for j in range(min(16, ny_all))))
        md.append("```\n")
    (out / "summary.md").write_text("\n".join(md) + "\n")
    write_json(out / "results.json", dict(
        data=str(args.data), hist=args.hist, masses=masses, lumi=args.lumi, zmax=zmax_all,
        window_correction=dict(slopes=corrections.get("corrected"), detail=window_detail, window=args.window),
        schemes=table, ranked_labels={str(N): lab.tolist() for N, lab in ranked.items()}))
    plot(ROOT, str(out / "signal_regions"), args.formats, B, S0, xe, ye, masses, table,
         (dymax, width, args.dphi_min), wide)
    print("\n".join(md[:2 + 4 + len(table)]))
    print(f"Written {out / 'summary.md'}, {out / 'results.json'}")


if __name__ == "__main__":
    main()
