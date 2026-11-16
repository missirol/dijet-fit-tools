# Categories in |Δy| and |Δφ|: S/√B study

`optimise_categories.py` compares ways of defining the analysis categories in the plane of
|Δy| and |Δφ| of the two leading jets, by the expected sensitivity Z = √(Σ S²/B) per signal
mass (S for σ = 1 nb, B = data yield; see the script's docstring and `--help`). It is a study,
not a step of the workflow: the categories of the analysis are still chosen with the regular
expressions of `prepare_inputs.py`, and new categories need their own m_jj histograms.

Run it inside the analysis environment, e.g. on the second input set of November 2026:

```
python studies/category_optimisation/optimise_categories.py \
  --data L1ScoutingSelection_Run2026D.root \
  --signal ZprimeToQQ_M200_TuneCP5_13p6TeV_pythia8.root ZprimeToQQ_M300_TuneCP5_13p6TeV_pythia8.root \
           ZprimeToQQ_M400_TuneCP5_13p6TeV_pythia8.root ZprimeToQQ_M500_TuneCP5_13p6TeV_pythia8.root \
  --hist L1SCaloJet_Std_Mjj150_J1J2_deltaRapAbs_deltaPhiAbs \
  --slice-hists 'L1SCaloJet_Std_Mjj150dPhi1p0dY(?P<lo>\d+p\d+)to(?P<hi>\d+p\d+)_J1J2_mass' \
  --exclude-slices 0p0to1p1 -o category_study
```

Outputs: `summary.md` (tables, regions of the ranked scheme), `results.json`, and the figure
`signal_regions.png` (mass-averaged s/b map with the boundaries; gain of each scheme).

## Results on that input set (with the mass-window correction)

| scheme | categories | Z / Z_current, M200 · M300 · M400 · M500 | Z / Z_max |
|---|---|---|---|
| current: 0.1 slices of |Δy| < 1.1 | 11 | 1 · 1 · 1 · 1 | 0.85–0.91 |
| one category |Δy| < 1.1 | 1 | 0.98 for all masses | |
| |Δy| edges 0, 0.5, 1, 1.5, 2, 2.5 | 5 | 1.08 · 1.09 · 1.09 · 1.11 | 0.95–0.98 |
| the same × back-to-back split | 10 | 1.09 · 1.10 · 1.12 · 1.15 | 0.98–0.99 |
| 2D bins ranked by s/b, 5 groups | 5 | 1.09 · 1.11 · 1.13 · 1.15 | 0.98–0.99 |

- The fine 0.1 slicing of |Δy| < 1.1 adds only 2% over a single category.
- Most of the gain comes from |Δy| between 1.1 and about 2.5.
- A separate back-to-back category (|Δφ| > π − π/40, the last |Δφ| bin) adds 2–5%.
- One mass-independent definition reaches 98–99% of each mass's own optimum.

## Limitations

- **B includes every m_jj.** The 2D histogram has no m_jj axis, so B counts data events at all
  masses, not the background under each signal peak.
  - Along |Δy|, the per-slice m_jj histograms give a correction for this (`--slice-hists`). It
    is at most 9% in S/√B for |Δy| < 1.1, and it is extrapolated linearly beyond that range.
  - Along |Δφ| there is no such check. The region of small |Δy| and |Δφ| below about 2.6 has
    high s/b in this proxy; that may be an artefact, and about 5% of Z_max comes from it.
- **The back-to-back cut is limited to the histogram's binning** (π/40 per |Δφ| bin).
- **The m_jj turn-on moves up at large |Δy|** (m_jj ≈ 2 pT cosh(Δy/2)). Since all categories
  share one fit range, check the m_jj spectra of new categories before adopting them.
- **S/√B ignores the m_jj shape.** It is a guide for comparing schemes, not a limit.
