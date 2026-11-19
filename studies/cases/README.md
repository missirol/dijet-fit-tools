# Analysis cases: how to configure and run them

An **analysis case** is a set of input files (data and signal) plus a choice of categories and
settings. The rules that decide whether a case works (and when it fails) are in
[docs/ANALYSIS_CONVENTIONS.md](../../docs/ANALYSIS_CONVENTIONS.md). This page explains how to run
one. The results of analysis cases are not kept in the repository.

## Inputs

- One data file and at least two signal files (the signal model is interpolated between the
  mass points), all with the same m_jj histograms, 1 GeV wide. The signal mass is read from the
  signal file name, which must contain `_M<mass>_`; each signal file has the histogram
  `eventsProcessed` used for the normalisation to σ = 1 nb.
- ROOT files are never committed; keep them outside the repository.
- Before running a case:
  - list the histograms of the data file and choose the regular expressions of the categories;
  - check from which m_jj value the histograms are filled (a cut on m_jj in their production):
    the lower edge of the fit range cannot be below it (`fit_background.py` refuses it);
  - if the data are identical to those of a case already run (e.g. the same histograms under
    another name), reuse its results instead of running again.

## Settings

The settings are the environment variables of `run_all.sh` (see its header):

| variable | meaning | default |
|---|---|---|
| `CATEGORIES` | regular expressions of the category histograms (first capture group = category name) | the 11 \|Δy\| slices of `run_all.sh` |
| `EXCLUDE` | regular expressions of histograms to leave out | the inclusive \|Δy\| slice |
| `XMIN`, `XMAX` | fit range [GeV]; the lower edge is raised in 5 GeV steps if needed | 150, 600 |
| `SCAN_MAX` | highest lower edge allowed [GeV] | 200 |
| `LUMI` | integrated luminosity [nb⁻¹] | 680000 |
| `LIMIT_STEP` | spacing of the mass grid of the expected limits [GeV] | 25 |

and those of the bias test in `studies/run_case.sh` (`BIAS_MU`, `BIAS_NTOYS`, `BIAS_GENERATE`,
`BIAS_MASSES`, `BIAS=0` to skip it; see its header).

## Running

```
CATEGORIES='(L1SCaloJet_Std_Opt1dRap\d+p\d+to\d+p\d+dPhi\d+p\d+to\d+p\d+)_J1J2_mass' EXCLUDE="" \
XMIN=170 XMAX=600 setup/envrun.sh studies/run_case.sh out/opt1 \
  DATA.root SIGNAL_M200_*.root SIGNAL_M300_*.root SIGNAL_M400_*.root SIGNAL_M500_*.root
```

runs the workflow (`run_all.sh`, steps 1–7 including the expected limits) and a short bias test
(100 toys per mass at μ = 1, generated with the default function of every category, fitted with
the envelope and with the generating function fixed). With about 10 categories this takes about
an hour on 4 CPUs; toys from every alternate function (`BIAS_GENERATE="alt1 alt2 alt3 alt4"`)
or more toys multiply the time.

The same command continues an interrupted run (finished parts are skipped, the bias test resumes
from its finished Combine jobs). The configuration is saved in `OUTDIR/case_config.txt`; a rerun
with other settings needs another output directory.

## What to look at

- `2_fit_background.log` / `background/background_summary.md`: the lower edge selected by the
  scan (and the categories that forced it up), the passing functions and the default per
  category. The case fails if a category has no passing function at the highest allowed edge.
- `signal/signal_summary.md`: the signal fits (every one must pass) and the masses kept.
- `limits/limits_summary.md`, `limits/expected_limits.png`: the expected limits.
- `bias/bias_summary.md`: the bias of the envelope and of the fixed-function fits (|bias| > 0.5 is
  reported with a warning), the failed toy fits, and how often the envelope picks the
  generating function.

When comparing cases, compare at the same masses, and keep in mind that a different selected fit
range changes the limits and can drop the 200 GeV mass point.
