# Analysis conventions

Decisions agreed for this analysis. The code implements them; change them only together with the
code, the tests and this file. How the workflow implements them is described in
[README.md](../README.md); the format of the final outputs in [OUTPUT_FORMAT.md](OUTPUT_FORMAT.md).

## Inputs and categories

- An **analysis case** is a set of input files (data and signal ROOT files) plus a choice of
  categories. Each category is an m_jj histogram; the categories of a case are exclusive data
  sets and are fitted independently.
- Categories are selected by regular expressions given to `prepare_inputs.py` (`--categories`,
  `--exclude`). No category names are hard-coded in the scripts; the default selection lives in
  `run_all.sh` only (currently the 11 |Δy| slices `L1SCaloJet_Std_Mjj150dPhi1p0dY<A>to<B>`, without
  the inclusive slice). Other categorisations are run by overriding `CATEGORIES`/`EXCLUDE` (see
  `studies/cases/`).
- Signal normalised to σ = 1 nb: × luminosity (nb⁻¹, default 680000) / integral of
  `eventsProcessed`.
- Yields always refer to the fit range.

## Fit range

- One fit range for all categories of a case, default 150–600 GeV. The range must stay
  configurable: no fit-range numbers in the scripts (defaults only in `run_all.sh`).
- If background fits fail, the lower edge is raised in 5 GeV steps, at most up to 200 GeV
  (`SCAN_MAX` in `run_all.sh`).
- The data must be filled from the lower edge up: the lower edge cannot be below a cut on m_jj
  applied to the input histograms (`fit_background.py` refuses it).
- Signal mass points closer than 10 GeV to either edge are dropped (M200 is dropped if the lower
  edge goes above 190 GeV).
- Rebinning is switched off (the workspace only supports 1 GeV bins); it stays available for
  studies (`--rebin-factors`), per category and logged.

## Background model

- A background fit is acceptable only with p > 5%, fit status 0, covariance quality 3 and no
  parameter at a limit. RooFit warnings do not by themselves fail a fit; there is no cut on pulls.
- Families: those of the analysis note's background table (dijet, moddijet, ua2, modexp,
  polypow, polyext), plus ua23 (not in the note).
- Per family (its orders), standard F-tests (5%) from the lowest acceptable order give the
  highest order. Default function = the family maximum with the largest p-value; up to 4
  alternates = the other acceptable orders up to their family maximum, ranked by p-value
  (RooMultiPdf: up to 5 functions).

## Signal model

- Signal fits: weighted templates with SumW2 errors by default; the MC-counts approach
  (`--weights counts`) stays available as an option.
- A signal fit fails if status ≠ 0, covariance quality ≠ 3, or a floating parameter ends at its
  limit (tail parameters that the procedure deliberately fixes at a limit do not count). Signal
  p-values are a quality figure, not a cut.
- Signal parameters vs mass: splines by default, polynomials via the command line. Each parameter
  is interpolated only through the mass points where it was floating in the fit (with fewer than
  two such points: through all points), and held constant beyond the outermost points used.

## When an analysis case fails

A case fails if any category has no passing background function at the highest allowed lower
edge, or if any signal fit fails in any category. Categories are never dropped: a failing case
needs another categorisation.

## Statistical tests (CMS Combine)

- One datacard per category; data, background model and signal are all taken from the workspace.
- All statistical tests use the combination of all categories of the case: background
  parameters, yields and function indices uncorrelated across categories; only the signal
  strength of a given signal is common to all of them.
- Expected limits: 95% CL upper limits on r, asymptotic CLs, blind (Asimov data of the
  background-only hypothesis), discrete profiling of the background function.
- Bias tests: toys with injected signal strengths per signal mass, generated from the default
  function of every category, fitted with background + signal, both with the envelope (discrete
  profiling) and with the generating function fixed. bias = (⟨r⟩ − μ) / std(r); |bias| > 0.5 is
  shown in the plots and reported with a warning, but does not fail a bias test.

## Plots

Every fit has a plot with a pull sub-pad, pull = (data − fit) / data error.
