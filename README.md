# Dijet resonance search: signal and background fits

Builds the inputs of the statistical analysis of a search for dijet resonances (Z'→qq) in CMS
data: parametric fits of the m_jj spectrum in 11 categories of the rapidity separation |Δy|
of the two jets, written to a RooWorkspace (format: [docs/OUTPUT_FORMAT.md](docs/OUTPUT_FORMAT.md)).

## Setup

ROOT 6.36.14 (PyROOT, RooFit), numpy, and CMS Combine (only for `RooMultiPdf` in the last
step). `setup/install_env.sh` creates such an environment from conda-forge.

## Workflow

```
./run_all.sh DATA.root ZprimeToQQ_M200_*.root ZprimeToQQ_M300_*.root ...
```

runs the five steps below (outputs in `output/`, logs `output/<N>_<step>.log`). Each step is
a separate script with `--help`; every step reads only the outputs of the previous ones, so
any step can be re-run on its own.

| step | script | does | main outputs |
|---|---|---|---|
| 1 | `scripts/prepare_inputs.py` | collects the m_jj histograms of the categories selected by regular expressions (`--categories`, `--exclude`; the category name is the first capture group); data unchanged, signal scaled to σ = 1 nb: × `--lumi` (nb⁻¹, default 680000) / integral of `eventsProcessed` | `inputs.root` with `CAT/data_obs`, `CAT/ZprimeQQ_<MASS>` |
| 2 | `scripts/fit_background.py` | fits candidate background functions to `data_obs` per category; selects the passing ones; with `--scan`, chooses the common fit range (and, if needed, per-category rebinning) | `background_fits.json`, `background_summary.md`, plots |
| 3 | `scripts/fit_signal.py` | fits each signal template with a double-sided Crystal Ball, in the same range/binning | `signal_fits.json`, `signal_summary.md`, plots |
| 4 | `scripts/parametrize_signal.py` | interpolates the signal parameters and yield in the mass (spline by default, `--default-model polN` / `--model PARAM=polN` for polynomials) | `signal_parametrisation.json`, plots |
| 5 | `scripts/build_workspace.py` | writes the workspace and checks that it reproduces the fit results | `workspace/dijet_workspace.root` (+ `.json` summary) |

Every fit produces a plot with the distribution, the fit and a pull sub-pad,
(data − fit)/σ_data. RooFit/Minuit warnings and errors are counted per fit and reported in
the plots and the summaries.

### Categories (step 1)

The scripts contain no category names. `prepare_inputs.py --categories REGEX [REGEX ...]`
selects the histograms whose whole name matches one of the patterns (and none of
`--exclude`); the first capture group of the matching pattern names the category. The
analysis defaults, set in `run_all.sh` (`CATEGORIES`, `EXCLUDE`), select the 11 |Δy| slices
`L1SCaloJet_Std_Mjj150dPhi1p0dY<A>to<B>_J1J2_mass` and leave out the inclusive `0p0to1p1`.
All later steps work on whatever categories `inputs.root` contains.

### Fit range and binning (step 2)

The fit range is the same for all categories and is always given on the command line
(`--xmin/--xmax`); the scripts have no built-in range. The analysis defaults live in
`run_all.sh` (150–600 GeV; scan limits 190 and 250 GeV), overridable with environment
variables. With `--scan` the lower edge is raised in 5 GeV steps (up to `--scan-max-unbinned`)
until every category has at least one passing function; then failing categories are rebinned
(2, 5, 10 bins merged), and only then is the lower edge raised further (up to `--scan-max`).
Rebinning is applied only where needed and is logged. Signal mass points closer than 10 GeV
to the lower edge are dropped (with the defaults, M200 is dropped above 190 GeV).

### Background functions and selection

Families `dijet` (CMS standard dijet function), `moddijet` ((1−x^{1/3}) variant) and `ua2`,
with 3, 4, 5 parameters (incl. normalisation); 6-parameter variants can be requested with
`--functions`. A function passes if the fit has status 0, covariance quality 3, no parameter
at a limit, goodness-of-fit p > 0.05, and its order is not above the order chosen by the Fisher
F-test in its family. All passing functions enter the `RooMultiPdf` of the category.

### Signal fits

Double-sided `RooCrystalBall` (common σ). The fit releases the parameters in stages
(core, then tail onsets α, then tail exponents n); a tail parameter that the template does not
determine (the tail lies outside the fit range, or n runs to its limit, i.e. the tail is
exponential) is fixed and recorded in `signal_fits.json`. The templates are weighted (MC counts
× L / N_processed); by default they are fitted as weighted histograms with SumW2-corrected
uncertainties. `--weights counts` fits the MC counts instead (Poisson errors) and scales the
yield back, which is equivalent for uniform weights.

## Inputs

ROOT input files are **never committed** to this repository (`*.root` is git-ignored).
