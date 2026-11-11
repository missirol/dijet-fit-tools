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
(2, 5, 10 bins merged, smallest factor that works), and only then is the lower edge raised
further (up to `--scan-max`), each category restarting from no rebinning at every step.
Rebinning is applied only where needed and is logged. Note that `build_workspace.py` refuses
categories with bins wider than 1 GeV, because Combine evaluates the pdfs at the bin centres
(see `docs/OUTPUT_FORMAT.md`). A category without any passing function stops the workflow,
unless `DROP_FAILING_CATEGORIES=1` (`fit_background.py --drop-failing-categories`) is set, in
which case it is left out and listed in `background_fits.json`. Signal mass points closer than
10 GeV to either edge of the fit range (or outside it) are dropped (with the defaults, M200 is
dropped above 190 GeV).

### Background functions and selection

Families `dijet` (CMS standard dijet function), `moddijet` ((1−x^{1/3}) variant) and `ua2`,
with 3, 4, 5, 6 parameters (incl. normalisation; `--functions` selects others). A fit is
acceptable if it has status 0, covariance quality 3, no parameter at a limit and goodness-of-fit
p > 0.05 (`--pvalue-threshold`); only acceptable fits are used. In each family, a sequence of
standard Fisher F-tests fixes the highest order: starting from the lowest acceptable order, the
next order is taken while it is significantly better (F-test p-value < 0.05,
`--ftest-threshold`) and its fit is acceptable. The `RooMultiPdf` of the category contains:

* the default: among the families' highest orders, the one with the largest p-value;
* up to 4 alternates (`--max-alternates`): the other acceptable functions, of order at most
  their family's highest order, ranked by p-value.

### Signal fits

Double-sided `RooCrystalBall` (common σ). The fit releases the parameters in stages
(core, then tail onsets α, then tail exponents n); a tail parameter that the template does not
determine (the tail lies outside the fit range, or n runs to its limit, i.e. the tail is
exponential) is fixed and recorded in `signal_fits.json`. The templates are weighted (MC counts
× L / N_processed); by default they are fitted as weighted histograms with SumW2-corrected
uncertainties; for the fit the weights are rescaled to a mean of 1 (and the yield scaled back),
because Minuit's convergence criterion depends on the overall scale of the likelihood.
`--weights counts` fits the MC counts instead (Poisson errors) and scales the yield back, which
is equivalent for uniform weights.

A signal fit is *failed* if it has status != 0, covariance quality != 3 or a parameter at its
limit (after the fixing rules); failed fits are kept in `signal_fits.json` but not used for
the parametrisation. The goodness-of-fit p-value is reported as a quality figure and does not
reject fits. In the parametrisation, each parameter is interpolated only through the mass points
where it was floating in the fit (values fixed by the procedure are not measurements), or
through all points if fewer than two of them float.

## Tests

Run inside the analysis environment, from the repository root:

```
pytest                                   # everything, fit range 150-600 GeV (~1 min)
pytest -m "not closure"                  # unit tests only (~10 s)
pytest --fit-xmin 160 --fit-xmax 650     # another fit range (inside 100-800 GeV)
```

The fit range used by all tests is a command-line option (`--fit-xmin`, `--fit-xmax`;
default 150-600 GeV, as in `run_all.sh`).

* **Unit tests** (`tests/test_*.py`) check the building blocks against exact answers: category
  selection, fit-range windows and rebinning, splines and their RooFit formulas, the background
  functions (formulas, decorrelated basis, starting values), the goodness-of-fit and F-test
  helpers, the Crystal Ball definition, and fits of Asimov histograms (contents equal to the
  expectation), which must return the generating parameters. One of them guards against the
  bias of bin-centre evaluation for wide bins.
* **Closure test** (`tests/test_closure.py`): `tests/synthetic.py` generates inputs in the real
  format from known functions (dijet background, Crystal Ball signal at 200-500 GeV, and an
  extra 350 GeV sample kept out of the fits), with settings that differ on purpose from the
  real ones (category names, histogram range, luminosity). `run_all.sh` runs on them with the
  configured fit range, and every step is compared with the truth: input normalisation and
  category selection, the fit range (must stay as configured), the fit of the generating
  background function, the signal fits (likelihood ratio against the truth), the interpolated
  signal model at 350 GeV, the workspace content, one plot per fit, and (if CMS Combine is
  installed) that Combine reads the workspace and that `-m` sets the mass.
* **Scan test** (also in `tests/test_closure.py`): a background with a turn-on just above the
  configured lower edge; `fit_background.py --scan` must raise the lower edge past it.

The synthetic inputs can also be generated on their own: `python tests/synthetic.py -o DIR
[--turnon-mid M]` prints the matching `run_all.sh` command.

## Inputs

ROOT input files are **never committed** to this repository (`*.root` is git-ignored).
