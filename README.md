# Dijet resonance search: signal and background fits

Builds the inputs of the statistical analysis of a search for dijet resonances (Z'→qq) in CMS
data: parametric fits of the m_jj spectrum in 11 categories of the rapidity separation |Δy|
of the two jets, written to a RooWorkspace (format: [docs/OUTPUT_FORMAT.md](docs/OUTPUT_FORMAT.md)).

## Setup

ROOT 6.36.14 (PyROOT, RooFit), numpy, and CMS Combine 11.1.0 (`RooMultiPdf` in step 5, the
datacards and statistical tests in steps 6-7). `setup/install_env.sh [PREFIX]` creates such
an environment from conda-forge (default `<repo>/env`, with micromamba in `<repo>/bin`; Linux
x86-64; it can be run from any directory). Then work in it in one of three ways (all commands
in this README assume one of them):

| command | does |
|---|---|
| `setup/shell.sh [PREFIX]` | opens a new bash shell (reading your `~/.bashrc`) with the environment active; leave it with `exit` |
| `source setup/activate.sh` | activates the environment in the current bash or zsh shell (it takes no arguments; another environment: `DIJET_ENV=PREFIX source setup/activate.sh`) |
| `setup/envrun.sh COMMAND ...` | runs one command in the environment, e.g. `setup/envrun.sh pytest` |

A non-default environment is given with `DIJET_ENV=PREFIX` (or as the argument of `shell.sh`).
`activate.sh` changes the current shell: besides activating the environment it sets
`MAMBA_ROOT_PREFIX` (if unset) and defines micromamba's shell function. All three first remove
any other ROOT installation set up in the calling shell (e.g. by `thisroot.sh`): `ROOTSYS`,
`PYTHONPATH`, `PYTHONHOME`, `LD_LIBRARY_PATH`, `DYLD_LIBRARY_PATH`, `LD_PRELOAD`,
`ROOT_INCLUDE_PATH`, `CPLUS_INCLUDE_PATH`, `C_INCLUDE_PATH`, `CLING_STANDARD_PCH`, and ROOT
`bin` directories in `PATH` (list in `setup/_clean_env.sh`); a line says what was removed
(settings of the analysis environment itself, when it is already active, are not reported).
Otherwise the environment's Python can load the other ROOT's libraries or bindings, and
`import ROOT` crashes.

## Workflow

```
./run_all.sh DATA.root ZprimeToQQ_M200_*.root ZprimeToQQ_M300_*.root ...
```

runs steps 1-6 below (outputs in `output/`, logs `output/<N>_<step>.log`). Each step is
a separate script with `--help`; every step reads only the outputs of the previous ones, so
any step can be re-run on its own. The bias test (step 7) is CPU-intensive and is run
separately (see [Bias test](#bias-test-step-7)).

| step | script | does | main outputs |
|---|---|---|---|
| 1 | `scripts/prepare_inputs.py` | collects the m_jj histograms of the categories selected by regular expressions (`--categories`, `--exclude`; the category name is the first capture group); data unchanged, signal scaled to σ = 1 nb: × `--lumi` (nb⁻¹, default 680000) / integral of `eventsProcessed` | `inputs.root` with `CAT/data_obs`, `CAT/ZprimeQQ_<MASS>` |
| 2 | `scripts/fit_background.py` | fits candidate background functions to `data_obs` per category; selects the passing ones; with `--scan`, chooses the common fit range (per-category rebinning only for studies, with `--rebin-factors`; `--rebin CAT=N` fixes it without a scan) | `background_fits.json`, `background_summary.md`, plots |
| 3 | `scripts/fit_signal.py` | fits each signal template with a double-sided Crystal Ball, in the same range/binning | `signal_fits.json`, `signal_summary.md`, plots |
| 4 | `scripts/parametrize_signal.py` | interpolates the signal parameters and yield in the mass (spline by default, `--default-model polN` / `--model PARAM=polN` for polynomials) | `signal_parametrisation.json`, plots |
| 5 | `scripts/build_workspace.py` | writes the workspace and checks that it reproduces the fit results | `workspace/dijet_workspace.root` (+ `.json` summary) |
| 6 | `scripts/make_datacards.py` | one Combine datacard per category (everything taken from the workspace), their combination, and `text2workspace.py`; checks the expected yields and the normalised shapes of the combined workspace against the source workspace | `datacards/datacard_CAT.txt`, `datacard_combined.{txt,root}`, `datacards.json` |
| 7 | `scripts/bias_test.py` | bias test of the signal strength with Combine toys, all categories combined | `bias/bias_results.json`, `bias_summary.md`, plots |

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
`run_all.sh` (150–600 GeV, scan up to a lower edge of 200 GeV, the highest allowed),
overridable with environment variables. With `--scan` the lower edge is raised in 5 GeV steps
(up to `--scan-max`) until every category has at least one passing function. Every step is
logged.

Rebinning is switched off by default: `build_workspace.py` refuses categories with bins wider
than 1 GeV, because the fits integrate the pdfs over such bins while Combine evaluates them at
the bin centres (see `docs/OUTPUT_FORMAT.md`).
For studies it can be switched on with `REBIN_FACTORS` (`--rebin-factors`, e.g. `2 5 10`): the
failing categories are then rebinned at `SCAN_MAX_UNBINNED` (default 190 GeV) with the smallest
factor that works, before the lower edge is raised further.

**An analysis case fails** (the workflow stops with an error) if a category still has no
passing background function at the highest allowed lower edge (`SCAN_MAX`, 200 GeV), or if any
signal fit fails in any category. Categories are never dropped: all categories of a case must
work, otherwise a different categorisation is needed. Signal mass points closer than
10 GeV to either edge of the fit range (or outside it) are dropped (with the defaults, M200 is
dropped above 190 GeV).

### Background functions and selection

The function families of the analysis note's background table, with x = m_jj/√s and the
number of parameters including the normalisation p0 (`--functions` selects a subset):

| family | parameters | form |
|---|---|---|
| `dijet` | 3–6 | p0 (1−x)^p1 / x^(p2 + p3 ln x + p4 ln²x + p5 ln³x) |
| `moddijet` | 3–6 | p0 (1−x^{1/3})^p1 / x^(p2 + p3 ln x + p4 ln²x + p5 ln³x) |
| `ua2` | 3–6 | p0 x^−p1 exp(−(p2 x + p3 x² + p4 x³ + p5 x⁴)) |
| `modexp` | 3–5 | p0 exp(p1 x^p2); p0 exp(p1 x^p2 + p1 (1−x)^p3); p0 exp(p1 x^p2 + p3 (1−x)^p4) |
| `polypow` | 3–5 | p0 (1 + p1 x + p2 x² + p3 x³)^−p_last |
| `polyext` | 5–6 | p0 (1−x)^p1 (1 + p4 x + p5 x²) x^−(p2 + p3 ln x) |

`ua23` (x^−p1 exp(−p2 x)) is not in the note's table, whose `ua2` family starts at 4
parameters; it is kept so that every family of this form starts at 3 parameters. The
`modexp`, `polypow` and `polyext` functions are not log-linear in their parameters: their
starting values come from a search over the parameters that enter non-linearly (the others
solved exactly at each point) followed by a full Poisson fit (`dijetfit/background.py`). A fit is
acceptable if it has status 0, covariance quality 3, no parameter at a limit and goodness-of-fit
p > 0.05 (`--pvalue-threshold`); only acceptable fits are used. A shape parameter that ends at
its limit has its range re-centred and is fitted again (up to 3 times) before it counts as at
its limit; the yield may vary between 0.5 and 1.5 times the events in the fit range. In each family, a sequence of
standard Fisher F-tests fixes the highest order: starting from the lowest acceptable order, the
next order is taken while it is significantly better (F-test p-value < 0.05,
`--ftest-threshold`) and its fit is acceptable. The `RooMultiPdf` of the category contains:

* the default: among the families' highest orders, the one with the largest p-value;
* up to 4 alternates (`--max-alternates`): the other acceptable functions, of order at most
  their family's highest order, ranked by p-value.

With `--alternates family-max`, only the highest order of each family is used (at most one
function per family). It was tried against the envelope bias of the signal strength and did
not reduce it (bias test of November 2026, results in the discussion of PR #2).

With 1 GeV bins the background and signal pdfs are evaluated at the bin centres and normalised
by their sum over the bin centres (not by their integral), in the fits and in the workspace:
the expected counts per bin then add up exactly to the yield, as the extended likelihood of
RooFit and Combine assumes. With the integral as normalisation they differ by O(10⁻⁵)
(background) and O(10⁻⁶) (signal); with ~2×10⁷ background events per category this pulled
the background parameters by O(1σ) and the signal strength in Combine by up to ~2σ (background)
and 0.3σ at 200 GeV (signal). The bin-sum normalised pdfs agree with exact bin integrals (the
model of Combine's RooParametricShapeBinPdf) to ~10⁻⁶ relative, ~10⁻⁴ of the statistical
uncertainty per bin.

`demo/bin_normalisation/` is a standalone workflow (artificial data with a known truth) that
shows how the integral normalisation biases the fits, and compares the bin-sum normalisation
with exact bin integration (see its README).

`studies/category_optimisation/` compares category definitions in |Δy| and |Δφ| of the two jets
by S/√B, from 2D histograms of the data and signal files (see its README); it is a study, not a
step of the workflow.

### Signal fits

Double-sided `RooCrystalBall` (common σ). The fit releases the parameters in stages
(core, then tail onsets α, then tail exponents n); a tail parameter that the template does not
determine (the tail lies outside the fit range, or n runs to its limit, i.e. the tail is
exponential) is fixed and recorded in `signal_fits.json`. The templates are weighted (MC counts
× L / N_processed); by default they are fitted as weighted histograms with SumW2-corrected
uncertainties; for the fit the weights are rescaled to a mean of 1 (and the yield scaled back),
because Minuit's convergence criterion depends on the overall scale of the likelihood.
`--weights counts` fits the MC counts instead (Poisson errors) and scales the yield back, which
is equivalent for uniform weights; for templates with non-uniform weights it falls back to the
default (the mode actually used is `weights_mode` of each fit in `signal_fits.json`). Two further modes exist for comparisons: `asymptotic`
(RooFit's AsymptoticError; RooFit warns that it is not designed for binned data) and `none`
(no error correction; uncertainties wrong by about 1/√w).

A signal fit is *failed* if it has status != 0, covariance quality != 3 or a parameter at its
limit (after the fixing rules). Every signal fit must pass in every category: a failed fit is
recorded in `signal_fits.json` (`usable: false`) and the workflow stops (the analysis case
fails). The goodness-of-fit p-value is reported as a quality figure and does not
reject fits. In the parametrisation, each parameter is interpolated only through the mass points
where it was floating in the fit (values fixed by the procedure are not measurements), or
through all points if fewer than two of them float.

## Bias test (step 7)

```
python3 scripts/bias_test.py output/datacards/datacards.json -o output/bias \
  [--masses 200 300 400 500] [--mu 0 0.5 1 5] [--ntoys 300] [--pull-toys 25 --pull-mu 1]
```

For every signal mass and injected signal strength μ (× σ = 1 nb), Combine generates binned toys
from the combined workspace (background: the default function of every category, with the
parameters of the fit to data; `--generate altK` uses the K-th alternates instead, `--generate family-NAME`
the highest order of family NAME in every category; a category with fewer alternates, or
without that family, uses its default) and fits them with signal + background, all categories
together (r common, background parameters uncorrelated):

* `envelope`: the background function is profiled (discrete profiling), as in the analysis;
* `fixed`: the background is the generating function (reference).

The bias is (⟨r_fit⟩ − μ) / std(r_fit) over the toys (uncertainty 1/√N). A bias above 0.5 in
absolute value (`--bias-threshold`) is shown in the summary plot (band), printed as a warning,
flagged in the tables (`above_threshold`) and listed at the end of `bias_summary.md`
(`points_above_threshold` in `bias_results.json`); it does not fail the bias test. For
reference: added in quadrature, a bias b increases the total uncertainty by √(1 + b²) − 1, i.e.
12% for b = 0.5 and 1% for b = 0.14 (a stricter threshold used in some analyses); a 68% interval
shifted by b covers the true value in 62.5% (b = 0.5) or 67.8% (b = 0.14) of cases. The expected uncertainty σ_r at each mass is measured first (Asimov, μ = 1,
`--robustFit`). With `--pull-toys K`, the first K toys of each point (taken from as many jobs
as needed) are also fitted with uncertainties (`--robustFit 1`, about a minute per toy) and the
pulls (r_fit − μ)/σ_r are reported.

Every job (`--toys-per-job` toys) has its own random seed. A toy fit counts as successful only if
it has a result in Combine's output tree, and the Combine log of every job is checked: the
numbers of successful (`n_reported_ok`) and failed fits (`n_failed_fit`) are stored per point,
together with toys lost in Combine jobs that failed (`n_lost_in_failed_jobs`; the run continues
and the failed jobs are listed in `failed_jobs`) and any other missing result
(`n_failed_unexplained`, which should be 0). `bias_results.json` also stores, per toy, r and,
for the envelope fits, the selected function index of every category (`envelope_index`, in the
order of `categories`).

An interrupted run (for example a machine that is shut down) is continued by running the same
command again with `--resume`: every finished Combine job leaves a marker (`done.json`) with its
settings, and jobs whose settings match are reused instead of being run again (also the Asimov
fits). The run then completes as if it had not been interrupted; `bias_results.json` records
how many jobs were reused (`resumed`).

## Tests

Run inside the analysis environment, from the repository root:

```
pytest                                   # everything, fit range 150-600 GeV (~8 min)
pytest -m "not closure and not slow"     # unit tests only (~1 min)
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
  signal model at 350 GeV, the workspace content, one plot per fit, and the Combine steps
  (below).
* **Combine tests** (also in `tests/test_closure.py`): the datacards of step 6; a Combine fit
  to its own Asimov data set (generating background function) returns the injected signal
  strength within 2% of its uncertainty, for no signal and 3σ (regression test for the
  normalisation of the background pdfs: with the old normalisation it was off by 1.2σ); `-m`
  sets the mass; the shape check of `make_datacards.py` catches a modified parameter;
  `bias_test.py` runs end to end (two generators, several jobs: every fit counted in the
  Combine log, unique seeds, pull toys from several jobs, fixed-function fits close) and keeps
  going when a Combine job fails.
* **Demo, study and setup tests** (`tests/test_demo.py`, `tests/test_studies.py`,
  `tests/test_setup.py`; marked `slow`): the bin-normalisation demo runs and shows the expected
  biases; the category study runs on synthetic histograms (no scheme above the upper bound);
  the setup scripts keep a ROOT installation of the calling shell out of the environment (with a
  fake one that fails if loaded).
* **Scan test** (also in `tests/test_closure.py`): a background with a turn-on just above the
  configured lower edge; `fit_background.py --scan` must raise the lower edge past it.

The synthetic inputs can also be generated on their own: `python tests/synthetic.py -o DIR
[--turnon-mid M]` prints the matching `run_all.sh` command.

## Inputs

ROOT input files are **never committed** to this repository (`*.root` is git-ignored).
