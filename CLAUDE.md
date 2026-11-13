# Notes for Claude Code

- Project: fits for a CMS dijet resonance search (Z'->qq), producing a RooWorkspace for Combine.
  See README.md (workflow) and docs/OUTPUT_FORMAT.md (final output format; keep it up to date).
- **Never commit ROOT files.** Input ROOT files are provided by the user via the prompt only.
- Do not commit or push without the user's explicit request.
- ROOT version: 6.36.14 (conda-forge), with CMS Combine 11.1.0 for RooMultiPdf. Run python
  inside the conda environment (setup/shell.sh, `source setup/activate.sh` or
  `setup/envrun.sh COMMAND`; these remove a ROOT set up in the calling shell); see setup/.
- Analysis conventions agreed with the user:
  - Categories are selected by regular expressions given to prepare_inputs.py (no category
    names hard-coded in the scripts; defaults in run_all.sh): currently the 11 slices
    `L1SCaloJet_Std_Mjj150dPhi1p0dY<A>to<B>` (inclusive 0p0to1p1 excluded), fitted independently.
  - Signal normalised to 1 nb: x lumi (nb^-1, default 680000) / integral of `eventsProcessed`.
  - One fit range for all categories, default 150-600 GeV. The range must stay configurable:
    no fit-range numbers hard-coded in the scripts (defaults only in run_all.sh).
    Raise the lower edge in 5 GeV steps if background fits fail, at most up to 200 GeV
    (SCAN_MAX in run_all.sh). Rebinning is switched off by
    default (the workspace only supports 1 GeV bins); it stays available for studies
    (--rebin-factors), per category and logged. Drop M200 if the lower edge goes above 190 GeV.
  - Background fits must reach p > 5% (plus status 0, covQual 3, no parameter at a limit);
    only such fits are used. Families: those of the analysis note's background table (dijet,
    moddijet, ua2, modexp, polypow, polyext; plus ua23, not in the note). Per family (its orders),
    standard F-tests (5%) from the lowest
    acceptable order give the highest order. Default = family maximum with the largest p-value;
    up to 4 alternates = other acceptable orders up to their family maximum, ranked by p-value
    (RooMultiPdf: up to 5 functions). RooFit warnings do not by themselves fail a fit; no cut
    on pulls.
  - An analysis case (a set of inputs and categories) fails if any category has no passing
    background function at the highest allowed lower edge, or if any signal fit fails in any
    category. Categories are never dropped; a failing case needs another categorisation.
  - Yields always refer to the fit range.
  - Signal fits: weighted templates with SumW2 errors by default; the MC-counts approach
    ("--weights counts") stays available as an option.
  - A signal fit is failed if status != 0, covQual != 3 or a floating parameter ends at its
    limit (tail parameters that the procedure deliberately fixes at a limit do not count); every
    signal fit must pass (see "analysis case fails"). Signal p-values are a quality figure, not a cut.
  - Signal parameters vs mass: splines by default, polynomials via command line; each parameter
    is interpolated only through the mass points where it was floating in the fit, and held
    constant beyond the outermost such points.
  - Combine: one datacard per category, everything (data, background model, signal) taken from
    the workspace. All statistical tests use the combination of the 11 categories (exclusive
    data sets): background parameters uncorrelated across categories, only the signal strength
    of a given signal common to all of them.
  - Bias tests: toys with injected signal strengths per signal mass (e.g. mu = 0.5, 1, 5),
    fitted with background + signal; check that the fits return the injected mu on average.
    Criterion |bias| < 0.5 sigma (bias = (<r> - mu) / std(r)).
- Every fit needs a plot with a pull sub-pad, pull = (data - fit) / data_error.
- Tests: run `pytest` (inside the environment) after changes; extend the unit tests and the
  closure test (tests/synthetic.py, tests/test_closure.py) together with new features.
