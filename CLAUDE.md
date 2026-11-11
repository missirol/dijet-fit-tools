# Notes for Claude Code

- Project: fits for a CMS dijet resonance search (Z'->qq), producing a RooWorkspace for Combine.
  See README.md (workflow) and docs/OUTPUT_FORMAT.md (final output format; keep it up to date).
- **Never commit ROOT files.** Input ROOT files are provided by the user via the prompt only.
- Do not commit or push without the user's explicit request.
- ROOT version: 6.36.14 (conda-forge), with CMS Combine 11.1.0 for RooMultiPdf. Run python
  inside the conda environment (e.g. `micromamba run -p ENV python ...`); see setup/.
- Analysis conventions agreed with the user:
  - Categories are selected by regular expressions given to prepare_inputs.py (no category
    names hard-coded in the scripts; defaults in run_all.sh): currently the 11 slices
    `L1SCaloJet_Std_Mjj150dPhi1p0dY<A>to<B>` (inclusive 0p0to1p1 excluded), fitted independently.
  - Signal normalised to 1 nb: x lumi (nb^-1, default 680000) / integral of `eventsProcessed`.
  - One fit range for all categories, default 150-600 GeV. The range must stay configurable:
    no fit-range numbers hard-coded in the scripts (defaults only in run_all.sh).
    Raise the lower edge in 5 GeV steps
    if background fits fail; rebin only if needed, per category, and log it. Drop M200 if the
    lower edge goes above 190 GeV.
  - Background fits must reach p > 5% (plus status 0, covQual 3, no parameter at a limit). Per
    function family, the acceptable function with the fewest parameters goes into the
    RooMultiPdf; the one with the highest p-value is the default. RooFit warnings do not by
    themselves fail a fit; no cut on pulls.
  - Yields always refer to the fit range.
  - Signal fits: weighted templates with SumW2 errors by default; the MC-counts approach
    ("--weights counts") stays available as an option.
  - A signal fit is failed if status != 0, covQual != 3 or a parameter at its limit; failed fits
    are not used for the parametrisation. Signal p-values are a quality figure, not a cut.
  - Signal parameters vs mass: splines by default, polynomials via command line; each parameter
    is interpolated only through the mass points where it was floating in the fit.
- Every fit needs a plot with a pull sub-pad, pull = (data - fit) / data_error.
- Tests: run `pytest` (inside the environment) after changes; extend the unit tests and the
  closure test (tests/synthetic.py, tests/test_closure.py) together with new features.
