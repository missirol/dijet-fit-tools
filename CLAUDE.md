# Notes for Claude Code

Fits for a CMS dijet resonance search (Z' → qq): background and signal fits of m_jj spectra per
category, a RooWorkspace and Combine datacards, expected limits and bias tests.

- [README.md](README.md): the workflow (`run_all.sh`, steps 1–8), the scripts, the tests.
- [docs/ANALYSIS_CONVENTIONS.md](docs/ANALYSIS_CONVENTIONS.md): the analysis decisions agreed with
  the user (fit range, background selection, failure rules, statistical tests). **Binding**:
  change them only when the user asks, together with the code, the tests and that file.
- [docs/OUTPUT_FORMAT.md](docs/OUTPUT_FORMAT.md): the final output format; keep it up to date.
- [studies/cases/README.md](studies/cases/README.md): how to configure and run an analysis case.
- The repository holds the code, the instructions to run it and the analysis criteria behind
  it, **not the results of analysis cases**: report results to the user (and as PR comments
  when asked), do not commit them.

## Repository and git

- Repository: `missirol/dijet-fit-tools`, default branch `main`.
- Do not commit or push without the user's explicit request. Never push to `main`: one branch
  per topic, started from the latest `main`, and a pull request to `main`. After a PR is merged,
  start the next topic on a new branch from the updated `main`.
- **Never commit ROOT files.** Input ROOT files are provided by the user via the prompt only.
- Commit authorship: the author (and committer) of every commit is the repository's owner, with
  their git credentials (name and email) as used in the repository history; Claude is credited
  only through the `Co-Authored-By:` trailer. Do not rely on the environment's git identity
  (cloud sessions are set up as Claude's own): take the owner's identity from the history, e.g.
  the author of the first commit of the default branch
  (`git log --reverse --format='%an <%ae>' origin/main | head -1`; in a shallow clone, where the
  history is cut, the author of the pull-request merge at the tip of `main`, i.e. a commit with
  committer `GitHub`: `git log -1 --format='%an <%ae> (committer %cn)' origin/main`), set it for
  the clone (`git config user.name ...`, `git config user.email ...`)
  or pass `--author`, and check with `git log --format='%an <%ae>'` before pushing.
- Pull requests: describe what changed and why, how it was tested, and the open items; post
  results of studies as PR comments (tables, not plots of the whole output).

## Environment and tests

- ROOT 6.36.14 (conda-forge) with CMS Combine 11.1.0 (RooMultiPdf, statistical tests). Create
  the environment with `setup/install_env.sh [PREFIX]` (several minutes); run python inside it
  with `setup/envrun.sh COMMAND`, `source setup/activate.sh` or `setup/shell.sh` (these remove a
  ROOT set up in the calling shell). With a PREFIX outside the repository, set
  `DIJET_ENV=PREFIX`.
- Tests: `setup/envrun.sh pytest` after every change (all tests, about 10 minutes);
  `pytest -m "not closure and not slow"` for the unit tests only (about 1 minute). Extend the
  unit tests and the closure test (`tests/synthetic.py`, `tests/test_closure.py`) together with
  new features.

## Running analysis cases

- Inputs: keep each set of uploaded files in its own directory outside the repository, renamed
  to their original names (without the prefix the upload adds; the signal file names contain
  `_M<mass>_`, from which the mass is read). Do not write the input file names in the
  repository. Before running, list the histograms, check
  where the m_jj histograms start being filled (XMIN must not be below it) and whether the data
  are identical to an earlier set (then reuse its results instead of rerunning).
- Run a case with `studies/run_case.sh` (workflow, expected limits, short bias test; resumable).
  Name cases as the user does (e.g. "Case 3", "Opt1").
- Time budget: about 4 hours per case. The short bias test (100 toys per mass at μ = 1, default
  generating function, envelope vs fixed fit) is the default; toys from every alternate function
  or with pulls multiply the run time and are run only when the user asks.
- Cloud sessions stop background jobs after 2 hours, and the machine has few CPUs (4 so far):
  run long jobs in the background with the longest timeout, rerun the same resumable command
  when one is stopped, schedule a check-in for when it should finish, and do not run two
  CPU-heavy studies at full parallelism at the same time.
- Report a case with: the fit range selected by the scan and why, the functions per category,
  the signal fits, the expected limits and the bias, compared with the earlier cases.

## Known pitfalls

- Combine `--freezeParameters rgx{pdfindex_.*}` crashes (`parseRegex`): list the parameter
  names explicitly (comma-separated; `datacards.json` has them).
- With the envelope, uncertainties on r need `--robustFit 1`; Hesse/Minos errors are wrong (the
  function index changes during the error scan).
- Expected limits with the envelope can differ by about 1% between identical runs.
- Combine prints stale "best fit" values after a failed fit; count fits from the output trees
  (as `bias_test.py` does), not from the printout.
- Input files have non-histogram keys (e.g. `TNamed`); check the class before reading a key as
  a directory or histogram.
- More categories do not by themselves give better limits: each category brings its own
  background parameters and envelope function choice.
