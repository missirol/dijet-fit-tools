#!/bin/bash
# Run one analysis case end to end, resumably: the workflow (run_all.sh: steps 1-7, including the
# expected limits) and then a short bias test (step 8). Run it inside the analysis environment:
#
#   setup/envrun.sh studies/run_case.sh OUTDIR DATA.root SIGNAL_M200_*.root SIGNAL_M300_*.root [...]
#
# Every part leaves a marker in OUTDIR/.run_case/ when it finishes; rerunning the same command
# after an interruption (e.g. a time limit on background jobs) skips the finished parts, and the
# bias test continues from its finished Combine jobs (bias_test.py --resume). The configuration
# (inputs and the settings below) is saved in OUTDIR/case_config.txt; a rerun with a different
# configuration in the same OUTDIR is refused (use another OUTDIR).
#
# Settings (environment variables):
#   run_all.sh settings, passed on unchanged: CATEGORIES, EXCLUDE, XMIN, XMAX, SCAN_MAX, LUMI,
#     LIMIT_STEP, FORMATS (default here: png), REBIN_FACTORS, SCAN_MAX_UNBINNED
#   BIAS               1 (default) to run the bias test, 0 to skip it
#   BIAS_MASSES        signal masses                    (default: all signal masses of the case)
#   BIAS_MU            injected signal strengths       (default: "1")
#   BIAS_NTOYS         toys per mass and signal strength (default: 100)
#   BIAS_GENERATE      generating functions            (default: "default", the default function
#                      of every category; "alt1 ... alt4" etc. multiply the run time)
#   BIAS_TOYS_PER_JOB  toys per Combine job             (default: 10; small jobs resume well)
#   BIAS_JOBS          Combine jobs in parallel         (default: number of CPUs)
#   BIAS_SEED          random seed                      (default: 2000)
#
# With the defaults, a case of about 10 categories takes about 1 hour on 4 CPUs, 20 categories
# about 2-3 hours. Outputs: those of run_all.sh in OUTDIR, the bias test in OUTDIR/bias, logs
# OUTDIR/<N>_<step>.log and OUTDIR/run_all.log.
set -euo pipefail

if [ $# -lt 4 ]; then
  echo "usage: $0 OUTDIR DATA.root SIGNAL.root SIGNAL.root [SIGNAL.root ...]" >&2
  exit 1
fi
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTDIR="$1"; shift
case "$OUTDIR" in /*) ;; *) OUTDIR="$PWD/$OUTDIR" ;; esac
INPUTS=()
for f in "$@"; do
  [ -f "$f" ] || { echo "ERROR: no such file: $f" >&2; exit 1; }
  INPUTS+=("$(cd "$(dirname "$f")" && pwd)/$(basename "$f")")
done

export FORMATS="${FORMATS:-png}"
BIAS="${BIAS:-1}"
BIAS_MASSES="${BIAS_MASSES:-}"
BIAS_MU="${BIAS_MU:-1}"
BIAS_NTOYS="${BIAS_NTOYS:-100}"
BIAS_GENERATE="${BIAS_GENERATE:-default}"
BIAS_TOYS_PER_JOB="${BIAS_TOYS_PER_JOB:-10}"
BIAS_JOBS="${BIAS_JOBS:-$(nproc 2>/dev/null || echo 1)}"
BIAS_SEED="${BIAS_SEED:-2000}"

mkdir -p "$OUTDIR/.run_case"
# configuration of this case (BIAS_JOBS only sets the parallelism, it does not change results)
config() {
  echo "inputs: ${INPUTS[*]}"
  for v in CATEGORIES EXCLUDE XMIN XMAX SCAN_MAX LUMI LIMIT_STEP FORMATS REBIN_FACTORS \
           SCAN_MAX_UNBINNED BIAS BIAS_MASSES BIAS_MU BIAS_NTOYS BIAS_GENERATE BIAS_TOYS_PER_JOB BIAS_SEED; do
    if [ -n "${!v+set}" ]; then echo "$v=${!v}"; else echo "$v (run_all.sh default)"; fi
  done
}
if [ -f "$OUTDIR/case_config.txt" ]; then
  if ! diff -q <(config) "$OUTDIR/case_config.txt" >/dev/null; then
    echo "ERROR: $OUTDIR was run with another configuration:" >&2
    diff <(config) "$OUTDIR/case_config.txt" >&2 || true
    echo "Use another OUTDIR (or remove this one)." >&2
    exit 1
  fi
else
  config > "$OUTDIR/case_config.txt"
fi

if [ -f "$OUTDIR/.run_case/workflow" ]; then
  echo "workflow (steps 1-7): done before, skipped"
else
  echo "workflow (steps 1-7): running, log $OUTDIR/run_all.log"
  if ! OUT="$OUTDIR" "$HERE/run_all.sh" "${INPUTS[@]}" > "$OUTDIR/run_all.log" 2>&1; then
    echo "ERROR: run_all.sh failed (the analysis case may fail; see $OUTDIR/run_all.log)" >&2
    exit 1
  fi
  touch "$OUTDIR/.run_case/workflow"
fi

if [ "$BIAS" = "0" ]; then
  echo "bias test (step 8): skipped (BIAS=0)"
elif [ -f "$OUTDIR/.run_case/bias" ]; then
  echo "bias test (step 8): done before, skipped"
else
  echo "bias test (step 8): running, log $OUTDIR/8_bias_test.log"
  read -r -a MU <<< "$BIAS_MU"
  read -r -a GEN <<< "$BIAS_GENERATE"
  read -r -a FMT <<< "$FORMATS"
  MASS_ARGS=()
  if [ -n "$BIAS_MASSES" ]; then read -r -a M <<< "$BIAS_MASSES"; MASS_ARGS=(--masses "${M[@]}"); fi
  if ! python3 "$HERE/scripts/bias_test.py" "$OUTDIR/datacards/datacards.json" -o "$OUTDIR/bias" \
      ${MASS_ARGS[@]+"${MASS_ARGS[@]}"} --mu "${MU[@]}" --ntoys "$BIAS_NTOYS" --generate "${GEN[@]}" --fits envelope fixed \
      --toys-per-job "$BIAS_TOYS_PER_JOB" --jobs "$BIAS_JOBS" --seed "$BIAS_SEED" \
      --formats "${FMT[@]}" --resume >> "$OUTDIR/8_bias_test.log" 2>&1; then
    echo "ERROR: bias_test.py failed (see $OUTDIR/8_bias_test.log)" >&2
    exit 1
  fi
  touch "$OUTDIR/.run_case/bias"
fi

echo "Case done: $OUTDIR"
echo "  fit range:      $(grep -h 'Selected fit range\|Fixed fit range' "$OUTDIR/2_fit_background.log" | tail -1)"
echo "  limits:         $OUTDIR/limits/limits_summary.md"
[ "$BIAS" = "0" ] || echo "  bias test:      $OUTDIR/bias/bias_summary.md"
