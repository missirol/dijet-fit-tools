#!/bin/bash
# Create a conda environment with ROOT 6.36.14 and CMS Combine 11.1.0 (from conda-forge),
# using micromamba (downloaded to <repo>/bin; Linux x86-64). Usage, from any directory:
#   setup/install_env.sh [PREFIX]   (default: <repo>/env)
#
# Afterwards, run the scripts inside the environment (ROOT's interpreter needs the
# environment's compiler settings), in one of three ways:
#   setup/shell.sh                  new shell with the environment active
#   source setup/activate.sh        activate it in the current shell (bash or zsh)
#   setup/envrun.sh COMMAND ...     run a single command in it
# All three first remove another ROOT installation set up in the calling shell (ROOTSYS,
# PYTHONPATH, LD_LIBRARY_PATH, ...; see setup/_clean_env.sh); with a non-default PREFIX,
# pass it to shell.sh, or set DIJET_ENV=PREFIX (activate.sh, envrun.sh).
#
# Combine is only needed for RooMultiPdf in the final workspace (build_workspace.py); without
# it, use build_workspace.py --no-multipdf.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"  # works from any directory
PREFIX="${1:-$REPO/env}"
case "$PREFIX" in /*) ;; *) PREFIX="$PWD/$PREFIX" ;; esac
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$REPO/.mamba}"
mkdir -p "$REPO/bin"
if [ ! -x "$REPO/bin/micromamba" ]; then
  curl -sSfL -o "$REPO/bin/micromamba" \
    https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
  chmod +x "$REPO/bin/micromamba"
fi
"$REPO/bin/micromamba" create -y -p "$PREFIX" -c conda-forge \
  "root=6.36.14" "cms-combine=11.1.0" "python=3.12" numpy pytest
DIJET_ENV="$PREFIX" "$REPO/setup/envrun.sh" python -c \
  "import ROOT; print('ROOT', ROOT.gROOT.GetVersion()); print('Combine:', ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit') >= 0)"
if [ "$PREFIX" = "$REPO/env" ]; then
  echo "Done. Use the environment with:  setup/shell.sh   (or: source setup/activate.sh)"
else
  echo "Done. Use the environment with:  setup/shell.sh $PREFIX   (or: DIJET_ENV=$PREFIX source setup/activate.sh)"
fi
