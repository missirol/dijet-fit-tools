#!/bin/bash
# Create a conda environment with ROOT 6.36.14 and CMS Combine 11.1.0 (from conda-forge),
# using micromamba. Usage:  setup/install_env.sh [PREFIX]   (default: ./env)
#
# Afterwards, run the scripts inside the environment (ROOT's interpreter needs the
# environment's compiler settings), in one of three ways:
#   setup/shell.sh                  new shell with the environment active
#   source setup/activate.sh        activate it in the current shell (bash or zsh)
#   setup/envrun.sh COMMAND ...     run a single command in it
# All three first remove another ROOT installation set up in the calling shell (ROOTSYS,
# PYTHONPATH, LD_LIBRARY_PATH, ...; see setup/_clean_env.sh); with a non-default PREFIX,
# pass it to shell.sh / activate.sh, or set DIJET_ENV=PREFIX.
#
# Combine is only needed for RooMultiPdf in the final workspace (build_workspace.py); without
# it, use build_workspace.py --no-multipdf.
set -euo pipefail
PREFIX="${1:-$PWD/env}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$PWD/.mamba}"
mkdir -p bin
if [ ! -x bin/micromamba ]; then
  curl -sSfL -o bin/micromamba \
    https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
  chmod +x bin/micromamba
fi
./bin/micromamba create -y -p "$PREFIX" -c conda-forge \
  "root=6.36.14" "cms-combine=11.1.0" "python=3.12" numpy pytest
DIJET_ENV="$PREFIX" setup/envrun.sh python -c \
  "import ROOT; print('ROOT', ROOT.gROOT.GetVersion()); print('Combine:', ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit') >= 0)"
echo "Done. Use the environment with:  setup/shell.sh${1:+ $PREFIX}   (or: source setup/activate.sh${1:+ $PREFIX})"
