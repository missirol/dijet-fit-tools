#!/bin/bash
# Create a conda environment with ROOT 6.36.14 and CMS Combine 11.1.0 (from conda-forge),
# using micromamba. Usage:  setup/install_env.sh [PREFIX]   (default: ./env)
#
# Afterwards, run the scripts inside the environment, e.g.
#   MAMBA_ROOT_PREFIX=$PWD/.mamba ./bin/micromamba run -p ./env ./run_all.sh ...
# (the environment must be activated / run through micromamba: ROOT's interpreter needs the
#  environment's compiler settings).
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
  "root=6.36.14" "cms-combine=11.1.0" "python=3.12" numpy
./bin/micromamba run -p "$PREFIX" python -c \
  "import ROOT; print('ROOT', ROOT.gROOT.GetVersion()); print('Combine:', ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit') >= 0)"
