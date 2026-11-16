#!/bin/bash
# Open a new interactive bash shell with the analysis environment active:
#
#   setup/shell.sh [PREFIX]      (PREFIX: default $DIJET_ENV, else <repo>/env)
#
# The new shell reads your ~/.bashrc as usual (aliases, prompt, ...) and then sources
# setup/activate.sh, which removes any other ROOT installation set up there and activates the
# environment. Leave it with `exit`.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export DIJET_SHELL_REPO="$REPO"
export DIJET_SHELL_PREFIX="${1:-}"
exec bash --rcfile "$REPO/setup/_shellrc.bash" -i
