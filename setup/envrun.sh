#!/bin/bash
# Run one command in the analysis environment, without activating it in the calling shell:
#
#   setup/envrun.sh COMMAND [ARGS...]      e.g.  setup/envrun.sh pytest
#
# The environment is <repo>/env, or $DIJET_ENV if set. Variables through which another ROOT
# installation in the calling shell would leak in (ROOTSYS, PYTHONPATH, LD_LIBRARY_PATH, ...;
# see setup/_clean_env.sh) are removed first. For interactive work use setup/shell.sh or
# `source setup/activate.sh` instead.
set -euo pipefail
if [ $# -eq 0 ]; then
  echo "usage: setup/envrun.sh COMMAND [ARGS...]" >&2
  exit 2
fi
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
. "$REPO/setup/_clean_env.sh"
_dijet_locate "$REPO"
_dijet_clean_env
exec "$DIJET_MICROMAMBA" run -p "$DIJET_PREFIX" "$@"
