# Activate the analysis environment in the current shell (bash or zsh). Source it, once per shell:
#
#   source setup/activate.sh                       (environment <repo>/env)
#   export DIJET_ENV=PREFIX; source setup/activate.sh   (another environment; exported, so that
#                                                      sourcing again and envrun.sh use it too)
#
# It takes no arguments: a sourced script sees the positional parameters of the caller (e.g. of
# a script that sources it), which must not be mistaken for an environment path.
#
# It first removes the variables through which another ROOT installation set up in this shell
# (ROOTSYS, PYTHONPATH, LD_LIBRARY_PATH, ...; see setup/_clean_env.sh) would leak into the
# environment, then activates the environment with micromamba. Afterwards python, pytest,
# ./run_all.sh, root, ... in this shell (and everything started from it) use the environment's
# ROOT and Combine. Sourcing it again is harmless. See also setup/shell.sh and setup/envrun.sh.

if [ -n "${ZSH_VERSION:-}" ]; then
  case "${ZSH_EVAL_CONTEXT:-}" in *:file*) ;; *)
    echo "setup/activate.sh must be sourced: source setup/activate.sh" >&2; exit 1 ;;
  esac
  _dijet_here="$(eval 'echo ${(%):-%x}')"
  _dijet_shell=zsh
else
  if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    echo "setup/activate.sh must be sourced: source setup/activate.sh" >&2
    exit 1
  fi
  _dijet_here="${BASH_SOURCE[0]}"
  _dijet_shell=bash
fi
_dijet_repo="$(cd "$(dirname "$_dijet_here")/.." && pwd)"
. "$_dijet_repo/setup/_clean_env.sh"

_dijet_activate() {
  _dijet_locate "$_dijet_repo" || return 1
  # leave a previously activated environment first (e.g. when sourced twice)
  if [ -n "${CONDA_PREFIX:-}" ] && command -v micromamba >/dev/null 2>&1; then
    micromamba deactivate >/dev/null 2>&1 || true
  fi
  _dijet_clean_env
  eval "$("$DIJET_MICROMAMBA" shell hook -s "$_dijet_shell")" || return 1
  micromamba activate "$DIJET_PREFIX" || return 1
  local found version
  found="$(command -v root-config 2>/dev/null || true)"
  if [ "$found" != "$DIJET_PREFIX/bin/root-config" ]; then
    echo "dijet setup: WARNING: root-config is '$found', not from $DIJET_PREFIX" >&2
    return 1
  fi
  version="$(root-config --version)"
  echo "dijet setup: environment $DIJET_PREFIX active (ROOT $version)" >&2
}

_dijet_activate
_dijet_status=$?
unset -f _dijet_activate
unset _dijet_here _dijet_shell _dijet_repo
return $_dijet_status 2>/dev/null || true
