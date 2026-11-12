# Shared helpers for setup/activate.sh, setup/envrun.sh and setup/shell.sh (bash and zsh).
# Not meant to be used directly.
#
# A ROOT installation set up in the calling shell (thisroot.sh, CMSSW, ...) leaks into the
# conda environment through the variables below: e.g. PYTHONPATH=$ROOTSYS/lib makes the
# environment's Python load the other ROOT's bindings, and `import ROOT` crashes. They are
# removed before the environment is used; nothing in the environment needs them.

DIJET_CLEAN_VARS="ROOTSYS PYTHONPATH PYTHONHOME LD_LIBRARY_PATH DYLD_LIBRARY_PATH LD_PRELOAD \
ROOT_INCLUDE_PATH CPLUS_INCLUDE_PATH C_INCLUDE_PATH CLING_STANDARD_PCH"

# Remove the variables above, and from PATH every ROOT bin directory ($ROOTSYS/bin, or any
# directory that contains thisroot.sh). Prints one line listing what was removed (to stderr).
_dijet_clean_env() {
  local removed="" v d newpath="" rootbin=""
  if [ -n "${ROOTSYS:-}" ]; then rootbin="${ROOTSYS%/}/bin"; fi
  for v in $(echo "$DIJET_CLEAN_VARS"); do
    if [ -n "$(eval "echo \"\${$v:-}\"")" ]; then
      removed="$removed $v"
    fi
    unset "$v"
  done
  local IFS=":"
  for d in $(echo "${PATH:-}"); do
    if [ -n "$d" ] && { [ "${d%/}" = "$rootbin" ] || [ -f "$d/thisroot.sh" ]; }; then
      removed="$removed PATH:$d"
      continue
    fi
    newpath="${newpath:+$newpath:}$d"
  done
  PATH="$newpath"
  export PATH
  if [ -n "$removed" ]; then
    echo "dijet setup: removed from the environment (other ROOT installation):$removed" >&2
  fi
}

# Locate micromamba and the environment; sets DIJET_MICROMAMBA, DIJET_PREFIX and
# MAMBA_ROOT_PREFIX. Usage: _dijet_locate REPO [PREFIX]
#   PREFIX: argument, else $DIJET_ENV, else REPO/env
#   micromamba: REPO/bin/micromamba (from install_env.sh), else $MAMBA_EXE, else on PATH
_dijet_locate() {
  local repo="$1"
  DIJET_PREFIX="${2:-${DIJET_ENV:-$repo/env}}"
  if [ -x "$repo/bin/micromamba" ]; then
    DIJET_MICROMAMBA="$repo/bin/micromamba"
  elif [ -n "${MAMBA_EXE:-}" ] && [ -x "$MAMBA_EXE" ]; then
    DIJET_MICROMAMBA="$MAMBA_EXE"
  else
    DIJET_MICROMAMBA="$(command -v micromamba 2>/dev/null || true)"
  fi
  if [ -z "$DIJET_MICROMAMBA" ]; then
    echo "dijet setup: micromamba not found (run setup/install_env.sh first)" >&2
    return 1
  fi
  if [ ! -d "$DIJET_PREFIX/conda-meta" ]; then
    echo "dijet setup: no conda environment at $DIJET_PREFIX (run setup/install_env.sh, or" \
         "give the environment path as argument or in DIJET_ENV)" >&2
    return 1
  fi
  DIJET_PREFIX="$(cd "$DIJET_PREFIX" && pwd)"
  export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$repo/.mamba}"
  export DIJET_ENV="$DIJET_PREFIX"
}
