# Startup file of the shell opened by setup/shell.sh (not meant to be used directly).
if [ -f "$HOME/.bashrc" ]; then
  . "$HOME/.bashrc"
fi
if [ -n "${DIJET_SHELL_PREFIX:-}" ]; then
  export DIJET_ENV="$DIJET_SHELL_PREFIX"
fi
if ! . "$DIJET_SHELL_REPO/setup/activate.sh"; then
  echo "dijet setup: WARNING: the analysis environment is NOT active in this shell (see above)" >&2
  PS1="[dijet env NOT active] ${PS1:-\$ }"
fi
unset DIJET_SHELL_REPO DIJET_SHELL_PREFIX
