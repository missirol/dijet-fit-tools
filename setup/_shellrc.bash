# Startup file of the shell opened by setup/shell.sh (not meant to be used directly).
if [ -f "$HOME/.bashrc" ]; then
  . "$HOME/.bashrc"
fi
if [ -n "${DIJET_SHELL_PREFIX:-}" ]; then
  . "$DIJET_SHELL_REPO/setup/activate.sh" "$DIJET_SHELL_PREFIX"
else
  . "$DIJET_SHELL_REPO/setup/activate.sh"
fi
unset DIJET_SHELL_REPO DIJET_SHELL_PREFIX
