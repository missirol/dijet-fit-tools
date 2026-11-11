"""Helpers for the CMS Combine steps: datacards, fit options, bias-test statistics."""
import math
import os
import shutil
import subprocess

from .common import DATA_NAME, SIGNAL_PROCESS

BKG_PROCESS = "bkg"
POI = "r"

# Options for fits with RooMultiPdf (discrete profiling), as recommended in the Combine
# documentation: start with Minuit strategy 0, freeze the parameters of the background functions
# not selected by the index, and in the discrete minimisation skip constant parameters,
# constraint terms and channels whose index is not being changed.
DISCRETE_PROFILING_OPTIONS = [
    "--cminDefaultMinimizerStrategy", "0",
    "--X-rtd", "MINIMIZER_freezeDisassociatedParams",
    "--X-rtd", "MINIMIZER_multiMin_hideConstants",
    "--X-rtd", "MINIMIZER_multiMin_maskConstraints",
    "--X-rtd", "MINIMIZER_multiMin_maskChannels=2",
]


def datacard_text(cat, ws_file, workspace="w", signal=SIGNAL_PROCESS, bkg=BKG_PROCESS,
                  data=DATA_NAME):
    """Datacard of one category. Signal and background are the pdfs <process>_<cat> of the
    workspace; their yields are the workspace functions <pdf>_norm (Combine multiplies them by
    the rate, here 1): the signal yield for sigma = 1 nb as a function of MH (x r), the floating
    background yield. The background parameters and the yield float freely; the index of the
    RooMultiPdf is declared as a discrete nuisance parameter."""
    w = max(len(cat), len(signal), len(bkg), len(data))

    def col(label, *xs):
        return (f"{label:<13}" + "".join(f"{x:<{w + 2}}" for x in xs)).rstrip()
    return "\n".join([
        f"# category {cat}",
        "imax 1  number of channels",
        "jmax 1  number of processes - 1",
        "kmax *  number of nuisance parameters",
        "-" * 40,
        f"shapes {signal:<{w}}  {cat}  {ws_file}  {workspace}:{signal}_{cat}",
        f"shapes {bkg:<{w}}  {cat}  {ws_file}  {workspace}:bkg_{cat}",
        f"shapes {data:<{w}}  {cat}  {ws_file}  {workspace}:{data}_{cat}",
        "-" * 40,
        col("bin", cat),
        col("observation", "-1"),
        "-" * 40,
        col("bin", cat, cat),
        col("process", signal, bkg),
        col("process", "0", "1"),
        col("rate", "1", "1"),
        "-" * 40,
        f"pdfindex_{cat}  discrete",
        "",
    ])


def pull(r_fit, r_true, err_lo, err_hi):
    """Pull of the fitted signal strength with asymmetric errors (Combine convention): the error
    on the side of the true value, i.e. the upper error if the fit is below the truth."""
    err = err_hi if r_fit < r_true else err_lo
    if not err or err <= 0 or not math.isfinite(err):
        return None
    return (r_fit - r_true) / err


def mean_and_error(values):
    """Mean, its standard error, standard deviation and its standard error (Gaussian approx.)."""
    n = len(values)
    if n < 2:
        return dict(n=n, mean=values[0] if n else float("nan"), mean_err=float("nan"),
                    std=float("nan"), std_err=float("nan"))
    mean = sum(values) / n
    std = math.sqrt(sum((v - mean) ** 2 for v in values) / (n - 1))
    return dict(n=n, mean=mean, mean_err=std / math.sqrt(n), std=std,
                std_err=std / math.sqrt(2 * (n - 1)))


def median(values):
    s = sorted(values)
    n = len(s)
    if not n:
        return float("nan")
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def require_tool(name):
    if shutil.which(name) is None:
        raise RuntimeError(f"'{name}' not found: CMS Combine must be installed and in the PATH")


def run(cmd, cwd, log):
    """Run a command in `cwd`, appending its output to `log`; raise on failure."""
    with open(log, "a") as f:
        f.write("$ " + " ".join(map(str, cmd)) + "\n")
        f.flush()
        p = subprocess.run(list(map(str, cmd)), cwd=cwd, stdout=f, stderr=subprocess.STDOUT,
                           env=dict(os.environ))
    if p.returncode:
        raise RuntimeError(f"command failed (exit code {p.returncode}), see {log}: "
                           + " ".join(map(str, cmd)))
