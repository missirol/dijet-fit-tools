"""Naming conventions, ROOT setup and small helpers shared by all workflow steps."""
import json
import os
import re
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path

SQRTS = 13600.0  # GeV

SIGNAL_PROCESS = "ZprimeQQ"
DATA_NAME = "data_obs"


def signal_hist_name(mass):
    return f"{SIGNAL_PROCESS}_{int(mass)}"


def select_categories(names, patterns, excludes=()):
    """Map analysis categories to input histogram names.

    A histogram is selected if its whole name matches one of the regular expressions
    `patterns` and none of `excludes`. The category name is the first capture group of the
    matching pattern, or the full histogram name if the pattern has no group.
    Returns ({category: histogram name}, {pattern: number of histograms it selected}).
    Raises ValueError if two histograms give the same category name."""
    pats = [re.compile(p) for p in patterns]
    excl = [re.compile(p) for p in excludes]
    cats, used = {}, {p: 0 for p in patterns}
    for name in names:
        if any(e.fullmatch(name) for e in excl):
            continue
        for text, pat in zip(patterns, pats):
            m = pat.fullmatch(name)
            if not m:
                continue
            cat = m.group(1) if pat.groups else name
            if cat in cats and cats[cat] != name:
                raise ValueError(f"histograms '{cats[cat]}' and '{name}' both give category '{cat}'")
            cats[cat] = name
            used[text] += 1
            break
    return cats, used


def category_labels(cats):
    """Short labels for plots: the category names without the prefix they all share (cut back
    to the last '_' in it, so that labels start at a word boundary)."""
    cats = list(cats)
    if len(cats) < 2:
        return {c: c for c in cats}
    prefix = os.path.commonprefix(cats)
    cut = prefix.rfind("_") + 1
    return {c: c[cut:] or c for c in cats}


def mass_from_filename(path):
    m = re.search(r"_M(\d+)_", Path(path).name)
    return int(m.group(1)) if m else None


def setup_root(load_combine=False):
    """Import ROOT in batch mode with a quiet default style. Returns the ROOT module."""
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.TH1.AddDirectory(False)
    ROOT.gErrorIgnoreLevel = ROOT.kWarning  # silence "Info in <TCanvas::Print>"
    # hide RooFit INFO/PROGRESS messages; warnings and errors are kept (and counted)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    if load_combine:
        if ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit") < 0:
            sys.exit("ERROR: cannot load libHiggsAnalysisCombinedLimit (CMS Combine), "
                     "which provides RooMultiPdf; see setup/install_env.sh, "
                     "or use --no-multipdf")
    return ROOT


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, payload):
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n")


@contextmanager
def capture_output():
    """Capture everything written to the process-level stdout/stderr (including C++ output
    from RooFit and Minuit) inside the block. Yields a dict whose 'text' is filled on exit."""
    out = {"text": ""}
    sys.stdout.flush()
    sys.stderr.flush()
    saved = os.dup(1), os.dup(2)
    with tempfile.TemporaryFile(mode="w+b") as tmp:
        os.dup2(tmp.fileno(), 1)
        os.dup2(tmp.fileno(), 2)
        try:
            yield out
        finally:
            try:
                import ROOT  # flush C stdio buffers before restoring the descriptors
                ROOT.gSystem.ProcessEvents()
                ROOT.gInterpreter.ProcessLine("fflush(stdout); fflush(stderr);")
            except Exception:
                pass
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(saved[0], 1)
            os.dup2(saved[1], 2)
            os.close(saved[0])
            os.close(saved[1])
            tmp.seek(0)
            out["text"] = tmp.read().decode(errors="replace")


def count_messages(text):
    """Number of warning- and error-like lines in captured RooFit/Minuit output."""
    nwarn = nerr = 0
    for line in text.splitlines():
        if _ERR_RE.search(line):
            nerr += 1
        elif _WARN_RE.search(line):
            nwarn += 1
    return nwarn, nerr


# RooFit tags messages with upper-case levels ("[#0] ERROR:", "[#1] WARNING:"); ROOT/Minuit
# use "Error in <...>" / "Warning in <...>"; Minuit2 prints "Minuit2:... Warning".
_ERR_RE = re.compile(r"\bERROR\b|Error in <|\bMinuit2?\w*.*\bError\b")
_WARN_RE = re.compile(r"\bWARNING\b|Warning in <|\bMinuit2?\w*.*\bWarning\b|\bWARNING:")


def message_summary(text, max_lines=20):
    """The distinct warning/error lines of a captured log (for reports)."""
    seen = []
    for line in text.splitlines():
        if (_ERR_RE.search(line) or _WARN_RE.search(line)) and line.strip() not in seen:
            seen.append(line.strip())
        if len(seen) >= max_lines:
            break
    return seen
