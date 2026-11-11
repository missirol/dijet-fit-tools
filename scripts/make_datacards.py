#!/usr/bin/env python3
"""Step 6: write the Combine datacards (one per category) and the combined workspace.

Reads the summary of the workspace written by build_workspace.py (dijet_workspace.json) and
writes, for every category CAT, a datacard taking everything from that workspace:
  data_obs  w:data_obs_CAT
  ZprimeQQ  w:ZprimeQQ_CAT   yield ZprimeQQ_CAT_norm (sigma = 1 nb, function of MH) x r
  bkg       w:bkg_CAT        RooMultiPdf, yield bkg_CAT_norm (floating)
  pdfindex_CAT discrete      the background function index is profiled (discrete profiling)
The signal process has the same name in all categories, so the signal strength r is common to
all of them; the background parameters and yields are separate objects per category, i.e.
uncorrelated. No systematic uncertainties are included.

The cards are then combined with combineCards.py (the category names are kept as bin names)
and converted with text2workspace.py into the workspace used by the statistical tests. The
combined workspace is checked: every category must be present with its index, and the expected
yields must equal r x ZprimeQQ_CAT_norm(MH) (signal) and bkg_CAT_norm (background).

Outputs (in --outdir):
  datacard_<CAT>.txt        one datacard per category
  datacard_combined.txt     combination of all categories
  datacard_combined.root    text2workspace.py output (POI r; MH set with combine -m)
  datacards.json            names used by the later steps (categories, indices, masses)
"""
import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dijetfit import combine as C  # noqa: E402
from dijetfit.common import SIGNAL_PROCESS, read_json, setup_root, write_json  # noqa: E402


def check_combined(ROOT, path, summary, cats, source):
    """Expected yields of the combined workspace against the source workspace."""
    f = ROOT.TFile.Open(str(path))
    w = f.Get("w")
    src_f = ROOT.TFile.Open(str(source))
    src = src_f.Get(summary.get("workspace", "w"))
    r, mh, mh_src = w.var(C.POI), w.var(summary["mass_variable"]), src.var(summary["mass_variable"])
    problems = []
    if not r or not mh:
        return [f"no variable '{C.POI}' or '{summary['mass_variable']}' in {path}"]
    masses = sorted({m for c in cats for m in summary["categories"][c]["signal"]["masses"]})
    for c in cats:
        if not w.cat(f"pdfindex_{c}"):
            problems.append(f"{c}: no pdfindex_{c}")
        sig = w.function(f"n_exp_final_bin{c}_proc_{SIGNAL_PROCESS}") or \
            w.function(f"n_exp_bin{c}_proc_{SIGNAL_PROCESS}")
        bkg = w.function(f"n_exp_final_bin{c}_proc_{C.BKG_PROCESS}") or \
            w.function(f"n_exp_bin{c}_proc_{C.BKG_PROCESS}")
        if not sig or not bkg:
            problems.append(f"{c}: expected-yield functions not found")
            continue
        for m in masses:
            for rv in (1.0, 2.5):
                mh.setVal(m)
                mh_src.setVal(m)
                r.setVal(rv)
                want = rv * src.function(f"{SIGNAL_PROCESS}_{c}_norm").getVal()
                if abs(sig.getVal() - want) > 1e-9 * max(1.0, abs(want)):
                    problems.append(f"{c}: signal yield {sig.getVal():g} != {want:g} "
                                    f"(MH {m}, r {rv})")
        want = src.var(f"bkg_{c}_norm").getVal()
        if abs(bkg.getVal() - want) > 1e-9 * want:
            problems.append(f"{c}: background yield {bkg.getVal():g} != {want:g}")
    r.setVal(1.0)
    return problems


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("workspace_summary", help="dijet_workspace.json from build_workspace.py")
    ap.add_argument("-o", "--outdir", default="output/datacards")
    ap.add_argument("--categories", nargs="*", help="subset of categories (default: all)")
    ap.add_argument("--no-text2workspace", action="store_true",
                    help="only write the datacards (no combined workspace, no checks)")
    args = ap.parse_args()

    summary = read_json(args.workspace_summary)
    if not summary.get("multipdf"):
        sys.exit("ERROR: the workspace has no RooMultiPdf (written with --no-multipdf)")
    ws_file = Path(summary["file"])
    if not ws_file.is_absolute():
        ws_file = (Path(args.workspace_summary).parent / ws_file).resolve()
    if not ws_file.exists():
        sys.exit(f"ERROR: workspace file {ws_file} not found")
    cats = args.categories or list(summary["categories"])
    unknown = [c for c in cats if c not in summary["categories"]]
    if unknown:
        sys.exit(f"ERROR: categories not in the workspace: {', '.join(unknown)}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rel_ws = os.path.relpath(ws_file, outdir.resolve())  # cards stay valid if moved together
    cards = {}
    for c in cats:
        card = outdir / f"datacard_{c}.txt"
        card.write_text(C.datacard_text(c, rel_ws, summary.get("workspace", "w")))
        cards[c] = card.name
    print(f"Written {len(cards)} datacards in {outdir}")

    meta = dict(workspace=str(ws_file), mass_variable=summary["mass_variable"],
                mass_range=summary["mass_range"], fit_range=summary["fit_range"], poi=C.POI,
                cards=cards, combined_card="datacard_combined.txt",
                combined_workspace=None, categories={})
    for c in cats:
        b = summary["categories"][c]["background"]
        meta["categories"][c] = dict(index=b["index"], default_index=b["default_index"],
                                     index_to_function=b["index_to_function"],
                                     masses=summary["categories"][c]["signal"]["masses"],
                                     n_data=summary["categories"][c]["n_data"])
    log = outdir / "combine_commands.log"
    log.write_text("")
    if not args.no_text2workspace:
        C.require_tool("combineCards.py")
        C.require_tool("text2workspace.py")
        import subprocess
        p = subprocess.run(["combineCards.py", *(f"{c}={cards[c]}" for c in cats)], cwd=outdir,
                           capture_output=True, text=True)
        if p.returncode:
            sys.exit(f"ERROR: combineCards.py failed:\n{p.stderr[-3000:]}")
        (outdir / "datacard_combined.txt").write_text(p.stdout)
        C.run(["text2workspace.py", "datacard_combined.txt", "-o", "datacard_combined.root"],
              outdir, log)
        meta["combined_workspace"] = "datacard_combined.root"
        ROOT = setup_root(load_combine=True)
        problems = check_combined(ROOT, outdir / "datacard_combined.root", summary, cats, ws_file)
        if problems:
            sys.exit("ERROR: the combined workspace does not reproduce the inputs:\n  "
                     + "\n  ".join(problems))
        print(f"Combined {len(cats)} categories: {outdir / 'datacard_combined.txt'}, "
              f"{outdir / 'datacard_combined.root'}")
        print("Checks on the combined workspace: all categories and indices present; expected "
              "yields = r x signal norm(MH) and the background norm, at every signal mass")
    write_json(outdir / "datacards.json", meta)
    print(f"Written {outdir / 'datacards.json'}")


if __name__ == "__main__":
    main()
