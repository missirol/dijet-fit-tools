#!/bin/bash
# Run the full workflow:  ./run_all.sh DATA.root SIGNAL_M200.root SIGNAL_M300.root ...
#
# This file holds the analysis configuration; the scripts themselves have no built-in fit range.
# Environment variables (optional):
#   OUT       output directory                         (default: output)
#   LUMI      integrated luminosity in nb^-1           (default: 680000)
#   XMIN      lower edge of the fit range [GeV]        (default: 150; scanned upwards if needed)
#   XMAX      upper edge of the fit range [GeV]        (default: 600)
#   SCAN_MAX_UNBINNED  largest lower edge tried before rebinning [GeV] (default: 190)
#   SCAN_MAX  largest lower edge tried at all [GeV]    (default: 250)
#   CATEGORIES  space-separated regular expressions selecting the category histograms
#             (whole-name match; first capture group = category name)
#             (default: the |dY| slices L1SCaloJet_Std_Mjj150dPhi1p0dY<A>to<B>_J1J2_mass)
#   EXCLUDE   space-separated regular expressions of histograms to leave out
#             (default: the inclusive slice 0p0to1p1)
#   FORMATS   plot formats                             (default: "png pdf")
set -euo pipefail

if [ $# -lt 3 ]; then
  echo "usage: $0 DATA.root SIGNAL.root [SIGNAL.root ...]" >&2
  exit 1
fi
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${OUT:-output}"
LUMI="${LUMI:-680000}"
XMIN="${XMIN:-150}"
XMAX="${XMAX:-600}"
SCAN_MAX_UNBINNED="${SCAN_MAX_UNBINNED:-190}"
SCAN_MAX="${SCAN_MAX:-250}"
CATEGORIES="${CATEGORIES:-(L1SCaloJet_Std_Mjj150dPhi1p0dY\d+p\d+to\d+p\d+)_J1J2_mass}"
EXCLUDE="${EXCLUDE:-L1SCaloJet_Std_Mjj150dPhi1p0dY0p0to1p1_J1J2_mass}"
FORMATS="${FORMATS:-png pdf}"
read -r -a CATEGORY_ARGS <<< "$CATEGORIES"
read -r -a EXCLUDE_ARGS <<< "$EXCLUDE"
DATA="$1"; shift

mkdir -p "$OUT"
set -x
python3 "$HERE/scripts/prepare_inputs.py" --data "$DATA" --signal "$@" --lumi "$LUMI" \
  --categories "${CATEGORY_ARGS[@]}" --exclude "${EXCLUDE_ARGS[@]}" -o "$OUT/inputs.root" | tee "$OUT/1_prepare_inputs.log"
python3 "$HERE/scripts/fit_background.py" "$OUT/inputs.root" -o "$OUT/background" \
  --xmin "$XMIN" --xmax "$XMAX" --scan --scan-max-unbinned "$SCAN_MAX_UNBINNED" \
  --scan-max "$SCAN_MAX" --formats $FORMATS | tee "$OUT/2_fit_background.log"
python3 "$HERE/scripts/fit_signal.py" "$OUT/inputs.root" \
  --background "$OUT/background/background_fits.json" -o "$OUT/signal" \
  --formats $FORMATS | tee "$OUT/3_fit_signal.log"
python3 "$HERE/scripts/parametrize_signal.py" "$OUT/signal/signal_fits.json" \
  --input "$OUT/inputs.root" -o "$OUT/signal_param" --formats $FORMATS \
  | tee "$OUT/4_parametrize_signal.log"
python3 "$HERE/scripts/build_workspace.py" --input "$OUT/inputs.root" \
  --background "$OUT/background/background_fits.json" \
  --signal "$OUT/signal_param/signal_parametrisation.json" \
  -o "$OUT/workspace/dijet_workspace.root" | tee "$OUT/5_build_workspace.log"
