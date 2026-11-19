# Analysis cases: configurations and results log

One entry per analysis case that was run: inputs, categories, settings, and the main results, so
that cases can be compared and rerun. ROOT files are never committed: the inputs are described by
their origin, and the user provides them again when a case is rerun.

Add an entry for every new case (or a rerun with different settings), with the commit of the
code it was run with. Newest entries at the bottom of each section.

## How to run a case

```
setup/envrun.sh studies/run_case.sh OUTDIR DATA.root SIGNAL_M200.root SIGNAL_M300.root ...
```

with the settings of the case as environment variables (see the header of
[`../run_case.sh`](../run_case.sh)). It runs the workflow (`run_all.sh`, steps 1–7 including the
expected limits) and a short bias test (100 toys per mass at μ = 1, generated with the default
function of every category, fitted with the envelope and with the generating function fixed).
It can be rerun after an interruption and continues where it stopped.

Example (Opt1 below):

```
CATEGORIES='(L1SCaloJet_Std_Opt1dRap\d+p\d+to\d+p\d+dPhi\d+p\d+to\d+p\d+)_J1J2_mass' EXCLUDE="" \
XMIN=170 XMAX=600 setup/envrun.sh studies/run_case.sh out/opt1 \
  L1ScoutingSelection_Run2026D.root ZprimeToQQ_M{200,300,400,500}_TuneCP5_13p6TeV_pythia8.root
```

## Input sets

All sets: `L1ScoutingSelection_Run2026D.root` (data) and `ZprimeToQQ_M{200,300,400,500}_TuneCP5_13p6TeV_pythia8.root`
(signal, 10⁶ processed events each), luminosity 680000 nb⁻¹. Each set was uploaded separately;
the histograms are 1 GeV wide.

| Set | Provided | m_jj histograms used | Filled from |
|---|---|---|---|
| 1 | first set (production case of the first PR) | `L1SCaloJet_Std_Mjj150dPhi1p0dY<A>to<B>`, `L1TJet_…` | 150 GeV |
| 2 | second set (four-case comparison) | as set 1 | 150 GeV |
| 3 | 2026-11-14 | `L1SCaloJet_Std_Mjj170dRap<A>to<B>dPhi<C>to<D>` (also `Cmb`) | 170 GeV |
| 4 | 2026-11-15 | `L1SCaloJet_Std_Opt0…`, `L1SCaloJet_Std_Opt1…` | 170 GeV |

Set 4 `Opt0` is identical to set 3 bin by bin in 150–700 GeV (its histograms extend to 720 GeV).

## Cases

Expected limits: median expected 95% CL upper limit on r (σ = 1 nb), blind, envelope (computed
before step 7 was added with the same Combine command as `scripts/expected_limits.py`). Bias at
μ = 1: envelope fit, the fixed-function fit in brackets (100 toys: ±0.10; 300 toys: ±0.058).

| Case | Set | Categories | Fit range (configured → selected) | Limits 200 / 300 / 400 / 500 GeV | Bias at μ = 1, 200 / 300 / 400 / 500 GeV | Code |
|---|---|---|---|---|---|---|
| Case 1 | 1 | 11 \|Δy\| slices 0–1.1 (`L1SCaloJet_Std_Mjj150dPhi1p0dY…`, inclusive excluded) | 150–600 → 150–600 | 1.04 / 0.373 / 0.153 / 0.102 | up to 0.76 in size (−, +, −, +), before the new background families | before d8d44f9 |
| Case 2 | 1 | as Case 1 with `L1TJet` | 150–600 → fails | – | – | d8d44f9 – 2546270 |
| Case 3 | 2 | as Case 1 | 150–600 → 170–600 | 1.36 / 0.261 / 0.119 / 0.074 | −0.42 / +0.25 / −0.21 / +0.39 (300 toys); −0.33 / +0.25 / −0.24 / +0.04 (200 toys, other seeds) | d8d44f9 – 2546270 |
| Case 4 | 2 | as Case 2 | 150–600 → fails | – | – | d8d44f9 – 2546270 |
| Case 5 | 3 | 20: \|Δy\| 0–2 in 4 slices × \|Δφ\| 2.10–2.80–2.90–3.00–3.05–3.20 (`Mjj170dRap…dPhi…`, `Std`) | 170–600 → 200–600 (M200 dropped) | – / 0.277 / 0.120 / 0.103 | see Opt0 | d8d44f9 – 2546270 |
| Opt0 | 4 | as Case 5 (`Opt0…`) | 170–600 → 200–600 (M200 dropped) | – / 0.277 / 0.120 / 0.103 | – / +0.22 (+0.04) / −0.25 (−0.03) / −0.08 (−0.06) | 2546270 |
| **Opt1** | 4 | 8: \|Δy\| 0–2 in 4 slices × \|Δφ\| 2.10–3.00, 3.00–3.20 (`Opt1…`) | 170–600 → **175–600** | **1.33 / 0.251 / 0.109 / 0.080** | −0.17 (+0.09) / +0.40 (+0.17) / −0.17 (+0.13) / +0.14 (−0.08) | 2546270 |
| Opt1-180 | 4 | as Opt1 | 180–600 → 185–600 | 1.48 / 0.270 / 0.106 / 0.088 | −0.10 (+0.19) / +0.49 (+0.11) / −0.46 (−0.05) / +0.29 (+0.07) | 2546270 |
| Opt1-580 | 4 | as Opt1 | 170–580 → 190–580 | 1.79 / 0.280 / 0.118 / 0.107 | +0.21 (+0.21) / −0.03 (−0.12) / −0.15 (−0.02) / −0.14 (−0.18) | 2546270 |

Cases 2 and 4 fail: the `L1TJet` spectra have structure from the quantised jet η and pT, and no
background function passes in any category up to 200 GeV.

## Notes per case

- **Case 1**: the production case of the first PR. Its bias test (with dijet, moddijet and ua2
  only) showed the envelope bias pattern; it was not rerun with the new families.
- **Case 3**: one category (dY0p5to0p6) has no passing function from 150 to 165 GeV. Pulls of
  the envelope fits (200 toys, `--robustFit`): means −0.36 / +0.28 / −0.25 / +0.08, widths 0.99
  / 1.04 / 1.10 / 1.10. The two bias runs disagree at 500 GeV (+0.39 ± 0.06 vs +0.04 ± 0.07),
  not understood. A run with toys from every alternate function was stopped at 72 of 192 jobs.
- **Case 5 / Opt0**: one category (dRap0p50to1p00dPhi3p00to3p05) has no passing function below
  200 GeV, so M200 is dropped. The categorisation gains in S/√B (+10% vs Case 3), but the many
  floating background parameters (order 5–6 defaults in 20 categories) and the envelope lose it:
  with the background shapes frozen the limits would be 0.111 (300 GeV) and 0.029 (500 GeV);
  the envelope costs +32% at 500 GeV. Spread of r, envelope over fixed: 1.21 / 1.38 / 1.51 at
  300 / 400 / 500 GeV.
- **Opt1**: the best case so far. Weakest category: dRap0p50to1p00dPhi3p00to3p20 (fails at
  170 GeV; only polypow5 passes at 175 GeV; several functions at p ≈ 0.03–0.047). Spread of r,
  envelope over fixed: 1.04 / 1.10 / 1.23 at 300 / 400 / 500 GeV. 1–6 envelope fits per 100 toys
  fail.
- **Opt1-180, Opt1-580**: raising the lower edge, or lowering the upper edge (which forces the
  lower edge up to 190 GeV), makes the limits worse at most masses; the 580–600 GeV region (CaloJet
  m_jj distortion from the cut on the uncorrected jet pT) does not cause the fit failures at the
  current p-value threshold.

## Open questions

- A lower edge below 170 GeV for Opt1 needs histograms filled from lower m_jj.
- The choice of the default categorisation (currently still Case 1's in `run_all.sh`).
