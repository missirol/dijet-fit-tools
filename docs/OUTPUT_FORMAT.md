# Output format: the workspace file

The final output of the workflow is **one ROOT file**
(default `output/workspace/dijet_workspace.root`, written by `scripts/build_workspace.py`)
containing **one `RooWorkspace` named `w`** with the signal and background models of all
categories. The naming follows the conventions of the CMS Combine tool, so that the
workspace can be used directly in a shape datacard (see the example below). Next to it,
`dijet_workspace.json` summarises the content (object names, binning, background function
index map, consistency checks).

## Categories

The categories are those selected in step 1 (`prepare_inputs.py --categories/--exclude`); each
category name `CAT` is the first capture group of the matching pattern. With the defaults of
`run_all.sh`, a category is a slice in the rapidity separation of the two jets, named after
the input histogram without the `_J1J2_mass` suffix, e.g.
`L1SCaloJet_Std_Mjj150dPhi1p0dY0p4to0p5` (0.4 < |Δy| < 0.5): 11 categories (0.0–0.1, …,
1.0–1.1; the inclusive 0.0–1.1 histogram is excluded). Categories are independent: no
parameter is shared between them, except the signal mass `MH`.

## Objects

| name | class | meaning |
|---|---|---|
| `MH` | `RooRealVar` (constant) | resonance mass hypothesis [GeV], range = mass points used (200–500); constant, so that it is a fixed hypothesis set by Combine's `-m` and not a fitted parameter |
| `mjj_CAT` | `RooRealVar` | observable m_jj [GeV]; range = fit range, binning = binning of the category |
| `data_obs_CAT` | `RooDataHist` | observed data in the fit range |
| `ZprimeQQ_CAT` | `RooCrystalBall` | signal shape: double-sided Crystal Ball with common sigma |
| `ZprimeQQ_CAT_mean`, `_sigma`, `_alphaL`, `_nL`, `_alphaR`, `_nR` | `RooFormulaVar(MH)` | signal shape parameters as functions of `MH` |
| `ZprimeQQ_CAT_norm` | `RooFormulaVar(MH)` | expected signal yield **inside the fit range** for σ(pp→Z'→qq) = **1 nb** and the luminosity given to `prepare_inputs.py` (default 680000 nb⁻¹) |
| `bkg_CAT` | `RooMultiPdf` (Combine) | background model: discrete choice among the accepted functions |
| `pdfindex_CAT` | `RooCategory` | index of the background function used by `bkg_CAT`; default (index 0) = the F-test-selected family maximum with the largest goodness-of-fit p-value |
| `bkg_CAT_<func>` | `RooGenericPdf` | each accepted background function (e.g. `bkg_CAT_dijet5`) |
| `<func>_p<k>_CAT` | `RooRealVar` | parameters of background function `<func>` (floating, set to the fitted values) |
| `bkg_CAT_norm` | `RooRealVar` | background yield in the fit range (floating; set to the fitted value of the default function) |

Notes:

* **Yields are always yields inside the fit range** (the range of `mjj_CAT`), for signal and
  background alike.
* **Signal normalisation**: the templates were scaled by L / N_processed (N_processed = integral
  of `eventsProcessed` in the signal file), i.e. to σ = 1 nb. A Combine signal strength `r`
  is therefore the cross section × branching fraction in nb.
* **Signal parametrisation**: by default a natural cubic spline through the values fitted at the
  mass points (polynomials can be chosen in `parametrize_signal.py`). Failed signal fits are not
  used. Each parameter is interpolated only through the mass points where it was floating in
  the fit (values fixed by the fitting procedure are not measurements), or through all points if
  fewer than two float; beyond its outermost points used, a parameter is held constant (e.g.
  `alphaR` when it was fixed at the highest mass because the right tail lies outside the fit
  range). The functions are clamped
  so that interpolation cannot produce an invalid shape: `alphaL`, `alphaR`, `nL`, `nR` to the
  ranges allowed in the signal fits (α in [0.1, 8], n in [0.5, 100]), `sigma` to at least 1 GeV,
  the yield to non-negative values; the mean is not clamped. `MH` is limited to the mass range
  covered by every category, so no parameter is extrapolated.
* **Background functions** (see `scripts/fit_background.py` and `dijetfit/background.py`):
  the families `dijet`, `moddijet` and `ua2` with 3–6 parameters (the number includes the
  normalisation). The formulas are written in a numerically decorrelated basis: the fitted
  parameters `<func>_p<k>_CAT` are linear combinations of the parameters of the published form,
  and have uncertainties of order 1. The function family is unchanged; the coefficients of the
  published form are in `background_fits.json` (`params_original`).
* **Background functions in `bkg_CAT`** (at most 5): a fit is acceptable if it has status 0,
  covariance quality 3, no parameter at a limit and goodness-of-fit p-value > 0.05. In each
  family, standard Fisher F-tests (threshold 0.05), starting from the lowest acceptable order,
  give the highest order. The default is the family maximum with the largest p-value; up to 4
  alternates are the other acceptable functions of order at most their family's maximum,
  ranked by p-value. `pdfindex_CAT` index 0 is the default, then the alternates in this order.
* The background pdfs carry their own numerical integrator settings (adaptive Gauss–Kronrod,
  relative precision 1e-10), which the fits need for these sample sizes (~2×10⁷ events/category).
* **Binning**: only 1 GeV bins are supported in the workspace (the binning is stored in
  `mjj_CAT`). Combine evaluates the pdfs at the bin centres; for wider bins (a category rebinned
  by the fit-range scan) the fits integrate the pdfs over each bin instead, and bin-centre
  evaluation would be biased by many standard deviations for these sample sizes. So
  `build_workspace.py` stops with an error for any category with bins wider than 1 GeV; the
  rebinned fit results remain available for inspection.
* `RooMultiPdf` is a Combine class: reading `bkg_CAT` needs the Combine library
  (`gSystem->Load("libHiggsAnalysisCombinedLimit")`). With `build_workspace.py --no-multipdf`
  the file contains everything except the `RooMultiPdf` objects.

## Example datacard (one category)

```
imax 1
jmax 1
kmax *
---------------
shapes ZprimeQQ  CAT  dijet_workspace.root  w:ZprimeQQ_CAT
shapes bkg       CAT  dijet_workspace.root  w:bkg_CAT
shapes data_obs  CAT  dijet_workspace.root  w:data_obs_CAT
---------------
bin          CAT
observation  -1
---------------
bin          CAT       CAT
process      ZprimeQQ  bkg
process      0         1
rate         1         1
---------------
pdfindex_CAT discrete
```

`rate 1` means that the normalisation is taken from `ZprimeQQ_CAT_norm` / `bkg_CAT_norm`.
The mass hypothesis is set with `-m`, e.g. `text2workspace.py card.txt -m 300`.
