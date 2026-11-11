# Demo: normalisation of pdfs evaluated at the bin centres

A standalone workflow, independent of the analysis inputs, that shows why the workflow
normalises its pdfs by their **sum over the bin centres** and not by their integral, and that
this is equivalent to exact bin integration (Combine's `RooParametricShapeBinPdf`).

```
python3 demo/bin_normalisation/run_demo.py all -o demo_output            # ~1 minute
python3 demo/bin_normalisation/run_demo.py all -o demo_output --toys 20  # + Poisson toys
```

Run it inside the analysis environment (ROOT, numpy). The steps can also be run one at a time:
`generate`, `fit`, `plot` (same `-o`). Options: `--masses` (200 300 400 500), `--mu` (0 1),
`--ncat` (11 identical, independent categories fitted together, as in the analysis), `--toys`,
`--seed`.

## What it does

**generate** writes `data.root` and `truth.json`. The truth is exact:

* background: the dijet function (1−x)^p1 / x^(p2 + p3 ln x), x = m / 13.6 TeV, with parameters
  that reproduce the background of one real |Δy| category (falling by ×1200 from 150 to 600 GeV),
  2.06×10⁷ events per category in 150–600 GeV, 1 GeV bins;
* signal: double-sided Crystal Ball with the parameters and yield (μ = 1) of the real signal fits;
* expected counts per bin = exact integrals over each bin (Simpson rule with 40 sub-intervals).

Per mass and μ: Asimov histograms (contents = expectation), and with `--toys` Poisson toys.
It also prints the midpoint mismatch of each shape, (Σ f(bin centre) × width) / ∫f − 1.

**fit** fits every data set (all categories together, r common, the background parameters of
each category floating, the signal shape fixed) with four models. They differ only in how the
expected count of a bin is computed from the pdf:

| model | background | signal |
|---|---|---|
| `OLD` | pdf(centre), normalised by its integral | pdf(centre), normalised by its integral (RooCrystalBall's analytical integral) |
| `BKG_ONLY` | pdf(centre), normalised by its **sum over the bin centres** | normalised by its integral |
| `NEW` (the workflow) | sum over the bin centres | sum over the bin centres |
| `EXACT` (= `RooParametricShapeBinPdf`) | exact bin integrals (Simpson formula: [f(c−½) + 4 f(c) + f(c+½)]/6) | exact bin integrals of the fixed shape |

The background function in the fit is the generating family (in the decorrelated basis of the
workflow), so the fit could reproduce the truth exactly: **any bias comes from the model of the
expected counts**. Outputs: `results.md` (the table of (r̂ − μ)/σ_r, and the per-bin predictions
minus the yield used in the likelihood), `results.json`.

**plot** writes `plots/`:

* `bias_vs_mass.png`: (r̂ − μ)/σ_r on Asimov data for every model, mass and μ;
* `prediction_differences_M<M>_mu<mu>.png`: the fitted predictions of each model minus the
  truth, in units of √N per bin; the dashed line is `OLD` evaluated at the true parameters,
  without refitting;
* `fit_M<M>_mu<mu>_<model>.png`: the fit of category 0, with the pull sub-pad.

## How the bias arises

RooFit and Combine fit binned data with an extended likelihood in which a pdf p(m; θ) is
evaluated at the bin centres c_i:

    −ln L = N − Σ_i n_i ln( N · p(c_i; θ) · Δ / Z(θ) )

where N is the expected total (the yield parameter) and Z(θ) the normalisation of the pdf.

* The term N is the expected total of the Poisson likelihood; the bin predictions are
  ν_i = N p(c_i) Δ / Z. These are consistent only if Σ_i ν_i = N, i.e. if **Z = Σ_i p(c_i) Δ**
  (the bin sum). Then −ln L is exactly the binned Poisson likelihood Σ_i (ν_i − n_i ln ν_i).
* With **Z = ∫ p** (the default) the bin predictions add up to N·S/I, with S/I = Σ p(c_i)Δ / ∫p
  = 1 + O(Δ² p'/p at the range edges), here ~10⁻⁶–10⁻⁵ (midpoint rule; the `generate` step
  prints it). The likelihood then contains the extra term N·(1 − S(θ)/I(θ)).
* That term is tiny, but its dependence on the shape parameters is multiplied by N ≈ 2×10⁷:
  its gradient is O(1) per standard deviation of a parameter. The fit therefore moves the shape
  parameters by O(1σ) to lower it, which distorts the predicted spectrum by up to ~0.5σ per bin
  (`prediction_differences`: compare `OLD` fitted with `OLD` at the true parameters, which is
  off by only ~0.02σ). The distortion is smooth and coherent in all categories, and the signal
  shape absorbs part of it: r̂ is biased by O(1σ), with a sign and size that depend on the mass.
* The same holds for the signal pdf: with its analytical integral, the signal predictions add up
  to r·N_sig·(1 + δ), δ ~ 10⁻⁶. This gives a term linear in r, of size r·N_sig·δ, which
  shifts r̂ by about (N_sig σ_r)·δ standard deviations. N_sig σ_r is the number of signal
  events corresponding to one standard deviation, ~4×10⁴ at 200 GeV (σ_r large because of the
  background), hence a bias of a few tenths of σ there (`BKG_ONLY`), negligible at higher masses.
* With the bin sum (`NEW`), the predictions are those of the pdf at the bin centres, which differ
  from the exact bin integrals by O(Δ² p''/p) ~ 10⁻⁴ relative: a smooth shape difference that
  the background parameters absorb. `NEW` and `EXACT` agree to ~10⁻⁴σ per bin and in r̂.

## Expected results (default settings, Asimov)

| M [GeV] | OLD | BKG_ONLY | NEW | EXACT |
|---|---|---|---|---|
| 200 | −4.5σ | +0.5σ | 0.00σ | 0.00σ |
| 300 | +1.0σ | +0.02σ | 0.00σ | 0.00σ |
| 400 | +1.2σ | 0.00σ | 0.00σ | 0.00σ |
| 500 | −0.5σ | 0.00σ | 0.00σ | 0.00σ |

(μ = 0; for μ = 1 similar, e.g. −3.5σ and +0.16σ at 200 GeV). The sizes differ from the real workspace (where `OLD` gave −1.7σ to
+2.2σ) because the truth here is a single steep 4-parameter function, but the pattern is the same.
