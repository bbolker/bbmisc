# Forstmeier & Schielzeth (2011) simulation reproduction

Reproduces and extends the type-I-error simulation from:

> Forstmeier, W. and Schielzeth, H. (2011). Cryptic multiple hypotheses
> testing in linear models: overestimated effect sizes and the winner's
> curse. *Behavioral Ecology and Sociobiology* 65:47-55.
> doi:10.1007/s00265-010-1038-5

The simulation generates null data (response independent of all
predictors), fits linear models with 1-6 two-level factors (and,
optionally, all their two-way interactions), and asks how often at least
one predictor comes out "significant" (p < 0.05) by chance alone -- for
three model-simplification scenarios and four significance criteria.

## Pipeline

Built with [`targets`](https://books.ropensci.org/targets/) as a
file-dependency ("Makefile in R") pipeline: run

```r
targets::tar_make(script = here::here("forstmeier_sim", "_targets.R"),
                   store = here::here("forstmeier_sim", "_targets"))
```

from anywhere in the repo. `output/corrections_results.rds` is the single
simulation output every figure reads from; it's cached and skipped on
rerun if it already exists -- delete it to force a genuine recompute
after changing the simulation code.

Single-threaded BLAS/OpenMP is enforced inside parallel workers
(`RhpcBLASctl`) to avoid oversubscribing cores.

## Simulation design

For N in {30, 50, 200}, m in 1-6 balanced, independently-randomised
two-level factors, with and without all two-way interactions (`k`
predictors total), each of 5000 replicates:

0. Predictors are coded as actual factors under the global sum-to-zero
   contrasts `options(contrasts = c("contr.sum", "contr.poly"))` (set
   inside `one_rep_corrections()` -- see its comment for why it can't just
   be set once in the calling script), so that in the presence of
   interactions each main-effect test is of the effect at the population
   mean of the other factors, rather than at an arbitrary reference level.
1. Fits the full model to null data (the **unselected** scenario).
2. Simplifies it with unrestricted backward `step()` (the **selected
   (all)** scenario) and with `step()` restricted to discarding only
   interactions, never main effects (the **selected (interax only)**
   scenario; see `step_interactions_only()`'s `scope = list(lower,
   upper)`).
3. Screens each of the three resulting models under four significance
   criteria: **uncorrected** (any p < 0.05), **Dunn-Sidak**, **Holm**, and
   **single-step max-\|T\|** (a Tukey-HSD-style correction using the
   model's actual multivariate-t null -- see `maxT_correction.R`); each
   yields both an experimentwise indicator (`sig_*`: is at least one term
   significant) and a per-parameter count (`n_sig_*`: how many terms are).

For the two *selected* scenarios, each of the above is computed two ways,
recorded in a `multcomp_k` column: **"minimal"** (calibrated on the
selected model's own surviving term count -- how a user who isn't
thinking about the search that produced it would do it) and
**"maximal"** (calibrated on the pre-selection full model's $k$, correlation
structure and residual df -- what Forstmeier & Schielzeth's Fig. 3b does,
applying the full model's Bonferroni threshold to the minimal model
rather than recomputing it from the minimal model's own smaller term
count). `multcomp_k` is `"none"` where the distinction doesn't apply: the
uncorrected criterion, and the unselected scenario (which has only one
calibration since fit and reference model coincide).

The plots additionally derive a fifth, uncorrected **per-parameter error
rate**, E[V]/k (V = how many of the original k full-model predictors are
significant), from the `n_sig_raw_*` counts divided by `k` --
`build_corrections_by_scenario_summary()` in
`plot_corrections_by_scenario_core.R`. Unlike the experimentwise
indicator, this is trivially ~alpha for the unselected scenario (each of
the k tests is independently exact) but inflated above alpha for the
selected scenarios: `step()`'s AIC-based retention and the final refit's
significance check draw on overlapping evidence about the same predictor,
so "retained" and "significant" are positively correlated rather than
independent.

## Code layout

- `R/simulate_fig1.R` -- data generation (`make_predictors`,
  `build_formula`), the sample-size grid (`sim_N_values`), and the two
  `step()`-based simplification strategies (`step_default()`,
  `step_interactions_only()`).
- `R/simulate_corrections.R` -- `screen_criteria()` (the four criteria,
  minimal/maximal calibration, both experimentwise and per-parameter
  counts), `one_rep_corrections()`, `run_condition_corrections()`.
- `R/maxT_correction.R` -- `maxT_crit()`: the single-step correction's
  critical value, via direct Monte Carlo simulation from the multivariate-t
  null rather than `mvtnorm::qmvt()`'s numerical integration (impractical
  at simulation scale for k ~ 20 correlated tests). See the file's roxygen
  docs for the full derivation and references.
- `R/run_corrections_core.R`, `R/run_corrections_simulation.R` -- the
  simulate-and-cache driver and its thin wrapper.
- `R/plot_fig1_core.R` (+ `make_plot.R`) -- 2-panel a/b reproduction of the
  paper's Fig. 1 layout, parametrized by which pair of scenario columns to
  plot.
- `R/plot_corrections_by_scenario_core.R` (+
  `make_corrections_by_scenario_plot.R`, `make_corrections_by_method_plot.R`,
  `make_corrections_by_N_plot.R`) -- three scenario x criterion
  comparisons built on one shared summary (`build_corrections_by_scenario_summary()`,
  which pivots out the `multcomp_k` column and derives the per-parameter
  criterion from `n_sig_raw_*`/`k`), all with-interactions cases only and
  all showing both the minimal and maximal calibration (linetype/shape): a
  1-column x n-row uncorrected panel (colour = N) beside an n-row x
  4-column facet_grid of the other four criteria (colour = N) combined
  with `patchwork`; a single facet_grid (rows = scenario, columns = N,
  colour = criterion); and its transpose (rows = N, columns = criterion,
  colour = scenario).
- `R/graphics_utils.R` -- `logit_breaks()` and the shared Okabe-Ito
  palette; every plot uses a logit y-scale.
- `smoke_test.R` -- correctness and timing checks, run standalone.

## Outputs (`output/`)

- `fig1.png` -- the paper's Fig. 1 layout (unselected vs. selected via
  unrestricted `step()`).
- `fig1_corrections.png` -- uncorrected vs. Dunn-Sidak vs. Holm vs.
  single-step max-\|T\|, unselected model only.
- `fig1_corrections_by_scenario.png` -- all three scenarios crossed with
  all five criteria (uncorrected, Dunn-Sidak, Holm, single-step max-\|T\|,
  per-parameter), uncorrected panel + other-criteria grid combined via
  `patchwork`, colour = N, linetype/shape = `multcomp_k`, with-interactions
  cases only.
- `fig1_corrections_by_method.png` -- the same comparison as a single
  facet_grid (rows = scenario, columns = N), colour = criterion,
  linetype/shape = `multcomp_k`, with-interactions cases only.
- `fig1_corrections_by_N.png` -- rows = N, columns = criterion, colour =
  scenario, linetype/shape = `multcomp_k`, with-interactions cases only.
- `corrections_results.rds` -- the cached simulation output underlying all
  of the above; includes both minimal and maximal k-counting and both
  experimentwise and per-parameter counts for the two selected scenarios
  (see Simulation design).

## Key findings so far

1. **The uncorrected/step()-simplified comparison reproduces the paper's
   Fig. 1**: the unselected full-model test tracks the independence
   formula alpha' = 1-(1-alpha)^k reasonably closely, while backward
   `step()` simplification inflates the type I error rate well beyond
   that, especially at low N relative to k (e.g. N=50 with six factors
   and their interactions, k=21).

2. **Dunn-Sidak correction does not show the same N=50 conservatism as
   the raw test.** In the uncorrected comparison, the N=50 with-
   interactions curve sits noticeably *below* both the N=200 curve and
   the theoretical alpha' line (correlation among the k test statistics --
   via both a shared residual-variance estimate and predictor collinearity
   -- makes "any significant" less likely than independence predicts, more
   so at low N). This conservatism mostly doesn't carry through to the
   Dunn-Sidak-corrected version: at k=21, N=50's empirical family-wise
   error rate is ~0.0496 -- indistinguishable from nominal 0.05. The
   degree of conservatism from correlated tests is threshold-dependent: a
   synthetic check isolating the shared-denominator effect alone gives a
   ~10% conservative bias at the raw alpha=0.05 threshold but only ~6% at
   the much smaller per-test threshold Dunn-Sidak uses for a 21-way
   correction. At the more extreme N=30/k=21 (df=8), this same mechanism
   is strong enough to show up clearly: Dunn-Sidak and Holm both swing
   noticeably conservative (down to ~0.03-0.04) at high m, while
   single-step max-\|T\| -- which uses the replicate's own exact
   correlation structure rather than an independence-based bound -- stays
   close to nominal throughout.

3. **Restricting `step()` to discard only interactions (never main
   effects) barely reduces the type I error inflation.** Across the
   `selected (all)` vs. `selected (interax only)` scenarios (see
   `fig1_corrections_by_scenario.png`), the restricted version is consistently slightly lower than the unrestricted
   one across every with-interactions condition (e.g. at N=50, k=21: 0.831
   vs. 0.844) but the gap is small (~1 percentage point at the most
   extreme condition) relative to the inflation itself. Most of the
   "cryptic multiple testing" problem comes from having many interaction
   terms available to prune in the first place, not from the additional
   freedom to also drop main effects.

4. **Calibrating a correction on the selected model's own surviving terms
   ("minimal" k) substantially undercorrects relative to calibrating on
   the pre-selection full model ("maximal" k, the paper's Fig. 3b
   approach).** At N=50, k=21 (selected-all scenario), minimal-k Dunn-
   Sidak/Holm/max-\|T\| all still show ~0.31-0.32 type I error, while
   maximal-k brings them down to ~0.18-0.19 -- both well above nominal
   0.05, but maximal is roughly half the inflation of minimal. The gap
   between the two calibrations grows with N/k (see
   `fig1_corrections_by_N.png`), and even the "correct" maximal-k
   calibration leaves substantial residual inflation at low N relative to
   k, matching the paper's own finding that P-value adjustment on the
   minimal model doesn't fully solve the problem when the initial full
   model was badly over-fit.
