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
three model-simplification scenarios, several significance criteria, and
(for the with-interactions conditions) both all predictors and main
effects only.

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

For the with-interactions conditions, `one_rep_corrections()` additionally
tracks the same breakdown restricted to just the m main-effect terms
(`is_main_effect_term()`; columns suffixed `_main_<scenario>`):

- **Unselected**: raw/uncorrected only (`sig_raw_main_unselected`), plus a
  cross-check (`sig_raw_mainonly_crosscheck`) against an independently fit
  main-effects-only model (no interaction terms at all) on the same data
  -- these needn't agree replicate-by-replicate (the random factor
  assignment isn't exactly orthogonal in any single finite draw) but
  should match in aggregate over many replicates. In practice the
  cross-check comes out systematically *higher* than the full model's
  own main-effect subset, growing with m (0.195 vs. 0.239 at N=30, m=6):
  the full interactive model burns residual df estimating the (truly
  null) interaction coefficients too, so its residual variance estimate
  is noisier and its critical value larger, making its main-effect tests
  more conservative than a leaner model's would be.
- **Selected (both scenarios)**: the full criteria x calibration cross.
  "Correct"/maximal calibration (`compute_calibration(fit_full, ref_subset
  = is_main_effect_term)`) uses k = m rather than k = m + choose(m,2),
  since these metrics were never intended to also test the interaction
  terms; "naive"/minimal calibration uses the count of the *selected*
  model's own surviving main effects (always m for the interactions-only
  restriction, since main effects are never dropped there; potentially
  fewer for unrestricted `step()`).

## Code layout

- `R/simulate_fig1.R` -- data generation (`make_predictors`,
  `build_formula`), the sample-size grid (`sim_N_values`), and the two
  `step()`-based simplification strategies (`step_default()`,
  `step_interactions_only()`).
- `R/simulate_corrections.R` -- `is_main_effect_term()`, `tidy_fit()`
  (shared tidy-and-filter helper), `compute_calibration()` (k, correlation
  matrix, and max-\|T\| critical value for a given reference model and
  optional term subset) and `screen_criteria()` (applies a calibration to
  a fit's terms, optionally restricted to a subset, under all four
  criteria) as separate steps -- letting the same calibration be reused
  across multiple term subsets without recomputing its Monte Carlo
  critical value redundantly. `screen_raw()` is a lightweight
  uncorrected-only variant (no calibration object needed). `one_rep_corrections()`
  ties these together per replicate (see its doc comment for exactly which
  calibration each scenario/subset combination uses); `run_condition_corrections()`
  runs replicates in parallel.
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
- `R/plot_corrections_by_scenario_core.R` -- `build_corrections_summary()`
  (the shared long-format summary builder, parametrized by which family of
  `sig_*` columns to pivot and what the derived per-parameter rate divides
  by) with two thin wrappers, `build_corrections_by_scenario_summary()`
  (all-effects, denominator k) and `build_main_effects_summary()`
  (main-effects-only, denominator m); `build_other_criteria_plot()` (the
  shared "rows = scenario, columns = criterion" facet_grid); and five plot
  functions built on these:
  - `make_corrections_by_scenario_plot.R` -- all three scenarios x all
    five criteria, uncorrected panel + other-criteria grid combined via
    `patchwork`, colour = N, linetype/shape = `multcomp_k`.
  - `make_corrections_by_method_plot.R` -- the same comparison as a single
    facet_grid (rows = scenario, columns = N), colour = criterion.
  - `make_corrections_by_N_plot.R` -- rows = N, columns = criterion,
    colour = scenario.
  - `make_corrections_by_scenario_minimal_plot.R` /
    `make_corrections_by_scenario_minimal_mainonly_plot.R` -- just
    single-step max-\|T\| and per-parameter (no patchwork needed), for the
    all-effects and main-effects-only metrics respectively, sharing a
    `make_minimal_scenario_plot()` body.

  All with-interactions cases only, all showing both the minimal and
  maximal calibration.
- `R/graphics_utils.R` -- `logit_breaks()` (breaks for a logit-scaled axis;
  `extend = TRUE` by default brackets the data range with the nearest
  candidate break outside it) and `scale_y_logit()` (builds the actual
  scale, computing breaks eagerly from the plotted data and widening the
  scale's own `limits` to match -- necessary since ggplot2 censors any
  break outside a scale's data-derived limits regardless of what the
  breaks function returns); `m_k_scale_x()` (the shared x-axis, labelling
  each tick "m\n(k)"); the shared Okabe-Ito palette. Every plot uses a
  logit y-scale via `scale_y_logit()`, except `make_corrections_plot.R`'s
  `facet_wrap(scales = "free_y")` plot, which still uses the plain
  `logit_breaks()` function directly (each panel needs its own
  independent range, incompatible with a single shared `limits=`).
- `smoke_test.R` -- correctness and timing checks, run standalone.
- `forstmeier_sig.qmd` -- a narrative Quarto report discussing and
  embedding the figures above, rendered to `forstmeier_sig.html`.

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
- `fig1_corrections_by_scenario_minimal.png` /
  `fig1_corrections_by_scenario_minimal_mainonly.png` -- just single-step
  max-\|T\| and per-parameter, no patchwork, for the all-effects and
  main-effects-only metrics respectively; the latter's "Unselected" row
  has no max-\|T\| line (that combination was never computed -- see
  Simulation design).
- `fig1_corrections_by_N.png` -- rows = N, columns = criterion, colour =
  scenario, linetype/shape = `multcomp_k`, with-interactions cases only.
- `corrections_results.rds` -- the cached simulation output underlying all
  of the above; includes both minimal and maximal k-counting, both
  experimentwise and per-parameter counts, and (for the with-interactions
  conditions) the main-effects-only breakdown and cross-check (see
  Simulation design).
- `forstmeier_sig.html` -- rendered from `forstmeier_sig.qmd`.

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
   error rate is ~0.042 -- close to nominal 0.05, unlike the raw test's
   more pronounced dip. The
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

5. **Restricting attention to just main effects doesn't rescue the
   correction.** Naively, calibrating "correctly" for main effects alone
   (k = m, the actual number of tests intended -- see
   `fig1_corrections_by_scenario_minimal_mainonly.png`) might be expected
   to land close to nominal, since m is so much smaller than k. It
   doesn't: at N=30, m=6 (selected-all scenario), main-effects-only
   single-step max-\|T\| is *higher* under "correct" (m-based) calibration
   than under "naive" (surviving-main-effects-based) calibration -- 0.41
   vs. 0.35 -- and both remain far above nominal, and even above the
   corresponding **all**-effects "correct" (k=21-based) rate of 0.32.
   `step()`'s search still explored the full k = m + choose(m,2)
   candidates even though only m of them are now being asked about, so
   calibrating for an exposure of m alone under-corrects for the search
   that actually happened; the per-parameter (uncorrected) rate, in
   contrast, comes out essentially the same whether measured over all
   effects or main effects only, since it isn't trying to correct for
   anything -- `step()`'s selection bias hits both term types at
   essentially the same per-test rate.
