# phyloslopes

Notes and benchmarks comparing ways to implement phylogenetic random-effect
terms in mixed models: RTMB with dense covariance, sparse tip-only
precision, sparse full-tree (root-dropped) precision, and edge-based
Z-matrix parameterizations; glmmTMB's `propto()`; and `phyr::pglmm_compare()`
as an independent check -- plus phylogenetic spline extensions (additive,
separable, tensor-product). Write-up: `phyloslopes.qmd`.

## Files

- `phyloslopes_utils.R` -- shared helpers: RTMB negative-log-likelihoods
  (`nllfun_*`), `phylo.to.Z()`, `drop_mrf_root()`, TMB-object accessor
  methods (`fixef.TMB`/`logLik.TMB`/`AIC.TMB`), and a `logLik.glmmTMB`/
  `AIC.glmmTMB` hack (registered via `registerS3method()`) that bypasses
  glmmTMB's NA-on-convergence-warning safety check so a flagged fit can
  still be compared with `bbmle::AICtab()`
- `phyloslopes_tests.R` -- fits every parameterization once on the real
  (49-species) data and checks they agree; saves the design-matrix/data
  ingredients (not the fits) to `phyloslopes_bench_setup.rds`
- `phyloslopes_bench1.R` / `phyloslopes_bench1_plot.R` -- microbenchmarks
  fitting time per method (reads `phyloslopes_bench_setup.rds`, doesn't
  refit); produces `phyloslopes_bench.png` and `phyloslopes_bench1.png`
- `phyloslopes_testinv.R` / `phyloslopes_testinv_plot.R` -- benchmarks
  `vcv()`/`solve()`/`MRFtools::mrf_penalty()` scaling across random
  subsamples (100-7500 tips) of the full fishtree phylogeny, parallelized
  via `future`/`furrr`; checkpoints to `phyloslopes_testinv.rds`
- `phyloslopes_tiny.R` -- clean/minimal 5-tip worked example: (1) how X,
  Z, and the penalty matrices (Sigma/Q) look under several phylogenetic
  random-effect parameterizations, and (2) the primary example -- fits
  and compares additive, separable, and tensor-product-smooth structures
  (phylo x low-k thin-plate spline) on one small simulated dataset; saves
  design matrices/basis objects/simulated data to `phyloslopes_tiny.rda`
- `phyloslopes_tiny_explore.R` -- archived record of the larger
  identifiability investigation `phyloslopes_tiny.R` grew out of: why the
  naive Kronecker-*sum* tensor-product construction has a genuine
  sigma_1/sigma_2 non-identifiability, and why the null/range-decomposed
  Kronecker-*product* construction avoids it
- `phyloslopes_tiny_predcovs.R` -- for each of the tiny example's
  structures (as estimated, plus additive/separable refit with their
  random-effects variance-component hyperparameters fixed at their true
  simulated values instead of estimated), predicts every species' curve
  value and slope at a new covariate value and the joint covariance
  matrix of those predictions; saves to `phyloslopes_tiny_predcovs.rda`
- `phyloslopes_linear.R` -- fits null/additive/independent/separable
  *linear* (log_bm) models on the real data (null/additive/independent
  via glmmTMB; separable via the RTMB `nllfun_sep` machinery, since
  `propto()` can't express its joint phylo x 2x2-trait Kronecker
  covariance); saves the fitted models to `phyloslopes_linear.rda`
- `phyloslopes_combo.R` -- the spline analog of `phyloslopes_linear.R`:
  fits null/additive/separable/tensor *smooth* (thin-plate spline in
  log_bm) models on the real data (null via glmmTMB; the rest via RTMB
  `nllfun_spline_*`, sharing the same machinery so they stay directly
  nested/comparable); saves to `phyloslopes_combo.rda`
- `phyloslopes.qmd` -- the write-up itself
- `README_tensor.qmd` -- standalone deep-dive on the tensor-product-smooth
  construction (`nllfun_spline_tensor`): Wood (2006, sec. 4.1.8)'s
  single-penalty-vs-multiple-term-penalty distinction, why the naive
  `te()`-extracted construction and the null/range-decomposed one differ,
  the derivation of the "hybrid" range-block fix (separate null-space
  scales, kept from the original construction, plus a Wood-style
  multiple-term penalty for the range block instead of a pure Kronecker
  product) now implemented in `nllfun_spline_tensor`, a validation of the
  `te()`-vs-RTMB equivalence against an independent worked example (Clark
  2024's blog post, `clark_2024_ex.R`), and the null block's optional
  `cor_null` correlation parameter (mapped to 0 by default -- see
  `nllfun_spline_tensor`'s header comment). Wired into `_targets.R`
  (`readme_tensor_html`, with `clark_2024_ex.R` as a file dependency).
- `clark_2024_ex.R` -- reproduces the worked example from Clark (2024)'s
  blog post verbatim (see `phyloslopes.bib`), used by `README_tensor.qmd`
  to validate the `te()` + `nllfun_tensor` construction against an
  independent, well-posed case (unlike the real 49-species data, where
  `gam()`'s own smoothing-parameter selection doesn't converge).

## Running things

`_targets.R` orchestrates most of the scripts above (`targets::tar_make()`).
Every sourced/rendered file has a small "source-tracking" sentinel target
(`utils_r`, `tiny_r`, `linear_r`, `combo_r`, `tests_r`, `bench1_r`,
`bench1_plot_r`, `qmd_source`, `readme_tensor_source`, `clark_r`, ...) that
just tracks that file's hash, so editing e.g. `phyloslopes_utils.R`
correctly invalidates every downstream target that sources it --
`targets`' static dependency analysis can't see *inside* a `source()`d
script's own nested `source()`/render() calls.

Two targets (`bench_setup_rds`: ~15 model refits + consistency checks;
`bench_rds`: a ~1000-fit microbenchmark) are expensive (~10+ min each) and
are marked `cue = "never"` so they don't rerun automatically -- but
`targets` still builds a `cue = "never"` target the *first* time (no
recorded output yet). After that, force a rerun by hand with
`targets::tar_invalidate()` once the fitting code actually changes.

`phyloslopes_testinv.R` / `phyloslopes_testinv_plot.R` are **not** wired
into `_targets.R`: `phyloslopes_testinv.R` is a multi-hour, 10-core
parallel job meant to be launched and checkpointed by hand (`nohup
Rscript phyloslopes_testinv.R &`), not auto-triggered by pipeline
dependency tracking, and its plot script isn't tracked either since its
input (`phyloslopes_testinv.rds`) isn't a managed target. Its size sweep
deliberately stops at 7500 tips -- 10000 is memory-prohibitive on this
machine. `phyloslopes.qmd` still displays its output
(`phyloslopes_testinv_plot.png`), which must already exist on disk before
rendering.

## `nllfun_*` functions (`phyloslopes_utils.R`)

All are RTMB negative-log-likelihoods with the same calling convention:
`getAll(params, chdat_x)` pulls in both the current parameter values and
whatever fixed design matrices/data the model needs (nothing is picked up
as a free variable from the calling environment), then each returns a
single penalized-negative-log-likelihood scalar for `MakeADFun()`. Listed
here are the 9 actively used ones, grouped by what they model; two more
(`nllfun_independent_slopes`, `nllfun_tensor`) are defined but excluded --
see the note at the end.

### Standalone phylogenetic covariance

Compared for fitting-time efficiency in `phyloslopes.qmd`'s "Cost of
fitting intercept-only phylogeny terms" section and `phyloslopes_bench1.R`
(alongside glmmTMB's `propto()`); all three fit the identical model (a
single phylogenetically-correlated random intercept), differing only in
how the phylogenetic covariance is represented:

- **`nllfun1`** -- dense covariance: `b ~ MVN(0, vcmat)` via `dmvnorm()`
  on the full `ntip x ntip` matrix.
- **`nllfun_prec`** -- sparse precision matrix via `dgmrf()`; used twice
  with different inputs -- once with the tip-only precision
  (`solve(vcmat)`, sparsified) and once with the all-internal-nodes
  precision (`MRFtools::mrf_penalty(..., internal_nodes = TRUE)` with the
  root dropped via `drop_mrf_root()`).
- **`nllfun_edge`** -- edge-based `Z` (`phylo.to.Z()`, tips-to-ancestral-edges
  mapping scaled by `sqrt(edge.length)`): `b` is plain iid `N(0, tau^2)`,
  no covariance/precision matrix needed at all, since the edge-weighted
  `Z` already induces the right tip covariance by construction.

### Linear x phylogenetic models

Combine phylogeny with a plain linear (`log_bm`) fixed effect, giving each
species a phylogenetically-correlated intercept and/or slope deviation.
Used in `phyloslopes.qmd`'s "linear effects" section and
`phyloslopes_linear.R`:

- **`nllfun_sep`** -- separable (Kronecker-product) phylo x intercept-slope
  covariance: `b` (ntip x 2) gets `Sigma = vcmat %x% Sigma_trait`, with
  `Sigma_trait` a free `unstructured(2)` covariance (two log-sds + a
  correlation), via `dseparable()`.
- **`nllfun_dense_slopes`** -- the same separable model, but built as one
  dense `2*ntip`-length multivariate normal (`kronecker(vcmat, Sigma2)`
  passed straight to `dmvnorm()`) rather than via `dseparable()` --
  brute-force ground truth for `nllfun_sep`/`nllfun_edge_slopes`.
- **`nllfun_edge_slopes`** -- the same separable covariance again, but
  realized edge-by-edge: each edge gets an iid 2-vector (intercept-
  innovation, slope-innovation) with shared `2x2` covariance `Sigma2`,
  and Brownian summation up the tree reproduces `Sigma2 %x% vcmat`
  exactly, without ever forming a `ntip x ntip` matrix.

### Spline x phylogenetic models

Combine phylogeny with a thin-plate-spline (TPS) smooth in `log_bm`
instead of a linear term. Used in `phyloslopes.qmd`'s "combinations of
thin-plate splines and phylogenetic effects" section and
`phyloslopes_combo.R`; see `README_tensor.qmd` for the full derivation:

- **`nllfun_spline_additive`** -- independent (summed, not crossed) iid
  spline-wiggle random effect plus a phylogenetically-correlated random
  intercept -- matches McGillycuddy et al.'s `propto()+s()` model exactly.
- **`nllfun_spline_separable`** -- additive's structure plus a second,
  phylogenetically-correlated copy of the wiggly spline coefficients
  (one extra diagonal scale across the wiggly dimensions, added on top of
  the still-present iid copy) -- nests additive as that scale -> 0.
- **`nllfun_spline_tensor`** -- a genuine tensor-product smooth of
  phylogeny x `log_bm`: the spline's own null space (constant + linear)
  and range space (wiggly) are split first (via `smooth2random()`, not a
  separate eigendecomposition), then each is crossed with phylogeny
  separately. The null block gets two independent per-direction scales
  (`logpsd_null`, plus an optional `cor_null` correlation, mapped to 0 by
  default) rather than one shared scale; the range block gets a Wood
  (2006, sec. 4.1.8)-style multiple-term penalty (phylo-shrinkage plus the
  TPS's own smoothness-eigenvalue shrinkage, each with its own scale)
  rather than a pure Kronecker product.

### Excluded: historical/superseded functions

- **`nllfun_independent_slopes`** -- independent (uncorrelated)
  phylogenetic intercept + phylogenetic slope penalties. Never called
  anywhere in the project; the "independent" linear model is instead fit
  via glmmTMB's `propto()` in `phyloslopes_linear.R`, which does the same
  job without a bespoke RTMB function.
- **`nllfun_tensor`** -- naive Kronecker-*sum* tensor-product random
  slopes (`Q = (1/sigma_1^2)*tp1 + (1/sigma_2^2)*tp2`, built from
  `tensor.prod.penalties()`). Only called from `phyloslopes_tiny_explore.R`
  (an archived investigation) and `README_tensor.qmd` (validating the
  naive `te()`-extraction approach before it's superseded by
  `nllfun_spline_tensor`'s null/range-decomposed construction) -- not used
  by the current pipeline. On a trivial (identity, no-null-space) trait
  penalty like `diag(2)`, this construction has a genuine, reproducible
  `sigma_1`/`sigma_2` non-identifiability (see `phyloslopes_tiny_explore.R`).
