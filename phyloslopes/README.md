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
