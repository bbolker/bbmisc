# phyloslopes

Notes and benchmarks comparing ways to implement phylogenetic random-effect
terms in mixed models: RTMB with dense covariance, sparse tip-only
precision, sparse full-tree (root-dropped) precision, and edge-based
Z-matrix parameterizations; glmmTMB's `propto()`; and `phyr::pglmm_compare()`
as an independent check -- plus phylogenetic spline extensions (additive,
separable, tensor-product). Write-up: `phyloslopes.qmd`.

## Files

- `phyloslopes_utils.R` -- shared helpers (RTMB negative-log-likelihoods,
  `phylo.to.Z()`, `drop_mrf_root()`, ...)
- `phyloslopes_tests.R` -- fits every parameterization once and checks they
  agree; saves the design-matrix/data ingredients (not the fits) to
  `phyloslopes_bench_setup.rds`
- `phyloslopes_bench1.R` / `phyloslopes_bench1_plot.R` -- microbenchmarks
  fitting time per method (reads `phyloslopes_bench_setup.rds`, doesn't
  refit); produces `phyloslopes_bench.png` and `phyloslopes_bench1.png`
- `phyloslopes_testinv.R` / `phyloslopes_testinv_plot.R` -- benchmarks
  `vcv()`/`solve()`/`MRFtools::mrf_penalty()` scaling across random
  subsamples (100-7500 tips) of the full fishtree phylogeny, parallelized
  via `future`/`furrr`; checkpoints to `phyloslopes_testinv.rds`
- `phyloslopes.qmd` -- the write-up itself

## Running things

`_targets.R` orchestrates the scripts above (`targets::tar_make()`). Two
targets are expensive (~10+ min: refitting ~15 models, then a 1000-fit
microbenchmark) and are marked `cue = "never"` so they don't rerun
automatically -- force a rerun with `targets::tar_invalidate()` when the
fitting code actually changes.

`phyloslopes_testinv.R` is **not** wired into `_targets.R`: it's a
multi-hour, 10-core parallel job meant to be launched and checkpointed by
hand (`nohup Rscript phyloslopes_testinv.R &`), not auto-triggered by
pipeline dependency tracking. Its size sweep deliberately stops at 7500
tips -- 10000 is memory-prohibitive on this machine.
