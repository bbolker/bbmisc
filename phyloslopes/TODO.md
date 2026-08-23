# TODO

* check PDF output for glitches
* consider other fits from  @garlandDoes1993, e.g. body mass vs hind limb, which has a more interesting pattern
* compare fitted spline covariance models? compare naive spline (no phylog) with phylog spline? (Predictions for this example will probably all be extremely similar to naive/no-phylog spline - v. small difference log-likelihood differences)
* LOOCV/cAIC?
* are there more interesting (larger) PGLMM examples to try?  What is the typical scale of phylogenies/models people are fitting these days, and how often would random slopes/tensor-product smooths be relevant?
* combine/integrate stuff from `phyloglmm`?
* downstream: can `MRFtools` do a better job maintaining sparsity throughout the `mrf_full`/`mrf_full_sparse` benchmark steps?
* try on larger examples where everything won't collapse to singular/simpler fits? evaluate performance on simulated data (time/model fit scaling?))
* explore the pattern of Hessian sparsity for `nllfun_edge`/`fit_edge`
  (edge-based `Z`, `phylo_to_Z()`) vs. `nllfun_prec`/`fit_prec_allnodes_noroot`
  (sparse `Q`, all internal nodes included but root dropped) vs.
  `nllfun1`/`fit_dense` (dense `Sigma = vcmat`) -- standard (R)TMB advice is
  that approaches mapping each observation to a *single* latent variable
  (definitely `fit_prec_allnodes_noroot`'s sparse Q; maybe `fit_dense`'s
  dense Sigma too, at least in the sense of not adding extra overlapping
  structure) should give a sparser/better-conditioned joint Hessian than one
  mapping each observation to several *overlapping* latent variables (the
  edge-based `Z`, where every observation touches every edge on its root
  path). We don't see that pattern in `phyloslopes_bench1.R`'s timings --
  maybe the example (49 tips) just isn't big enough for the difference to
  show up yet?
* the ML vs. REML distinction (`gam(..., method=)` vs. our RTMB fits'
  implicit ML, since only `b` is integrated out via Laplace, not `beta` --
  see the `phyloslopes_combo.R` `te()`-vs-`gam()` comparison) may matter
  more for "new-style" random effects (*sensu* Hodges) -- where the
  random-effect variances themselves (e.g. spline wiggliness) are the
  focal quantity of interest -- than for "old-style" random effects, where
  they're usually a nuisance parameter to be integrated out/ignored
* Discovered via `phyloslopes_tiny.R`: `MRFtools::mrf_penalty(tree, "brownian",
  internal_nodes = FALSE)` does **not** equal `solve(vcv(tree))`, even though
  both present themselves as "the tips-only phylogenetic precision matrix."
   * `internal_nodes = FALSE` Schur-complements the internal nodes (including
     the root) out of the rank-deficient joint precision. That implicitly
     gives the root an *improper* (flat/infinite-variance) prior, which
     induces nonzero covariance between tips in different top-level clades --
     `vcv()`/`solve(vcmat)` has those as an exact zero.
   * The tips-only precision route actually validated/used elsewhere in this
     project (`phyloslopes_tests.R`'s `Qprec_tip`, the qmd's `phylo_prec`
     chunk) is always `Matrix(solve(vcmat), sparse = TRUE)` -- never
     `mrf_penalty(..., internal_nodes = FALSE)`.
   * The *correct* MRFtools-based route (checked numerically in
     `phyloslopes_tiny.R`) is `internal_nodes = TRUE` + `drop_mrf_root()`: inverting
     that joint precision and taking the tip x tip block reproduces `vcmat`
     exactly.
   * Add a note about this to the "Cost of computing Sigma / Q"
     section of phyloslopes.qmd (~L68-86) -- it currently presents
     `mrf_tips`/`mrf_tips_sparse` timings without flagging that they measure a
     *different model* from the `solve(vcmat)`-based fits used everywhere
     else in the doc, not just a faster way to compute the same matrix.
