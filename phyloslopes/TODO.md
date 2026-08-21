# TODO

* clean up this file!
* keep working on `phyloslopes`
   * clean up/integrate/throw out junk
   * prettify tables (check PDF output for glitches)
* downstream: can `MRFtools` do a better job maintaining sparsity throughout the `mrf_full`/`mrf_full_sparse` benchmark steps?
* `phylo.to.Z` FIXME in `phyloslopes_utils.R`: rename to `phylo_to_Z`?
* ~~see if we can do the tensor-product construction via `te()` and pull out the stuff we need, rather than doing the eigendecomposition ourselves? Is this problem somehow specific to what we're trying to do here? (How does it differ from @clarkPhylogeneticSmoothing2024 ?)~~ done -- see `README_tensor.qmd`. Yes, straightforwardly (`nllfun_tensor`, reusing `smooth.construct()`'s raw `$X`/`$S` directly, no decomposition needed -- `S1` alone keeps the summed `Q` full rank regardless of `S2`'s deficiency). Not specific to our problem: it's Wood (2006, sec 4.1.8)'s standard "multiple term penalty" recipe, the same one mgcv's own `te()` uses by default. It turned out to still lose to a *hybrid* construction -- keep `nllfun_spline_tensor`'s separate null-space scales, but rebuild its range block the multiple-term way too -- which is now implemented (`Qr_phylo`/`Qr_smooth` in `nllfun_spline_tensor`, built via `smooth2random()`'s `trans.D` rather than a separate `eigen()` call).
* figure out how to implement the "single penalty" approach (what we call `separable`) for the case where one or both components have a genuinely *rank-deficient* penalty -- i.e. a `separable`-style model built directly from a raw (rank-deficient) marginal penalty, unlike `nllfun_sep`/`nllfun_dense_slopes` (trait covariance always full rank by construction) or `nllfun_spline_separable` (uses the already-pre-whitened `Xr`, not the raw penalty). Probably easiest via `smooth2random()` (its `Xf`/`Xr`/`trans.D`, as now used for `nllfun_spline_tensor`'s range block), passing the resulting full-rank/reduced penalty components into something like the current `nllfun_sep`, but with an added `Xf` (null-space fixed-effect) term.
* can we do the TP+phylo tensor smooth with `internal_nodes = TRUE` (with root dropped)? Does that improve efficiency as in the intercept-only case?
* try on larger examples where everything won't collapse to singular/simpler fits?
* understand Schur-complement vs `solve(vcmat)` differences (see below)
* explore the pattern of Hessian sparsity for `nllfun_edge`/`fit_edge`
  (edge-based `Z`, `phylo.to.Z()`) vs. `nllfun_prec`/`fit_prec_allnodes_noroot`
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
* decide what if anything needs to be kept from last item below
* the ML vs. REML distinction (`gam(..., method=)` vs. our RTMB fits'
  implicit ML, since only `b` is integrated out via Laplace, not `beta` --
  see the `phyloslopes_combo.R` `te()`-vs-`gam()` comparison) may matter
  more for "new-style" random effects (*sensu* Hodges) -- where the
  random-effect variances themselves (e.g. spline wiggliness) are the
  focal quantity of interest -- than for "old-style" random effects, where
  they're usually a nuisance parameter to be integrated out/ignored
* ~~investigate a model more complex than the current "hybrid" tensor
  (`nllfun_spline_tensor`), which gives the phylo-intercept and
  phylo-linear-slope null-space directions *independent* (diagonal,
  uncorrelated) penalties (`logpsd_null`, length 2) -- what about letting
  them *correlate*, the way `nllfun_sep`/`nllfun_dense_slopes` give
  intercept/slope a full unstructured 2x2 covariance instead of
  independent scales?~~ done -- see `README_tensor.qmd`'s "null block"
  section. No new function needed: `nllfun_spline_tensor`'s null block now
  takes a `cor_null` parameter (`kronecker(phylomat, us2$corr(cor_null))`
  in place of the old hard-coded `diag(Kn)`), and the old
  `nllfun_tensor_nullspace` was removed entirely -- it, and this
  generalization, are both just `nllfun_sep` (+ `map=`) reused, since the
  null-space design *is* `model.matrix(~1+log_bm)` up to a linear
  reparameterization. As guessed, not worth it on the real 49-species
  example: freeing the correlation gains a little logLik but lands on a
  correlation boundary (rho ~ -1), so `cor_null`/`corval` stay mapped to 0
  by default everywhere
* "null" (as in "null-space") for the phylo-intercept/phylo-slope terms in
  `nllfun_spline_tensor` is potentially misleading terminology: those
  terms are in the null space of the *thin-plate spline's* wiggliness
  penalty, but they are NOT unpenalized overall -- they still carry a
  phylogenetic penalty (`logpsd_null`). "Null" could easily be misread as
  "fixed effect/no penalization," which is wrong here (the phylogenetic
  component itself has no null space of its own, given proper handling
  like dropping the root)

## 2. Document the MRFtools tips-only vs. solve(vcmat) precision mismatch

Discovered via `phyloslopes_tiny.R`: `MRFtools::mrf_penalty(tree, "brownian",
internal_nodes = FALSE)` does **not** equal `solve(vcv(tree))`, even though
both present themselves as "the tips-only phylogenetic precision matrix."

- `internal_nodes = FALSE` Schur-complements the internal nodes (including
  the root) out of the rank-deficient joint precision. That implicitly
  gives the root an *improper* (flat/infinite-variance) prior, which
  induces nonzero covariance between tips in different top-level clades --
  `vcv()`/`solve(vcmat)` has those as an exact zero.
- The tips-only precision route actually validated/used elsewhere in this
  project (`phyloslopes_tests.R`'s `Qprec_tip`, the qmd's `phylo_prec`
  chunk) is always `Matrix(solve(vcmat), sparse = TRUE)` -- never
  `mrf_penalty(..., internal_nodes = FALSE)`.
- The *correct* MRFtools-based route (checked numerically in
  `phyloslopes_tiny.R`) is `internal_nodes = TRUE` + `drop_mrf_root()`: inverting
  that joint precision and taking the tip x tip block reproduces `vcmat`
  exactly.
- Add a note about this to the "Cost of computing $\bSigma$ / $\bQ$"
  section of phyloslopes.qmd (~L68-86) -- it currently presents
  `mrf_tips`/`mrf_tips_sparse` timings without flagging that they measure a
  *different model* from the `solve(vcmat)`-based fits used everywhere
  else in the doc, not just a faster way to compute the same matrix.
