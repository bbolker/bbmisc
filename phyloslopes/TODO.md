# TODO

## 1. Clean up phyloslopes.qmd text

- `### independent` (~L251) is an empty stub section header -- either fill
  in or remove.
- FIXME at ~L241: add an *independent* intercept/slope variation example
  (`propto(1|species, vcmat) + propto(0 + log_bm | species, vcmat)`, or
  `propto(1||species, vcmat)`) -- watch for `||` silently getting converted
  to a `diag` structure instead of two independent terms.
- `## junk` (~L468) is an empty section header at the end of the doc --
  remove or populate.
- Reconcile the existing "## To do" list (~L447) with this file so
  there's one place to check, not two.
- FIXME at ~L86: can `MRFtools` do a better job maintaining sparsity
  throughout the `mrf_full`/`mrf_full_sparse` benchmark steps?
- FIXME at ~L302 (just added): `tensor.prod.model.matrix()` /
  `tensor.prod.penalties()` argument-order mismatch in the
  "### tensor-product" section (`tp` built trait-outer/phylo-inner, `tX`
  built phylo-outer/trait-inner) -- fix the order if this ever gets wired
  into a real `nllfun_tensor` fit.
- `phylo.to.Z` FIXME in `phyloslopes_utils.R`: rename to `phylo_to_Z`?

## 2. Document the MRFtools tips-only vs. solve(vcmat) precision mismatch

Discovered via `phylo_tiny.R`: `MRFtools::mrf_penalty(tree, "brownian",
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
  `phylo_tiny.R`) is `internal_nodes = TRUE` + `drop_mrf_root()`: inverting
  that joint precision and taking the tip x tip block reproduces `vcmat`
  exactly.
- Add a note about this to the "Cost of computing $\bSigma$ / $\bQ$"
  section of phyloslopes.qmd (~L68-86) -- it currently presents
  `mrf_tips`/`mrf_tips_sparse` timings without flagging that they measure a
  *different model* from the `solve(vcmat)`-based fits used everywhere
  else in the doc, not just a faster way to compute the same matrix.

## 3. Continue the tiny (5-tip) example (phylo_tiny.R)

- Currently builds/prints/compares Z and Q under several parameterizations
  and simulates one dataset under the Kronecker-sum tensor-product model
  (`sigma_1 = sigma_2 = 2`, `sigma_resid = 0.1`, `beta0 = 0`, `beta1 = 1`)
  -- no actual model fit yet.
- Wire up an `nllfun_tensor` (RTMB), fit it to the simulated data, and
  check whether `sigma_1`/`sigma_2 >> sigma_resid` avoids the
  saddle-point/identifiability problem noted for this model on the real
  running-speed data (phyloslopes.qmd ~L324).
- Compare recovered `beta0`/`beta1`/`sigma_1`/`sigma_2`/`sigma_resid`
  against the true simulated values.
- Consider also fitting the *separable* (Kronecker-product) model to data
  simulated under the tensor-product (Kronecker-sum) model, to see how
  badly misspecified it is -- mirrors the qmd's existing
  separable-vs-tensor-product discussion, but small enough to check exactly
  (5x5 / 10x10 matrices, not 49-tip ones).
