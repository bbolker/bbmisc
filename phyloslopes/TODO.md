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
* decide what if anything needs to be kept from last item below
* the ML vs. REML distinction (`gam(..., method=)` vs. our RTMB fits'
  implicit ML, since only `b` is integrated out via Laplace, not `beta` --
  see the `phyloslopes_combo.R` `te()`-vs-`gam()` comparison) may matter
  more for "new-style" random effects (*sensu* Hodges) -- where the
  random-effect variances themselves (e.g. spline wiggliness) are the
  focal quantity of interest -- than for "old-style" random effects, where
  they're usually a nuisance parameter to be integrated out/ignored
* investigate a model more complex than the current "hybrid" tensor
  (`nllfun_spline_tensor`), which gives the phylo-intercept and
  phylo-linear-slope null-space directions *independent* (diagonal,
  uncorrelated) penalties (`logpsd_null`, length 2) -- what about letting
  them *correlate*, the way `nllfun_sep`/`nllfun_dense_slopes` give
  intercept/slope a full unstructured 2x2 covariance instead of
  independent scales? Probably not worth it for the real 49-species
  example (already pushing the data's limits), but worth having in the
  model space for completeness
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

## 3. ~~Continue the tiny (5-tip) example~~ done -- see findings below

`nllfun_tensor` was wired up (`phyloslopes_utils.R`) and fit to data
simulated under the Kronecker-sum tensor-product model, along with an
extended investigation (sections (a)-(e)) into *why* it doesn't recover the
true parameters. That investigation, and the naive Kronecker-sum
construction itself, now live only in `phyloslopes_tiny_explore.R` -- it's
not worth using (see below). The tiny worked example it grew out of,
`phyloslopes_tiny.R`, has since been repurposed: part 1 is still a
clean/minimal illustration of how `X`, `Z`, and the penalty matrices
(Sigma/Q) look for small examples under different phylogenetic
random-effect parameterizations; part 2 (the primary content now) fits and
compares additive, separable, and tensor-product (null/range-decomposed,
*not* naive Kronecker-sum) combinations of phylogenetic variation with a
low-k (k=4) thin-plate spline on one small simulated data set --
`nllfun_spline_additive`/`_separable`/`_tensor`, the same three structures
phyloslopes.qmd's "phylogenetic splines" section fits to the real
49-species data, scaled down to something fully inspectable. Findings from
the investigation that produced this design:

- The naive `tensor.prod.penalties()` + single shared `dgmrf()` construction
  (`nllfun_tensor`) has a genuine, reproducible non-identifiability between
  its two variance-component scales (`sigma_1`, the phylo-structured
  component; `sigma_2`, the "flat-ridge" component) -- one of the two always
  runs off toward a boundary (confirmed via multi-start: 15-40 diverse
  starting points all converge to the identical optimum, so it's the true
  global optimum, not a bad start). This holds regardless of `ntip` (checked
  at 5 and 10), regardless of `nobs` vs. `n` random effects (checked with
  `nrep` replicates/tip so `sigma_resid` isn't just collapsing to the
  interpolation boundary), and across repeated simulated datasets (8/10
  independent seeds at `ntip = 10` showed the divergence; the
  correctly-specified tensor fit never once beat a deliberately misspecified
  separable fit's logLik, despite having fewer parameters).
- It is **not** a fixed-effect/random-effect (`X`/`tX`) aliasing problem --
  refitting with `X` restricted to intercept-only (no fixed slope competing
  with the random tensor-slope term) does not fix it.
- The actual cause: `Strait = diag(2)` (the trait-side "penalty" for a plain
  intercept+slope model) is a trivial, structureless, full-rank identity
  matrix -- indistinguishable from a generic flat ridge. Summing a
  phylo-structured term (`kron(Sphylo, I_2)`) with a term that's already
  just `I` gives the optimizer no way to tell "shrink toward phylogenetic
  similarity" apart from "shrink generically," so one scale trades off
  against the other without cost.
- This does **not** generalize to a genuine tensor-product smooth of
  phylogeny x a thin-plate spline, even using the same naive
  `tensor.prod.penalties()` + shared-`dgmrf()` construction with the *raw*
  (rank-deficient) TPS penalty: because the TPS penalty has real internal
  (null vs. range space) structure, unlike `diag(2)`, both variance
  components land at finite, stable values (confirmed via profiling +
  multi-start). The qmd's actual `nllfun_spline_tensor` (eigendecompose
  null/range space first, give each its own separate Kronecker-*product*
  term rather than summing into one shared precision) is also well-posed.
- Applying that same "give the (here, trivial) null space its own
  Kronecker-product term, no shared sum" recipe back to the plain
  intercept+slope case (`nllfun_tensor_nullspace`) is also well-identified
  on the same data that broke `nllfun_tensor` -- it's structurally the
  separable-model family (diagonal, not fully unstructured, trait
  covariance).
- General lesson: the pathology is specific to *naive Kronecker-sum with a
  trivial (identity-like, no-null-space) trait penalty*, not a generic
  property of phylo x trait tensor-product models. Flagged in
  phyloslopes.qmd's tensor-product section (~L299-332) with a pointer here.
