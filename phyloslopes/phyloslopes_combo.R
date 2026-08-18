## phyloslopes_combo.R -- fit and compare four models for how phylogeny
## enters a *smooth* (thin-plate spline in log_bm, not just linear) trend,
## on the real (49-species) carniherbi49 data:
##   - null:      no log_bm effect at all, just a phylo-correlated intercept
##   - additive:  fixed intercept+linear trend + iid spline wiggle +
##                phylo-correlated intercept (nllfun_spline_additive)
##   - separable: additive plus a phylogenetically-correlated copy of the
##                wiggle coefficients on top of the (still-present) iid
##                one (nllfun_spline_separable)
##   - tensor:    genuine tensor-product smooth of phylogeny x log_bm,
##                built by eigendecomposing the spline's own null (constant
##                + linear) / range (wiggly) space first, then crossing
##                each with phylogeny separately (nllfun_spline_tensor)
##
## Analogous to phyloslopes_linear.R (which does the same null/additive/
## separable comparison for a *linear*, not spline, log_bm effect) but not
## a parallel structure: there's no "independent" analog here (a phylo
## intercept + non-phylogenetic random *wiggle* is just... additive), and
## unlike phyloslopes_linear.R, additive/separable/tensor are all fit via
## RTMB here (not glmmTMB) even though additive alone could be (as
## phyloslopes.qmd's fit_spline_glmmTMB validates) -- keeping all three on
## the same nllfun_spline_* family, sharing the same b_spline/b_phylo
## machinery, is what makes them directly nested/comparable via
## bbmle::AICtab() the way separable properly nests additive.
##
## Setup and fitting code below is adapted directly from phyloslopes.qmd's
## "phylogenetic splines" section (spline_setup/spline_additive/
## spline_separable/spline_tensor_setup/spline_tensor chunks) -- see that
## section for the full derivation/discussion; only `null` is new here.
##
## Run with: Rscript phyloslopes_combo.R

suppressMessages({
  library(ade4)
  library(ape)
  library(dplyr)
  library(glmmTMB)
  library(RTMB)
  library(reformulas)
  library(mgcv)
  library(Matrix)
  library(tibble)
})

source("phyloslopes_utils.R")

## -- data (mirrors phyloslopes_tests.R / phyloslopes.qmd's setup chunk) -----
data("carniherbi49", package = "ade4")
chtree <- read.tree(text = carniherbi49$tre2)
chtree$tip.label <- fixup_species(chtree$tip.label)
vcmat <- vcv(chtree)
chdat <- carniherbi49$tab2 |>
  mutate(log_bm = log(bodymass),
         log_rs = log(runningspeed),
         g = 1  ## dummy for grouping variable
         ) |>
  tibble::rownames_to_column("species") |>
  mutate(across(species, fixup_species)) |>
  mutate(across(species, \(x) factor(x, levels = rownames(vcmat))))
stopifnot(identical(sort(rownames(vcmat)), sort(as.character(chdat$species))))

## -- null: no log_bm effect at all, just a phylo-correlated intercept -------
fit_null <- glmmTMB(log_rs ~ 1 + propto(0 + species | g, vcmat), data = chdat)

## -- spline basis setup (phyloslopes.qmd's spline_setup chunk) -------------
## absorb.cons=TRUE matches glmmTMB's own s() handling -- essential for an
## exact match to fit_spline_glmmTMB (phyloslopes.qmd), not just a
## convenience
sm <- smoothCon(s(log_bm, k = 10), data = chdat, absorb.cons = TRUE)
sm2ran <- smooth2random(sm[[1]], "", type = 2)
Xf <- sm2ran$Xf              ## fixed (linear-in-log_bm) design, 49 x 1
Xfull <- cbind(1, Xf)        ## + intercept column, so beta = c(intercept, linear coef)
Xr <- sm2ran$rand$Xr         ## wiggly random-effect design, 49 x 8, iid-equivalent
Kw <- ncol(Xr)

## cross the wiggly basis with the phylo indicator for the separable model
## below (species-outer/wiggly-inner, matching tensor.prod.model.matrix()
## convention -- Zphylo needs transposing first, same fac2sparse gotcha as
## always)
Zphylo <- fac2sparse(chdat$species)
Xr_joint <- tensor.prod.model.matrix(list(as(t(Zphylo), "dgCMatrix"), as(Xr, "dgCMatrix")))

## -- additive ----------------------------------------------------------
p0_spadd <- list(beta = rep(0, 2), b_spline = rep(0, Kw), b_phylo = rep(0, nrow(vcmat)),
                 logsd = 0, logsd_f = 0, logpsd = 0)
chdat_x <- c(chdat, lst(Xfull, Xr, Zphylo, vcmat))
obj_spadd <- MakeADFun(nllfun_spline_additive, p0_spadd, silent = TRUE,
                       random = c("b_spline", "b_phylo"))
fit_additive <- TMBfit(obj_spadd)

## -- separable -----------------------------------------------------------
## adds a phylogenetically-correlated copy of the wiggly coefficients on
## top of additive's baseline (non-phylogenetic) wiggle -- a single extra
## scale (logpsd_f), diagonal across the Kw wiggly dimensions, ADDED TO
## (not substituted for) b_spline, so the model properly nests additive
## (recoverable exactly as logpsd_f -> -Inf)
p0_spsep <- modifyList(p0_spadd, list(b_wiggly = rep(0, nrow(vcmat)*Kw), logpsd_f = -10))
chdat_x <- c(chdat, lst(Xfull, Xr, Xr_joint, Zphylo, Kw, vcmat))
obj_spsep <- MakeADFun(nllfun_spline_separable, p0_spsep, silent = TRUE,
                       random = c("b_spline", "b_wiggly", "b_phylo"))
fit_separable <- TMBfit(obj_spsep)

## -- tensor-product smooth -------------------------------------------------
## needs a *fresh* (non-absorb.cons) smoothCon() call: sm above already has
## the constant direction folded away (a 9-dim basis with a 1-dim null
## space), which would leave only one combined null-space scale instead of
## two separately-shrinkable ones
X <- model.matrix(~ log_bm, data = chdat)
sm_full <- smoothCon(s(log_bm, k = 10), data = chdat)  ## no absorb.cons -- want the full null space
S_spline <- sm_full[[1]]$S[[1]]
ev <- eigen(S_spline, symmetric = TRUE)
null_idx <- which(ev$values < 1e-8 * max(ev$values))
range_idx <- which(ev$values >= 1e-8 * max(ev$values))
Xf_null <- sm_full[[1]]$X %*% ev$vectors[, null_idx, drop = FALSE]                 ## 49 x 2: [linear, constant]
Xr_range <- sm_full[[1]]$X %*% ev$vectors[, range_idx, drop = FALSE] %*%
  diag(1/sqrt(ev$values[range_idx]), nrow = length(range_idx))                    ## 49 x 8, Wood's sqrt(D) reparam
Kn <- ncol(Xf_null); Kr <- ncol(Xr_range)

tZphylo <- as(t(Zphylo), "dgCMatrix")
Xnull_joint <- tensor.prod.model.matrix(list(tZphylo, as(Xf_null, "dgCMatrix")))   ## 49 x 98
Xrange_joint <- tensor.prod.model.matrix(list(tZphylo, as(Xr_range, "dgCMatrix"))) ## 49 x 392

## the two null-space directions (constant and linear-in-log_bm) need
## SEPARATE scales -- a single shared scale can't shrink away the
## (unsupported) phylo-slope variation without also killing the
## (well-supported) phylo-intercept variation
p0_sptensor <- list(beta = rep(0, 2), b_null = rep(0, nrow(vcmat)*Kn), b_range = rep(0, nrow(vcmat)*Kr),
                    logsd = 0, logpsd_null = rep(0, Kn), logpsd_range = 0)
chdat_x <- c(chdat, lst(X, Xnull_joint, Xrange_joint, Kr, vcmat))
obj_sptensor <- MakeADFun(nllfun_spline_tensor, p0_sptensor, silent = TRUE,
                          random = c("b_null", "b_range"))
fit_tensor <- TMBfit(obj_sptensor)

## -- store as a named list ---------------------------------------------
phyloslopes_combo_models <- list(
  null = fit_null,
  additive = fit_additive,
  separable = fit_separable,
  tensor = fit_tensor
)

## -- quick comparison summary --------------------------------------------
## fixef()/logLik()/AIC() dispatch uniformly across the glmmTMB fit and the
## RTMB TMBfit objects (fixef.TMB/logLik.TMB/AIC.TMB, phyloslopes_utils.R)
cat("=== phyloslopes_combo_models: fixed effects + logLik/AIC ===\n")
for (nm in names(phyloslopes_combo_models)) {
  fit <- phyloslopes_combo_models[[nm]]
  fx <- if (inherits(fit, "glmmTMB")) fixef(fit)$cond else fixef(fit)
  cat(sprintf("\n--- %s ---\n", nm))
  print(round(fx, 4))
  cat(sprintf("logLik = %.4f, AIC = %.4f\n", as.numeric(logLik(fit)), AIC(fit)))
}

## -- save the fitted models themselves (mirrors phyloslopes_linear.R) ------
## NB: both glmmTMB and RTMB TMBfit objects wrap a TMB ADFun with a C++
## pointer that goes stale on a save()/load() round-trip into a fresh R
## session -- basic accessors used above (fixef/logLik/AIC/summary,
## already-computed at fit time) keep working on a reloaded object, but
## anything that needs to re-evaluate the objective (obj$fn()/obj$gr(),
## sdreport(), update(), refitting) will not.
save(phyloslopes_combo_models, chtree, vcmat, chdat, file = "phyloslopes_combo.rda")
