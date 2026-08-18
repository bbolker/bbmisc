## phyloslopes_linear.R -- fit and compare four *linear* (log_bm) models for
## how phylogeny enters a species-level random intercept and/or slope, all
## on the real (49-species, one obs/species) carniherbi49 data:
##   - null:        no log_bm effect at all, just a phylo-correlated intercept
##   - additive:    fixed log_bm slope + phylo-correlated intercept
##   - independent: phylo-correlated intercept + a SEPARATE, non-
##                  phylogenetic (iid) random slope for log_bm -- two
##                  independent random-effect terms, not a joint/correlated
##                  structure
##   - separable:   intercept AND slope jointly phylo-correlated, via a full
##                  phylo x 2x2-trait Kronecker covariance
##
## null/additive/independent are plain glmmTMB fits -- propto() natively
## handles a phylogenetically-correlated single random-effect term across
## species levels. separable is NOT expressible via propto(): tested
## directly, `propto(0 + log_bm + species | g, vcmat)` errors ("matrix is
## not the correct dimensions"), because propto() replaces one random
## effect's correlation structure across its grouping-factor *levels*, not
## a crossed species x trait block with a full correlated covariance. So
## separable reuses the existing RTMB nllfun_sep machinery
## (phyloslopes_utils.R) instead -- exactly the model fit as fit_sep in
## phyloslopes_tests.R.
##
## NB: independent does not converge cleanly on this dataset ("Model
## convergence problem; non-positive-definite Hessian matrix", logLik NA)
## -- with exactly one observation per species and no replication, the
## 49-dimensional phylo-correlated intercept can already interpolate the
## data on its own, so stacking an *additional* random-effect term (the
## iid slope) on top just gives the optimizer room to drive the residual
## variance to a zero boundary, destabilizing the Hessian (the same
## overparameterization/interpolation pathology documented at length in
## phyloslopes_tiny_explore.R). Kept in the list anyway, unconverged
## warnings and all -- they're recoverable from the stored glmmTMB object
## itself (e.g. via its $sdr, $fit$convergence, or diagnose()).
##
## Run with: Rscript phyloslopes_linear.R

suppressMessages({
  library(ade4)
  library(ape)
  library(dplyr)
  library(glmmTMB)
  library(RTMB)
  library(reformulas)
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

## -- additive: fixed log_bm slope + phylo-correlated intercept -------------
fit_additive <- glmmTMB(log_rs ~ log_bm + propto(0 + species | g, vcmat), data = chdat)

## -- independent: phylo-correlated intercept + a SEPARATE, non-
## phylogenetic (iid) random slope for log_bm -- see header re: convergence
fit_independent <- glmmTMB(log_rs ~ log_bm + (0 + log_bm | species) +
                             propto(0 + species | g, vcmat), data = chdat)

## -- separable: intercept AND slope jointly phylo-correlated, via a full
## phylo x 2x2-trait Kronecker covariance (dseparable(), nllfun_sep in
## phyloslopes_utils.R) -- not expressible through propto() (see header)
X <- model.matrix(~ log_bm, data = chdat)
us2 <- unstructured(2)
rt <- mkReTrms(findbars(~ (1 + log_bm | species)), fr = chdat, calc.lambdat = FALSE)
p2 <- list(beta = rep(0, 2), b = matrix(0, nrow = nrow(vcmat), ncol = 2),
          logrsd = rep(0, 2), logsdres = 0, corval = 0)
chdat_x <- c(chdat, lst(X, phylomat = vcmat, scale = 1, Z = t(rt$Zt), us2))
fit_separable <- TMBfit(MakeADFun(nllfun_sep, p2, silent = TRUE, random = "b"))

## -- store as a named list ---------------------------------------------
phyloslopes_linear_models <- list(
  null = fit_null,
  additive = fit_additive,
  independent = fit_independent,
  separable = fit_separable
)

## -- quick comparison summary --------------------------------------------
## fixef()/logLik()/AIC() dispatch uniformly across the glmmTMB fits and
## the RTMB TMBfit object (fixef.TMB/logLik.TMB/AIC.TMB, phyloslopes_utils.R)
cat("=== phyloslopes_linear_models: fixed effects + logLik/AIC ===\n")
for (nm in names(phyloslopes_linear_models)) {
  fit <- phyloslopes_linear_models[[nm]]
  fx <- if (inherits(fit, "glmmTMB")) fixef(fit)$cond else fixef(fit)
  cat(sprintf("\n--- %s ---\n", nm))
  print(round(fx, 4))
  cat(sprintf("logLik = %.4f, AIC = %.4f\n", as.numeric(logLik(fit)), AIC(fit)))
}

## -- save the fitted models themselves, as requested ------------------------
## NB: both glmmTMB and RTMB TMBfit objects wrap a TMB ADFun with a C++
## pointer that goes stale on a save()/load() round-trip into a fresh R
## session -- basic accessors used above (fixef/logLik/AIC/summary,
## already-computed at fit time) keep working on a reloaded object, but
## anything that needs to re-evaluate the objective (obj$fn()/obj$gr(),
## sdreport(), update(), refitting) will not.
save(phyloslopes_linear_models, chtree, vcmat, chdat, file = "phyloslopes_linear.rda")
