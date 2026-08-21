## phyloslopes_tiny_predcovs.R -- for each of the three tiny-example
## structures (additive, separable, tensor-product; see phyloslopes_tiny.R)
## as normally estimated, plus additive/separable refit with their
## random-effects variance-component hyperparameters held fixed at their
## true simulated values instead of estimated (so they can't collapse to
## a degenerate boundary the way separable's logpsd_f does when
## estimated), predict every species' curve VALUE and SLOPE at a new
## covariate value (x0 = 3) and get the joint 10x10 (5 values + 5 slopes)
## covariance matrix of those predictions.
##
## Approach: predicted value/slope are *linear* combinations of the full
## parameter vector (beta, b...) for all three models, so no delta-method
## approximation beyond that linearity is needed --
##   1. get the (fixed+random) basis contribution at x0 by evaluating the
##      already-fit smooth object at x0 via mgcv::PredictMat() (same basis/
##      knots/reparameterization as fitting -- not re-derived), and its
##      derivative by central finite difference of that basis evaluation
##      (exact up to O(h^2) truncation error; the basis functions are smooth
##      in x, so this is not a piecewise-linear/derivative-doesn't-exist
##      situation)
##   2. assemble a 10 x npar matrix A mapping the full parameter vector to
##      the 10 predicted quantities
##   3. get the joint (fixed+random) parameter covariance via
##      sdreport(obj, getJointPrecision = TRUE), invert, and compute
##      A %*% Cov %*% t(A)
## Run with: Rscript phyloslopes_tiny_predcovs.R

library(ape)
library(Matrix)
library(mgcv)      ## PredictMat()
library(MRFtools)
library(MASS)
library(RTMB)
library(tibble)    ## lst()

source("phyloslopes_utils.R")

## -- load the tiny-example setup/data from phyloslopes_tiny.R ---------------
## (design matrices, simulated data, and the smoothCon() results themselves,
## so PredictMat() below uses exactly the basis/knots/reparameterization the
## models were fit with, not a re-derived approximation of it)
load("phyloslopes_tiny.rda")
ntip <- length(chtree$tip.label)
Kw <- ncol(Xr)
Xf <- Xfull[, -1, drop = FALSE]  ## Xfull = cbind(1, Xf); only Xfull was saved
## tensor range-block penalty components (pure functions of data, not
## parameters, so cheaply recomputed here rather than saved -- see
## nllfun_spline_tensor / phyloslopes_tiny.R)
Sphylo <- solve(vcmat)
Qr_phylo <- kronecker(Sphylo, diag(Kr))
Qr_smooth <- kronecker(diag(ntip), diag(d_range))

## -- refit the three structures -----------------------------------------
## (TMB/RTMB ADFun objects hold C++ pointers and don't survive a save()/
## load() round-trip, so they're always refit here from the loaded design
## matrices/data rather than loaded directly)
## each chdat_x snapshot is also saved under a model-specific name
## (chdat_x_add, chdat_x_sep, ...): RTMB's report() re-evaluates the
## objective's R closure fresh on every call rather than replaying a
## frozen tape, so it needs chdat_x (looked up by nllfun_spline_*() via
## ordinary lexical scoping from .GlobalEnv) to still match *that* model
## at report()-time -- not just at the MakeADFun() call that built it.
## The sanity-check loop below restores the right snapshot before each
## model's report() call rather than relying on chdat_x still happening
## to hold the last-assigned value.
chdat_x <- lst(log_rs = sim_dat$y, Xfull, Xr, Zphylo, vcmat)
chdat_x_add <- chdat_x
p0_add <- list(beta = rep(0, 2), b_spline = rep(0, Kw), b_phylo = rep(0, ntip),
               logsd = 0, logsd_f = 0, logpsd = 0)
obj_add <- MakeADFun(nllfun_spline_additive, p0_add, silent = TRUE,
                     random = c("b_spline", "b_phylo"))
fit_add <- TMBfit(obj_add)

chdat_x <- lst(log_rs = sim_dat$y, Xfull, Xr, Xr_joint, Zphylo, Kw, vcmat)
chdat_x_sep <- chdat_x
p0_sep <- modifyList(p0_add, list(b_wiggly = rep(0, ntip * Kw), logpsd_f = -10))
obj_sep <- MakeADFun(nllfun_spline_separable, p0_sep, silent = TRUE,
                     random = c("b_spline", "b_wiggly", "b_phylo"))
fit_sep <- TMBfit(obj_sep)

## cor_null (fixed at 0 via map=) is nllfun_spline_tensor's optional
## null-space-direction correlation -- see its header comment
us2 <- unstructured(2)
chdat_x <- lst(log_rs = sim_dat$y, X = Xfull, Xnull_joint, Xrange_joint, Qr_phylo, Qr_smooth, vcmat, us2)
chdat_x_tensor <- chdat_x
p0_tensor <- list(beta = rep(0, 2), b_null = rep(0, ntip * Kn), b_range = rep(0, ntip * Kr),
                  logsd = 0, logpsd_null = rep(0, Kn), cor_null = 0,
                  logsigma1_range = 0, logsigma2_range = 0)
obj_tensor <- MakeADFun(nllfun_spline_tensor, p0_tensor, silent = TRUE,
                        random = c("b_null", "b_range"),
                        map = list(cor_null = factor(NA)))
fit_tensor <- TMBfit(obj_tensor)

## -- additive/separable refit with the random-effects variance-component
## hyperparameters FIXED at their true simulated values (sd_f, sd_wiggly,
## sd_phylo -- see phyloslopes_tiny.R), via TMB's map = (rather than let
## them be estimated and potentially collapse to a boundary, as logpsd_f
## does for the MLE separable fit above). beta and logsd (residual --
## not a random-effects covariance parameter) are still freely estimated,
## and the latent b_* are still profiled out by the inner Laplace
## approximation as usual; only the variance-component map = factor(NA)
## entries are held fixed. tensor has no clean true-value mapping (its
## null/range basis doesn't correspond to any of the DGP's three latent
## components) so it's not included here -- see phyloslopes.qmd/session
## discussion.
## same chdat_x content additive/separable already used (only p0/map
## differ), so just restore those snapshots rather than rebuilding
p0_add_true <- modifyList(p0_add, list(logsd_f = log(sd_f), logpsd = log(sd_phylo)))
map_add_true <- list(logsd_f = factor(NA), logpsd = factor(NA))
chdat_x <- chdat_x_add
obj_add_true <- MakeADFun(nllfun_spline_additive, p0_add_true, silent = TRUE,
                          random = c("b_spline", "b_phylo"), map = map_add_true)
fit_add_true <- TMBfit(obj_add_true)

p0_sep_true <- modifyList(p0_sep, list(logsd_f = log(sd_f), logpsd_f = log(sd_wiggly),
                                        logpsd = log(sd_phylo)))
map_sep_true <- list(logsd_f = factor(NA), logpsd_f = factor(NA), logpsd = factor(NA))
chdat_x <- chdat_x_sep
obj_sep_true <- MakeADFun(nllfun_spline_separable, p0_sep_true, silent = TRUE,
                          random = c("b_spline", "b_wiggly", "b_phylo"), map = map_sep_true)
fit_sep_true <- TMBfit(obj_sep_true)

## -- prediction-basis machinery ----------------------------------------
x0 <- 3            ## covariate value to predict/differentiate at
h <- 1e-4           ## central-finite-difference step for the basis derivative

## evaluate a fitted smooth's basis (and its derivative) at a new x
basis_deriv <- function(smobj, x0, h) {
  Xp <- PredictMat(smobj, data.frame(x_rep = x0 + h))
  Xm <- PredictMat(smobj, data.frame(x_rep = x0 - h))
  X0 <- PredictMat(smobj, data.frame(x_rep = x0))
  list(val = X0, deriv = (Xp - Xm) / (2 * h))
}

## additive/separable: Xf/Xr (from smooth2random(), used when fitting) are a
## linear reparameterization of sm[[1]]$X (the absorb.cons-constrained
## basis); recover that reparameterization matrix by (exact) least squares
## against the already-known Xf/Xr, then apply it to new-x evaluations
Tmat <- solve(crossprod(sm[[1]]$X), crossprod(sm[[1]]$X, cbind(Xf, Xr)))
get_XfXr <- function(x0, h) {
  bd <- basis_deriv(sm[[1]], x0, h)
  val <- bd$val %*% Tmat; der <- bd$deriv %*% Tmat
  list(Xf = val[, 1, drop = FALSE], Xr = val[, -1, drop = FALSE],
       Xf_d = der[, 1, drop = FALSE], Xr_d = der[, -1, drop = FALSE])
}

## tensor: Xf_null/Xr_range (from smooth2random(), used when fitting) are a
## linear reparameterization of sm_full[[1]]$X (the raw, unconstrained
## basis) -- recover that reparameterization matrix by (exact) least
## squares against the already-known Xf_null/Xr_range, same trick as
## get_XfXr()/Tmat above, rather than redoing the eigendecomposition
Tmat_nr <- solve(crossprod(sm_full[[1]]$X), crossprod(sm_full[[1]]$X, cbind(Xf_null, Xr_range)))
get_null_range <- function(x0, h) {
  bd <- basis_deriv(sm_full[[1]], x0, h)
  val <- bd$val %*% Tmat_nr; der <- bd$deriv %*% Tmat_nr
  list(Xf_null = val[, seq_len(Kn), drop = FALSE],
       Xr_range = val[, Kn + seq_len(Kr), drop = FALSE],
       Xf_null_d = der[, seq_len(Kn), drop = FALSE],
       Xr_range_d = der[, Kn + seq_len(Kr), drop = FALSE])
}

## -- sanity check: predicted value at each species' OWN x must reproduce
## REPORT(mu) exactly (up to floating point) -- validates the basis-recovery
## machinery above before trusting its extrapolation to x0
get_parlist <- function(obj) obj$env$parList(par = obj$env$last.par.best)
pl_add <- get_parlist(obj_add); pl_sep <- get_parlist(obj_sep)
pl_tensor <- get_parlist(obj_tensor)
pl_add_true <- get_parlist(obj_add_true); pl_sep_true <- get_parlist(obj_sep_true)

## additive and separable share the same value formula (beta + wiggle
## shape + phylo intercept); separable adds a species-specific wiggly
## deviation on top when wig_idx is supplied. chdat_x_i is restored into
## the *global* chdat_x (<<-) right before obj$report(): nllfun_spline_*()
## looks chdat_x up by ordinary lexical scoping from .GlobalEnv, and
## report() re-evaluates that closure fresh rather than replaying a
## frozen tape, so it needs the matching snapshot live at call time.
check_additive_style <- function(obj, fit, pl, bi, i, chdat_x_i, wig_idx = NULL) {
  chdat_x <<- chdat_x_i
  mu_report <- unique(obj$report()$mu[species_rep == chtree$tip.label[i]])
  mu_pred <- c(1, bi$Xf) %*% fit$fit$par[1:2] + bi$Xr %*% pl$b_spline + pl$b_phylo[i]
  if (!is.null(wig_idx)) mu_pred <- mu_pred + bi$Xr %*% pl$b_wiggly[wig_idx]
  stopifnot(isTRUE(all.equal(c(mu_pred), mu_report, tolerance = 1e-6)))
}
check_tensor_style <- function(obj, fit, pl, bi, bi_nr, i, null_idx_i, range_idx_i, chdat_x_i) {
  chdat_x <<- chdat_x_i
  mu_report <- unique(obj$report()$mu[species_rep == chtree$tip.label[i]])
  mu_pred <- c(1, bi$Xf) %*% fit$fit$par[1:2] +
    bi_nr$Xf_null %*% pl$b_null[null_idx_i] + bi_nr$Xr_range %*% pl$b_range[range_idx_i]
  stopifnot(isTRUE(all.equal(c(mu_pred), mu_report, tolerance = 1e-6)))
}

for (i in seq_len(ntip)) {
  bi <- get_XfXr(x[i], h)
  wig_idx <- ((i - 1) * Kw + 1):(i * Kw)
  check_additive_style(obj_add, fit_add, pl_add, bi, i, chdat_x_add)
  check_additive_style(obj_add_true, fit_add_true, pl_add_true, bi, i, chdat_x_add)
  check_additive_style(obj_sep, fit_sep, pl_sep, bi, i, chdat_x_sep, wig_idx)
  check_additive_style(obj_sep_true, fit_sep_true, pl_sep_true, bi, i, chdat_x_sep, wig_idx)

  bi_nr <- get_null_range(x[i], h)
  null_idx_i <- ((i - 1) * Kn + 1):(i * Kn); range_idx_i <- ((i - 1) * Kr + 1):(i * Kr)
  check_tensor_style(obj_tensor, fit_tensor, pl_tensor, bi, bi_nr, i, null_idx_i, range_idx_i, chdat_x_tensor)
}
cat("Prediction-basis sanity checks passed for all 5 models (all", ntip, "species).\n")

## -- build the 10 x npar "A" matrix (5 values + 5 slopes, at x0) for one
## model's parameter vector, given named blocks of the value/derivative
## contribution -- a block is either a single row (shared across all
## species, e.g. beta/b_spline) or a list of ntip species-specific rows
## (e.g. b_phylo, b_wiggly, b_null, b_range)
make_A <- function(pnames, blocks_val, blocks_deriv) {
  npar <- length(pnames)
  A_val <- matrix(0, ntip, npar); A_slope <- matrix(0, ntip, npar)
  for (nm in names(blocks_val)) {
    idx <- which(pnames == nm)
    v <- blocks_val[[nm]]; d <- blocks_deriv[[nm]]
    if (is.list(v)) {
      chunk <- length(idx) / ntip
      for (i in seq_len(ntip)) {
        ii <- idx[((i - 1) * chunk + 1):(i * chunk)]
        A_val[i, ii] <- v[[i]]; A_slope[i, ii] <- d[[i]]
      }
    } else {
      for (i in seq_len(ntip)) { A_val[i, idx] <- v; A_slope[i, idx] <- d }
    }
  }
  rbind(A_val, A_slope)
}

b0_XfXr <- get_XfXr(x0, h)
b0_nullrange <- get_null_range(x0, h)
rowlabs <- c(paste0("value_", chtree$tip.label), paste0("slope_", chtree$tip.label))

A_add <- make_A(names(obj_add$env$last.par.best),
  list(beta = c(1, b0_XfXr$Xf), b_spline = b0_XfXr$Xr,
       b_phylo = lapply(seq_len(ntip), function(i) 1)),
  list(beta = c(0, b0_XfXr$Xf_d), b_spline = b0_XfXr$Xr_d,
       b_phylo = lapply(seq_len(ntip), function(i) 0)))

A_sep <- make_A(names(obj_sep$env$last.par.best),
  list(beta = c(1, b0_XfXr$Xf), b_spline = b0_XfXr$Xr,
       b_wiggly = lapply(seq_len(ntip), function(i) b0_XfXr$Xr),
       b_phylo = lapply(seq_len(ntip), function(i) 1)),
  list(beta = c(0, b0_XfXr$Xf_d), b_spline = b0_XfXr$Xr_d,
       b_wiggly = lapply(seq_len(ntip), function(i) b0_XfXr$Xr_d),
       b_phylo = lapply(seq_len(ntip), function(i) 0)))

A_tensor <- make_A(names(obj_tensor$env$last.par.best),
  list(beta = c(1, b0_XfXr$Xf),
       b_null = lapply(seq_len(ntip), function(i) b0_nullrange$Xf_null),
       b_range = lapply(seq_len(ntip), function(i) b0_nullrange$Xr_range)),
  list(beta = c(0, b0_XfXr$Xf_d),
       b_null = lapply(seq_len(ntip), function(i) b0_nullrange$Xf_null_d),
       b_range = lapply(seq_len(ntip), function(i) b0_nullrange$Xr_range_d)))

## additive/true, separable/true: identical block specification to
## A_add/A_sep (beta, b_spline, b_phylo[, b_wiggly] all still present and
## unmapped -- only the variance-component hyperparameters were fixed),
## just re-indexed against the *_true objects' own (shorter, since the
## mapped-out hyperparameters drop out of last.par.best entirely) parameter
## vectors
A_add_true <- make_A(names(obj_add_true$env$last.par.best),
  list(beta = c(1, b0_XfXr$Xf), b_spline = b0_XfXr$Xr,
       b_phylo = lapply(seq_len(ntip), function(i) 1)),
  list(beta = c(0, b0_XfXr$Xf_d), b_spline = b0_XfXr$Xr_d,
       b_phylo = lapply(seq_len(ntip), function(i) 0)))

A_sep_true <- make_A(names(obj_sep_true$env$last.par.best),
  list(beta = c(1, b0_XfXr$Xf), b_spline = b0_XfXr$Xr,
       b_wiggly = lapply(seq_len(ntip), function(i) b0_XfXr$Xr),
       b_phylo = lapply(seq_len(ntip), function(i) 1)),
  list(beta = c(0, b0_XfXr$Xf_d), b_spline = b0_XfXr$Xr_d,
       b_wiggly = lapply(seq_len(ntip), function(i) b0_XfXr$Xr_d),
       b_phylo = lapply(seq_len(ntip), function(i) 0)))

## -- point estimates + joint covariance, all 5 models ------------------
## additive/true and separable/true (random-effects variance-component
## hyperparameters fixed at their true simulated values -- see the
## obj_add_true/obj_sep_true construction above) come last, after the
## three normally-estimated fits
predcovs <- list()
## chdat_x restored (<<-) before sdreport() for the same reason as in the
## sanity-check loop above -- it re-evaluates obj$fn()/obj$gr() fresh, so
## needs the matching data snapshot live, not whatever chdat_x last
## happened to hold
model_defs <- list(additive = list(obj = obj_add, A = A_add, chdat_x = chdat_x_add),
                    separable = list(obj = obj_sep, A = A_sep, chdat_x = chdat_x_sep),
                    tensor = list(obj = obj_tensor, A = A_tensor, chdat_x = chdat_x_tensor),
                    `additive/true` = list(obj = obj_add_true, A = A_add_true, chdat_x = chdat_x_add),
                    `separable/true` = list(obj = obj_sep_true, A = A_sep_true, chdat_x = chdat_x_sep))
for (nm in names(model_defs)) {
  obj_i <- model_defs[[nm]]$obj
  A_i <- model_defs[[nm]]$A
  chdat_x <- model_defs[[nm]]$chdat_x
  sdr <- sdreport(obj_i, getJointPrecision = TRUE)
  Cov_joint <- as.matrix(solve(sdr$jointPrecision))
  point <- drop(A_i %*% obj_i$env$last.par.best)
  Cov_pred <- A_i %*% Cov_joint %*% t(A_i)
  names(point) <- rowlabs; dimnames(Cov_pred) <- list(rowlabs, rowlabs)

  cat(sprintf("\n=== %s: value/slope at x=%g, all %d species ===\n", nm, x0, ntip))
  print(round(point, 4))
  cat(sprintf("\n=== %s: SEs ===\n", nm))
  print(round(sqrt(diag(Cov_pred)), 4))

  predcovs[[nm]] <- list(point = point, cov = Cov_pred)
}

## -- save (plain-data results only -- point estimates and covariance
## matrices, never the obj/fit TMB objects themselves, which hold C++
## pointers that don't survive a save()/load() round-trip) -----------------
save(x0, predcovs, file = "phyloslopes_tiny_predcovs.rda")
