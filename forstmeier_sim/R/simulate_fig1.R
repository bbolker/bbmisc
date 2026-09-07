## Simulation functions for reproducing Fig. 1 of Forstmeier & Schielzeth (2011)
## Behav Ecol Sociobiol 65:47-55, doi:10.1007/s00265-010-1038-5

## Sample sizes used throughout: 200/50 are Fig. 1's two cases; 30 is the
## extra low-N/k case from the paper's Fig. 3b.
sim_N_values <- c(30, 50, 200)

## Balanced, independently-randomised two-level factors, one per column.
## Balance is exact per factor; approximate independence across factors
## follows from randomising each factor's assignment separately. Actual
## factors (rather than a numeric -0.5/0.5 recoding) under the global
## sum-to-zero contrasts set in one_rep_corrections() so that, in the
## presence of interactions, each main-effect test is of the effect at the
## population mean of the other factors, not at a reference level.
make_predictors <- function(N, m) {
  half <- N %/% 2
  rest <- N - 2 * half
  lvls <- c(rep("lo", half), rep("hi", half + rest))
  X <- replicate(m, factor(sample(lvls), levels = c("lo", "hi")), simplify = FALSE)
  names(X) <- paste0("X", seq_len(m))
  as.data.frame(X)
}

## RHS-only predictor spec, shared by build_formula() (two-sided) and the
## scope formulas used by step_interactions_only() (one-sided).
predictor_terms <- function(m, interactions) {
  preds <- paste0("X", seq_len(m))
  if (interactions) sprintf("(%s)^2", paste(preds, collapse = " + ")) else paste(preds, collapse = " + ")
}

build_formula <- function(m, interactions) {
  stats::as.formula(paste("y ~", predictor_terms(m, interactions)))
}

## number of predictors (including two-way interactions) for m factors
n_predictors <- function(m, interactions) {
  if (interactions) m + choose(m, 2) else m
}

any_significant <- function(fit, alpha = 0.05) {
  tt <- broom::tidy(fit)
  tt <- tt[tt$term != "(Intercept)", , drop = FALSE]
  if (nrow(tt) == 0) return(FALSE)
  any(tt$p.value < alpha, na.rm = TRUE)
}

## -- model-simplification strategies -----------------------------------
##
## Each takes (fit_full, m, interactions) and returns the simplified fit.
## Swappable via one_rep_corrections()'s simplify_fns argument (see
## simulate_corrections.R), so the simulation machinery doesn't need to
## know which variant is in use.

## Unrestricted backward AIC elimination (Forstmeier & Schielzeth's original
## design): step() on a model built from a two-way-interaction terms object
## respects marginality on its own (won't drop a main effect while its
## interaction remains), so no explicit scope is needed.
step_default <- function(fit_full, m, interactions) {
  step(fit_full, trace = 0)
}

## Restricted backward AIC elimination that may only discard interactions,
## never main effects: scope$lower pins the main-effects-only model as the
## floor, so step()'s drop candidates are always a subset of the interaction
## terms. With no interactions in the model to begin with, this is a no-op.
step_interactions_only <- function(fit_full, m, interactions) {
  if (!interactions) return(fit_full)
  lower <- stats::as.formula(paste("~", predictor_terms(m, FALSE)))
  upper <- stats::as.formula(paste("~", predictor_terms(m, TRUE)))
  step(fit_full, scope = list(lower = lower, upper = upper), direction = "backward", trace = 0)
}
