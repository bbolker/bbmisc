#' Fast Monte Carlo single-step max-|T| critical value
#'
#' Computes the critical value for the "single-step" simultaneous-inference
#' correction (the multiple-comparisons analogue of Tukey's HSD, generalised
#' from pairwise mean contrasts to arbitrary correlated regression
#' coefficients) by simulating directly from the multivariate-t null
#' distribution of the test statistics, rather than by numerically inverting
#' its CDF.
#'
#' @details
#' # Background
#' For a linear model with k estimated coefficients, the k marginal
#' t-statistics are *individually* exact under their null, but they are not
#' independent: they share one residual-variance estimate (so all k are
#' driven by the same random denominator) and, unless the design is
#' orthogonal, their numerators are correlated too. The joint null
#' distribution of the standardized coefficients is exactly multivariate-t
#' with correlation matrix `R = cov2cor(vcov(fit))` and `df =
#' fit$df.residual`. The single-step correction rejects if `max(abs(T)) >
#' c`, where `c` is the `(1 - alpha)` quantile of that multivariate-t's
#' max-of-|T| distribution -- this exactly controls the family-wise error
#' rate (not just bounds it, unlike Bonferroni/Holm), because it uses the
#' actual joint null rather than a dependence-agnostic bound.
#'
#' # Why not `mvtnorm::qmvt()` directly
#' `mvtnorm::qmvt()` computes exactly this quantity (and is what
#' `multcomp::glht(test = adjusted("single-step"))` uses internally), but
#' its root-finding search over repeated `pmvt()` numerical-integration
#' calls costs several seconds per call at k ~ 20 dimensions regardless of
#' the requested precision -- the cost is the number of root-finding
#' iterations, not the integration budget -- which is impractical at
#' simulation scale (thousands of replicates).
#'
#' # This function's approach
#' Since we don't need `qmvt()`'s generality, we replace it with direct
#' simulation: draw `nsim` vectors from the multivariate-t(`R`, `df`) null
#' via [mvtnorm::rmvt()], take `max(abs(.))` of each draw, and use the
#' empirical `level` quantile of those maxima as the critical value -- the
#' standard "simulate, then take the order statistic" way of estimating a
#' quantile of a known distribution. This matches `qmvt()`'s answer closely
#' at a small fraction of the cost, making a per-replicate-exact correction
#' feasible across a full simulation grid.
#'
#' # Provenance
#' This specific shortcut (parametric simulation from the fitted null,
#' rather than numerical integration or data resampling) does not appear
#' to be packaged as a multiple-*comparisons* correction in any CRAN or
#' Bioconductor package as of this writing. The mathematically identical
#' technique -- simulate from the fitted (multivariate normal/t)
#' distribution of the coefficients, take the maximum absolute
#' standardized deviate per draw, use the empirical quantile as a
#' simultaneous critical value -- is, however, well established for
#' *simultaneous confidence bands* on GAM smooths (Marra & Wood 2012),
#' and is implemented in `gratia::confint.gam(..., type = "simultaneous")`
#' and discussed publicly in Simpson (2016). The two standard packaged
#' alternatives for multiple-comparisons correction proper are (a) the
#' exact numerical-integration route in `mvtnorm`/`multcomp` (slow here,
#' as above), and (b) nonparametric permutation-based maxT procedures
#' (Westfall & Young 1993), implemented in Bioconductor's
#' `multtest::mt.maxT` and CRAN's `permuco::compute_maxT` -- these
#' resample the *data* to build an empirical null rather than simulating
#' from a fitted parametric one, so they are a related but distinct
#' mechanism from the one used here.
#'
#' @param R a correlation matrix (k x k) for the k test statistics, e.g.
#'   `cov2cor(vcov(fit)[-1, -1, drop = FALSE])` with the intercept dropped.
#' @param df residual degrees of freedom shared by all k test statistics.
#' @param nsim number of Monte Carlo draws used to estimate the quantile.
#'   5000 gives a critical value stable to about +/-0.01 in these examples;
#'   increase for a smoother/more reproducible threshold, decrease for speed.
#' @param level confidence level (1 - family-wise alpha); default 0.95.
#'
#' @return A single numeric critical value `c` such that, under the joint
#'   null, `P(max(abs(T)) > c) approx 1 - level`.
#'
#' @references
#' Genz, A. and Bretz, F. (2009). *Computation of Multivariate Normal and
#' t Probabilities*. Lecture Notes in Statistics 195, Springer.
#'
#' Hothorn, T., Bretz, F. and Westfall, P. (2008). Simultaneous inference
#' in general parametric models. *Biometrical Journal*, 50(3), 346-363.
#' (the `multcomp` package and its single-step method)
#'
#' Westfall, P.H. and Young, S.S. (1993). *Resampling-Based Multiple
#' Testing: Examples and Methods for p-Value Adjustment*. Wiley.
#' (permutation-based maxT; see also Bioconductor's `multtest::mt.maxT`
#' and CRAN's `permuco::compute_maxT`)
#'
#' Marra, G. and Wood, S.N. (2012). Coverage properties of confidence
#' intervals for generalized additive model components. *Scandinavian
#' Journal of Statistics*, 39(1), 53-74. (simulation-based simultaneous
#' bands via the same simulate/max/quantile mechanic, for smooths rather
#' than hypothesis tests)
#'
#' Simpson, G. (2016). Simultaneous intervals for smooths, revisited.
#' \url{https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/}
#'
#' @examples
#' \dontrun{
#' fit <- lm(y ~ (X1 + X2 + X3)^2, data = dat)
#' R <- cov2cor(vcov(fit)[-1, -1, drop = FALSE])
#' crit <- maxT_crit(R, df = fit$df.residual)
#' tt <- broom::tidy(fit)
#' any(abs(tt$statistic[tt$term != "(Intercept)"]) > crit)
#' }
maxT_crit <- function(R, df, nsim = 5000, level = 0.95) {
  k <- nrow(R)
  if (is.null(k) || k <= 1) {
    ## single test: no correction needed, single-step reduces to the
    ## ordinary two-sided t critical value
    return(stats::qt(1 - (1 - level) / 2, df))
  }
  Tsim <- mvtnorm::rmvt(nsim, sigma = R, df = df)
  as.numeric(stats::quantile(apply(abs(Tsim), 1, max), level))
}

#' Single-step-corrected significance test for a fitted linear model
#'
#' Convenience wrapper around [maxT_crit()]: fits the single-step
#' correction to a specific `lm` object (using its actual coefficient
#' correlation matrix and residual df) and reports whether any non-intercept
#' predictor is significant after correction. Mirrors the uncorrected
#' `any_significant()` helper in `simulate_fig1.R`, so the two can be
#' compared directly on the same fitted model.
#'
#' @param fit a fitted `lm` object.
#' @param alpha family-wise significance level (default 0.05).
#' @param nsim number of Monte Carlo draws passed to [maxT_crit()].
#'
#' @return `TRUE`/`FALSE`: whether any predictor's t-statistic exceeds the
#'   single-step critical value.
maxT_significant <- function(fit, alpha = 0.05, nsim = 5000) {
  V <- stats::vcov(fit)[-1, -1, drop = FALSE]  ## drop intercept
  R <- stats::cov2cor(V)
  crit <- maxT_crit(R, df = fit$df.residual, nsim = nsim, level = 1 - alpha)
  tt <- broom::tidy(fit)
  tvals <- tt$statistic[tt$term != "(Intercept)"]
  any(abs(tvals) > crit)
}
