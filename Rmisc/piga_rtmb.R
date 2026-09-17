library(RTMB)

## RTMB implementation of the Poisson-inverse gamma (PIGA) log-likelihood,
## see dpiga.R for the distribution, references, and the base-R version.
##
## RTMB overloads R's standard math primitives (lgamma, log, sqrt, besselK,
## ...) for its AD "advector" type, so the objective below is written as the
## same formula as dpiga(), just against `data$k`/`mu`/`phi` computed from
## the free parameters. RTMB's besselK() *is* differentiable with respect to
## both its argument and its order (confirmed against numDeriv::grad here),
## which is notable because the source paper's own EM algorithm avoids
## needing that derivative in closed form, using numerical differentiation
## (numDeriv::grad on ghyp::Egig) instead, precisely because it isn't
## generally available/reliable. RTMB's automatic differentiation sidesteps
## the entire EM construction and lets a joint (mu, phi) Newton-type
## optimizer be used directly on the exact PIGA log-likelihood, at least for
## the "distribution" (no-covariates) case set up here.
##
## However, exactly as documented in dpiga.R, RTMB's besselK() inherits the
## same crash/hang hazard as base R's for extreme (order, argument) pairs
## (confirmed: it segfaults on order = 1e10), so the fitted parameters must
## again be kept off a numerically dangerous region -- here by giving
## nlminb() box constraints on the working (log-scale) parameters.

piga_nll <- function(parameters) {
  getAll(parameters)
  mu  <- exp(log_mu)
  phi <- exp(log_phi)

  nu  <- data$k - phi - 1
  arg <- 2 * sqrt(mu * phi)

  logK <- log(besselK(arg, nu, expon.scaled = TRUE)) - arg

  loglik <- log(2) - lgamma(data$k + 1) +
    0.5 * (data$k + phi + 1) * log(mu * phi) -
    lgamma(phi + 1) +
    logK

  ADREPORT(mu)
  ADREPORT(phi)
  -sum(loglik)
}

if (FALSE) {

  source("dpiga.R")

  set.seed(1)
  mu_true <- 3; phi_true <- 2.5
  lambda <- 1 / rgamma(20000, shape = phi_true + 1, rate = phi_true)
  k <- rpois(20000, lambda * mu_true)

  data <- list(k = k)
  parameters <- list(log_mu = 0, log_phi = 0)

  ## piga_nll() refers to `data` by lexical scope, as is idiomatic for RTMB
  f <- function(parameters) piga_nll(parameters)
  environment(f) <- environment()

  obj <- MakeADFun(f, parameters, silent = TRUE)
  fit <- nlminb(obj$par, obj$fn, obj$gr, lower = c(-6, -6), upper = c(6, 6))
  stopifnot(fit$convergence == 0)

  sdr <- sdreport(obj)
  print(summary(sdr, select = "all"))

  ## cross-check against the plain-R optim + dpiga() fit
  nll_r <- function(par) -sum(dpiga(k, exp(par[1]), exp(par[2]), log = TRUE))
  fit_r <- optim(c(0, 0), nll_r, method = "L-BFGS-B",
                 lower = c(-6, -6), upper = c(6, 6))
  stopifnot(all.equal(fit$par, fit_r$par, tolerance = 1e-6))
  stopifnot(all.equal(obj$fn(fit$par), fit_r$value, tolerance = 1e-6))
}
