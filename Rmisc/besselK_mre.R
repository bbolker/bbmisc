## Minimal reproducible examples: besselK() segfaults / hangs / OOM-aborts
## for large orders `nu`, with no sanity bound on its magnitude.
##
## Found while implementing the Poisson-inverse gamma (PIGA) log-likelihood
## (see dpiga.R, piga_rtmb.R), whose density needs besselK() with an order
## that depends on data and a fitted dispersion parameter -- exactly the
## setup where an ordinary MLE fit can wander into this.
##
## WARNING: several calls below deliberately crash or hang the R process.
## Run this script with `Rscript besselK_mre.R` (one crash just ends that
## process) rather than sourcing it inside a session you care about, and
## run one example block at a time.
##
## Environment this was found on:
##   R version 4.6.1 (2026-06-24), x86_64-pc-linux-gnu, Ubuntu 24.04.5 LTS
##
## Suggested filing: https://bugs.r-project.org/ (component likely "Basic
## math functions" / nmath), or raise on r-devel first per usual etiquette.
## The same underlying routine is reachable (and crashes the same way) via
## RTMB::besselK() inside RTMB::MakeADFun() objectives -- see the note at
## the bottom -- which may warrant a separate report against
## https://github.com/kaskr/RTMB.

## --- 1. Minimal, direct crash -----------------------------------------
## besselK(x, nu, expon.scaled = TRUE) segfaults for some finite orders
## with no sanity bound. Not recoverable with tryCatch(); the R session
## just aborts.

besselK(1000, 3e9, expon.scaled = TRUE)
## *** caught segfault ***
## address 0x734d5d1d0040, cause 'memory not mapped'

## Non-finite order also segfaults:
# besselK(Inf, -Inf, expon.scaled = TRUE)

## --- 2. Magnitude sweep (informational; not all lines are safe to run) --
## besselK(1000, nu, expon.scaled = TRUE) for increasing |nu|:
##
##   nu       result
##   -------  --------------------------------------------------
##   1e3-1e6  returns instantly
##   1e7      returns Inf, ~0.27s
##   1e8      returns Inf, ~2.3s
##   1e9      returns Inf, ~23s
##   3e9      *** segfault ***
##   5e9      *** segfault ***
##   8e9      *** segfault ***
##   1e10     Error: cannot allocate vector of size 74.5 Gb  (recoverable)
##
## The non-monotonic transition -- segfault at 3e9/5e9/8e9, but a *graceful*
## allocation-failure error at 1e10 -- suggests this isn't simply "runs out
## of memory once nu is large enough": something behaves inconsistently in
## that range, consistent with (but not confirmed against the Fortran/C
## source) a 32-bit signed integer overflow near 2^31 (~2.1e9) in an
## internal work-array size/index computation.

## --- 3. Realistic trigger: an entirely ordinary MLE fit -----------------
## No extreme values are chosen by hand anywhere here; optim() reaches the
## crashing region on its own during a standard two-parameter fit.

## log-density of the Poisson-inverse gamma distribution
## (Tzougas 2020, Risks 8(3):97, Eq. 3); see dpiga.R for the safe version.
dpiga_unsafe <- function(x, mu, phi, log = FALSE) {
  nu   <- x - phi - 1
  arg  <- 2 * sqrt(mu * phi)
  logK <- log(besselK(arg, nu, expon.scaled = TRUE)) - arg
  logdens <- log(2) - lgamma(x + 1) + 0.5 * (x + phi + 1) * log(mu * phi) -
    2 * lgamma(phi + 1) + logK
  if (log) logdens else exp(logdens)
}

set.seed(1)
mu_true <- 3; phi_true <- 2.5
lambda <- 1 / rgamma(20000, shape = phi_true + 1, rate = phi_true)
k <- rpois(20000, lambda * mu_true)

nll <- function(par) {
  mu <- exp(par[1]); phi <- exp(par[2])
  -sum(dpiga_unsafe(k, mu, phi, log = TRUE))
}

fit <- optim(c(0, 0), nll, method = "BFGS")   # sane start, no bounds

## *** caught segfault ***
## address 0x5fe60d6933a8, cause 'memory not mapped'
##
## Traceback:
##  1: besselK(arg, nu, expon.scaled = TRUE)
##  2: dpiga_unsafe(k, mu, phi, log = TRUE)
##  3: fn(par, ...)
##  4: (function (par) fn(par, ...))(c(2291.76964906897, 3433.88951605988))
##  5: optim(c(0, 0), nll, method = "BFGS")
## An irrecoverable exception occurred. R is aborting now ...
##
## BFGS's finite-difference line search overshoots to
## par = c(2291.8, 3433.9), i.e. mu = exp(2291.8), phi = exp(3433.9) (both
## effectively Inf in double precision), giving nu = -Inf, which crashes
## inside besselK().

## --- 4. Note on RTMB -----------------------------------------------------
## RTMB's AD-enabled besselK() shares the same hazard for large finite
## order:
##
##   library(RTMB)
##   f <- function(par) besselK(par[1], par[2])
##   F <- MakeADFun(f, c(2, 1.5), silent = TRUE)
##   F$fn(c(1e10, -1e10))
##   # *** caught segfault ***
##
## Unlike base R's optim(..., method = "BFGS") above, an RTMB fit using
## exact analytic gradients (obj$gr) did not overshoot into this region in
## the cases tried (sane start, and a deliberately poor start); only
## derivative-free Nelder-Mead from a poor start reproduced a severe (>120s,
## not observed to complete) slowdown consistent with the same O(|nu|) cost,
## short of an actual crash.
