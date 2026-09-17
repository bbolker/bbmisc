## Poisson-inverse gamma (PIGA) distribution
##
## Tzougas, G. (2020) "EM Estimation for the Poisson-Inverse Gamma Regression
## Model with Varying Dispersion: An Application to Insurance Ratemaking",
## Risks 8(3):97, Eqs. (1)-(3). https://doi.org/10.3390/risks8030097
##
## Mixture representation:
##   k | lambda ~ Poisson(lambda * mu)
##   lambda     ~ InverseGamma(shape = phi + 1, scale = phi)   [E(lambda) = 1]
##
## Unconditional pmf (paper's Eq. 3):
##   P(k) = (2/k!) * (mu*phi)^{(k+phi+1)/2} / Gamma(phi+1) * K_{k-phi-1}(2*sqrt(mu*phi))
## where K_nu(.) is the modified Bessel function of the second kind, order nu
## (`?besselK`).
##
## Moments (paper's Eqs. 7-8): E(k) = mu, Var(k) = mu + mu^2/(phi-1), phi > 1.
## (The pmf itself is defined, and typically fit, for any phi > 0.)

##' Density of the Poisson-inverse gamma (PIGA) distribution
##'
##' @param x vector of (non-negative integer) counts
##' @param mu mean parameter, mu > 0
##' @param phi dispersion parameter, phi > 0 (variance is finite only for phi > 1)
##' @param log if TRUE, return the log-density
##' @return vector of (log-)densities, recycled to the length of the longest
##'   of \code{x}, \code{mu}, \code{phi}
dpiga <- function(x, mu, phi, log = FALSE) {
  n   <- max(length(x), length(mu), length(phi))
  x   <- rep_len(x,   n)
  mu  <- rep_len(mu,  n)
  phi <- rep_len(phi, n)

  nu  <- x - phi - 1
  arg <- 2 * sqrt(mu * phi)

  logdens <- rep(NaN, n)

  ## valid parameter region
  valid <- mu > 0 & phi > 0 & !is.na(x)
  ## support: non-negative integers
  intx  <- valid & x == round(x) & x >= 0
  logdens[valid & !intx] <- -Inf

  ## Guard against calling besselK() with extreme arguments: base R's
  ## besselK() (and RTMB's AD version, see piga_rtmb.R) can *segfault the
  ## whole R session* on a non-finite order or argument, and its cost is
  ## O(|nu|) with no internal sanity bound -- even a large but finite order
  ## (~1e9) can hang for tens of seconds or attempt a multi-gigabyte
  ## allocation. This is a real hazard during unconstrained ML optimization,
  ## where an optimizer can transiently propose astronomically large mu/phi.
  ## Realistic dispersion parameters and claim counts never approach this
  ## range, so capping |nu| and arg at 1e4 costs nothing in practice.
  BESSEL_CAP <- 1e4
  finite_ok <- intx & is.finite(arg) & is.finite(nu) &
      abs(nu) < BESSEL_CAP & arg < BESSEL_CAP
  if (any(intx & !finite_ok)) {
    warning("mu/phi out of a safe numerical range for besselK(); ",
            "returning NaN instead of risking a hang or crash")
  }
  logdens[intx & !finite_ok] <- NaN

  if (any(finite_ok)) {
    a <- arg[finite_ok]; v <- nu[finite_ok]
    xk <- x[finite_ok]; m <- mu[finite_ok]; p <- phi[finite_ok]

    ## log(K_nu(a)), computed via the exponentially-scaled Bessel K for
    ## numerical stability: besselK(a, v, expon.scaled = TRUE) = K_v(a)*exp(a)
    logK <- log(besselK(a, v, expon.scaled = TRUE)) - a

    logdens[finite_ok] <- log(2) - lgamma(xk + 1) +
      0.5 * (xk + p + 1) * log(m * p) -
      lgamma(p + 1) +
      logK
  }

  if (any(x != round(x) | x < 0, na.rm = TRUE)) {
    warning("non-integer or negative x: density is zero")
  }

  if (log) logdens else exp(logdens)
}

if (FALSE) {

  ## --- sanity check against brute-force numerical integration -----------
  brute <- function(k, mu, phi) {
    integrand <- function(lam) {
      dpois(k, lambda = lam * mu) *
        (phi^(phi + 1) / gamma(phi + 1)) * lam^(-phi - 2) * exp(-phi / lam)
    }
    integrate(integrand, lower = 0, upper = Inf, rel.tol = 1e-10)$value
  }
  pars <- expand.grid(k = 0:6, mu = c(0.5, 2, 5), phi = c(1.5, 3, 8))
  pars$analytic <- with(pars, dpiga(k, mu, phi))
  pars$brute    <- mapply(brute, pars$k, pars$mu, pars$phi)
  stopifnot(max(abs(pars$analytic - pars$brute)) < 1e-8)

  ## sums to (approximately) 1 over the support
  stopifnot(abs(sum(dpiga(0:150, 2, 1.5)) - 1) < 1e-4)

  ## --- ML fit on simulated data ------------------------------------------
  set.seed(1)
  mu_true <- 3; phi_true <- 2.5
  lambda <- 1 / rgamma(20000, shape = phi_true + 1, rate = phi_true)
  k <- rpois(20000, lambda * mu_true)

  nll <- function(par) {
    mu <- exp(par[1]); phi <- exp(par[2])
    -sum(dpiga(k, mu, phi, log = TRUE))
  }
  ## bound the search region well within the besselK() safety cap above
  ## (see the comment on BESSEL_CAP): this is standard practice for any
  ## likelihood with a numerically fragile special function, independent
  ## of the specific bug in besselK() noted here.
  fit <- optim(c(0, 0), nll, method = "L-BFGS-B",
               lower = c(-6, -6), upper = c(6, 6), hessian = TRUE)
  mu_hat  <- exp(fit$par[1])
  phi_hat <- exp(fit$par[2])
  se <- sqrt(diag(solve(fit$hessian))) * c(mu_hat, phi_hat)  # delta method
  cat("mu_hat =", mu_hat, "(se", se[1], "), phi_hat =", phi_hat,
      "(se", se[2], ")\n")
}
