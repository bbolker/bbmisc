
## thermal performance curves
## basic version from Rezende
tpcfun <- function(T, C, Q10, Tth, d) {
  b <- pmin(1, 1-d*(T-Tth)^2)
  b*C*exp(T/10*log(Q10))
}

## differentiable version
tpcfun_sp <- function(temp, C, Q10, Tth, d, eps = 1e-3) {
  b <- 1-squareplus(d*(temp-Tth)^2, eps)
  b*C*exp(temp/10*log(Q10))
}

## further refinement: clamp at 0!
## FIXME: allow passing log_Q10 directly? (at present we
##  exponentiate to get from fitting parameter → model parameter,
##  then take the log again to put it into the formula)
tpcfun_sp_c <- function(temp, C, Q10, Tth, d, eps = 1e-3) {
  b <- squareplus(1-squareplus(d*(temp-Tth)^2, eps))
  b*C*exp(temp/10*log(Q10))
}

squareplus <- function(x, eps = 1e-3) {
  (x + sqrt(x^2 + eps^2))/2
}

prior_ranges <- list(log_Q10 = c(0.2, 2),
                     Tth = c(15, 25),
                     log_d = c(-7.5, -3.5)
                     )

## tried to protect insulate from 'prior_ranges' def above ... couldn't figure it out,
##  so using explicit `use_priors` in tmbdata ...
tpc_nll <-
  local(
    function(pars) {
      getAll(pars, tmbdata)
      mu <- tpcfun_sp_c(Tvec,
                        C = exp(log_C),
                  Q10 = exp(log_Q10),
                  Tth = Tth,
                  d = exp(log_d))
  prob <- mu/(1+mu)
  REPORT(prob)
  ADREPORT(prob)
  if (plotdebug) {
    lines(Tvec, prob, col = (debug_ctr %% 8) + 1)
    debug_ctr <<- debug_ctr + 1
    if (!is.na(debug_sleep)) Sys.sleep(debug_sleep)
  }
  disp <- exp(log_disp)
  shape1 <- prob*disp
  shape2 <- (1-prob)*disp
  ## -sum(lgamma(shape1) + lgamma(shape2) - lgamma(shape1 + shape2) +
  ## (shape1-1)*log(y)+(shape2-1)*log(1-y))
  res <- -sum(dbeta(y, shape1 = prob*disp, shape2 = (1-prob)*disp, log = TRUE))
  if (use_priors) {
    ## does iterating in this way actually work?
    for (p in names(prior_ranges)) {
      res <- res - dnorm(pars[[p]],
                         mean = mean(prior_ranges[[p]]),
                         sd = diff(prior_ranges[[p]])/4,
                         log = TRUE)
    }
  }
      return(res)
    })


mk_tpc_nll <- local(
  function(pars, tmbdata, silent = TRUE, ...) {
    MakeADFun(tpc_nll, pars, silent = silent, ...)
  })


uni_prior_samples <- function(nm, vals, breadth = 4, n = 1000) {
  sd <- diff(vals)/breadth
  m <- mean(vals)
  res <- rnorm(n, m, sd)
  if (grepl("^log_", nm)) res <- exp(res)
  res
}

draw_prior_samples <- function(prior_ranges, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  Map(\(nm, vals) uni_prior_samples(nm, vals, ...),
      names(prior_ranges), prior_ranges) |>
    do.call(what = cbind)
}


tpc_prob <- function(pars = true_pars, fn = tpcfun_sp_c, Tvec = 15:30,
                     invlink = function(mu) mu/(1+mu)) {
  ## FIXME:: do.call() instead?
  for (nm in names(pars)) {
    if (grepl("^log_", nm)) {
      newnm <- gsub("^log_", "", nm)
      pars[[newnm]] <- exp(pars[[nm]])
      nm <- newnm
    }
    assign(nm, pars[[nm]])
  }
  mu <- fn(Tvec, C = C, Q10 = Q10, Tth = Tth, d = d)
  prob <- invlink(mu)
  return(prob)
}
