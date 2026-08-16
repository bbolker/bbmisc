## Utility functions for phyloslopes.qmd

## -- plotting / TMB helpers ------------------------------------------------

ifun <- function(x) image(Matrix(x), sub = "",
                          xlab = "", ylab = "")

TMBfit <- function(obj, optimizer = nlminb, response = NULL) {
  fit <- optimizer(obj$par, obj$fn, obj$gr)
  ret <- list(fit = fit, obj = obj)
  if (!is.null(response)) attr(ret, "response") <- response
  class(ret) <- c("TMB", "list")
  ret
}
AIC.TMB <- function(x) { 2*x$fit$objective + 2*npar(x) }
logLik.TMB <- function(x) { -1*x$fit$objective }
fixef.TMB <- function(x) { x$obj$env$parList()$beta }
npar <- function(x, ...) {
  UseMethod("npar")
}
npar.TMB <- function(x) { length(x$fit$par) }
## these assume we have REPORT(mu) and REPORT(resid) and ADREPORT(mu)
predict.TMB <- function(x, se.fit = FALSE) {
  if (!se.fit) x$obj$report()$mu
  sdr <- RTMB::sdreport(x$obj)  ## need the RTMB version
  sdr2 <- summary(sdr, select = "report")
  list(fit = sdr2[,"Estimate"],
       se.fit = sdr2[,"Std. Error"])
}

## -- phylogenetic Z-matrix construction -------------------------------------

##' construct a Z matrix from a phylogeny
##' @name phylo_machinery
##' @param r a \code{phylo} object
##' @param stand standardize edge lengths by determinant of phylogenetic covariance matrix
##' @export
##' @importFrom Matrix t diag
##' @importFrom ape vcv
## FIXME:: rename to phylo_to_Z ?
phylo.to.Z <- function(r, stand = FALSE) {
  ntip <- length(r$tip.label)
  Zid <- Matrix::Matrix(0.0, ncol = length(r$edge.length), nrow = ntip)
  nodes <- (ntip + 1):max(r$edge)
  root <- nodes[!(nodes %in% r$edge[, 2])]
  for (i in 1:ntip) {
    cn <- i ## current node
    while (cn != root) {
      ce <- which(r$edge[, 2] == cn) ## find current edge
      Zid[i, ce] <- 1 ## set Zid to 1
      cn <- r$edge[ce, 1] ## find previous node
    }
  }
  tZid <- t(Zid)
  Z <- t(sqrt(r$edge.length) * tZid)
  if (stand) {
    V <- ape::vcv(r)
    ## V <- V/max(V)
    sig <- exp(as.numeric(determinant(V)["modulus"]) / ntip)
    ## sig <- det(V)^(1/ntip)
    ## all.equal(Z/sqrt(sig),
    ##  t(sqrt(r$edge.length / sig) * tZid), tolerance=2e-16)
    Z <- Z/sqrt(sig)
  }
  rownames(Z) <- r$tip.label
  colnames(Z) <- 1:length(r$edge.length)
  return(Z)
}

## -- data cleanup helpers ---------------------------------------------------

gsub2 <- function(x, pattern, replacement) {
  gsub(pattern, replacement, x)
}

## typos/mismatches between data and tip labels
## (fix upstream?)
fixup_species <- function(x) {
  x |>
    gsub2("\\.", "_") |>
    gsub2("nelsoni$", "") |>
    gsub2("cinereorenteus", "cinereoargenteus") |>
    gsub2("hemonius", "hemionus") |>
    gsub2("thomsonii", "thompsonii") |>
    gsub2("Odicoileus", "Odocoileus") |>
    gsub2("_$", "")
}

## -- RTMB penalty-function factories ----------------------------------------

mk_f_phylo <- function(vcmat, scale) {
  function(x) dmvnorm(x, 0, Sigma = vcmat, scale = scale, log = TRUE) ## scale=1 for identifiability
}
## need to pass params
mk_f_cov <- function(corval, logrsd) {
    function(x) dmvnorm(x, rep(0,2),
                        Sigma = us2$corr(corval), scale = exp(logrsd),
                        log = TRUE)
}

## -- RTMB negative log-likelihoods for alternative parameterizations -------
## all rely on `chdat_x` (data + extras, via getAll) and on X/Z being set
## up in the calling environment

## dense covariance ("propto"-equivalent) parameterization
nllfun1 <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(X %*% beta + Z %*% b)
  REPORT(mu)
  ADREPORT(mu)
  resid <- log_rs - mu
  REPORT(resid)
  pen <- -dmvnorm(b, 0, Sigma = vcmat, scale = exp(logpsd), log = TRUE)
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsd), log = TRUE))
  lik + pen
}

## sparse precision-matrix parameterization
nllfun_prec <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(X %*% beta + Z %*% b)
  REPORT(mu)
  ADREPORT(mu)
  resid <- log_rs - mu
  REPORT(resid)
  pen <- -dgmrf(b, 0, Q = phyloprec, scale = exp(-logpsd), log = TRUE)
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsd), log = TRUE))
  lik + pen
}

## separable (phylogenetic x intercept-slope covariance) parameterization
## modular: should be able to handle various separable cor structures
nllfun_sep <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(X %*% beta + Z %*% c(b))
  REPORT(mu)
  ADREPORT(mu)
  resid <- log_rs - mu
  REPORT(resid)
  pen <- -dseparable(mk_f_phylo(phylomat, scale), mk_f_cov(corval, logrsd))(b)
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsdres), log = TRUE))
  lik + pen
}
