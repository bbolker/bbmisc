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

##' drop the root row/column from a whole-tree (tips + internal nodes) MRF
##' precision matrix, removing the Brownian-motion "no fixed origin"
##' rank deficiency
##' @param tree a \code{phylo} object
##' @param Q the whole-tree penalty matrix, e.g. from
##'   \code{MRFtools::mrf_penalty(tree, internal_nodes = TRUE)}
##' @param sparse convert the result to a sparse (\code{CsparseMatrix}) matrix?
##' @export
drop_mrf_root <- function(tree, Q, sparse = TRUE) {
  ntip <- length(tree$tip.label)
  ## mrf_penalty() orders [tips][internal nodes], with dimnames = tip labels
  ## then "N<node id>"; reorder to [root][other internal nodes][tips] *by
  ## name*, not position, so this doesn't silently misalign if mrf_penalty()'s
  ## internal ordering convention ever changes. Root identified the same way
  ## as in phylo.to.Z(): the internal node that never appears as an edge's
  ## child
  internal_ids <- (ntip + 1):(ntip + tree$Nnode)
  root_id <- internal_ids[!(internal_ids %in% tree$edge[, 2])]
  stopifnot(length(root_id) == 1)
  internal_names <- paste0("N", internal_ids)
  root_name <- paste0("N", root_id)
  ord_names <- c(root_name, setdiff(internal_names, root_name), tree$tip.label)
  stopifnot(all(ord_names %in% rownames(Q)))
  Q_full <- Q[ord_names, ord_names]
  if (sparse) Q_full <- as(Q_full, "CsparseMatrix")
  Q_full[-1, -1]
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
## us2 is an explicit argument (not looked up in the calling nllfun_*'s
## local frame) because mk_f_cov is a separate function -- getAll(params,
## chdat_x) inside nllfun_sep only binds names into *that* function's own
## execution frame, not into mk_f_cov's (mk_f_cov's free variables resolve
## via its own lexical closure, i.e. wherever it was *defined*, not wherever
## it's *called* from). Callers must pull `us2` out of chdat_x themselves
## (via their own getAll()) and pass it in explicitly.
mk_f_cov <- function(us2, corval, logrsd) {
    function(x) dmvnorm(x, rep(0,2),
                        Sigma = us2$corr(corval), scale = exp(logrsd),
                        log = TRUE)
}

## -- RTMB negative log-likelihoods for alternative parameterizations -------
## all rely on `chdat_x` (via getAll) for EVERYTHING they reference beyond
## `params` -- data (log_rs, ...), design/incidence matrices (X, Z and their
## analogues Zdense/KR/Zphylo/Xfull/Xr/...), the phylogenetic covariance
## (vcmat), precision matrices (phyloprec, tp1/tp2), dimension constants
## (Kw, Kr), and helper objects (us2). Nothing is read from the calling
## environment -- every nllfun_* call site must build a self-contained
## chdat_x with all of these, or the fit will error (missing binding) rather
## than silently pick up a stale/wrong-shape object left over from a
## previous fit in the same session (see phyloslopes_tiny_explore.R for a case where the
## latter actually happened, with `vcmat`/`Kr`)

## dense covariance ("propto"-equivalent) parameterization
## @knitr nllfun1
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

## independent phylogenetic intercepts and slopes: extends nllfun1 by adding
## a separate phylogenetic penalty for species-level slope deviations.
## Z (nobs x ntip) maps species intercepts to observations; Z_slope (nobs x
## ntip) maps species slopes to observations, typically as a weighted
## species indicator (diag(covariate) %*% Z or more generally Khatri-Rao
## product). Both b_int and b_slope are penalized independently via vcmat
nllfun_independent_slopes <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(X %*% beta + Z %*% b_int + Z_slope %*% b_slope)
  REPORT(mu)
  ADREPORT(mu)
  resid <- log_rs - mu
  REPORT(resid)
  pen_int <- -dmvnorm(b_int, 0, Sigma = vcmat, scale = exp(logpsd_int), log = TRUE)
  pen_slope <- -dmvnorm(b_slope, 0, Sigma = vcmat, scale = exp(logpsd_slope), log = TRUE)
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsd), log = TRUE))
  lik + pen_int + pen_slope
}

## @knitr nllfun_prec
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

#' @knitr nllfun_edge
## edge-based Z-matrix parameterization: Z (from phylo.to.Z()) maps tips to
## their ancestral edges, scaled by sqrt(edge.length), so b ~ iid N(0, tau^2)
## gives Cov(Z %*% b) = tau^2 * vcv(tree) -- no covariance/precision matrix
## needed, just a plain iid Gaussian prior on b
nllfun_edge <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(X %*% beta + Z %*% b)
  REPORT(mu)
  ADREPORT(mu)
  resid <- log_rs - mu
  REPORT(resid)
  pen <- -sum(dnorm(b, 0, sd = exp(logpsd), log = TRUE))
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsd), log = TRUE))
  lik + pen
}

#' @knitr nllfun_sep
## separable (phylogenetic x intercept-slope covariance) parameterization
## modular: should be able to handle various separable cor structures
## NB: Z (= t(rt$Zt) from mkReTrms) has columns ordered [species-outer,
## trait-inner] (each species' intercept/slope columns adjacent); b is
## ntip x 2 (species x trait), so it must be flattened the same way via
## c(t(b)), *not* c(b) (which would flatten trait-outer/species-inner and
## silently permute b relative to Z's columns)
nllfun_sep <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(X %*% beta + Z %*% c(t(b)))
  REPORT(mu)
  ADREPORT(mu)
  resid <- log_rs - mu
  REPORT(resid)
  pen <- -dseparable(mk_f_phylo(phylomat, scale), mk_f_cov(us2, corval, logrsd))(b)
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsdres), log = TRUE))
  lik + pen
}

## edge-based random-slopes: KR (from KhatriRao(t(phylo.to.Z(tree)), t(J)))
## maps each observation to 2 columns per edge (intercept-innovation,
## slope-innovation), so b (length 2*nedge) is naturally 96 iid-across-edges
## 2-blocks; Sigma2 gives their (shared, per-edge) 2x2 covariance. Brownian
## summation up the tree induces exactly the same species-level covariance
## structure as the separable/Kronecker model above (Sigma2 %x% vcmat) --
## see nllfun_dense_slopes for the direct (dense, ground-truth) equivalent
nllfun_edge_slopes <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(X %*% beta + KR %*% b)
  REPORT(mu)
  ADREPORT(mu)
  resid <- log_rs - mu
  REPORT(resid)
  b_mat <- matrix(b, ncol = 2, byrow = TRUE)
  D <- exp(logrsd)
  Sigma2 <- (matrix(D, ncol = 1) %*% matrix(D, nrow = 1)) * us2$corr(corval)
  pen <- -sum(dmvnorm(b_mat, rep(0, 2), Sigma = Sigma2, log = TRUE))
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsd), log = TRUE))
  lik + pen
}

## dense-Kronecker random-slopes: brute-force ground truth for the
## separable/edge-based models above. b (length 2*ntip) ~ N(0, Sigma2 %x%
## vcmat) directly, with Zdense (= t(rt$Zt), species-outer/trait-inner)
## mapping b to each observation's fitted intercept+slope contribution
nllfun_dense_slopes <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(X %*% beta + Zdense %*% b)
  REPORT(mu)
  ADREPORT(mu)
  resid <- log_rs - mu
  REPORT(resid)
  D <- exp(logrsd)
  Sigma2 <- (matrix(D, ncol = 1) %*% matrix(D, nrow = 1)) * us2$corr(corval)
  Sigma_full <- kronecker(vcmat, Sigma2)
  pen <- -dmvnorm(b, rep(0, length(b)), Sigma = Sigma_full, log = TRUE)
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsd), log = TRUE))
  lik + pen
}

## tensor-product (Kronecker-*sum*) random-slopes: b (length 2*ntip, ordered
## [species-outer, intercept/slope-inner] to match tX from
## tensor.prod.model.matrix(list(Z_species, X)), same layout nllfun_dense_slopes
## uses for its Kronecker-*product* Zdense/b) ~ N(0, Q_tensor^-1) with
## Q_tensor = (1/sigma_1^2)*tp1 + (1/sigma_2^2)*tp2. This is a Kronecker sum
## of two fixed penalty components (tp1 = phylo-structured, tp2 = flat-ridge,
## both from tensor.prod.penalties()), not a Kronecker product, so it can't
## be written via dgmrf()'s scalar `scale` argument on one fixed Q the way
## nllfun_prec/nllfun_spline_* are -- Q_tensor has to be rebuilt from tp1/tp2
## at every evaluation using the current sigma_1/sigma_2
nllfun_tensor <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(X %*% beta + tX %*% b)
  REPORT(mu)
  ADREPORT(mu)
  resid <- log_rs - mu
  REPORT(resid)
  Q_tensor <- (1/exp(logsigma1)^2) * tp1 + (1/exp(logsigma2)^2) * tp2
  pen <- -dgmrf(b, 0, Q = Q_tensor, log = TRUE)
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsd_resid), log = TRUE))
  lik + pen
}

## "tensor-product-smooth" construction (see nllfun_spline_tensor below)
## applied to a plain linear (intercept+slope) trait basis, which has no
## genuine null/range split of its own -- the whole 2-dim trait space *is*
## the (unpenalized) null space, nothing wiggly to separate out. So the
## "eigendecompose into null/range, give each its own Kronecker PRODUCT"
## recipe collapses to a single Kronecker-product term with no range
## component at all: structurally the separable-model family
## (nllfun_dense_slopes/nllfun_sep), just with independent (not correlated)
## per-trait-direction scales instead of a full 2x2 covariance. Verified
## (phyloslopes_tiny_explore.R) that this avoids the sigma_1/sigma_2 non-identifiability
## nllfun_tensor's Kronecker *sum* has for a trivial (identity, no-null-
## space) trait penalty like diag(2)
nllfun_tensor_nullspace <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(X %*% beta + tX %*% b)
  REPORT(mu)
  ADREPORT(mu)
  resid <- log_rs - mu
  REPORT(resid)
  Sigma_full <- kronecker(vcmat, diag(exp(2*logpsd_null), nrow = length(logpsd_null)))
  pen <- -dmvnorm(b, rep(0, length(b)), Sigma = Sigma_full, log = TRUE)
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsd_resid), log = TRUE))
  lik + pen
}

## -- RTMB negative log-likelihoods for phylogenetic + spline models ---------
## (McGillycuddy et al.'s propto+s() additive model, and separable/tensor-
## product-smooth extensions of it). All rely on chdat_x supplying the
## precomputed (fixed, not rebuilt per-iteration) model matrices, `vcmat`,
## and any dimension constants (Kw, Kr) -- nothing from the calling
## environment.

## additive: independent (i) baseline iid spline-wiggle random effect and
## (ii) phylogenetic random intercept, exactly matching the propto+s()
## model in McGillycuddy et al. (main.tex). Xf/Xr come from
## smooth2random(smoothCon(s(x), ..., absorb.cons=TRUE)[[1]], "", type=2).
## beta = c(intercept, linear-in-x coefficient), Xfull = cbind(1, Xf), so
## that fixef.TMB() (which reads params$beta) works as it does elsewhere
nllfun_spline_additive <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(Xfull %*% beta + Xr %*% b_spline + Zphylo %*% b_phylo)
  REPORT(mu)
  ADREPORT(mu)
  resid <- log_rs - mu
  REPORT(resid)
  pen_spline <- -sum(dnorm(b_spline, 0, sd = exp(logsd_f), log = TRUE))
  pen_phylo <- -dmvnorm(b_phylo, 0, Sigma = vcmat, scale = exp(logpsd), log = TRUE)
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsd), log = TRUE))
  lik + pen_spline + pen_phylo
}

## separable: additive model's structure, plus a phylogenetically-correlated
## copy of the spline-wiggle coefficients (single extra scale, diagonal --
## not a full covariance -- across the Kw wiggly dimensions). Properly
## *nests* the additive model as logpsd_f -> -Inf (needs the baseline
## b_spline term above, or it loses that flexibility entirely and can only
## collapse all the way to the intercept-only model)
nllfun_spline_separable <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(Xfull %*% beta + Xr %*% b_spline + Xr_joint %*% b_wiggly + Zphylo %*% b_phylo)
  REPORT(mu)
  ADREPORT(mu)
  resid <- log_rs - mu
  REPORT(resid)
  pen_spline <- -sum(dnorm(b_spline, 0, sd = exp(logsd_f), log = TRUE))
  ## Kw (wiggly dimension) comes from chdat_x (set alongside Xr, before
  ## chdat_x is built) -- not recomputed here, so it can't silently diverge
  ## from the Xr actually used to build Xr_joint
  Sigma_wiggly <- exp(2*logpsd_f) * diag(Kw)
  Sigma_full <- kronecker(vcmat, Sigma_wiggly)
  pen_wiggly <- -dmvnorm(b_wiggly, rep(0, length(b_wiggly)), Sigma = Sigma_full, log = TRUE)
  pen_phylo <- -dmvnorm(b_phylo, 0, Sigma = vcmat, scale = exp(logpsd), log = TRUE)
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsd), log = TRUE))
  lik + pen_spline + pen_wiggly + pen_phylo
}

## tensor-product smooth: phylo x [null-space of spline] (constant + linear,
## phylo-varying, diagonal/uncorrelated with *separate* scales for the two
## directions) + phylo x [range-space of spline] (wiggly, phylo-varying,
## single scale). Both terms use vcmat directly (always full rank), never a
## rank-deficient precision matrix -- the spline's own null space is
## resolved via eigendecomposition of its *marginal* penalty before
## crossing with phylo (see phyloslopes.qmd), not by feeding mgcv's raw
## (rank-deficient) tensor penalty straight into dgmrf()/dmvnorm(), which
## gives NaN. NB: needs *two* separate null-space scales (logpsd_null,
## length 2) -- a single shared scale can't shrink away the (unsupported)
## phylo-slope variation without also killing the (well-supported)
## phylo-intercept variation, giving a much worse fit
nllfun_spline_tensor <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(X %*% beta + Xnull_joint %*% b_null + Xrange_joint %*% b_range)
  REPORT(mu)
  ADREPORT(mu)
  resid <- log_rs - mu
  REPORT(resid)
  ## Kr (wiggly range-space dimension) comes from chdat_x (set alongside
  ## Xr_range, before chdat_x is built) -- not recomputed here, so it can't
  ## silently diverge from the basis actually used to build Xrange_joint
  Sigma_full_null <- kronecker(vcmat, diag(exp(2*logpsd_null), nrow = length(logpsd_null)))
  Sigma_full_range <- kronecker(vcmat, exp(2*logpsd_range) * diag(Kr))
  pen_null <- -dmvnorm(b_null, rep(0, length(b_null)), Sigma = Sigma_full_null, log = TRUE)
  pen_range <- -dmvnorm(b_range, rep(0, length(b_range)), Sigma = Sigma_full_range, log = TRUE)
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsd), log = TRUE))
  lik + pen_null + pen_range
}
