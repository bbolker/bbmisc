## Utility functions for phyloslopes.qmd

## -- plotting / TMB helpers ------------------------------------------------

## help("image-methods")
ifun <- function(x, main = "") image(Matrix(x), sub = "",
                          xlab = "", ylab = "", main = main)

## coerce a (possibly AD-typed) dense/sparse matrix into something
## RTMB::dgmrf() accepts as Q: dgmrf() needs an object with a proper "Dim"
## slot (a Matrix-package S4 class, or RTMB's own "adsparse"), which neither
## a plain base "matrix" nor a bare AD-typed dense matrix ("advector" with a
## dim attribute) has -- passing either directly gives an opaque
## "length-0 dimension vector" error from inside dgmrf()
to_Q <- function(P) {
  if (inherits(P, "advector")) as(as(P, "sparseMatrix"), "adsparse")
  else Matrix::Matrix(P, sparse = TRUE)
}

TMBfit <- function(obj, optimizer = nlminb, response = NULL) {
  fit <- optimizer(obj$par, obj$fn, obj$gr)
  ret <- list(fit = fit, obj = obj)
  if (!is.null(response)) attr(ret, "response") <- response
  class(ret) <- c("TMB", "list")
  ret
}
AIC.TMB <- function(x) { 2*x$fit$objective + 2*npar(x) }
## proper logLik-classed return (df = number of estimated fixed-effect
## parameters), not just a bare number -- needed for anything that relies
## on the standard logLik() contract (bbmle::AICtab(), anova(), the
## generic AIC()/BIC() fallback, ...), and still behaves as a plain number
## anywhere it's just used arithmetically/printed
logLik.TMB <- function(x) {
  val <- -1*x$fit$objective
  attr(val, "df") <- npar(x)
  class(val) <- "logLik"
  val
}
fixef.TMB <- function(x) { x$obj$env$parList()$beta }
npar <- function(x, ...) {
  UseMethod("npar")
}

## HACK, per glmmTMB's own troubleshooting vignette
## (https://cran.r-project.org/web/packages/glmmTMB/vignettes/troubleshooting.html):
## glmmTMB deliberately sets logLik() to NA for models with convergence
## problems (e.g. a non-positive-definite Hessian), to keep users from
## silently including a bad fit in a model comparison. This OVERRIDES that
## safety check globally for glmmTMB objects (once phyloslopes_utils.R is
## sourced) by reading the raw nlminb objective directly off the fit
## (fit$fit$objective), exactly as the vignette suggests -- needed so
## bbmle::AICtab() (and anything else calling logLik() generically) can
## compare phyloslopes_linear_models' `independent` fit alongside the
## others. For a normally-converged fit this returns the identical value
## glmmTMB's own logLik() would (both are just -objective at the nlminb
## optimum); it only changes behavior for the NA-flagged case. Use with
## the same caution the vignette gives: a flagged model may not actually
## be at its true optimum, so treat its logLik/AIC as approximate, not a
## substitute for fixing the underlying fit.
logLik.glmmTMB <- function(object, ...) {
  val <- -object$fit$objective
  attr(val, "df") <- length(object$fit$par)
  ## nrow(object$frame), not stats::nobs(object): nobs() would do its own
  ## S3 dispatch to nobs.glmmTMB, which is only registered if the glmmTMB
  ## package is actually attached (library(glmmTMB)) in the calling
  ## session -- object$frame is a plain data.frame already stored on the
  ## object, so this works even when phyloslopes_utils.R is sourced (and
  ## a saved phyloslopes_linear.rda loaded) without glmmTMB attached
  attr(val, "nobs") <- nrow(object$frame)
  class(val) <- "logLik"
  val
}
## same hack, same caveat -- AIC.glmmTMB computes AIC from the NA-gated
## internal logLik rather than by calling the (now-overridden) logLik()
## generic, so it still propagates NA for a flagged fit unless this is
## *also* overridden; bbmle::AICtab() calls AIC() directly when available,
## so both overrides are needed for it to include a flagged model
AIC.glmmTMB <- function(object, ..., k = 2) {
  2*object$fit$objective + k*length(object$fit$par)
}
## logLik/AIC become S4 generics once stats4 is loaded (e.g. via
## library(bbmle)); their "ANY" fallback method's UseMethod() dispatch
## resolves S3 methods starting from the *calling package's own sealed
## namespace*, not from .GlobalEnv -- so a plain top-level function
## definition above is only visible to dispatch triggered directly from
## the global environment (e.g. typing logLik(fit) at the console/in a
## sourced script), not to dispatch triggered from *inside* another
## package's internal code (e.g. bbmle::ICtab()'s own
## sapply(L, function(x) logLik(x)), which resolved to NA even after the
## overrides above -- confirmed by testing bbmle::AICtab(...,
## logLik = TRUE) directly). registerS3method() inserts these into R's
## actual S3 method table, which every namespace's dispatch consults
## regardless of where the call originates.
registerS3method("logLik", "glmmTMB", logLik.glmmTMB)
registerS3method("AIC", "glmmTMB", AIC.glmmTMB)
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
phylo_to_Z <- function(r, stand = FALSE) {
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
  ## as in phylo_to_Z(): the internal node that never appears as an edge's
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
  resid <- log_rs - mu
  pen <- -dmvnorm(b, 0, Sigma = vcmat, scale = exp(logpsd), log = TRUE)
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsd), log = TRUE))
  lik + pen
}

## @knitr add_reports
## overly fancy machinery  for adding REPORT/ADREPORT statements to an existing objective function
## could add 'resid <- log_rs - mu' as well but we might want to be clever/fancy about substituting
##  the name of the response variable (i.e. not hard-coding log_rs)
## (the main point is so that the code included in the body of the report can be cleaner)
add_reports <- function(f) {
  b <- as.list(body(f))
  n <- length(b)                     ## includes the leading `{`
  new_stmts <- list(quote(REPORT(resid)), quote(REPORT(mu)), quote(ADREPORT(mu)))
  b <- append(b, new_stmts, after = n - 1)   ## insert right before the final statement
  body(f) <- as.call(b)
  f
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
## edge-based Z-matrix parameterization: Z (from phylo_to_Z()) maps tips to
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

## @knitr nllfun_sep_comment

## NB for code below: Z (= t(rt$Zt) from mkReTrms) has columns ordered [species-outer,
## trait-inner] (each species' intercept/slope columns adjacent); b is
## ntip x 2 (species x trait), so it must be flattened the same way via
## c(t(b)), *not* c(b) (which would flatten trait-outer/species-inner and
## silently permute b relative to Z's columns)

## @knitr nllfun_sep
## separable (phylogenetic x intercept-slope covariance) parameterization
## modular: should be able to handle various separable cor structures
nllfun_sep <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(X %*% beta + Z %*% c(t(b)))
  resid <- log_rs - mu
  pen <- -dseparable(mk_f_phylo(phylomat, scale), mk_f_cov(us2, corval, logrsd))(b)
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsdres), log = TRUE))
  lik + pen
}

## @knitr nllfun_edge
## edge-based random-slopes: KR (from KhatriRao(t(phylo_to_Z(tree)), t(J)))
## maps each observation to 2 columns per edge (intercept-innovation,
## slope-innovation), so b (length 2*nedge) is naturally 96 iid-across-edges
## 2-blocks; Sigma2 gives their (shared, per-edge) 2x2 covariance. Brownian
## summation up the tree induces exactly the same species-level covariance
## structure as the separable/Kronecker model above (Sigma2 %x% vcmat) --
## see nllfun_dense_slopes for the direct (dense, ground-truth) equivalent
nllfun_edge_slopes <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(X %*% beta + KR %*% b)
  resid <- log_rs - mu
  b_mat <- matrix(b, ncol = 2, byrow = TRUE)
  Sigma2 <- us2$corr(corval)
  pen <- -sum(dmvnorm(b_mat, Sigma = Sigma2, scale = exp(logrsd), log = TRUE))
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsd), log = TRUE))
  lik + pen
}

## dense-Kronecker random-slopes: brute-force ground truth for the
## separable/edge-based models above. b (length 2*ntip) ~ N(0, vcmat %x%
## Sigma2) directly, with Zdense (= t(rt$Zt), species-outer/trait-inner)
## mapping b to each observation's fitted intercept+slope contribution
nllfun_dense_slopes <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(X %*% beta + Zdense %*% b)
  REPORT(mu)
  ADREPORT(mu)
  resid <- log_rs - mu
  REPORT(resid)
  Sigma_full <- kronecker(vcmat, us2$corr(corval))
  pen <- -dmvnorm(b, Sigma = Sigma_full, scale = rep(exp(logrsd), nrow(vcmat)), log = TRUE)
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
## (nllfun_dense_slopes/nllfun_sep) -- which is why this no longer needs its
## own function. It used to live here as `nllfun_tensor_nullspace`, with
## independent (not correlated) per-trait-direction scales via
## `diag(exp(2*logpsd_null))`; `nllfun_sep` already IS that model plus a
## correlation parameter (`kronecker(phylomat, us2$corr(corval))`, `scale =
## exp(logrsd)`), so the old function is just `nllfun_sep` with
## `map = list(corval = factor(NA))` (and, for a single shared scale instead
## of two, additionally `logrsd = factor(c(1, 1))`) -- confirmed to reproduce
## it exactly (see README_tensor.qmd's "Isolating the mechanism" section).
## Verified (phyloslopes_tiny_explore.R) that giving the trait dimensions
## their own scales (however parameterized) avoids the sigma_1/sigma_2
## non-identifiability `nllfun_tensor`'s Kronecker *sum* has for a trivial
## (identity, no-null-space) trait penalty like diag(2).

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
# @knitr nllfun_spline_separable
nllfun_spline_separable <- function(params) {
  getAll(params, chdat_x)
  mu <- drop(Xfull %*% beta + Xr %*% b_spline + Xr_joint %*% b_wiggly + Zphylo %*% b_phylo)
  resid <- log_rs - mu
  Sigma_wiggly <- diag(Kw)
  Sigma_full <- kronecker(vcmat, Sigma_wiggly)
  pen_spline <- -sum(dnorm(b_spline, 0, sd = exp(logsd_f), log = TRUE))
  pen_phylo <- -dmvnorm(b_phylo, Sigma = vcmat, scale = exp(logpsd), log = TRUE)
  pen_wiggly <- -dmvnorm(b_wiggly, Sigma = Sigma_full, scale = exp(logpsd_f), log = TRUE)
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsd), log = TRUE))
  lik + pen_spline + pen_wiggly + pen_phylo
}

# @knitr tensor-comments
## tensor-product smooth: phylo x [null-space of spline] (constant + linear,
## phylo-varying, diagonal/uncorrelated with *separate* scales for the two
## directions) + phylo x [range-space of spline] (wiggly, phylo-varying,
## MULTIPLE-TERM Kronecker-*sum* penalty: phylo-shrinkage (Sphylo x I_Kr)
## plus the TPS's own smoothness-eigenvalue shrinkage (I_ntip x
## diag(d_range)), each with its own scale). The spline's own null space is
## resolved via eigendecomposition of its *marginal* penalty before
## crossing with phylo (see phyloslopes.qmd/README_tensor.qmd), not by
## feeding mgcv's raw (rank-deficient) tensor penalty straight into
## dgmrf()/dmvnorm(), which gives NaN.
##
## Two design choices, both justified at length in README_tensor.qmd (which
## also explains why RTMB doesn't choke on the rank-deficient TPS penalty
## S_spline in the first place -- only ever used inside an always-full-rank
## sum, never alone):
## - the null block needs *two* separate scales (logpsd_null, length 2): a
##   single shared scale can't shrink away the (unsupported) phylo-slope
##   variation without also killing the (well-supported) phylo-intercept
##   variation -- worth ~4.4 nats on the real 49-species data;
## - the range block needs the multiple-term (not pure-Kronecker-*product*)
##   recipe: per Wood (2006, sec 4.1.8), a pure product `Sphylo x I_Kr`
##   alone gives the optimizer only one knob (phylo-shrinkage), with no way
##   to separately trade off smoothness against fit the way GAM smoothing
##   normally does -- worth another ~0.5 nats, and beats even the fully
##   naive te()-extracted construction (which gets multiple-term range but
##   loses the separate null-space scales, since S2's null block is exactly
##   zero there, forcing one shared null scale no matter how the rest of
##   the model is built)
##
## The null block's `cor_null` parameter (new -- README_tensor.qmd, "null
## block" section) generalizes the two null-space scales' independence: it's
## the same `kronecker(phylomat, us2$corr(corval))` recipe `nllfun_sep` uses
## for the linear intercept+slope model, transplanted onto the spline's
## reparameterized null-space basis (`Xnull_joint`/`Xf_null`, which *is*
## `model.matrix(~1+x)` up to a linear reparameterization -- so this is
## really the same phylogenetic prior on an intercept+slope pair either way).
## `map = list(cor_null = factor(NA))` fixes it at 0 (independent, the
## behavior above), reproducing the old hard-coded `diag(Kn)` construction
## exactly; freeing it lets the phylo-intercept and phylo-slope null-space
## directions correlate (TODO item 2) -- worth ~0.56 nats on the real data,
## but lands near a correlation boundary (rho ~ -0.996), so treat that gain
## as suggestive, not a reason to make it the default

## @knitr nllfun_spline_tensor
nllfun_spline_tensor <- function(params) {
  getAll(params, chdat_x)
  ## b_null is an ntip x Kn matrix (natural shape for dseparable()),
  ## so it must be flattened via c(t(b_null)) -- not c(b_null) -- to match
  ## Xnull_joint's [species-outer, null-dim-inner] column order
  mu <- drop(X %*% beta + Xnull_joint %*% c(t(b_null)) + Xrange_joint %*% b_range)
  resid <- log_rs - mu
  ## null block: the same phylo x unstructured(Kn) covariance nllfun_sep puts
  ## on a plain intercept+slope pair (see README_tensor.qmd's "null block"
  ## section)
  pen_null <- -dseparable(mk_f_phylo(vcmat, 1), mk_f_cov(us2, cor_null, logpsd_null))(b_null)
  ## range block
  Q_range <- exp(-2*logsigma1_range) * Qr_phylo + exp(-2*logsigma2_range) * Qr_smooth
  pen_range <- -dgmrf(b_range, Q = to_Q(Q_range), log = TRUE)
  lik <- -sum(dnorm(log_rs, mean = mu, sd = exp(logsd), log = TRUE))
  lik + pen_null + pen_range
}
