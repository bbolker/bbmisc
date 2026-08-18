## phyloslopes_tiny.R -- clean/minimal tiny (5-tip) worked example. Two
## parts: (1) how X, Z, and the penalty matrices (Sigma/Q) look for several
## phylogenetic random-effect parameterizations (phyloslopes_utils.R /
## phyloslopes.qmd); (2) the primary example -- fit and compare additive,
## separable, and tensor-product combinations of phylogenetic variation with
## a low-k (k=4) thin-plate spline on the same small simulated data set. For
## the (much larger) fitting/identifiability investigation this grew out of
## -- including why the naive Kronecker-sum tensor-product construction is
## avoided here -- see phyloslopes_tiny_explore.R and TODO.md item 3.

library(ape)
library(Matrix)
library(mgcv)      ## smoothCon(), tensor.prod.model.matrix(), tensor.prod.penalties()
library(MRFtools)
library(MASS)      ## mvrnorm()

source("phyloslopes_utils.R")  ## phylo.to.Z(), drop_mrf_root()

## -- tiny deterministic tree ------------------------------------------------
## (((1,2),3),(4,5)), all edges length 1 -- not ultrametric, but simpler
## numbers and doesn't affect the comparison of parameterizations
chtree <- read.tree(text = "(((1:1,2:1):1,3:1):1,(4:1,5:1):1);")
ntip <- length(chtree$tip.label)
stopifnot(ntip == 5)

## -- covariance / precision under several parameterizations -----------------

vcmat <- vcv(chtree)                          ## dense Sigma (tips only)
Q_dense <- solve(vcmat)                       ## dense Q, hand-inverted from Sigma

Q_tips <- Matrix(Q_dense, sparse = TRUE)
## sparse tips-only precision. This is literally solve(vcmat) in sparse
## storage -- matching `Qprec_tip` in phyloslopes_tests.R / the `phylo_prec`
## chunk in phyloslopes.qmd, which is the "sparse tip-only precision"
## parameterization actually used for fitting in this project.

Q_full_raw <- mrf_penalty(chtree, "brownian", internal_nodes = TRUE)
Q_noroot <- drop_mrf_root(chtree, Q_full_raw)
## whole-tree (tips + internal nodes) joint precision, root row/col dropped.
## Z only has nonzero entries in the tip columns (internal-node b's never
## enter mu), so the observation-level covariance this implies is the tip
## x tip block of solve(Q_noroot) -- checked below to equal vcmat exactly.

Z_edge <- phylo.to.Z(chtree)                  ## tip x edge Z (sqrt(edge.length) entries)
## Z_edge %*% t(Z_edge) reproduces vcmat exactly (checked below)

Z_species <- Diagonal(ntip, x = 1)
dimnames(Z_species) <- list(chtree$tip.label, chtree$tip.label)
## one observation per species -> plain species-indicator Z is just I_ntip

## sanity checks
Sigma_noroot <- solve(as.matrix(Q_noroot))
stopifnot(all.equal(Sigma_noroot[chtree$tip.label, chtree$tip.label], vcmat,
                     check.attributes = FALSE))
stopifnot(all.equal(as.matrix(Z_edge %*% t(Z_edge)), vcmat, check.attributes = FALSE))

cat("=== vcmat (dense Sigma) ===\n"); print(round(vcmat, 3))
cat("\n=== Q_dense / Q_tips (solve(vcmat), sparse tips-only precision) ===\n"); print(round(Q_dense, 3))
cat("\n=== Q_noroot (MRFtools::mrf_penalty, full tree, root dropped) ===\n"); print(round(as.matrix(Q_noroot), 3))
cat("\n=== Z_edge (tip x edge) ===\n"); print(round(as.matrix(Z_edge), 3))

## NB: MRFtools::mrf_penalty(chtree, "brownian", internal_nodes = FALSE)
## computes a *different* tips-only precision than solve(vcmat) -- it Schur-
## complements the internal nodes (including the root) out of the
## rank-deficient joint precision, which implicitly gives the root an
## improper (flat/infinite-variance) prior rather than fixing it at a known
## value the way ape::vcv()/solve(vcmat) does. That induces nonzero
## covariance even between tips in different top-level clades (e.g. tips 1
## and 4 here), which vcmat has as an exact zero -- not a bug in MRFtools or
## in Q_dense, just a different model. The project's own code never
## actually uses this combination for fitting (the validated tips-only
## precision route is always solve(vcmat), as above); kept here only to
## document the discrepancy.
Q_tips_mrf_schur <- mrf_penalty(chtree, "brownian", internal_nodes = FALSE)
cat("\n=== Q_tips_mrf_schur (MRFtools tips-only Schur complement -- NOT solve(vcmat), see comment) ===\n")
print(round(as.matrix(Q_tips_mrf_schur), 3))

## =============================================================================
## Primary example: fit and compare additive, separable, and tensor-product
## combinations of phylogenetic variation with a low-k (k=4) thin-plate
## spline on the same small simulated data set -- the same three structures
## phyloslopes.qmd's "phylogenetic splines" section fits to the real
## 49-species data (nllfun_spline_additive/_separable/_tensor,
## phyloslopes_utils.R), scaled down to something fully inspectable. The
## naive Kronecker-*sum* tensor-product construction (tensor.prod.penalties()
## fed into one shared dgmrf()) is deliberately NOT used here -- it has a
## genuine, reproducible non-identifiability whenever the trait-side penalty
## is trivial, and even where it doesn't (this spline case) it's just the
## less-clean way to build the same tensor-product-smooth model; see
## phyloslopes_tiny_explore.R for that whole investigation.
## =============================================================================

library(RTMB)
library(tibble)    ## lst(), for chdat_x construction without name = name duplication

## -- k=4 thin-plate-spline bases + phylogenetic incidence matrix -----------
## Replicate nrep observations/tip (nobs = nrep*ntip = 20) so the random
## effects below aren't hopelessly overparameterized relative to the data --
## a single obs/tip collapses everything to an interpolating, non-identified
## fit (see phyloslopes_tiny_explore.R for that failure mode).
nrep <- 4
x <- as.numeric(chtree$tip.label)             ## == 1:5 in tip order
species_rep <- rep(chtree$tip.label, each = nrep)
x_rep <- rep(x, each = nrep)
## nobs x ntip species-incidence matrix (fac2sparse() returns the transposed
## "levels x observations" convention, so transpose it back)
Zphylo <- t(fac2sparse(factor(species_rep, levels = chtree$tip.label)))

## additive/separable: absorb.cons=TRUE (matches glmmTMB's own s() handling)
## splits the k=4 (3 after the sum-to-zero constraint) basis into a fixed
## linear-in-x column (Xf) and Kw wiggly random-effect columns (Xr)
sm <- smoothCon(s(x_rep, k = 4), data = data.frame(x_rep = x_rep), absorb.cons = TRUE)
sm2ran <- smooth2random(sm[[1]], "", type = 2)
Xf <- sm2ran$Xf                               ## fixed (linear-in-x) design, nobs x 1
Xfull <- cbind(1, Xf)                         ## + intercept, beta = c(intercept, linear coef)
Xr <- sm2ran$rand$Xr                          ## wiggly random-effect design, nobs x Kw
Kw <- ncol(Xr)
## cross the wiggly basis with the phylo indicator for the separable model
## (species-outer/wiggly-inner, matching tensor.prod.model.matrix() convention)
Xr_joint <- tensor.prod.model.matrix(list(as(Zphylo, "dgCMatrix"), as(Xr, "dgCMatrix")))

## tensor-product smooth: needs the *full* (non-absorb.cons) null space --
## eigendecompose the spline's own marginal penalty into null (constant +
## linear, unpenalized, Kn = 2) and range (wiggly, penalized, Kr = 2) space
## first, rather than feeding the raw rank-deficient penalty into a shared
## Kronecker-sum precision (see phyloslopes_tiny_explore.R for why not)
sm_full <- smoothCon(s(x_rep, k = 4), data = data.frame(x_rep = x_rep))
S_spline <- sm_full[[1]]$S[[1]]               ## 4x4 penalty, rank-deficient (2-dim null space)
Xspline <- sm_full[[1]]$X                     ## nobs x 4 spline basis
ev <- eigen(S_spline, symmetric = TRUE)
null_idx <- which(ev$values < 1e-8 * max(ev$values))
range_idx <- which(ev$values >= 1e-8 * max(ev$values))
Xf_null <- Xspline %*% ev$vectors[, null_idx, drop = FALSE]
Xr_range <- Xspline %*% ev$vectors[, range_idx, drop = FALSE] %*%
  diag(1/sqrt(ev$values[range_idx]), nrow = length(range_idx))   ## Wood's sqrt(D) reparameterization
Kn <- ncol(Xf_null); Kr <- ncol(Xr_range)
stopifnot(Kn == 2, Kr == 2)
Xnull_joint <- tensor.prod.model.matrix(list(as(Zphylo, "dgCMatrix"), as(Xf_null, "dgCMatrix")))
Xrange_joint <- tensor.prod.model.matrix(list(as(Zphylo, "dgCMatrix"), as(Xr_range, "dgCMatrix")))

cat(sprintf("\n=== Xfull (fixed-effect design, %d x %d: intercept + linear-in-x) ===\n",
            nrow(Xfull), ncol(Xfull)))
print(round(Xfull, 3))
cat(sprintf("\n=== Xr (wiggly spline basis, %d x %d) ===\n", nrow(Xr), ncol(Xr)))
print(round(Xr, 3))
cat(sprintf("\n=== Xr_joint (phylo x wiggly tensor-product matrix, %d x %d) ===\n",
            nrow(Xr_joint), ncol(Xr_joint)))
print(round(as.matrix(Xr_joint), 3))
cat(sprintf("\n=== Xnull_joint (phylo x null-space tensor-product matrix, %d x %d) ===\n",
            nrow(Xnull_joint), ncol(Xnull_joint)))
print(round(as.matrix(Xnull_joint), 3))
cat(sprintf("\n=== Xrange_joint (phylo x range-space tensor-product matrix, %d x %d) ===\n",
            nrow(Xrange_joint), ncol(Xrange_joint)))
print(round(as.matrix(Xrange_joint), 3))

## -- simulate one data set under the separable model -------------------------
## mu = intercept + linear-in-x + (shared/average wiggle shape) +
## (phylogenetically-correlated per-species DEVIATION from that shape) +
## phylogenetic intercept. Separable is now the correctly-specified model;
## additive is its nested special case (sigma_wiggly -> 0) and is therefore
## misspecified whenever sigma_wiggly is truly nonzero, as here.
## Tensor-product remains a genuinely different (non-nested) family -- see
## phyloslopes.qmd's tensor-product section for why summing Kronecker
## penalties can't reproduce a Kronecker product, and vice versa.
set.seed(101)
beta0 <- 0; beta1 <- 1
sd_f <- 0.3          ## SD of the shared (average) spline wiggle shape
sd_wiggly <- 0.5     ## SD of the phylogenetically-correlated per-species wiggle deviation
sd_phylo <- 1        ## SD of the phylogenetic intercept
sigma_resid <- 0.1

b_spline_true <- rnorm(Kw, sd = sd_f)
b_wiggly_true <- drop(MASS::mvrnorm(1, mu = rep(0, ntip * Kw),
                                     Sigma = sd_wiggly^2 * kronecker(vcmat, diag(Kw))))
b_phylo_true <- drop(MASS::mvrnorm(1, mu = rep(0, ntip), Sigma = sd_phylo^2 * vcmat))
mu <- drop(Xfull %*% c(beta0, beta1) + Xr %*% b_spline_true + Xr_joint %*% b_wiggly_true +
             Zphylo %*% b_phylo_true)
y <- mu + rnorm(length(mu), sd = sigma_resid)

sim_dat <- data.frame(species = species_rep, x = x_rep, mu = mu, y = y)
cat(sprintf("\n=== simulated data (%d reps/tip, %d obs total) ===\n", nrep, nrow(sim_dat)))
print(sim_dat)

## -- fit all three structures to the same simulated data --------------------
chdat_x <- lst(log_rs = y, Xfull, Xr, Zphylo, vcmat)
p0_add <- list(beta = rep(0, 2), b_spline = rep(0, Kw), b_phylo = rep(0, ntip),
               logsd = 0, logsd_f = 0, logpsd = 0)
obj_add <- MakeADFun(nllfun_spline_additive, p0_add, silent = TRUE,
                     random = c("b_spline", "b_phylo"))
fit_add <- TMBfit(obj_add)

chdat_x <- lst(log_rs = y, Xfull, Xr, Xr_joint, Zphylo, Kw, vcmat)
p0_sep <- modifyList(p0_add, list(b_wiggly = rep(0, ntip * Kw), logpsd_f = -10))
obj_sep <- MakeADFun(nllfun_spline_separable, p0_sep, silent = TRUE,
                     random = c("b_spline", "b_wiggly", "b_phylo"))
fit_sep <- TMBfit(obj_sep)

chdat_x <- lst(log_rs = y, X = Xfull, Xnull_joint, Xrange_joint, Kr, vcmat)
p0_tensor <- list(beta = rep(0, 2), b_null = rep(0, ntip * Kn), b_range = rep(0, ntip * Kr),
                  logsd = 0, logpsd_null = rep(0, Kn), logpsd_range = 0)
obj_tensor <- MakeADFun(nllfun_spline_tensor, p0_tensor, silent = TRUE,
                        random = c("b_null", "b_range"))
fit_tensor <- TMBfit(obj_tensor)

cat("\n=== fixed effects (intercept, linear-in-x) across structures ===\n")
print(rbind(true = c(beta0, beta1),
            additive = fixef.TMB(fit_add),
            separable = fixef.TMB(fit_sep),
            tensor = fixef.TMB(fit_tensor)))
cat("\n=== logLik across structures ===\n")
print(c(additive = logLik.TMB(fit_add), separable = logLik.TMB(fit_sep),
        tensor = logLik.TMB(fit_tensor)))
## Result: even though sigma_wiggly is genuinely nonzero here (0.5, larger
## than sigma_f = 0.3), separable's extra phylogenetically-correlated wiggle
## scale (logpsd_f) *still* collapses to its floor, exactly reproducing the
## additive fit (fixef and logLik match) -- and this holds across a range of
## true sigma_wiggly values from 0.3 up to 2 (checked interactively, not
## repeated here). So the earlier finding (separable always collapsing to
## additive) isn't just separable correctly recognizing an absent signal --
## at this tree's size (ntip = 5, nrep = 4 => 20 obs), there simply isn't
## enough information to identify a phylogenetically-correlated deviation
## from the shared wiggle shape (10 extra latent dimensions, b_wiggly,
## governed by one extra scale parameter) at all, regardless of its true
## magnitude. A genuine power limitation of this tiny example, not a
## property of the separable model itself. Tensor-product is NOT nested in
## either (different Kronecker family -- see phyloslopes.qmd's
## tensor-product section) and its fit relative to additive/separable is
## data-dependent, as before.

## -- save everything for downstream use/fitting -----------------------------
save(chtree, vcmat, Q_dense, Q_tips, Q_full_raw, Q_noroot, Q_tips_mrf_schur,
     Z_edge, Z_species, nrep, x, x_rep, species_rep, Zphylo,
     Xfull, Xr, Xr_joint, Xspline, S_spline, Xf_null, Xr_range, Kn, Kr,
     Xnull_joint, Xrange_joint, sim_dat,
     beta0, beta1, sd_f, sd_wiggly, sd_phylo, sigma_resid,
     b_spline_true, b_wiggly_true, b_phylo_true,
     file = "phyloslopes_tiny.rda")
