## phylo_tiny.R -- tiny 5-tip worked example: compare Z/Q under several
## phylogenetic random-effect parameterizations (see phyloslopes_utils.R /
## phyloslopes.qmd), and simulate data under the Kronecker-sum
## tensor-product (phylo x intercept-slope) model from phyloslopes.qmd's
## "tensor-product" section.

library(ape)
library(Matrix)
library(mgcv)      ## tensor.prod.model.matrix(), tensor.prod.penalties()
library(MRFtools)
library(MASS)       ## mvrnorm()

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

## -- tensor-product (phylo x intercept-slope) machinery ---------------------
## x = 1:5 (species 1..5, tip-label order); intercept + slope fixed-effect design
x <- as.numeric(chtree$tip.label)             ## == 1:5 in tip order
Xmat <- cbind(1, x)                           ## ntip x 2 fixed-effect design

## NB: tensor.prod.model.matrix()/tensor.prod.penalties() both use the
## convention "first list element = outer block, second = inner" (checked
## directly against small test matrices). The qmd's own tensor-product
## snippet builds tX from list(Z, X) (phylo outer, trait inner) but tp from
## list(diag(2), solve(vcmat)) (trait outer, phylo inner) -- reversed, so
## its tX/tp columns don't actually line up. Use the matched order here.
tX <- tensor.prod.model.matrix(list(as(as(Z_species, "generalMatrix"), "CsparseMatrix"),
                                     as(Xmat, "dgCMatrix")))
## columns ordered [species-outer, intercept/slope-inner]: sp1-int, sp1-slope,
## sp2-int, sp2-slope, ...

Sphylo <- Q_dense
Strait <- diag(2)
tp <- tensor.prod.penalties(list(Sphylo, Strait))
## tp[[1]] = kron(Sphylo, I_2): phylogenetic smoothness, shared across intercept & slope
## tp[[2]] = kron(I_ntip, I_2) = I_{2*ntip}: flat ridge, no phylo structure

cat(sprintf("\n=== tX (tensor-product model matrix, %d x %d) ===\n", nrow(tX), ncol(tX)))
print(round(as.matrix(tX), 3))
cat("\n=== tp[[1]] (phylo-structured penalty component) ===\n")
print(round(as.matrix(tp[[1]]), 3))
cat("\n=== tp[[2]] (flat-ridge penalty component) ===\n")
print(round(as.matrix(tp[[2]]), 3))

## -- simulate data under the Kronecker-sum tensor-product model -------------
## Q_tensor = (1/sigma_1^2) * tp[[1]]  +  (1/sigma_2^2) * tp[[2]]
## b ~ N(0, solve(Q_tensor)); sigma_1 = sigma_2 = 2 >> sigma_resid = 0.1 so
## both components dominate the residual noise (stress test for the
## identifiability problems noted in phyloslopes.qmd for this model fit to
## the real running-speed data)

set.seed(101)
sigma_1 <- 2         ## SD scale on the phylo-structured component
sigma_2 <- 2         ## SD scale on the flat-ridge component
sigma_resid <- 0.1
beta0 <- 0
beta1 <- 1

Q_tensor <- (1/sigma_1^2) * tp[[1]] + (1/sigma_2^2) * tp[[2]]
Sigma_tensor <- solve(Q_tensor)
b <- drop(MASS::mvrnorm(1, mu = rep(0, ncol(tX)), Sigma = as.matrix(Sigma_tensor)))

mu <- drop(Xmat %*% c(beta0, beta1) + tX %*% b)
y <- mu + rnorm(ntip, sd = sigma_resid)

sim_dat <- data.frame(species = chtree$tip.label, x = x, mu = mu, y = y)
cat("\n=== simulated data ===\n"); print(sim_dat)

## -- save everything for downstream use/fitting -----------------------------
save(chtree, vcmat, Q_dense, Q_tips, Q_full_raw, Q_noroot, Q_tips_mrf_schur,
     Z_edge, Z_species, Xmat, tX, tp, Q_tensor, b, sim_dat,
     sigma_1, sigma_2, sigma_resid, beta0, beta1,
     file = "phylo_tiny.rda")
