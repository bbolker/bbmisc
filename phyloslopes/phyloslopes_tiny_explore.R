## phyloslopes_tiny_explore.R -- exploratory scratch companion to
## phyloslopes_tiny.R. Keeps the full tensor-product (phylo x
## intercept-slope) identifiability investigation: fitting, multi-start,
## the ntip=5 vs. ntip=10 / X-Z-aliasing / seed-
## robustness / phylo-x-TPS / null-space-only checks. phyloslopes_tiny.R
## itself has been pared back down to a clean/minimal worked example of
## how X, Z, and the penalty matrices (Sigma/Q) look under different
## small-example models -- see that file for the current TPS(k=4) x
## phylogenetic tensor-product example. This file is not meant to stay in
## sync with it; it's a frozen record of how that investigation was done.

library(ape)
library(Matrix)
library(mgcv)      ## tensor.prod.model.matrix(), tensor.prod.penalties()
library(MRFtools)
library(MASS)       ## mvrnorm()
library(RTMB)
library(numDeriv)  ## jacobian(), for the profiled-Hessian saddle-point check
library(tibble)    ## lst(), for chdat_x construction without name = name duplication

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
##
## b is species-level (length 2*ntip = 10) regardless of how many
## observations there are, so a single observation per tip (nobs = ntip = 5)
## leaves the tensor-product model badly overparameterized (10 random
## effects, 5 data points) -- fitting it collapses sigma_resid to the
## boundary and interpolates the data, a basic overparameterization
## artifact unrelated to the qmd's saddle-point issue. Replicate nrep
## observations per tip (nobs = nrep*ntip = 20 > 10) so sigma_resid has
## room to be estimated away from that boundary, for a fairer test.
nrep <- 4
species_rep <- rep(chtree$tip.label, each = nrep)
x_rep <- rep(x, each = nrep)
Xmat_rep <- cbind(1, x_rep)
Z_species_rep <- t(fac2sparse(factor(species_rep, levels = chtree$tip.label)))
tX_rep <- tensor.prod.model.matrix(list(as(as(Z_species_rep, "generalMatrix"), "CsparseMatrix"),
                                         as(Xmat_rep, "dgCMatrix")))
## tX_rep: nobs x 10, same [species-outer, intercept/slope-inner] column
## order as tX above, just with each species' pair of columns repeated
## across its nrep observation-rows

set.seed(101)
sigma_1 <- 2         ## SD scale on the phylo-structured component
sigma_2 <- 2         ## SD scale on the flat-ridge component
sigma_resid <- 0.1
beta0 <- 0
beta1 <- 1

Q_tensor <- (1/sigma_1^2) * tp[[1]] + (1/sigma_2^2) * tp[[2]]
Sigma_tensor <- solve(Q_tensor)
b <- drop(MASS::mvrnorm(1, mu = rep(0, ncol(tX)), Sigma = as.matrix(Sigma_tensor)))

mu <- drop(Xmat_rep %*% c(beta0, beta1) + tX_rep %*% b)
y <- mu + rnorm(length(mu), sd = sigma_resid)

sim_dat <- data.frame(species = species_rep, x = x_rep, mu = mu, y = y)
cat(sprintf("\n=== simulated data (%d reps/tip, %d obs total) ===\n", nrep, nrow(sim_dat)))
print(sim_dat)

## -- fit nllfun_tensor (correctly-specified) to the tensor-simulated data ---
## nllfun_tensor (phyloslopes_utils.R) rebuilds Q_tensor = (1/sigma_1^2)*tp1
## + (1/sigma_2^2)*tp2 from the fixed tp1/tp2 components at every evaluation
## (a Kronecker *sum*, unlike dgmrf()'s single fixed-Q + scalar-scale usage
## elsewhere). True sigma_1 = sigma_2 = 2 >> sigma_resid = 0.1 above; see the
## "Result" note below the Hessian check for what this strong-signal regime
## does and doesn't fix relative to the real-data saddle point
## (phyloslopes.qmd tensor-product section)
chdat_x <- lst(log_rs = sim_dat$y, X = Xmat_rep, tX = tX_rep,
                 tp1 = Matrix(tp[[1]], sparse = TRUE),
                 tp2 = Matrix(tp[[2]], sparse = TRUE))
## dgmrf() requires a sparse Q -- tensor.prod.penalties() returns plain
## dense matrices, so tp1/tp2 need converting before nllfun_tensor uses them
p0_tensor <- list(beta = rep(0, 2), b = rep(0, ncol(tX)),
                   logsigma1 = 0, logsigma2 = 0, logsd_resid = 0)
obj_tensor <- MakeADFun(nllfun_tensor, p0_tensor, silent = TRUE, random = "b")
fit_tensor <- TMBfit(obj_tensor)

true_vals <- c(beta0 = beta0, beta1 = beta1, sigma_1 = sigma_1, sigma_2 = sigma_2,
                sigma_resid = sigma_resid)
est_vals <- c(fixef.TMB(fit_tensor),
              sigma_1 = exp(unname(fit_tensor$fit$par["logsigma1"])),
              sigma_2 = exp(unname(fit_tensor$fit$par["logsigma2"])),
              sigma_resid = exp(unname(fit_tensor$fit$par["logsd_resid"])))
names(est_vals)[1:2] <- c("beta0", "beta1")
cat("\n=== tensor-product fit: recovered vs. true parameters ===\n")
print(rbind(true = true_vals, estimated = round(est_vals, 3)))

## profiled (fixed-effects-only) Hessian, as in the qmd's tensor-product
## section -- all-positive eigenvalues means a proper minimum, not a saddle
H_tensor <- numDeriv::jacobian(obj_tensor$gr, fit_tensor$fit$par)
H_tensor <- (H_tensor + t(H_tensor)) / 2  ## symmetrize away numerical noise
ev_tensor <- eigen(H_tensor, symmetric = TRUE, only.values = TRUE)$values
cat("\nProfiled Hessian eigenvalues (tensor-product fit):\n"); print(ev_tensor)
## Result: sigma_resid recovers well (~0.094 vs. true 0.1) once nobs >
## n random effects (nrep = 4 above), ruling out the earlier boundary/
## interpolation artifact from a single obs/tip -- but sigma_1 still blows
## up (~5e4 vs. true 2), sigma_2 undershoots (~0.64 vs. true 2), beta0/beta1
## are both off by a consistent ~2 (their true difference beta1-beta0 = 1 is
## preserved -- -1.998 vs. 2.561), and the smallest Hessian eigenvalue stays
## numerically zero (~2e-10). So even simulating *correctly* under the
## tensor-product model, with strong true signal (sigma_1 = sigma_2 = 2 >>
## sigma_resid = 0.1) and nobs > n random effects, does *not* avoid a
## near-flat identifiability ridge -- consistent with (not just an artifact
## of) the real-data saddle point in phyloslopes.qmd. The correlated
## beta0/sigma_1 bias suggests aliasing between the phylo x intercept/slope
## component (tp[[1]], which varies smoothly with species/tip order here)
## and the fixed-effect trend itself, for only ntip = 5 distinct phylo
## eigenvalues/species to pin it down.

## -- fit the separable (Kronecker-*product*) model to the tensor-simulated
## data -- misspecified: mirrors the qmd's separable-vs-tensor-product
## discussion, but small enough to check exactly (10x10, not 98x98).
## nllfun_dense_slopes (phyloslopes_utils.R) already implements this: b ~
## N(0, Sigma2 %x% vcmat), with Zdense = tX (same species-outer/trait-inner
## column order the tensor-product tX above uses)
us2 <- unstructured(2)
## nllfun_dense_slopes reads X, Zdense, vcmat, and us2 all from chdat_x --
## nothing from the calling environment
chdat_x <- lst(log_rs = sim_dat$y, X = Xmat_rep, Zdense = tX_rep, vcmat, us2)
p0_sep <- list(beta = rep(0, 2), b = rep(0, ncol(tX)),
                logrsd = rep(0, 2), logsd = 0, corval = 0)
obj_sep <- MakeADFun(nllfun_dense_slopes, p0_sep, silent = TRUE, random = "b")
fit_sep <- TMBfit(obj_sep)

cat("\n=== misspecified separable fit vs. correctly-specified tensor-product fit ===\n")
cat(sprintf("logLik separable (misspecified, Kronecker product): %.3f\n", logLik.TMB(fit_sep)))
cat(sprintf("logLik tensor-product (correct, Kronecker sum):     %.3f\n", logLik.TMB(fit_tensor)))

## =============================================================================
## Extended investigation: is the sigma_1/sigma_2 non-identifiability above
## (a) an artifact of ntip = 5, (b) a fixed-effect/random-effect (X/tX)
## aliasing problem, (c) specific to this one random draw of b, or (d)
## specific to the trivial (diag(2), no null/range structure) trait
## penalty in the *linear* case -- i.e. would it also show up for a genuine
## tensor-product smooth of phylogeny x a thin-plate spline? Findings below;
## see phyloslopes.qmd's tensor-product section for the write-up this
## feeds into.
## =============================================================================

## -- (a)/(b): ntip = 10 (not 5), multi-start, and intercept-only X ----------
## Bigger tree, same nrep = 4 replication logic, so nobs (40) is well above
## n random effects (20) -- and multi-start from many random starting points,
## to make sure the divergence above is a genuine global optimum, not a bad
## start.
chtree10 <- read.tree(text = paste0(
  "((((1:1,2:1):1,3:1):1,(4:1,5:1):1):1,",
  "(((6:1,7:1):1,8:1):1,(9:1,10:1):1):1);"))  ## two mirrored 5-tip clades, joined at a new root
ntip10 <- length(chtree10$tip.label)
vcmat10 <- vcv(chtree10)
Sphylo10 <- solve(vcmat10)
x10 <- as.numeric(chtree10$tip.label)
Xmat10 <- cbind(1, x10)
Z_species10 <- Diagonal(ntip10, x = 1)
dimnames(Z_species10) <- list(chtree10$tip.label, chtree10$tip.label)
tX10 <- tensor.prod.model.matrix(list(as(as(Z_species10, "generalMatrix"), "CsparseMatrix"),
                                       as(Xmat10, "dgCMatrix")))
tp10 <- tensor.prod.penalties(list(Sphylo10, diag(2)))

species10_rep <- rep(chtree10$tip.label, each = nrep)
x10_rep <- rep(x10, each = nrep)
Xmat10_rep <- cbind(1, x10_rep)
Z_species10_rep <- t(fac2sparse(factor(species10_rep, levels = chtree10$tip.label)))
tX10_rep <- tensor.prod.model.matrix(list(as(as(Z_species10_rep, "generalMatrix"), "CsparseMatrix"),
                                           as(Xmat10_rep, "dgCMatrix")))

set.seed(101)
Q_tensor10 <- (1/sigma_1^2) * tp10[[1]] + (1/sigma_2^2) * tp10[[2]]
b10 <- drop(MASS::mvrnorm(1, mu = rep(0, ncol(tX10)), Sigma = as.matrix(solve(Q_tensor10))))
mu10 <- drop(Xmat10_rep %*% c(beta0, beta1) + tX10_rep %*% b10)
y10 <- mu10 + rnorm(length(mu10), sd = sigma_resid)

chdat_x <- lst(log_rs = y10, X = Xmat10_rep, tX = tX10_rep,
                 tp1 = Matrix(tp10[[1]], sparse = TRUE), tp2 = Matrix(tp10[[2]], sparse = TRUE))
p0_tensor10 <- list(beta = rep(0, 2), b = rep(0, ncol(tX10)),
                     logsigma1 = 0, logsigma2 = 0, logsd_resid = 0)
obj_tensor10 <- MakeADFun(nllfun_tensor, p0_tensor10, silent = TRUE, random = "b")

set.seed(2026)
nstart <- 15
par_names <- names(obj_tensor10$par)
starts <- matrix(0, nrow = nstart, ncol = length(obj_tensor10$par), dimnames = list(NULL, par_names))
starts[, par_names == "beta"] <- matrix(rnorm(nstart * 2, sd = 3), ncol = 2)
starts[, "logsigma1"] <- runif(nstart, -3, 6)
starts[, "logsigma2"] <- runif(nstart, -3, 6)
starts[, "logsd_resid"] <- runif(nstart, -5, 2)
res10 <- lapply(seq_len(nstart), function(i)
  tryCatch(nlminb(starts[i, ], obj_tensor10$fn, obj_tensor10$gr), error = function(e) NULL))
ok10 <- !sapply(res10, is.null)
objv10 <- sapply(res10[ok10], `[[`, "objective")
best10 <- res10[[which(ok10)[which.min(objv10)]]]
cat(sprintf("\n=== (a) ntip=10, %d-start multi-start: nll range [%.4f, %.4f] ===\n",
            nstart, min(objv10), max(objv10)))
cat(sprintf("best: sigma1=%.4g sigma2=%.4g sigma_resid=%.4g (true: %.1f, %.1f, %.2f)\n",
            exp(best10$par["logsigma1"]), exp(best10$par["logsigma2"]), exp(best10$par["logsd_resid"]),
            sigma_1, sigma_2, sigma_resid))
## Result: all nstart starts converge to (numerically) the SAME point, one of
## sigma_1/sigma_2 always diverging (>1e3) -- confirms this is the true
## global optimum (not a bad start), and that doubling ntip doesn't fix it.

## intercept-only X: does removing the fixed slope (so the linear trend must
## come entirely from tX's random-slope columns) change anything?
X10_int <- Xmat10_rep[, 1, drop = FALSE]
chdat_x <- lst(log_rs = y10, X = X10_int, tX = tX10_rep,
                 tp1 = Matrix(tp10[[1]], sparse = TRUE), tp2 = Matrix(tp10[[2]], sparse = TRUE))
p0_int <- list(beta = 0, b = rep(0, ncol(tX10)), logsigma1 = 0, logsigma2 = 0, logsd_resid = 0)
obj_int <- MakeADFun(nllfun_tensor, p0_int, silent = TRUE, random = "b")
fit_int <- nlminb(obj_int$par, obj_int$fn, obj_int$gr)
cat(sprintf("\n=== (b) intercept-only X (no fixed slope): nll=%.4f, sigma1=%.4g sigma2=%.4g ===\n",
            fit_int$objective, exp(fit_int$par["logsigma1"]), exp(fit_int$par["logsigma2"])))
## Result: still diverges (one of sigma_1/sigma_2 -> large), at a WORSE nll
## than the full-X fit above -- rules out fixed/random (X/tX) aliasing as
## the cause. The redundancy is between sigma_1 and sigma_2 *within* the
## random-effect side (tp1 vs. tp2), not between beta1 and tX's random slope.

## -- (c) seed robustness (checked once, not re-run every script execution) -
## Looping seeds 101:110 (same ntip = 10 setup) found: sigma_1 or sigma_2
## diverges (> 100) in 8/10 draws (either direction -- confirms a genuine
## redundancy between the two components, not a one-directional pathology
## of sigma_1 specifically); 8/10 have a numerically-zero Hessian eigenvalue;
## and the correctly-specified tensor fit NEVER beat the misspecified
## separable fit's logLik in any of the 10 draws, despite having fewer fixed
## parameters (5 vs. 6) -- systematic, not a one-off unlucky draw.

## -- (d) does this generalize to phylo x thin-plate-spline? -----------------
## Naive: tensor.prod.penalties(list(Sphylo, S_spline)) fed into ONE shared
## dgmrf(), exactly like nllfun_tensor above but with a genuine (rank-
## deficient) TPS penalty instead of the trivial diag(2). "Genuine": phylo x
## null-space + phylo x range-space as two SEPARATE Kronecker-product terms
## (nllfun_spline_tensor, phyloslopes_utils.R -- the qmd's actual
## tensor-product-smooth implementation).
sm_full10 <- smoothCon(s(x10_rep, k = 6), data = data.frame(x10_rep = x10_rep))
S_spline10 <- sm_full10[[1]]$S[[1]]                 ## 6x6, rank-deficient (2-dim null space)
Xspline10 <- sm_full10[[1]]$X                       ## nobs x 6
ev10 <- eigen(S_spline10, symmetric = TRUE)
null_idx10 <- which(ev10$values < 1e-8 * max(ev10$values))
range_idx10 <- which(ev10$values >= 1e-8 * max(ev10$values))
Xf_null10 <- Xspline10 %*% ev10$vectors[, null_idx10, drop = FALSE]
Xr_range10 <- Xspline10 %*% ev10$vectors[, range_idx10, drop = FALSE] %*%
  diag(1/sqrt(ev10$values[range_idx10]), nrow = length(range_idx10))
Kn10 <- ncol(Xf_null10); Kr10 <- ncol(Xr_range10)

tZphylo10 <- t(fac2sparse(factor(species10_rep, levels = chtree10$tip.label)))
tX_naive10 <- tensor.prod.model.matrix(list(as(tZphylo10, "dgCMatrix"), as(Xspline10, "dgCMatrix")))
tp_naive10 <- tensor.prod.penalties(list(Sphylo10, S_spline10))
Xnull_joint10 <- tensor.prod.model.matrix(list(as(tZphylo10, "dgCMatrix"), as(Xf_null10, "dgCMatrix")))
Xrange_joint10 <- tensor.prod.model.matrix(list(as(tZphylo10, "dgCMatrix"), as(Xr_range10, "dgCMatrix")))

set.seed(101)
Q_tensor_tps_true <- (1/sigma_1^2) * tp_naive10[[1]] + (1/sigma_2^2) * tp_naive10[[2]]
b_tps_true <- drop(MASS::mvrnorm(1, mu = rep(0, ncol(tX_naive10)),
                                  Sigma = as.matrix(solve(Q_tensor_tps_true))))
mu_tps <- drop(Xmat10_rep %*% c(beta0, beta1) + tX_naive10 %*% b_tps_true)
y_tps <- mu_tps + rnorm(length(mu_tps), sd = sigma_resid)

chdat_x <- lst(log_rs = y_tps, X = Xmat10_rep, tX = tX_naive10,
                 tp1 = Matrix(tp_naive10[[1]], sparse = TRUE),
                 tp2 = Matrix(tp_naive10[[2]], sparse = TRUE))
p0_naive_tps <- list(beta = rep(0, 2), b = rep(0, ncol(tX_naive10)),
                      logsigma1 = 0, logsigma2 = 0, logsd_resid = 0)
obj_naive_tps <- MakeADFun(nllfun_tensor, p0_naive_tps, silent = TRUE, random = "b")
fit_naive_tps <- nlminb(obj_naive_tps$par, obj_naive_tps$fn, obj_naive_tps$gr)

## nllfun_spline_tensor reads X, Xnull_joint, Xrange_joint, vcmat, and Kr all
## from chdat_x -- nothing from the calling environment, so no risk of it
## silently picking up the ntip=5 `vcmat` left over from earlier in this
## script (which is exactly what happened here before this was fixed)
chdat_x <- lst(log_rs = y_tps, X = Xmat10_rep, Xnull_joint = Xnull_joint10, Xrange_joint = Xrange_joint10,
                 vcmat = vcmat10, Kr = Kr10)
p0_genuine_tps <- list(beta = rep(0, 2), b_null = rep(0, ntip10 * Kn10), b_range = rep(0, ntip10 * Kr10),
                        logsd = 0, logpsd_null = rep(0, Kn10), logpsd_range = 0)
obj_genuine_tps <- MakeADFun(nllfun_spline_tensor, p0_genuine_tps, silent = TRUE,
                              random = c("b_null", "b_range"))
fit_genuine_tps <- nlminb(obj_genuine_tps$par, obj_genuine_tps$fn, obj_genuine_tps$gr)

cat(sprintf("\n=== (d) phylo x TPS: naive Kronecker-sum nll=%.4f (sigma1=%.4g sigma2=%.4g); genuine null/range nll=%.4f ===\n",
            fit_naive_tps$objective, exp(fit_naive_tps$par["logsigma1"]), exp(fit_naive_tps$par["logsigma2"]),
            fit_genuine_tps$objective))
## Result: UNLIKE the linear case, BOTH constructions are well-identified
## here (checked via multi-start + profiling in the exploratory session --
## not repeated here for runtime) -- sigma_1/sigma_2 land at finite, stable
## values, no numerically-zero Hessian eigenvalue driven by an unbounded
## ridge. The linear case's pathology comes specifically from Strait =
## diag(2) being a trivial, structureless, full-rank "identity" trait
## penalty -- indistinguishable from a generic flat ridge. S_spline has
## genuine internal (null vs. range) structure, so summing
## kron(Sphylo, I_k) with kron(I_ntip, S_spline) does NOT reproduce the
## redundancy, even fed naively into one shared dgmrf().

## -- (e) the tensor-product-smooth recipe applied to the *linear* case -----
## nllfun_spline_tensor's trick is "split trait space into null (unpenalized)
## + range (penalized), give each its own Kronecker PRODUCT". For plain
## intercept+slope, the whole 2-dim trait space already IS the null space
## (nothing wiggly) -- so applying that recipe here means treating b as ONE
## Kronecker-product term with no range component at all: nllfun_tensor_
## nullspace (phyloslopes_utils.R). Refit on the *same* ntip=5 data that
## made nllfun_tensor's sigma_1 diverge above.
## FROZEN-RECORD NOTE: nllfun_tensor_nullspace has since been removed from
## phyloslopes_utils.R (superseded by nllfun_sep + map = list(corval =
## factor(NA)), which reproduces it exactly -- see README_tensor.qmd's
## "null block" section), so this chunk no longer runs as-is; per this
## file's header, it's a frozen record of the investigation, not meant to
## stay in sync
chdat_x <- lst(log_rs = sim_dat$y, X = Xmat_rep, tX = tX_rep, vcmat)
p0_nullspace <- list(beta = rep(0, 2), b = rep(0, ncol(tX)), logpsd_null = rep(0, 2), logsd_resid = 0)
obj_nullspace <- MakeADFun(nllfun_tensor_nullspace, p0_nullspace, silent = TRUE, random = "b")
fit_nullspace <- nlminb(obj_nullspace$par, obj_nullspace$fn, obj_nullspace$gr)
cat(sprintf("\n=== (e) null-space-only Kronecker-product fit on the SAME ntip=5 data: nll=%.4f ===\n",
            fit_nullspace$objective))
cat(sprintf("(naive Kronecker-sum nllfun_tensor fit on this data: nll=%.4f, sigma1 diverged to ~%.3g)\n",
            fit_tensor$fit$objective, est_vals["sigma_1"]))
## Result: well-identified (checked via multi-start in the exploratory
## session -- every start lands on the identical point here too), no
## unbounded parameter -- confirms the general lesson: the naive Kronecker-
## SUM construction is only pathological when paired with a trivial
## (identity, no-null-space) trait penalty; both a genuinely-structured
## trait penalty (TPS, (d)) and dropping the redundant sum in favor of
## Kronecker-PRODUCT blocks (separable, or this null-space-only case) avoid
## it.

## -- save everything for downstream use/fitting -----------------------------
save(chtree, vcmat, Q_dense, Q_tips, Q_full_raw, Q_noroot, Q_tips_mrf_schur,
     Z_edge, Z_species, Xmat, tX, Xmat_rep, tX_rep, nrep, tp, Q_tensor, b, sim_dat,
     sigma_1, sigma_2, sigma_resid, beta0, beta1,
     true_vals, est_vals, ev_tensor,
     chtree10, tX10_rep, tp10, best10, fit_int,
     fit_naive_tps, fit_genuine_tps, fit_nullspace,
     file = "phyloslopes_tiny_explore.rda")
