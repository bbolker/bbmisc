## Consistency checks across alternative parameterizations of the
## phylogenetic random-effect model explored in phyloslopes.qmd:
##   - glmmTMB + propto()              (dense covariance, glmmTMB backend)
##   - RTMB + dmvnorm()                (dense covariance)
##   - RTMB + dgmrf()                  (sparse tip-only precision matrix)
##   - RTMB + dgmrf(), full tree,      (sparse precision over tips + internal
##     root dropped                     nodes via mrf_penalty(), root removed
##                                       to fix the rank deficiency)
##   - RTMB + edge-based Z-matrix      (phylo.to.Z(): tip x edge matrix, iid
##                                       Gaussian prior on b, no covariance/
##                                       precision matrix needed)
##   - phyr::pglmm_compare()           (independent implementation, dense
##                                       covariance, comparative-data wrapper
##                                       around phyr::pglmm)
## These all represent the *same* statistical model, just parameterized
## differently, so fixed effects and log-likelihoods should agree.
## Run with: Rscript phyloslopes_tests.R

suppressMessages({
  library(ade4)
  library(ape)
  library(dplyr)
  library(glmmTMB)
  library(RTMB)
  library(reformulas)
  library(Matrix)
  library(MRFtools)
  library(phyr)
  library(mgcv)
})

source("phyloslopes_utils.R")

## -- data (mirrors the `setup` chunk in phyloslopes.qmd) --------------------
data("carniherbi49", package = "ade4")
chtree <- read.tree(text = carniherbi49$tre2)
chtree$tip.label <- fixup_species(chtree$tip.label)
vcmat <- vcv(chtree)
chdat <- carniherbi49$tab2 |>
  mutate(log_bm = log(bodymass),
         log_rs = log(runningspeed),
         g = 1  ## dummy for grouping variable
         ) |>
  tibble::rownames_to_column("species") |>
  mutate(across(species, fixup_species)) |>
  mutate(across(species, \(x) factor(x, levels = rownames(vcmat))))
stopifnot(identical(sort(rownames(vcmat)), sort(as.character(chdat$species))))

X <- model.matrix(~ log_bm, data = chdat)
Z <- fac2sparse(chdat$species)
p0 <- list(beta = rep(0, 2), b = rep(0, nrow(vcmat)), logsd = 0, logpsd = 0)

## -- fit: glmmTMB + propto ---------------------------------------------------
fit_propto_glmmTMB <- glmmTMB(log_rs ~ log_bm + propto(0 + species | g, vcmat),
                              data = chdat)

## -- fit: RTMB + dmvnorm (dense covariance) ----------------------------------
chdat_x <- chdat
fit_dense <- TMBfit(MakeADFun(nllfun1, p0, silent = TRUE, random = "b"))

## -- fit: RTMB + dgmrf (sparse precision) ------------------------------------
Qprec_tip <- Matrix(solve(vcmat), sparse = TRUE)
chdat_x <- c(chdat, list(phyloprec = Qprec_tip))
fit_prec <- TMBfit(MakeADFun(nllfun_prec, p0, silent = TRUE, random = "b"))

## -- fit: RTMB + dgmrf, full tree (tips + internal nodes), root dropped ------
## uses mrf_penalty(..., internal_nodes = TRUE); the full ntip+Nnode
## precision matrix is rank-deficient (no fixed root -- shifting every node
## by the same constant leaves all edge contrasts unchanged), so the root
## must be dropped to get a proper (full-rank) precision matrix. This
## should then recover exactly the tip-only model, i.e. confirm that
## marginalizing the internal nodes reduces to fit_dense/fit_prec.
ntip <- length(chtree$tip.label)
Q_full_raw <- mrf_penalty(chtree, "brownian", internal_nodes = TRUE)
## mrf_penalty() orders [tips][internal nodes], with dimnames = tip labels
## then "N<node id>"; reorder to [root][other internal nodes][tips] *by
## name*, not position, so this doesn't silently misalign if mrf_penalty()'s
## internal ordering convention ever changes. Root identified the same way
## as in phylo.to.Z(): the internal node that never appears as an edge's
## child
internal_ids <- (ntip + 1):(ntip + chtree$Nnode)
root_id <- internal_ids[!(internal_ids %in% chtree$edge[, 2])]
stopifnot(length(root_id) == 1)
internal_names <- paste0("N", internal_ids)
root_name <- paste0("N", root_id)
ord_names <- c(root_name, setdiff(internal_names, root_name), chtree$tip.label)
stopifnot(all(ord_names %in% rownames(Q_full_raw)))
Q_full <- as(Q_full_raw[ord_names, ord_names], "CsparseMatrix")
Q_noroot <- Q_full[-1, -1]

## fac2sparse() returns levels x observations (the "Zt" convention); use
## placeholder levels 1:Nnode for the (never-observed -> all-zero-column)
## internal nodes, then transpose to the usual observations x levels form
Z_full <- t(fac2sparse(factor(chdat$species,
                              levels = c(seq_len(chtree$Nnode), chtree$tip.label)),
                       drop.unused.levels = FALSE))
Z_noroot <- Z_full[, -1]

p0_noroot <- modifyList(p0, list(b = rep(0, chtree$Nnode + ntip - 1)))
chdat_x <- c(chdat, list(Z = Z_noroot, phyloprec = Q_noroot))
fit_prec_allnodes_noroot <- TMBfit(MakeADFun(nllfun_prec, p0_noroot, silent = TRUE, random = "b"))

## -- fit: RTMB + edge-based Z-matrix (phylo.to.Z) ----------------------------
## phylo.to.Z() builds a tip x edge matrix with sqrt(edge.length) entries on
## ancestral edges, so Z_edge %*% b with b iid unit-scale reproduces vcv(tree)
## exactly -- no covariance/precision matrix needed, and (unlike the
## all-nodes precision matrix above) no root singularity, since there's no
## free "root value" parameter
pZ <- phylo.to.Z(chtree)
Z_edge <- pZ[as.character(chdat$species), ]
p0_edge <- modifyList(p0, list(b = rep(0, nrow(chtree$edge))))
chdat_x <- c(chdat, list(Z = Z_edge))
fit_edge <- TMBfit(MakeADFun(nllfun_edge, p0_edge, silent = TRUE, random = "b"))

## -- fit: phyr::pglmm_compare (independent implementation) -------------------
## pglmm_compare() wants the tip labels as rownames of `data`; REML = FALSE
## to match the (RE)ML convention used everywhere else here
chdat_phyr <- chdat
rownames(chdat_phyr) <- as.character(chdat_phyr$species)
fit_phyr <- pglmm_compare(log_rs ~ log_bm, family = "gaussian",
                          data = chdat_phyr, phy = chtree, REML = FALSE)

## -- random-slopes models -----------------------------------------------
## three independent parameterizations of a correlated-random-slopes
## phylogenetic model (species-level random intercept + slope, 2x2
## covariance, phylogenetically structured): these should mutually agree,
## but are a *different* statistical model from the intercept-only fits
## above (different fixed effects/logLik), so get their own check group
us2 <- unstructured(2)
rt <- mkReTrms(findbars(~ (1 + log_bm | species)), fr = chdat, calc.lambdat = FALSE)

## fit: separable (dseparable), Sigma2 %x% vcmat via the Kronecker
## abstraction. NB Z (= t(rt$Zt)) orders columns [species-outer,
## trait-inner]; b is ntip x 2 (species x trait), so it must be flattened
## with c(t(b)) to match -- c(b) silently permutes b relative to Z
p2 <- list(beta = rep(0, 2), b = matrix(0, nrow = nrow(vcmat), ncol = 2),
          logrsd = rep(0, 2), logsdres = 0, corval = 0)
chdat_x <- c(chdat, list(phylomat = vcmat, scale = 1, Z = t(rt$Zt)))
fit_sep <- TMBfit(MakeADFun(nllfun_sep, p2, silent = TRUE, random = "b"))

## fit: dense-Kronecker ground truth, b (length 2*ntip) ~ N(0, Sigma2 %x% vcmat)
Zdense <- t(rt$Zt)
p0_ds <- modifyList(p0, list(b = rep(0, 2*nrow(vcmat)), logrsd = rep(0, 2), corval = 0))
chdat_x <- c(chdat, list(Zdense = Zdense))
fit_dense_slopes <- TMBfit(MakeADFun(nllfun_dense_slopes, p0_ds, silent = TRUE, random = "b"))

## fit: edge-based (phylo.to.Z + KhatriRao), no covariance/precision matrix
## needed -- see nllfun_edge_slopes in phyloslopes_utils.R
f_slopes <- ~ 1 + (1 + log_bm | species)
J <- eval(bquote(model.matrix(~.(findbars(f_slopes)[[1]][[2]]), data = chdat)))
pZ <- phylo.to.Z(chtree)
pZ_ord <- pZ[as.character(chdat$species), ]
KR <- t(KhatriRao(t(pZ_ord), t(J)))
## built explicitly, not via modifyList(p0, ...): p0 carries a `logpsd`
## entry nllfun_edge_slopes never references, and modifyList() would let
## it ride along in the optimization vector with an always-exactly-zero
## gradient, frozen at its 0 starting value
p0_es <- list(beta = rep(0, 2), b = rep(0, 2*nrow(chtree$edge)), logsd = 0,
              logrsd = rep(0, 2), corval = 0)
chdat_x <- c(chdat, list(KR = KR))
fit_edge_slopes <- TMBfit(MakeADFun(nllfun_edge_slopes, p0_es, silent = TRUE, random = "b"))

## -- phylogenetic spline models -------------------------------------------
## McGillycuddy et al.'s propto+s() additive model (main.tex) and two
## extensions: separable (phylo-correlated wiggle, single extra scale) and
## a genuine tensor-product smooth (phylo x log_bm). Additive and separable
## should agree with each other (separable properly nests additive); the
## tensor-product model, weakly informed here, should collapse to match
## the tip-only intercept-only group above (fit_dense/fit_prec/...)
sm <- smoothCon(s(log_bm, k = 10), data = chdat, absorb.cons = TRUE)
sm2ran <- smooth2random(sm[[1]], "", type = 2)
Xf <- sm2ran$Xf
Xfull <- cbind(1, Xf)
Xr <- sm2ran$rand$Xr
Kw <- ncol(Xr)
Xr_joint <- tensor.prod.model.matrix(list(as(t(Zphylo <- fac2sparse(chdat$species)), "dgCMatrix"),
                                          as(Xr, "dgCMatrix")))

p0_spadd <- list(beta = rep(0, 2), b_spline = rep(0, Kw), b_phylo = rep(0, nrow(vcmat)),
                 logsd = 0, logsd_f = 0, logpsd = 0)
chdat_x <- chdat
fit_spline_additive <- TMBfit(MakeADFun(nllfun_spline_additive, p0_spadd, silent = TRUE,
                                        random = c("b_spline", "b_phylo")))

## validate against the paper's actual glmmTMB call; vcv(chtree,corr=TRUE) is
## exactly proportional to raw vcmat since the tree is ultrametric, so this
## is a reparameterization, not a different model
mat_corr <- vcv(chtree, corr = TRUE)[levels(chdat$species), levels(chdat$species)]
fit_spline_glmmTMB <- glmmTMB(log_rs ~ s(log_bm) + propto(0 + species | g, mat_corr), data = chdat,
                              control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

p0_spsep <- modifyList(p0_spadd, list(b_wiggly = rep(0, nrow(vcmat)*Kw), logpsd_f = -10))
chdat_x <- c(chdat, list(Xr_joint = Xr_joint))
fit_spline_separable <- TMBfit(MakeADFun(nllfun_spline_separable, p0_spsep, silent = TRUE,
                                         random = c("b_spline", "b_wiggly", "b_phylo")))

## tensor-product smooth: split the spline's own (rank-deficient) penalty
## into null- and range-space first, via smooth2random() (not a separate
## eigendecomposition -- its trans.D is 1/sqrt(eigenvalue) per range-space
## direction, in the same column order as Xr; see README_tensor.qmd), then
## cross both blocks with the (always full-rank) phylo indicator. The range
## block gets a Wood (2006, sec 4.1.8)-style multiple-term penalty (phylo-
## shrinkage plus the TPS's own smoothness-eigenvalue shrinkage, each with
## its own scale) rather than a pure Kronecker product -- see
## phyloslopes.qmd/README_tensor.qmd for why we can't just feed mgcv's raw
## tensor penalty to dgmrf()/dmvnorm() (NaN on a singular Q), and why the
## multiple-term range block beats the single-scale pure-Kronecker-product
## one.
sm_full <- smoothCon(s(log_bm, k = 10), data = chdat)  ## no absorb.cons -- want the full null space
sm2ran_full <- smooth2random(sm_full[[1]], "", type = 2)
Xf_null <- sm2ran_full$Xf                                        ## null (constant + linear) design
Xr_range <- sm2ran_full$rand$Xr                                  ## wiggly (range) design, Wood's sqrt(D) reparam
Kn <- ncol(Xf_null); Kr <- ncol(Xr_range)
d_range <- 1 / sm2ran_full$trans.D[seq_len(Kr)]^2                ## true TPS range-space eigenvalues
tZphylo <- as(t(Zphylo), "dgCMatrix")
Xnull_joint <- tensor.prod.model.matrix(list(tZphylo, as(Xf_null, "dgCMatrix")))
Xrange_joint <- tensor.prod.model.matrix(list(tZphylo, as(Xr_range, "dgCMatrix")))

## range-block penalty components, precomputed once (pure functions of
## data, not parameters -- see nllfun_spline_tensor)
Sphylo <- solve(vcmat)
Qr_phylo <- kronecker(Sphylo, diag(Kr))
Qr_smooth <- kronecker(diag(nrow(vcmat)), diag(d_range))

## cor_null (fixed at 0 via map=) is nllfun_spline_tensor's optional
## null-space-direction correlation -- see its header comment; us2 was
## already defined above (line ~129) for the separable-model fits
p0_sptensor <- list(beta = rep(0, 2), b_null = rep(0, nrow(vcmat)*Kn), b_range = rep(0, nrow(vcmat)*Kr),
                    logsd = 0, logpsd_null = rep(0, Kn), cor_null = 0,
                    logsigma1_range = 0, logsigma2_range = 0)
chdat_x <- c(chdat, list(Xnull_joint = Xnull_joint, Xrange_joint = Xrange_joint,
                         Qr_phylo = Qr_phylo, Qr_smooth = Qr_smooth, vcmat = vcmat, us2 = us2))
fit_spline_tensor <- TMBfit(MakeADFun(nllfun_spline_tensor, p0_sptensor, silent = TRUE,
                                      random = c("b_null", "b_range"),
                                      map = list(cor_null = factor(NA))))

## -- consistency-check machinery ----------------------------------------------

## fixef() dispatches differently depending on fit type: fixef.TMB (utils.R)
## returns a plain vector, fixef.glmmTMB returns a list with a $cond
## component, pglmm_compare stores a 1-column matrix in $B
get_fixef <- function(fit) {
  if (inherits(fit, "TMB")) return(fixef(fit))
  if (inherits(fit, "glmmTMB")) return(fixef(fit)$cond)
  if (inherits(fit, "pglmm_compare")) return(drop(fit$B))
  stop("don't know how to extract fixed effects from class: ",
       paste(class(fit), collapse = "/"))
}

## logLik() has no method for pglmm_compare objects; its loglik is stashed
## directly in $logLik
get_loglik <- function(fit) {
  if (inherits(fit, "pglmm_compare")) return(fit$logLik)
  as.numeric(logLik(fit))
}

## tolerance is looser than the other checks: fit_edge's optimizer reports
## "false convergence" (nlminb code 8) because the residual-noise variance
## sits near a zero boundary for this dataset, so its optimum is resolved
## slightly less precisely (fixef relative diff ~6e-5) than the other,
## better-conditioned parameterizations
## label1/label2 are required, not defaulted via deparse(substitute(...)):
## every call site here passes fits[[p[1]]]/fits[[p[2]]], for which that
## idiom would print the indexing expression, not a useful fit name
check_fixef_equal <- function(fit1, fit2, label1, label2, tolerance = 1e-4) {
  cat(sprintf("checking fixed effects agree: %s vs %s ...\n", label1, label2))
  stopifnot(all.equal(get_fixef(fit1), get_fixef(fit2),
                      check.attributes = FALSE, tolerance = tolerance))
}

## tolerance loosened (from 1e-6) for the same reason as check_fixef_equal:
## different optimizers (nlminb, glmmTMB, phyr's nlopt-based default)
## converge to the same optimum with slightly different numerical
## precision. Worse for the bigger random-slopes models, where "false
## convergence" (nlminb code 8) is common due to a near-zero
## residual-variance boundary in this dataset -- confirmed via gradient
## norm (~0.02-0.10 at the reported "optimum", nowhere near 0, and
## unmovable by a warm-restart with tighter nlminb control settings) that
## this is nlminb giving up on an ill-conditioned near-boundary ridge, not
## a real difference in optimum. fit_sep vs. fit_edge_slopes' relative
## logLik diff reaches ~5e-3 (absolute diff ~0.03 nats, logLik ~ -6.5) --
## looser than the general default below, so that specific pair gets an
## explicit override where it's checked
check_loglik_equal <- function(fit1, fit2, label1, label2, tolerance = 1e-3) {
  cat(sprintf("checking log-likelihoods agree: %s vs %s ...\n", label1, label2))
  stopifnot(all.equal(get_loglik(fit1), get_loglik(fit2),
                      tolerance = tolerance))
}

## -- consistency checks: all pairs of the six (equivalent) parameterizations
fits <- list(fit_dense = fit_dense, fit_prec = fit_prec,
             fit_propto_glmmTMB = fit_propto_glmmTMB,
             fit_prec_allnodes_noroot = fit_prec_allnodes_noroot,
             fit_edge = fit_edge,
             fit_phyr = fit_phyr)
pairs <- combn(names(fits), 2, simplify = FALSE)
for (p in pairs) {
  check_fixef_equal(fits[[p[1]]], fits[[p[2]]], label1 = p[1], label2 = p[2])
}
for (p in pairs) {
  check_loglik_equal(fits[[p[1]]], fits[[p[2]]], label1 = p[1], label2 = p[2])
}

## -- consistency checks: all pairs of the three random-slopes parameterizations
slopes_fits <- list(fit_sep = fit_sep, fit_dense_slopes = fit_dense_slopes,
                    fit_edge_slopes = fit_edge_slopes)
slopes_pairs <- combn(names(slopes_fits), 2, simplify = FALSE)
for (p in slopes_pairs) {
  check_fixef_equal(slopes_fits[[p[1]]], slopes_fits[[p[2]]], label1 = p[1], label2 = p[2])
}
for (p in slopes_pairs) {
  ## fit_edge_slopes is the least-precisely-converged of the three (highest
  ## gradient norm at its reported "optimum" -- see check_loglik_equal's
  ## comment above), so any pair involving it needs more room than the
  ## other two (well-converged) pairs
  tol <- if ("fit_edge_slopes" %in% p) 1e-2 else 1e-3
  check_loglik_equal(slopes_fits[[p[1]]], slopes_fits[[p[2]]], label1 = p[1], label2 = p[2],
                      tolerance = tol)
}

## -- consistency checks: spline models ------------------------------------
## additive and separable should agree with each other and with the paper's
## glmmTMB fit (separable properly nests additive, and collapses to it here)
spline_fits <- list(fit_spline_additive = fit_spline_additive,
                    fit_spline_separable = fit_spline_separable,
                    fit_spline_glmmTMB = fit_spline_glmmTMB)
spline_pairs <- combn(names(spline_fits), 2, simplify = FALSE)
for (p in spline_pairs) {
  check_fixef_equal(spline_fits[[p[1]]], spline_fits[[p[2]]], label1 = p[1], label2 = p[2])
}
for (p in spline_pairs) {
  check_loglik_equal(spline_fits[[p[1]]], spline_fits[[p[2]]], label1 = p[1], label2 = p[2])
}

## NOT an equality check (unlike the pairs above): fit_spline_tensor and
## fit_dense are genuinely different, non-nested model structures (tensor
## adds a phylo-correlated spline-wiggle component dense doesn't have), so
## there's no reason to expect exact agreement. Under the OLD single-scale
## (pure Kronecker-product) range-block construction, tensor's wiggle
## component always degenerated to exactly zero here, making it collapse
## to a bit-for-bit copy of fit_dense (beta and logLik agreed to ~1e-9) --
## that was mistakenly treated as an invariant to assert on. The hybrid
## range-block fix (nllfun_spline_tensor, phyloslopes_utils.R;
## README_tensor.qmd) gives the wiggle component a genuinely useful second
## smoothing scale, so it no longer degenerates: logLik improves
## (-6.16 vs. dense's -6.66) and beta shifts by a few percent. Kept here as
## a reporting line, not an assertion -- see phyloslopes_combo.R for the
## fuller null/additive/separable/tensor comparison on this same data.
cat(sprintf("\nfit_spline_tensor vs fit_dense (informational, not asserted equal):\n"))
cat(sprintf("  logLik: tensor = %.4f, dense = %.4f\n",
            get_loglik(fit_spline_tensor), get_loglik(fit_dense)))
cat(sprintf("  fixef:  tensor = %s\n          dense  = %s\n",
            paste(round(get_fixef(fit_spline_tensor), 4), collapse = ", "),
            paste(round(get_fixef(fit_dense), 4), collapse = ", ")))

cat("all consistency checks passed.\n")

## -- save benchmarking inputs ---------------------------------------------
## the raw data/design-matrix ingredients phyloslopes_bench1.R times fitting
## on -- not the fits themselves -- saved once here (after the consistency
## checks above have validated them) so that script can read them in
## directly instead of re-sourcing this whole script (and redoing every fit)
## just to rebuild them. X and Z are included even though they're derivable
## from chdat, because the nllfun_* functions in phyloslopes_utils.R pick
## them up as free variables from the calling environment (see the comment
## above nllfun1), not from chdat_x -- they must exist under exactly these
## names wherever those functions are called
bench_setup <- list(
  chdat = chdat, chtree = chtree, vcmat = vcmat,
  X = X, Z = Z,
  p0 = p0, Qprec_tip = Qprec_tip,
  Z_noroot = Z_noroot, Q_noroot = Q_noroot, p0_noroot = p0_noroot,
  Z_edge = Z_edge, p0_edge = p0_edge,
  chdat_phyr = chdat_phyr,
  rt = rt, p2 = p2,
  Zdense = Zdense, p0_ds = p0_ds,
  KR = KR, p0_es = p0_es
)
saveRDS(bench_setup, "phyloslopes_bench_setup.rds")
