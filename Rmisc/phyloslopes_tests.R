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
p0_es <- modifyList(p0, list(b = rep(0, 2*nrow(chtree$edge)), logrsd = rep(0, 2), corval = 0))
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

## tensor-product smooth: eigendecompose the spline's own (rank-deficient)
## penalty directly, cross both null- and range-space with the (always
## full-rank) phylo indicator -- see phyloslopes.qmd for why we can't just
## feed mgcv's raw tensor penalty to dgmrf()/dmvnorm() (NaN on a singular Q)
sm_full <- smoothCon(s(log_bm, k = 10), data = chdat)  ## no absorb.cons -- want the full null space
S_spline <- sm_full[[1]]$S[[1]]
ev <- eigen(S_spline, symmetric = TRUE)
null_idx <- which(ev$values < 1e-8 * max(ev$values))
range_idx <- which(ev$values >= 1e-8 * max(ev$values))
Xf_null <- sm_full[[1]]$X %*% ev$vectors[, null_idx, drop = FALSE]
Xr_range <- sm_full[[1]]$X %*% ev$vectors[, range_idx, drop = FALSE] %*%
  diag(1/sqrt(ev$values[range_idx]), nrow = length(range_idx))
Kn <- ncol(Xf_null); Kr <- ncol(Xr_range)
tZphylo <- as(t(Zphylo), "dgCMatrix")
Xnull_joint <- tensor.prod.model.matrix(list(tZphylo, as(Xf_null, "dgCMatrix")))
Xrange_joint <- tensor.prod.model.matrix(list(tZphylo, as(Xr_range, "dgCMatrix")))

p0_sptensor <- list(beta = rep(0, 2), b_null = rep(0, nrow(vcmat)*Kn), b_range = rep(0, nrow(vcmat)*Kr),
                    logsd = 0, logpsd_null = rep(0, Kn), logpsd_range = 0)
chdat_x <- c(chdat, list(Xnull_joint = Xnull_joint, Xrange_joint = Xrange_joint))
fit_spline_tensor <- TMBfit(MakeADFun(nllfun_spline_tensor, p0_sptensor, silent = TRUE,
                                      random = c("b_null", "b_range")))

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
## residual-variance boundary in this dataset -- relative diffs there can
## reach ~1e-4, hence the looser tolerance for that group
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
  check_loglik_equal(slopes_fits[[p[1]]], slopes_fits[[p[2]]], label1 = p[1], label2 = p[2])
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

## the tensor-product smooth is weakly informed here and should collapse
## back to the tip-only intercept-only group -- a cross-group check
check_fixef_equal(fit_spline_tensor, fit_dense, label1 = "fit_spline_tensor", label2 = "fit_dense")
check_loglik_equal(fit_spline_tensor, fit_dense, label1 = "fit_spline_tensor", label2 = "fit_dense")

cat("all consistency checks passed.\n")

## -- benchmark (optional) -----------------------------------------------
## set do_bench <- TRUE before sourcing/running this script to time just
## the fitting procedure (MakeADFun + optimization, or glmmTMB()) for each
## method -- not the data/matrix setup above. Results saved to
## phylo_bench.rds.
if (!exists("do_bench")) do_bench <- FALSE

if (do_bench) {
  library(microbenchmark)
  library(ggplot2)

  ## glmmTMB() bundles R-level formula/data preprocessing together with
  ## MakeADFun() + optimization; to see how much of its relative slowness is
  ## processing overhead vs. actual fitting, pre-process once with
  ## doFit=FALSE (excluded from the timed comparison, like the other
  ## pre-built matrices above), then time only fitTMB(doOptim=FALSE)
  ## (glmmTMB's MakeADFun() call) + a plain nlminb() optimization -- this
  ## reproduces glmmTMB_propto's fixef/logLik exactly (verified separately)
  pre_glmmTMB <- glmmTMB(log_rs ~ log_bm + propto(0 + species | g, vcmat),
                        data = chdat, doFit = FALSE)

  bench <- microbenchmark(
    glmmTMB_propto = glmmTMB(log_rs ~ log_bm + propto(0 + species | g, vcmat),
                             data = chdat),
    glmmTMB_obj_only = {
      obj <- fitTMB(pre_glmmTMB, doOptim = FALSE)
      nlminb(obj$par, obj$fn, obj$gr)
    },
    RTMB_dense = {
      chdat_x <- chdat
      TMBfit(MakeADFun(nllfun1, p0, silent = TRUE, random = "b"))
    },
    RTMB_prec_tip = {
      chdat_x <- c(chdat, list(phyloprec = Qprec_tip))
      TMBfit(MakeADFun(nllfun_prec, p0, silent = TRUE, random = "b"))
    },
    RTMB_prec_allnodes_noroot = {
      chdat_x <- c(chdat, list(Z = Z_noroot, phyloprec = Q_noroot))
      TMBfit(MakeADFun(nllfun_prec, p0_noroot, silent = TRUE, random = "b"))
    },
    RTMB_edge = {
      chdat_x <- c(chdat, list(Z = Z_edge))
      TMBfit(MakeADFun(nllfun_edge, p0_edge, silent = TRUE, random = "b"))
    },
    phyr_compare = pglmm_compare(log_rs ~ log_bm, family = "gaussian",
                                 data = chdat_phyr, phy = chtree, REML = FALSE),
    times = 100
  )

  ## microbenchmark objects support rbind(), so a new method can be timed
  ## on its own and appended to a previously-saved result without redoing
  ## the (slower) full run, e.g.:
  ##   old <- readRDS("phylo_bench.rds")
  ##   new <- microbenchmark(newmethod = ..., times = 100)
  ##   bench <- rbind(old, new)

  saveRDS(bench, "phylo_bench.rds")
  print(bench)

  ## random-slopes methods are a different modeling task (more parameters,
  ## bigger random-effects vector) than the intercept-only methods above,
  ## so they get their own benchmark/plot rather than being mixed in
  bench_slopes <- microbenchmark(
    sep_dseparable = {
      chdat_x <- c(chdat, list(phylomat = vcmat, scale = 1, Z = t(rt$Zt)))
      TMBfit(MakeADFun(nllfun_sep, p2, silent = TRUE, random = "b"))
    },
    dense_Kronecker = {
      chdat_x <- c(chdat, list(Zdense = Zdense))
      TMBfit(MakeADFun(nllfun_dense_slopes, p0_ds, silent = TRUE, random = "b"))
    },
    edge_KhatriRao = {
      chdat_x <- c(chdat, list(KR = KR))
      TMBfit(MakeADFun(nllfun_edge_slopes, p0_es, silent = TRUE, random = "b"))
    },
    times = 100
  )

  saveRDS(bench_slopes, "phylo_bench_slopes.rds")
  print(bench_slopes)

  ## combined plot: hand-built (not autoplot()) violin plot, times on x,
  ## methods on y, faceted top/bottom by model class with free x scales
  ## (the two classes differ by ~1-2 orders of magnitude in fitting time).
  ## phyr_compare is excluded here (its timings are so much faster than
  ## everything else that it squashes the other violins) but stays in
  ## `bench`/phylo_bench.rds -- see its median/5%/95% printed below
  theme_set(theme_bw())
  bench_df_intercept <- as.data.frame(bench)
  bench_df_intercept <- droplevels(bench_df_intercept[bench_df_intercept$expr != "phyr_compare", ])
  bench_df <- rbind(
    transform(bench_df_intercept, model_type = "random-intercept only"),
    transform(as.data.frame(bench_slopes), model_type = "random-slopes (+ intercept)")
  )
  bench_df$time_ms <- bench_df$time / 1e6
  bench_df$model_type <- factor(bench_df$model_type,
                                levels = c("random-slopes (+ intercept)",
                                           "random-intercept only"))

  bench_plot <- ggplot(bench_df, aes(x = time_ms, y = expr)) +
    geom_violin(fill = "gray") +
    scale_x_log10() +
    facet_wrap(~ model_type, ncol = 1, scales = "free") +
    labs(x = "time (ms)", y = NULL,
        title = "phyloslopes: fitting time by method")
  ggsave("phyloslopes_bench.png", bench_plot, width = 8, height = 7)
}
