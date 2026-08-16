## Consistency checks across alternative parameterizations of the
## phylogenetic random-effect model explored in phyloslopes.qmd:
##   - glmmTMB + propto()              (dense covariance, glmmTMB backend)
##   - RTMB + dmvnorm()                (dense covariance)
##   - RTMB + dgmrf()                  (sparse tip-only precision matrix)
##   - RTMB + dgmrf(), full tree,      (sparse precision over tips + internal
##     root dropped                     nodes via mrf_penalty(), root removed
##                                       to fix the rank deficiency)
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
## mrf_penalty() orders [tips][internal nodes]; reorder to [internal][tips]
## to match the factor levels used to build Z_full below, so that the root
## (the first internal node) ends up at position 1
ord <- c(ntip + seq_len(chtree$Nnode), seq_len(ntip))
Q_full <- as(Q_full_raw[ord, ord], "CsparseMatrix")
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

## -- consistency-check machinery ----------------------------------------------

## fixef() dispatches differently depending on fit type: fixef.TMB (utils.R)
## returns a plain vector, fixef.glmmTMB returns a list with a $cond component
get_fixef <- function(fit) {
  if (inherits(fit, "TMB")) return(fixef(fit))
  if (inherits(fit, "glmmTMB")) return(fixef(fit)$cond)
  stop("don't know how to extract fixed effects from class: ",
       paste(class(fit), collapse = "/"))
}

check_fixef_equal <- function(fit1, fit2, tolerance = 1e-6,
                              label1 = deparse(substitute(fit1)),
                              label2 = deparse(substitute(fit2))) {
  cat(sprintf("checking fixed effects agree: %s vs %s ...\n", label1, label2))
  stopifnot(all.equal(get_fixef(fit1), get_fixef(fit2),
                      check.attributes = FALSE, tolerance = tolerance))
}

check_loglik_equal <- function(fit1, fit2, tolerance = 1e-6,
                               label1 = deparse(substitute(fit1)),
                               label2 = deparse(substitute(fit2))) {
  cat(sprintf("checking log-likelihoods agree: %s vs %s ...\n", label1, label2))
  stopifnot(all.equal(as.numeric(logLik(fit1)), as.numeric(logLik(fit2)),
                      tolerance = tolerance))
}

## -- consistency checks: all pairs of the four (equivalent) parameterizations
fits <- list(fit_dense = fit_dense, fit_prec = fit_prec,
             fit_propto_glmmTMB = fit_propto_glmmTMB,
             fit_prec_allnodes_noroot = fit_prec_allnodes_noroot)
pairs <- combn(names(fits), 2, simplify = FALSE)
for (p in pairs) {
  check_fixef_equal(fits[[p[1]]], fits[[p[2]]], label1 = p[1], label2 = p[2])
}
for (p in pairs) {
  check_loglik_equal(fits[[p[1]]], fits[[p[2]]], label1 = p[1], label2 = p[2])
}

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

  bench <- microbenchmark(
    glmmTMB_propto = glmmTMB(log_rs ~ log_bm + propto(0 + species | g, vcmat),
                             data = chdat),
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
    times = 20
  )

  saveRDS(bench, "phylo_bench.rds")
  print(bench)

  bench_plot <- autoplot(bench) + theme_bw() +
    ggtitle("phyloslopes: fitting time by method")
  ggsave("phyloslopes_bench.png", bench_plot, width = 8, height = 5)
}
