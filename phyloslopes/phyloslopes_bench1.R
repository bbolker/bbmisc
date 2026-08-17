## Benchmark of model-fitting time (MakeADFun + optimization, or glmmTMB())
## across the alternative parameterizations checked for consistency in
## phyloslopes_tests.R -- not the data/matrix setup, just the fit itself.
## Reads the design matrices/data built (and validated) there from
## phyloslopes_bench_setup.rds, rather than re-sourcing that whole script
## (and redoing every fit) just to rebuild them. Run
## `Rscript phyloslopes_tests.R` first if that file doesn't exist yet.
## Results saved to phylo_bench.rds / phylo_bench_slopes.rds; see
## phyloslopes_bench1_plot.R for the accompanying plot.
## Run with: Rscript phyloslopes_bench1.R

suppressMessages({
  library(glmmTMB)
  library(RTMB)
  library(Matrix)
  library(phyr)
  library(microbenchmark)
})
source("phyloslopes_utils.R")

list2env(readRDS("phyloslopes_bench_setup.rds"), envir = globalenv())

## nllfun_sep/nllfun_dense_slopes/nllfun_edge_slopes pick up `us2` as a free
## variable (via mk_f_cov() in phyloslopes_utils.R), not from chdat_x; it's
## a fixed, data-independent glmmTMB helper, so just rebuild it here rather
## than serializing it
us2 <- unstructured(2)

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
