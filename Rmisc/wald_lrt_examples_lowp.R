## Exact coverage and median CIs for the two-group logistic regression
## example in wald_lrt_examples.qmd, but with a much lower baseline
## probability (p0 = 0.01 rather than 0.2). Produces the analogues of
## the last two figures in the qmd (median estimates/CI limits and
## one- and two-sided coverage).
##
## The CIs for every possible outcome (y0, y1) don't depend on p0 (it
## only enters through the binomial weights), so they are the same as
## in the qmd; they're saved in an .rds file so that the weighting step
## can be rerun cheaply for other baselines.

library(glmmTMB)
library(car)
suppressPackageStartupMessages(library(tidyverse))
library(futurize)
library(patchwork)

## avoid overcommitting threads
library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)

## adjust this to something that makes sense on your machine (memory
## as well as cores), or set the WALD_WORKERS environment variable
nworkers <- as.integer(Sys.getenv("WALD_WORKERS",
                                  min(28, parallel::detectCores() - 1)))
plan(multicore, workers = nworkers)

theme_set(theme_bw())
zmargin <- theme(panel.spacing = grid::unit(0, "pt"))
## Okabe-Ito minus black and yellow
oi_cols <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
scale_colour_discrete <- function(...) scale_colour_manual(..., values = oi_cols)

source("wald_lrt_funs.R")

p0 <- 0.01
nvec <- c(20, 50, 100)

## estimates and CIs for all outcomes
cis_file <- "wald_lrt_exact_cis.rds"
if (file.exists(cis_file)) {
  cis_2g <- readRDS(cis_file)
} else {
  cis_2g <- map(nvec, exact_cis_n)
  saveRDS(cis_2g, cis_file)
}

cov_2g <- map_dfr(cis_2g, exact_cov_n, p0 = p0)
cis_2g <- bind_rows(cis_2g)

## the unpenalized estimate is NaN when both groups have zero successes;
## wmedian() drops these outcomes, which is harmless for p0 = 0.2 but
## not necessarily here. Probability of y0 = y1 = 0 at each beta:
p_both0 <- expand_grid(n = nvec, beta = 3:6) |>
  mutate(prob = dbinom(0, n, p0) * dbinom(0, n, plogis(qlogis(p0) + beta)))
print(p_both0, n = Inf)

## median estimates and CI limits
med_2g <- map_dfr(3:6, exact_medians, cis = cis_2g, p0 = p0)
p_med <- plot_medians(med_2g)
ggsave("wald_lrt_lowp_medians.png", p_med, width = 9, height = 5)

## coverage (no direct penalty labels, so a colour legend instead)
p_cov <- plot_coverage(cov_2g)
ggsave("wald_lrt_lowp_coverage.png", p_cov, width = 10, height = 7)
