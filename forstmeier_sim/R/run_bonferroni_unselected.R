## Writes output/bonferroni_unselected.rds: true-Bonferroni (p < alpha/k)
## experimentwise indicators and per-replicate significant-term counts for
## the unselected full model, over the same
## (N, m, interactions) grid as corrections_results.rds. A stopgap for
## fig1_corrections.png until corrections_results.rds is regenerated with
## its own sig_bonferroni_* columns (see TODO.md).

library(here)
library(furrr)
library(future)
library(dplyr)
library(tidyr)
library(purrr)

source(here::here("forstmeier_sim", "R", "simulate_fig1.R"))
source(here::here("forstmeier_sim", "R", "simulate_corrections.R"))

set.seed(20261111)

RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

nrep <- 5000
n_workers <- min(24, max(1, parallel::detectCores() - 1, na.rm = TRUE))
future::plan(future::multisession, workers = n_workers)

results <- tidyr::expand_grid(
  N = sim_N_values,
  m = 1:6,
  interactions = c(FALSE, TRUE)
) |>
  dplyr::mutate(
    k = purrr::map2_dbl(m, interactions, n_predictors),
    sims = purrr::pmap(list(N, m, interactions), function(N, m, interactions) {
      run_condition_bonferroni_unselected(N, m, interactions, nrep)
    })
  ) |>
  tidyr::unnest(sims)

future::plan(future::sequential)

saveRDS(results, here::here("forstmeier_sim", "output", "bonferroni_unselected.rds"))
