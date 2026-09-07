## Smoke test / timing check for the corrections simulation grid. Verifies
## correctness on a couple of extreme conditions, then times the grid at a
## small nrep, scaled up, across worker counts.

library(here)
library(furrr)
library(future)
library(broom)
library(dplyr)
library(tidyr)
library(purrr)
library(mvtnorm)

source(here::here("forstmeier_sim", "R", "simulate_fig1.R"))
source(here::here("forstmeier_sim", "R", "maxT_correction.R"))
source(here::here("forstmeier_sim", "R", "simulate_corrections.R"))

RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

simplify_fns <- list(step = step_default, step_restricted = step_interactions_only)

## ---- correctness checks ----
set.seed(1)
print(one_rep_corrections(200, 1, FALSE, simplify_fns = simplify_fns))
set.seed(1)
print(one_rep_corrections(50, 6, TRUE, simplify_fns = simplify_fns))

## ---- timing scan across worker counts ----
conditions <- expand_grid(
  N = sim_N_values,
  m = 1:6,
  interactions = c(FALSE, TRUE)
)
nrep_smoke <- 150
full_nrep <- 5000

for (n_workers in c(8, 16, 24, 32)) {
  plan(multisession, workers = n_workers)
  ## warm up worker pool so startup cost doesn't distort the timing
  invisible(furrr::future_map(1:n_workers, ~ Sys.getpid(),
                               .options = furrr::furrr_options(seed = TRUE)))

  set.seed(20260906)
  t <- system.time({
    results <- conditions |>
      mutate(
        sims = pmap(list(N, m, interactions), function(N, m, interactions) {
          run_condition_corrections(N, m, interactions, nrep_smoke, simplify_fns = simplify_fns)
        })
      ) |>
      unnest(sims)
  })
  plan(sequential)

  est_sec <- t[["elapsed"]] * (full_nrep / nrep_smoke)
  cat(sprintf(
    "workers=%2d  elapsed=%.2fs for %d reps  -> est. full run (nrep=%d): %.1f min\n",
    n_workers, t[["elapsed"]], nrow(conditions) * nrep_smoke, full_nrep, est_sec / 60
  ))
}
