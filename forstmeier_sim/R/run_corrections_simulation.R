## Runs the multiple-comparisons-correction simulation grid and writes
## output/corrections_results.rds. See run_corrections_core.R for the
## shared simulation logic.
##
## With the environment variable FORSTMEIER_TOUCH_SIM set (`make
## touch-sim`), keeps the existing results file instead, so that targets
## records it as up to date -- the equivalent of `touch`, for code changes
## known not to affect the results. Use with care: it is exactly how stale
## results can survive real changes to the simulation.
##
## The environment variable FORSTMEIER_WORKERS, if set, overrides the
## default worker count (still capped at one less than the core count).

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
source(here::here("forstmeier_sim", "R", "run_corrections_core.R"))

results_path <- here::here("forstmeier_sim", "output", "corrections_results.rds")

if (nzchar(Sys.getenv("FORSTMEIER_TOUCH_SIM"))) {
  if (!file.exists(results_path)) stop("FORSTMEIER_TOUCH_SIM is set but ", results_path, " doesn't exist")
  message("FORSTMEIER_TOUCH_SIM set: keeping existing ", basename(results_path), " without rerunning")
} else {
  env_workers <- Sys.getenv("FORSTMEIER_WORKERS")
  if (nzchar(env_workers)) {
    n_workers <- min(as.integer(env_workers),
                     max(1, parallel::detectCores() - 1, na.rm = TRUE))
    run_corrections_simulation(results_path, n_workers = n_workers)
  } else {
    run_corrections_simulation(results_path)
  }
}
