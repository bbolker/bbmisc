## Runs the multiple-comparisons-correction simulation grid and writes
## output/corrections_results.rds. See run_corrections_core.R for the
## shared simulate-and-cache logic.

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

run_corrections_simulation(
  here::here("forstmeier_sim", "output", "corrections_results.rds")
)
