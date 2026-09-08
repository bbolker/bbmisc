## Reads output/corrections_results.rds and writes
## output/fig1_corrections_by_scenario_minimal.png: single-step max-|T| and
## per-parameter error rates only (rows = model-simplification scenario,
## columns = criterion), no patchwork. See plot_corrections_by_scenario_core.R.

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

source(here::here("forstmeier_sim", "R", "plot_corrections_by_scenario_core.R"))

make_corrections_by_scenario_minimal_plot(
  here::here("forstmeier_sim", "output", "corrections_results.rds"),
  here::here("forstmeier_sim", "output", "fig1_corrections_by_scenario_minimal.png")
)
