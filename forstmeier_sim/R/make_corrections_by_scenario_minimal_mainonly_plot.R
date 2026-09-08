## Reads output/corrections_results.rds and writes
## output/fig1_corrections_by_scenario_minimal_mainonly.png: like
## fig1_corrections_by_scenario_minimal.png, but using the main-effects-only
## error rates (with-interactions rows only). See
## plot_corrections_by_scenario_core.R.

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

source(here::here("forstmeier_sim", "R", "plot_corrections_by_scenario_core.R"))

make_corrections_by_scenario_minimal_mainonly_plot(
  here::here("forstmeier_sim", "output", "corrections_results.rds"),
  here::here("forstmeier_sim", "output", "fig1_corrections_by_scenario_minimal_mainonly.png")
)
