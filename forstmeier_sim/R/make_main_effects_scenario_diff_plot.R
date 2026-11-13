## Reads output/corrections_results.rds and writes
## output/fig1_mainonly_scenario_diff.png: paired differences in the
## main-effects-only max-|T| and per-parameter type I error rates between
## the two selection scenarios (all terms vs. interactions only).

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

source(here::here("forstmeier_sim", "R", "plot_corrections_by_scenario_core.R"))

make_main_effects_scenario_diff_plot(
  here::here("forstmeier_sim", "output", "corrections_results.rds"),
  here::here("forstmeier_sim", "output", "fig1_mainonly_scenario_diff.png")
)
