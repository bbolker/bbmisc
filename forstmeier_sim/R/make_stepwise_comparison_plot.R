## Reads output/corrections_results.rds and writes
## output/fig1_stepwise_comparison.png: a facet_grid comparison across
## model-simplification methods (rows), excluding the multiple-comparisons-
## corrected cases (those are compared separately in fig1_corrections.png
## and fig1_corrections_by_scenario.png).

library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

source(here::here("forstmeier_sim", "R", "plot_stepwise_comparison_core.R"))

results_path <- here::here("forstmeier_sim", "output", "corrections_results.rds")

specs <- list(
  "(a) Without model simplification" = list(path = results_path, column = "sig_raw_unselected"),
  "(b) With model simplification (unrestricted step())" = list(path = results_path, column = "sig_raw_step"),
  "(c) With model simplification (interactions only)" = list(path = results_path, column = "sig_raw_step_restricted")
)

make_stepwise_comparison_plot(
  specs,
  here::here("forstmeier_sim", "output", "fig1_stepwise_comparison.png")
)
