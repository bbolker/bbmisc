## Same data as make_stepwise_comparison_plot.R, but with facet/aesthetic
## roles swapped: facets are (N, interactions) and colour is the
## simplification method, so the three methods can be compared directly
## against each other within a fixed condition. Writes
## output/fig1_stepwise_comparison_by_condition.png.

library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

source(here::here("forstmeier_sim", "R", "plot_stepwise_comparison_core.R"))

results_path <- here::here("forstmeier_sim", "output", "corrections_results.rds")

specs <- list(
  "Unselected" = list(path = results_path, column = "sig_raw_unselected"),
  "Selected (unrestricted step())" = list(path = results_path, column = "sig_raw_step"),
  "Selected (interactions only)" = list(path = results_path, column = "sig_raw_step_restricted")
)

make_stepwise_comparison_plot_by_condition(
  specs,
  here::here("forstmeier_sim", "output", "fig1_stepwise_comparison_by_condition.png")
)
