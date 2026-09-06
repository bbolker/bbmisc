## Reads output/corrections_results.rds and writes
## output/fig1_corrections_by_scenario.png: uncorrected type I error rate
## (1 column x 3 rows) next to Dunn-Sidak/Holm/max-|T| (3 columns x 3 rows),
## rows = model-simplification scenario, combined side by side with
## patchwork. See plot_corrections_by_scenario_core.R for why this is two
## ggplot objects rather than one facet_grid.

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

source(here::here("forstmeier_sim", "R", "plot_corrections_by_scenario_core.R"))

make_corrections_by_scenario_plot(
  here::here("forstmeier_sim", "output", "corrections_results.rds"),
  here::here("forstmeier_sim", "output", "fig1_corrections_by_scenario.png")
)
