## Reads output/corrections_results.rds and writes
## output/fig1_corrections_by_N.png: rows = N, columns = criterion
## (uncorrected, Dunn-Sidak, Holm, single-step max-|T|), colour/shape =
## model-simplification scenario (unselected/selected/selected-
## interactions-only) -- with-interactions cases only.

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

source(here::here("forstmeier_sim", "R", "plot_corrections_by_scenario_core.R"))

make_corrections_by_N_plot(
  here::here("forstmeier_sim", "output", "corrections_results.rds"),
  here::here("forstmeier_sim", "output", "fig1_corrections_by_N.png")
)
