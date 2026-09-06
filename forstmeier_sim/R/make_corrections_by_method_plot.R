## Reads output/corrections_results.rds and writes
## output/fig1_corrections_by_method.png: all four significance criteria
## (uncorrected, Dunn-Sidak, Holm, single-step max-|T|) as colour/shape,
## rows = model-simplification scenario, columns = N -- with-interactions
## cases only, on one shared y-scale (a single facet_grid, no patchwork).

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

source(here::here("forstmeier_sim", "R", "plot_corrections_by_scenario_core.R"))

make_corrections_by_method_plot(
  here::here("forstmeier_sim", "output", "corrections_results.rds"),
  here::here("forstmeier_sim", "output", "fig1_corrections_by_method.png")
)
