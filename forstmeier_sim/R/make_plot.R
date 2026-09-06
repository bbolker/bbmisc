## Reads output/corrections_results.rds and writes output/fig1.png
## (reproduction of Fig. 1 in Forstmeier & Schielzeth 2011, unrestricted
## step() variant). See plot_fig1_core.R for the shared plot-building logic.

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

source(here::here("forstmeier_sim", "R", "plot_fig1_core.R"))

make_fig1_plot(
  here::here("forstmeier_sim", "output", "corrections_results.rds"),
  col_unselected = "sig_raw_unselected",
  col_selected = "sig_raw_step",
  here::here("forstmeier_sim", "output", "fig1.png")
)
