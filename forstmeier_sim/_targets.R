## targets pipeline for the Forstmeier & Schielzeth simulation reproduction.
##
## The simulation/plotting code is written as plain scripts rather than
## pure functions, so this pipeline is built from `format = "file"`
## targets: each script and each artifact it produces is tracked as a
## file, and targets re-sources a script whenever any file it depends on
## (declared by referencing the upstream target name) changes. This is the
## "Makefile in R" pattern.
##
## output/corrections_results.rds is the single simulation output that
## every figure reads from: it screens all three model-simplification
## scenarios (unselected, unrestricted step(), interactions-only step())
## under all four significance criteria (uncorrected, Dunn-Sidak, Holm,
## single-step max-|T|). Shared logic lives in a handful of "core" files
## sourced by multiple thin wrapper scripts:
##   simulate_fig1.R          -- data generation, the two swappable
##                               step()-simplification strategies
##   simulate_corrections.R   -- one_rep_corrections()/run_condition_corrections()
##   maxT_correction.R        -- single-step max-|T| critical value
##   run_corrections_core.R   -- run_corrections_simulation(): the simulate-and-cache logic
##   plot_fig1_core.R         -- make_fig1_plot(): parametrized 2-panel a/b reproduction
##   plot_stepwise_comparison_core.R -- n-row/n-condition scenario comparisons
##   plot_corrections_by_scenario_core.R -- patchwork scenario x correction-method figure

library(targets)

list(
  tar_target(sim_functions_file, here::here("forstmeier_sim", "R", "simulate_fig1.R"),
             format = "file"),

  tar_target(maxT_functions_file, here::here("forstmeier_sim", "R", "maxT_correction.R"),
             format = "file"),

  tar_target(corrections_functions_file, here::here("forstmeier_sim", "R", "simulate_corrections.R"),
             format = "file"),

  tar_target(run_corrections_core_file, here::here("forstmeier_sim", "R", "run_corrections_core.R"),
             format = "file"),

  tar_target(graphics_utils_file, here::here("forstmeier_sim", "R", "graphics_utils.R"),
             format = "file"),

  tar_target(run_corrections_script, here::here("forstmeier_sim", "R", "run_corrections_simulation.R"),
             format = "file"),

  tar_target(corrections_results_file, {
    sim_functions_file          ## dependency only: invalidates this target if edited
    maxT_functions_file         ## dependency only
    corrections_functions_file  ## dependency only
    run_corrections_core_file   ## dependency only
    source(run_corrections_script)
    here::here("forstmeier_sim", "output", "corrections_results.rds")
  }, format = "file"),

  ## -- Fig. 1 reproduction: 2-panel a/b layout, one figure per --
  ## -- model-simplification variant --

  tar_target(plot_fig1_core_file, here::here("forstmeier_sim", "R", "plot_fig1_core.R"),
             format = "file"),

  tar_target(make_plot_script, here::here("forstmeier_sim", "R", "make_plot.R"),
             format = "file"),

  tar_target(fig1_png_file, {
    corrections_results_file  ## dependency only
    plot_fig1_core_file         ## dependency only
    graphics_utils_file          ## dependency only
    source(make_plot_script)
    here::here("forstmeier_sim", "output", "fig1.png")
  }, format = "file"),

  tar_target(make_plot_restricted_script,
             here::here("forstmeier_sim", "R", "make_plot_restricted.R"),
             format = "file"),

  tar_target(fig1_restricted_png_file, {
    corrections_results_file  ## dependency only
    plot_fig1_core_file         ## dependency only
    graphics_utils_file          ## dependency only
    source(make_plot_restricted_script)
    here::here("forstmeier_sim", "output", "fig1_restricted.png")
  }, format = "file"),

  ## -- multiple-comparisons-correction comparison, unselected model only --

  tar_target(make_corrections_plot_script, here::here("forstmeier_sim", "R", "make_corrections_plot.R"),
             format = "file"),

  tar_target(fig1_corrections_png_file, {
    corrections_results_file  ## dependency only
    graphics_utils_file        ## dependency only
    source(make_corrections_plot_script)
    here::here("forstmeier_sim", "output", "fig1_corrections.png")
  }, format = "file"),

  ## -- stepwise-simplification-method comparison (unselected / --
  ## -- unrestricted step() / interactions-only step()), no    --
  ## -- multiple-comparisons correction --

  tar_target(plot_stepwise_comparison_core_file,
             here::here("forstmeier_sim", "R", "plot_stepwise_comparison_core.R"),
             format = "file"),

  tar_target(make_stepwise_comparison_plot_script,
             here::here("forstmeier_sim", "R", "make_stepwise_comparison_plot.R"),
             format = "file"),

  tar_target(fig1_stepwise_comparison_png_file, {
    corrections_results_file            ## dependency only
    plot_stepwise_comparison_core_file  ## dependency only
    graphics_utils_file                  ## dependency only
    source(make_stepwise_comparison_plot_script)
    here::here("forstmeier_sim", "output", "fig1_stepwise_comparison.png")
  }, format = "file"),

  tar_target(make_stepwise_comparison_plot_by_condition_script,
             here::here("forstmeier_sim", "R", "make_stepwise_comparison_plot_by_condition.R"),
             format = "file"),

  tar_target(fig1_stepwise_comparison_by_condition_png_file, {
    corrections_results_file            ## dependency only
    plot_stepwise_comparison_core_file  ## dependency only
    graphics_utils_file                  ## dependency only
    source(make_stepwise_comparison_plot_by_condition_script)
    here::here("forstmeier_sim", "output", "fig1_stepwise_comparison_by_condition.png")
  }, format = "file"),

  ## -- scenario x correction-method comparison, combined via patchwork --

  tar_target(plot_corrections_by_scenario_core_file,
             here::here("forstmeier_sim", "R", "plot_corrections_by_scenario_core.R"),
             format = "file"),

  tar_target(make_corrections_by_scenario_plot_script,
             here::here("forstmeier_sim", "R", "make_corrections_by_scenario_plot.R"),
             format = "file"),

  tar_target(fig1_corrections_by_scenario_png_file, {
    corrections_results_file                ## dependency only
    plot_corrections_by_scenario_core_file  ## dependency only
    graphics_utils_file                      ## dependency only
    source(make_corrections_by_scenario_plot_script)
    here::here("forstmeier_sim", "output", "fig1_corrections_by_scenario.png")
  }, format = "file"),

  tar_target(make_corrections_by_method_plot_script,
             here::here("forstmeier_sim", "R", "make_corrections_by_method_plot.R"),
             format = "file"),

  tar_target(fig1_corrections_by_method_png_file, {
    corrections_results_file                ## dependency only
    plot_corrections_by_scenario_core_file  ## dependency only
    graphics_utils_file                      ## dependency only
    source(make_corrections_by_method_plot_script)
    here::here("forstmeier_sim", "output", "fig1_corrections_by_method.png")
  }, format = "file"),

  tar_target(make_corrections_by_N_plot_script,
             here::here("forstmeier_sim", "R", "make_corrections_by_N_plot.R"),
             format = "file"),

  tar_target(fig1_corrections_by_N_png_file, {
    corrections_results_file                ## dependency only
    plot_corrections_by_scenario_core_file  ## dependency only
    graphics_utils_file                      ## dependency only
    source(make_corrections_by_N_plot_script)
    here::here("forstmeier_sim", "output", "fig1_corrections_by_N.png")
  }, format = "file")
)
