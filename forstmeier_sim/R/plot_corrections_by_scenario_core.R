## Combines corrections_results.rds's three model-simplification scenarios
## (unselected, unrestricted step(), interactions-only step()) with its four
## significance criteria (uncorrected, Dunn-Sidak, Holm, single-step
## max-|T|) into one figure: a 1-column x n-row uncorrected-only panel next
## to an n-row x 3-column corrected-methods facet_grid, combined side by
## side with patchwork. Built as two separate ggplot objects rather than
## one facet_grid, since the uncorrected panel needs the theoretical
## alpha' reference curve while the corrected panels need the flat
## alpha=0.05 line -- mixing those into one facet_grid call would need
## per-column conditional layers.

library(ggplot2)
library(scales)
library(patchwork)

source(here::here("forstmeier_sim", "R", "graphics_utils.R"))

scenario_labels <- c(
  unselected = "Unselected",
  step = "Selected\n(all)",
  step_restricted = "Selected\n(interax only)"
)

criterion_labels <- c(
  raw = "Uncorrected",
  bonf = "Dunn–Šidák",
  holm = "Holm",
  maxT = "Single-step max|T|"
)

## Long-format summary: one row per (scenario, criterion, multcomp_k, N, m,
## interactions) with the proportion and its Gaussian-approximation 95% CI.
## multcomp_k distinguishes which model's k a correction was calibrated
## on: "minimal" (the selected model's own surviving term count -- the
## naive approach), "maximal" (the pre-selection full model's k -- the
## correct approach), or "none" where the distinction doesn't apply (the
## uncorrected criterion, or the unselected scenario, which has only one
## calibration since fit and ref coincide). Deliberately not left as NA:
## interaction(scenario, multcomp_k) -- used for grouping in
## make_corrections_by_N_plot() -- returns NA for *every* row when
## multcomp_k is NA, regardless of scenario, which collapses all scenarios
## into one undifferentiated group and produces stray lines connecting
## points across scenarios.
build_corrections_by_scenario_summary <- function(results_path) {
  readRDS(results_path) |>
    tidyr::pivot_longer(
      cols = dplyr::starts_with("sig_"),
      names_to = c("criterion", "scenario", "multcomp_k"),
      names_pattern = "^sig_(raw|bonf|holm|maxT)_(unselected|step_restricted|step)(?:_(naive|correct))?$",
      values_to = "sig"
    ) |>
    dplyr::mutate(multcomp_k = dplyr::recode_values(multcomp_k,
      from = c("naive", "correct"), to = c("minimal", "maximal"), default = "none")) |>
    dplyr::group_by(scenario, criterion, multcomp_k, N, m, interactions, k) |>
    dplyr::summarise(n = dplyr::n(), prop = mean(sig), .groups = "drop") |>
    dplyr::mutate(
      scenario = factor(scenario_labels[scenario], levels = scenario_labels),
      criterion = factor(criterion_labels[criterion], levels = criterion_labels),
      se = sqrt(prop * (1 - prop) / n),
      ci_lo = pmax(0, prop - 1.96 * se),
      ci_hi = pmin(1, prop + 1.96 * se),
      group = interaction(
        ifelse(interactions, "With 2-way interactions", "Without interactions"),
        paste0("N = ", N),
        sep = ", "
      )
    )
}

## Zero spacing between facet panels (project graphics preference), used
## throughout this file since none of these plots have tick labels at risk
## of overlapping between adjacent panels.
zmargin <- theme(panel.spacing = grid::unit(0, "pt"))

## Colour = N, linetype/shape = multcomp_k -- with-interactions cases only.
## Unlike make_corrections_by_method_plot(), this keeps both the minimal and
## maximal calibrations (distinguished by linetype/shape) rather than
## defaulting to "maximal" via filter_default_calibration().
make_corrections_by_scenario_plot <- function(results_path, png_path) {
  summ <- build_corrections_by_scenario_summary(results_path) |>
    dplyr::filter(interactions) |>
    dplyr::mutate(multcomp_k = factor(multcomp_k, levels = c("minimal", "maximal", "none")))
  theme_set(theme_bw())

  ## -- left: uncorrected only, 1 column x n rows --
  left_data <- summ |> dplyr::filter(criterion == criterion_labels[["raw"]])

  grey_lines <- left_data |>
    dplyr::distinct(scenario, m, k) |>
    dplyr::mutate(alpha_prime = 1 - (1 - 0.05)^k)

  p_left <- ggplot(left_data, aes(m, prop, colour = factor(N), linetype = multcomp_k, shape = multcomp_k)) +
    ## ymin = -Inf rather than 0: 0 is not representable on a logit scale
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.05,
             fill = "grey80", alpha = 0.5) +
    geom_line(data = grey_lines, aes(x = m, y = alpha_prime), inherit.aes = FALSE,
              colour = "grey50", linewidth = 2, alpha = 0.5) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = factor(N), group = N),
                colour = NA, alpha = 0.3) +
    geom_line() +
    geom_point() +
    facet_grid(rows = vars(scenario), switch = "y") +
    scale_colour_manual(values = okabe_ito, name = "N") +
    scale_fill_manual(values = okabe_ito, guide = "none") +
    scale_linetype_manual(values = c(minimal = "22", maximal = "solid", none = "solid"), name = "multcomp k") +
    scale_shape_manual(values = c(minimal = 17, maximal = 16, none = 1), name = "multcomp k") +
    scale_x_continuous(name = "Number of explanatory variables", breaks = 1:6) +
    scale_y_continuous(trans = "logit", breaks = logit_breaks()) +
    ylab("Proportion of models with type I errors") +
    zmargin +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0))

  ## -- right: corrected methods, n rows x 3 columns --
  right_data <- summ |> dplyr::filter(criterion != criterion_labels[["raw"]])

  nominal_ref <- right_data |> dplyr::distinct(scenario, criterion)

  p_right <- ggplot(right_data, aes(m, prop, colour = factor(N), linetype = multcomp_k, shape = multcomp_k)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.05,
             fill = "grey80", alpha = 0.5) +
    geom_hline(data = nominal_ref, aes(yintercept = 0.05), inherit.aes = FALSE,
               colour = "grey40", linetype = "dashed", linewidth = 1) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = factor(N), group = interaction(N, multcomp_k)),
                colour = NA, alpha = 0.3) +
    geom_line() +
    geom_point() +
    facet_grid(rows = vars(scenario), cols = vars(criterion)) +
    scale_colour_manual(values = okabe_ito, name = "N") +
    scale_fill_manual(values = okabe_ito, guide = "none") +
    scale_linetype_manual(values = c(minimal = "22", maximal = "solid", none = "solid"), name = "multcomp k") +
    scale_shape_manual(values = c(minimal = 17, maximal = 16, none = 1), name = "multcomp k") +
    scale_x_continuous(name = "Number of explanatory variables", breaks = 1:6) +
    scale_y_continuous(trans = "logit", breaks = logit_breaks()) +
    ylab(NULL) +
    zmargin +
    theme(strip.text.y = element_blank())

  fig <- (p_left + p_right) +
    plot_layout(widths = c(1, 3), guides = "collect") &
    theme(legend.position = "right")

  ggsave(png_path, fig, width = 15, height = 3.2 * length(scenario_labels), dpi = 150)
}

## All four criteria (including uncorrected) as colour, both multcomp_k
## calibrations as linetype/shape, so a single facet_grid still suffices (no
## patchwork): rows = scenario (matching the row layout of the other
## figures), columns = N (with-interactions cases only). Mixing the
## uncorrected and corrected criteria onto one shared y-scale means the
## "Selected" rows' corrected lines get visually compressed near 0.05 next
## to the much larger uncorrected values -- kept simple here rather than
## reaching for free_y or a secondary axis.
make_corrections_by_method_plot <- function(results_path, png_path) {
  summ <- build_corrections_by_scenario_summary(results_path) |>
    dplyr::filter(interactions) |>
    dplyr::mutate(
      N_label = factor(paste0("N = ", N), levels = paste0("N = ", sort(unique(N)))),
      multcomp_k = factor(multcomp_k, levels = c("minimal", "maximal", "none"))
    )

  theme_set(theme_bw())

  fig <- ggplot(summ, aes(m, prop, colour = criterion, linetype = multcomp_k, shape = multcomp_k)) +
    geom_hline(yintercept = 0.05, colour = "grey40", linetype = "dashed", linewidth = 1) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = criterion, group = interaction(criterion, multcomp_k)),
                colour = NA, alpha = 0.3) +
    geom_line() +
    geom_point() +
    facet_grid(rows = vars(scenario), cols = vars(N_label)) +
    scale_colour_manual(values = okabe_ito, name = NULL) +
    scale_fill_manual(values = okabe_ito, guide = "none") +
    scale_linetype_manual(values = c(minimal = "22", maximal = "solid", none = "solid"), name = "multcomp k") +
    scale_shape_manual(values = c(minimal = 17, maximal = 16, none = 1), name = "multcomp k") +
    scale_x_continuous(name = "Number of explanatory variables", breaks = 1:6) +
    scale_y_continuous(trans = "logit", breaks = logit_breaks()) +
    ylab("Proportion of models with type I errors") +
    zmargin

  ggsave(png_path, fig, width = 10, height = 8, dpi = 150)
}

## Rows = N, columns = criterion (uncorrected first, then the three
## corrections), colour = model-simplification scenario, linetype/shape =
## multcomp_k ("minimal"/"maximal" calibration, or "none" where the
## distinction doesn't apply -- uncorrected, or the unselected scenario) --
## with-interactions cases only. Unlike the other two plots in this file,
## this one keeps both calibrations rather than defaulting to "maximal".
make_corrections_by_N_plot <- function(results_path, png_path) {
  summ <- build_corrections_by_scenario_summary(results_path) |>
    dplyr::filter(interactions) |>
    dplyr::mutate(
      N_label = factor(paste0("N = ", N), levels = paste0("N = ", sort(unique(N)))),
      multcomp_k = factor(multcomp_k, levels = c("minimal", "maximal", "none"))
    )

  theme_set(theme_bw())

  fig <- ggplot(summ, aes(m, prop, colour = scenario, linetype = multcomp_k, shape = multcomp_k)) +
    geom_hline(yintercept = 0.05, colour = "grey40", linetype = "dashed", linewidth = 1) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = scenario, group = interaction(scenario, multcomp_k)),
                colour = NA, alpha = 0.3) +
    geom_line() +
    geom_point() +
    facet_grid(rows = vars(N_label), cols = vars(criterion)) +
    scale_colour_manual(values = okabe_ito, name = NULL) +
    scale_fill_manual(values = okabe_ito, guide = "none") +
    scale_linetype_manual(values = c(maximal = "solid", minimal = "22", none = "solid"), name = "multcomp k") +
    scale_shape_manual(values = c(maximal = 16, minimal = 17, none = 1), name = "multcomp k") +
    scale_x_continuous(name = "Number of explanatory variables", breaks = 1:6) +
    scale_y_continuous(trans = "logit", breaks = logit_breaks()) +
    ylab("Proportion of models with type I errors") +
    zmargin +
    theme(strip.text.y = element_text(angle = 0))

  ggsave(png_path, fig, width = 10, height = 8, dpi = 150)
}
