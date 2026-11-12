## Combines corrections_results.rds's three model-simplification scenarios
## (unselected, unrestricted step(), interactions-only step()) with its five
## per-model summary criteria (uncorrected, Dunn-Sidak, Holm, single-step
## max-|T|, and the uncorrected per-parameter rate) into one figure: a
## 1-column x n-row uncorrected-only panel next to an n-row x 4-column
## "other criteria" facet_grid, combined side by side with patchwork. Built
## as two separate ggplot objects rather than one facet_grid, since the
## uncorrected panel needs the theoretical alpha' reference curve while the
## other panels need the flat alpha=0.05 line -- mixing those into one
## facet_grid call would need per-column conditional layers.

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
  dunnsidak = "Dunn–Šidák",
  holm = "Holm",
  maxT = "Single-step max|T|",
  perparam = "Per-parameter"
)

## Builds the long-format summary shared by every by-scenario/method/N
## plot: one row per (scenario, criterion, multcomp_k, N, m, interactions)
## with the proportion and its Gaussian-approximation 95% CI. `variant`
## selects which family of sig_* columns to pivot -- "" (the default) for
## the all-effects metrics used throughout this project's main figures, or
## "main_" for the main-effects-only metrics simulate_corrections.R also
## writes (sig_*_main_*, with-interactions rows only); `denom_col` ("k" or
## "m") is what the derived per-parameter rate divides
## n_sig_raw_<variant><scenario> by. See build_corrections_by_scenario_summary()/
## build_main_effects_summary() below for the two variants actually used.
##
## The column-matching regex is used for both selecting which columns to
## pivot and parsing their names, so a column outside the requested variant
## (e.g. a main-effects column when variant = "", or the differently-shaped
## sig_raw_mainonly_crosscheck column always) is simply excluded from the
## pivot rather than silently producing an NA criterion/scenario that later
## collapses into a meaningless combined group (tidyr::pivot_longer()
## doesn't error on a names_pattern that fails to match a selected column).
##
## multcomp_k distinguishes which model's k a correction was calibrated
## on: "minimal" (the selected model's own surviving term count -- the
## naive approach), "maximal" (the pre-selection full model's k -- the
## correct approach), or "none" where the distinction doesn't apply (the
## uncorrected and per-parameter criteria, or the unselected scenario, which
## has only one calibration since fit and ref coincide). Deliberately not
## left as NA: interaction(scenario, multcomp_k) -- used for grouping in
## make_corrections_by_N_plot() -- returns NA for *every* row when
## multcomp_k is NA, regardless of scenario, which collapses all scenarios
## into one undifferentiated group and produces stray lines connecting
## points across scenarios.
##
## The "perparam" criterion is the per-parameter error rate E[V]/denom
## (V = how many of the relevant term set are significant, uncorrected --
## see n_sig_raw_* in simulate_corrections.R), as opposed to the "raw"
## criterion's experimentwise indicator (is at least one significant). Its
## per-replicate "sig" value is already a rate in [0, 1] rather than a 0/1
## indicator, so grouped mean(sig) still gives E[V]/denom, but the shared
## Gaussian CI below (derived for a mean of 0/1 indicators) is only an
## approximation to its true sampling variance.
build_corrections_summary <- function(results_path, variant = "", denom_col = "k") {
  col_pattern <- paste0(
    "^sig_(raw|dunnsidak|holm|maxT|perparam)_", variant,
    "(unselected|step_restricted|step)(?:_(naive|correct))?$"
  )
  scenario_names <- names(scenario_labels)

  raw <- readRDS(results_path)
  perparam <- stats::setNames(
    lapply(scenario_names, function(nm) raw[[paste0("n_sig_raw_", variant, nm)]] / raw[[denom_col]]),
    paste0("sig_perparam_", variant, scenario_names)
  )

  dplyr::bind_cols(raw, perparam) |>
    tidyr::pivot_longer(
      cols = dplyr::matches(col_pattern),
      names_to = c("criterion", "scenario", "multcomp_k"),
      names_pattern = col_pattern,
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

## All-effects metrics (the project's main figures) -- see
## build_corrections_summary().
build_corrections_by_scenario_summary <- function(results_path) {
  build_corrections_summary(results_path, variant = "", denom_col = "k")
}

## Main-effects-only metrics, with-interactions rows only (see
## simulate_corrections.R) -- see build_corrections_summary().
build_main_effects_summary <- function(results_path) {
  build_corrections_summary(results_path, variant = "main_", denom_col = "m")
}

## Zero spacing between facet panels (project graphics preference), used
## throughout this file since none of these plots have tick labels at risk
## of overlapping between adjacent panels.
zmargin <- theme(panel.spacing = grid::unit(0, "pt"))

## Facet layout for the by-scenario plots: rows = scenario and columns =
## criterion by default, or the transpose with flip = TRUE. Row strips go
## on the left, outside the axis.
scenario_facets <- function(flip) {
  if (flip) {
    facet_grid(rows = vars(criterion), cols = vars(scenario), switch = "y")
  } else {
    facet_grid(rows = vars(scenario), cols = vars(criterion), switch = "y")
  }
}

## Theme for the scenario strips: horizontal text when they are row strips
## (flip = FALSE), and blanked when show = FALSE because a neighbouring
## patchwork panel already labels them.
scenario_strip_theme <- function(flip, show) {
  if (flip) {
    if (show) theme() else theme(strip.text.x = element_blank())
  } else {
    if (show) {
      theme(strip.text.y.left = element_text(angle = 0))
    } else {
      theme(strip.text.y.left = element_blank())
    }
  }
}

## Scenario x criterion grid, colour = N, linetype/shape = multcomp_k --
## the "everything but uncorrected" facet_grid shared by
## make_corrections_by_scenario_plot() (as its second patchwork panel,
## scenario-strip text hidden since the uncorrected panel already labels
## scenarios, criteria = dunnsidak/holm/maxT/perparam) and
## make_minimal_scenario_plot() (standalone, scenario strips shown,
## criteria = maxT/perparam only).
build_other_criteria_plot <- function(data, show_scenario_strip, ylab_text, flip = FALSE) {
  fig <- ggplot(data, aes(m, prop, colour = factor(N), linetype = multcomp_k, shape = multcomp_k)) +
    ## No `data =`: yintercept is a fixed parameter rather than an
    ## aesthetic mapping, so geom_hline() draws in every panel the facet
    ## grid renders, even one with no data of its own.
    geom_hline(yintercept = 0.05, colour = "grey40", linetype = "dashed", linewidth = 1) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = factor(N), group = interaction(N, multcomp_k)),
                colour = NA, alpha = 0.3) +
    geom_line() +
    geom_point() +
    ## switch = "y" + strip.placement = "outside" regardless of
    ## show_scenario_strip, even when strip text ends up blanked: when
    ## this plot is patchworked next to another with row strips on the
    ## outside (make_corrections_by_scenario_plot()'s uncorrected panel),
    ## structuring both plots' strips the same way (just blanking text on
    ## one side) avoids a stray gap patchwork's cross-plot alignment
    ## otherwise inserts between the bottom axis ticks and panel border.
    scenario_facets(flip) +
    scale_colour_manual(values = okabe_ito, name = "N") +
    scale_fill_manual(values = okabe_ito, guide = "none") +
    scale_linetype_manual(values = c(minimal = "22", maximal = "solid", none = "solid"), name = "multcomp k") +
    scale_shape_manual(values = c(minimal = 17, maximal = 16, none = 1), name = "multcomp k") +
    m_k_scale_x() +
    scale_y_logit(c(data$ci_lo, data$ci_hi)) +
    ylab(ylab_text) +
    zmargin +
    theme(strip.placement = "outside") +
    scenario_strip_theme(flip, show_scenario_strip)
}

## Colour = N, linetype/shape = multcomp_k -- with-interactions cases only.
## Unlike make_corrections_by_method_plot(), this keeps both the minimal and
## maximal calibrations (distinguished by linetype/shape) rather than
## defaulting to "maximal" via filter_default_calibration().
## flip = FALSE: scenarios as rows, the uncorrected panel as a left-hand
## column. flip = TRUE: scenarios as columns, the uncorrected panel as a
## top row above the other criteria.
make_corrections_by_scenario_plot <- function(results_path, png_path, flip = FALSE) {
  summ <- build_corrections_by_scenario_summary(results_path) |>
    dplyr::filter(interactions) |>
    dplyr::mutate(multcomp_k = factor(multcomp_k, levels = c("minimal", "maximal", "none")))
  theme_set(theme_bw())

  ## -- uncorrected only: 1 column x n rows (or 1 row x n columns if flipped) --
  left_data <- summ |> dplyr::filter(criterion == criterion_labels[["raw"]])

  grey_lines <- left_data |>
    dplyr::distinct(scenario, m, k) |>
    dplyr::mutate(alpha_prime = 1 - (1 - 0.05)^k)

  p_left <- ggplot(left_data, aes(m, prop, colour = factor(N), linetype = multcomp_k, shape = multcomp_k)) +
    geom_line(data = grey_lines, aes(x = m, y = alpha_prime), inherit.aes = FALSE,
              colour = "grey50", linewidth = 2, alpha = 0.5) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = factor(N), group = N),
                colour = NA, alpha = 0.3) +
    geom_line() +
    geom_point() +
    ## the criterion facet has the single level "Uncorrected": it labels the
    ## row when flipped, and is dropped otherwise
    (if (flip) scenario_facets(flip) else facet_grid(rows = vars(scenario), switch = "y")) +
    scale_colour_manual(values = okabe_ito, name = "N") +
    scale_fill_manual(values = okabe_ito, guide = "none") +
    ## guide = "none": left_data's multcomp_k is always "none" (the
    ## uncorrected criterion), so this scale's own legend would just be a
    ## redundant single-key duplicate of p_right's three-key one once
    ## combined via guides = "collect" below.
    scale_linetype_manual(values = c(minimal = "22", maximal = "solid", none = "solid"), guide = "none") +
    scale_shape_manual(values = c(minimal = 17, maximal = 16, none = 1), guide = "none") +
    m_k_scale_x() +
    ## Includes grey_lines$alpha_prime (the theoretical reference curve,
    ## plotted from separate data) alongside the ribbon extent: since
    ## scale_y_logit() sets the scale's limits explicitly, anything plotted
    ## but left out of this range would otherwise get silently clipped.
    scale_y_logit(c(left_data$ci_lo, left_data$ci_hi, grey_lines$alpha_prime)) +
    zmargin +
    theme(strip.placement = "outside") +
    scenario_strip_theme(flip, show = TRUE)

  ## -- other criteria: n rows x 4 columns (or 4 rows x n columns if flipped) --
  right_data <- summ |> dplyr::filter(criterion != criterion_labels[["raw"]])
  p_right <- build_other_criteria_plot(right_data, show_scenario_strip = FALSE,
                                       ylab_text = NULL, flip = flip)

  ylab_text <- "Proportion with type I errors"
  n_scen <- length(scenario_labels)
  n_other <- nlevels(droplevels(right_data$criterion))
  if (flip) {
    ## top row loses its x axis (the panels below share it); the y title
    ## goes on the taller bottom block
    p_left <- p_left + ylab(NULL) +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    p_right <- p_right + ylab(ylab_text)
    fig <- (p_left / p_right) + plot_layout(heights = c(1, n_other), guides = "collect")
    width <- 3.5 * n_scen + 3
    height <- 2.6 * (n_other + 1) + 1
  } else {
    p_left <- p_left + ylab(ylab_text)
    fig <- (p_left + p_right) + plot_layout(widths = c(1, n_other), guides = "collect")
    width <- 15
    height <- 3.2 * n_scen
  }
  fig <- fig & theme(legend.position = "right")

  ggsave(png_path, fig, width = width, height = height, dpi = 150)
}

## Shared body of the minimal by-scenario plots: just max-|T| and
## per-parameter, no uncorrected panel and so no patchwork needed (a single
## facet_grid suffices), reusing the same build_other_criteria_plot() the
## full plot uses for its other-criteria panel. Scenario strips are shown
## here since there's no uncorrected panel to supply them. flip = TRUE puts
## criteria in rows and scenarios in columns.
## `summary_fn` selects which family of metrics to plot -- all-effects
## (build_corrections_by_scenario_summary()) or main-effects-only
## (build_main_effects_summary()); `effects_label` names that family in
## the y-axis title, since "all effects" vs. "main effects only" isn't
## otherwise distinguishable between the two plots this builds.
make_minimal_scenario_plot <- function(results_path, png_path, summary_fn, effects_label,
                                       flip = FALSE) {
  summ <- summary_fn(results_path) |>
    dplyr::filter(interactions) |>
    dplyr::mutate(multcomp_k = factor(multcomp_k, levels = c("minimal", "maximal", "none")))
  theme_set(theme_bw())

  minimal_data <- summ |> dplyr::filter(criterion %in% criterion_labels[c("maxT", "perparam")])
  fig <- build_other_criteria_plot(minimal_data, show_scenario_strip = TRUE,
                                    ylab_text = paste0("Proportion with type I errors (", effects_label, ")"),
                                    flip = flip)

  n_scen <- length(scenario_labels)
  if (flip) {
    ggsave(png_path, fig, width = 3.5 * n_scen + 2.5, height = 3.2 * 2 + 1.2, dpi = 150)
  } else {
    ggsave(png_path, fig, width = 8, height = 3.2 * n_scen, dpi = 150)
  }
}

make_corrections_by_scenario_minimal_plot <- function(results_path, png_path, flip = FALSE) {
  make_minimal_scenario_plot(results_path, png_path, build_corrections_by_scenario_summary,
                              effects_label = "all effects", flip = flip)
}

## Same layout as make_corrections_by_scenario_minimal_plot(), but using
## the main-effects-only metrics (see build_main_effects_summary()).
make_corrections_by_scenario_minimal_mainonly_plot <- function(results_path, png_path, flip = FALSE) {
  make_minimal_scenario_plot(results_path, png_path, build_main_effects_summary,
                              effects_label = "main effects only", flip = flip)
}

## All five criteria (including uncorrected and per-parameter) as colour,
## both multcomp_k calibrations as linetype/shape, so a single facet_grid
## still suffices (no patchwork): rows = scenario (matching the row layout
## of the other figures), columns = N (with-interactions cases only).
## Mixing the uncorrected/per-parameter and corrected criteria onto one
## shared y-scale means the "Selected" rows' corrected lines get visually
## compressed near 0.05 next to the much larger uncorrected values -- kept
## simple here rather than reaching for free_y or a secondary axis.
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
    m_k_scale_x() +
    scale_y_logit(c(summ$ci_lo, summ$ci_hi)) +
    ylab("Proportion with type I errors") +
    zmargin

  ggsave(png_path, fig, width = 10, height = 8, dpi = 150)
}

## Rows = N, columns = criterion (uncorrected first, then the three
## corrections, then per-parameter), colour = model-simplification scenario,
## linetype/shape =
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
    m_k_scale_x() +
    scale_y_logit(c(summ$ci_lo, summ$ci_hi)) +
    ylab("Proportion with type I errors") +
    zmargin +
    theme(strip.text.y = element_text(angle = 0))

  ggsave(png_path, fig, width = 10, height = 8, dpi = 150)
}
