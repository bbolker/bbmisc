## Shared plotting logic for the Fig. 1 reproduction (uncorrected vs.
## step()-simplified type I error rate). Reused by make_plot.R (unrestricted
## step()) and make_plot_restricted.R (interactions-only-removal step()).
## Reads the "raw" (uncorrected) columns out of corrections_results.rds --
## col_unselected/col_selected name which pair of scenario columns to plot.

library(ggplot2)
library(scales)

source(here::here("forstmeier_sim", "R", "graphics_utils.R"))

make_fig1_plot <- function(results_path, col_unselected, col_selected, png_path) {
  results <- readRDS(results_path)

  panel_labels <- c(
    unselected = "(a) Without model simplification",
    selected   = "(b) With model simplification"
  )

  summ_raw <- results |>
    dplyr::group_by(N, m, interactions, k) |>
    dplyr::summarise(n = dplyr::n(),
                      prop_unselected = mean(.data[[col_unselected]]),
                      prop_selected = mean(.data[[col_selected]]),
                      .groups = "drop")

  ## Gaussian-approximation 95% CI for a binomial proportion: p +/- 1.96*se
  long_panel <- function(data, prop_col, panel_name) {
    data |>
      dplyr::transmute(N, m, interactions, k, panel = panel_name, prop = .data[[prop_col]], n) |>
      dplyr::mutate(
        se = sqrt(prop * (1 - prop) / n),
        ci_lo = pmax(0, prop - 1.96 * se),
        ci_hi = pmin(1, prop + 1.96 * se)
      )
  }

  summ <- dplyr::bind_rows(
    long_panel(summ_raw, "prop_unselected", panel_labels[["unselected"]]),
    long_panel(summ_raw, "prop_selected", panel_labels[["selected"]])
  ) |>
    dplyr::mutate(
      panel = factor(panel, levels = panel_labels),
      group = interaction(
        ifelse(interactions, "With 2-way interactions", "Without interactions"),
        paste0("N = ", N),
        sep = ", "
      )
    )

  ## theoretical expectation alpha' = 1 - (1-alpha)^k, one line per
  ## interaction scenario, same for both panels
  grey_lines <- summ |>
    dplyr::distinct(m, interactions, k) |>
    dplyr::mutate(alpha_prime = 1 - (1 - 0.05)^k) |>
    tidyr::crossing(panel = factor(panel_labels, levels = panel_labels))

  theme_set(theme_bw())

  fig <- ggplot(summ, aes(m, prop, colour = group, linetype = factor(N))) +
    ## ymin = -Inf rather than 0: 0 is not representable on a logit scale
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.05,
             fill = "grey80", alpha = 0.5) +
    geom_line(data = grey_lines, aes(x = m, y = alpha_prime, group = interactions),
              inherit.aes = FALSE, colour = "grey50", linewidth = 2, alpha = 0.5) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = group, group = interaction(group, N)),
                colour = NA, alpha = 0.3, show.legend = FALSE) +
    geom_line() +
    geom_point() +
    facet_wrap(~panel) +
    scale_colour_manual(values = okabe_ito, name = NULL) +
    scale_fill_manual(values = okabe_ito, guide = "none") +
    scale_linetype_discrete(name = "N") +
    scale_x_continuous(name = "Number of explanatory variables", breaks = 1:6) +
    scale_y_continuous(trans = "logit", breaks = logit_breaks()) +
    ylab("Proportion of models with type I errors")

  ggsave(png_path, fig, width = 10, height = 5, dpi = 150)
}
