## Combines one or more Fig.1-style results .rds files (each with columns
## N, m, interactions, k, sig_full, sig_step) into a comparison across
## model-simplification methods. Two orientations share one summary builder:
## make_stepwise_comparison_plot() facets by method (rows) with (N,
## interactions) as colour/linetype; make_stepwise_comparison_plot_by_condition()
## facets by (N, interactions) with method as colour, for comparing the
## simplification methods directly against each other within a fixed
## condition.
##
## `specs` is a named list, in the row/colour order desired; each name
## becomes a label, and each element is list(path = <rds path>, column =
## "sig_full" or "sig_step"). E.g. two specs sharing one path can pull out
## the unselected (sig_full) and selected (sig_step) rows from the same
## simulation run; a third spec from a different .rds adds another
## simplification variant (e.g. interactions-only-removal step()).

library(ggplot2)
library(scales)

source(here::here("forstmeier_sim", "R", "graphics_utils.R"))

## Long-format summary shared by both plot orientations: one row per
## (method, N, m, interactions) with the proportion and its Gaussian-
## approximation 95% CI.
build_stepwise_summary <- function(specs) {
  ## purrr::imap() calls .f(value, name, ...), i.e. (spec, label)
  long_one <- function(spec, label) {
    readRDS(spec$path) |>
      dplyr::group_by(N, m, interactions, k) |>
      dplyr::summarise(n = dplyr::n(), prop = mean(.data[[spec$column]]), .groups = "drop") |>
      dplyr::mutate(
        method = label,
        se = sqrt(prop * (1 - prop) / n),
        ci_lo = pmax(0, prop - 1.96 * se),
        ci_hi = pmin(1, prop + 1.96 * se)
      )
  }

  purrr::imap(specs, long_one) |>
    dplyr::bind_rows() |>
    dplyr::mutate(method = factor(method, levels = names(specs)))
}

## Rows = method, colour/linetype = (N, interactions) -- no free scale
## needed, since every row uses the same (uncorrected) significance
## criterion and so shares a comparable value range.
make_stepwise_comparison_plot <- function(specs, png_path) {
  summ <- build_stepwise_summary(specs) |>
    dplyr::mutate(
      group = interaction(
        ifelse(interactions, "With 2-way interactions", "Without interactions"),
        paste0("N = ", N),
        sep = ", "
      )
    )

  ## theoretical expectation alpha' = 1 - (1-alpha)^k, one line per
  ## interaction scenario, repeated in every row
  grey_lines <- summ |>
    dplyr::distinct(method, m, interactions, k) |>
    dplyr::mutate(alpha_prime = 1 - (1 - 0.05)^k)

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
    facet_grid(rows = vars(method)) +
    scale_colour_manual(values = okabe_ito, name = NULL) +
    scale_fill_manual(values = okabe_ito, guide = "none") +
    scale_linetype_discrete(name = "N") +
    scale_x_continuous(name = "Number of explanatory variables", breaks = 1:6) +
    scale_y_continuous(trans = "logit", breaks = logit_breaks()) +
    ylab("Proportion of models with type I errors")

  ggsave(png_path, fig, width = 7, height = 3.2 * length(specs), dpi = 150)
}

## Rows = interactions, columns = N, colour = method -- lets the three
## simplification methods be compared directly within a fixed condition.
make_stepwise_comparison_plot_by_condition <- function(specs, png_path) {
  summ <- build_stepwise_summary(specs) |>
    dplyr::mutate(
      interactions_label = factor(interactions, levels = c(FALSE, TRUE),
                                   labels = c("Without interactions", "With 2-way interactions")),
      N_label = factor(paste0("N = ", N), levels = paste0("N = ", sort(unique(N))))
    )

  grey_lines <- summ |>
    dplyr::distinct(interactions_label, N_label, m, k) |>
    dplyr::mutate(alpha_prime = 1 - (1 - 0.05)^k)

  theme_set(theme_bw())

  fig <- ggplot(summ, aes(m, prop, colour = method)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.05,
             fill = "grey80", alpha = 0.5) +
    geom_line(data = grey_lines, aes(x = m, y = alpha_prime),
              inherit.aes = FALSE, colour = "grey50", linewidth = 2, alpha = 0.5) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = method),
                colour = NA, alpha = 0.3, show.legend = FALSE) +
    geom_line() +
    geom_point() +
    facet_grid(rows = vars(interactions_label), cols = vars(N_label)) +
    scale_colour_manual(values = okabe_ito[seq_along(specs)], name = NULL) +
    scale_fill_manual(values = okabe_ito[seq_along(specs)], guide = "none") +
    scale_x_continuous(name = "Number of explanatory variables", breaks = 1:6) +
    scale_y_continuous(trans = "logit", breaks = logit_breaks()) +
    ylab("Proportion of models with type I errors")

  ggsave(png_path, fig, width = 8, height = 6, dpi = 150)
}
