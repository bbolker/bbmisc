## Reads output/corrections_results.rds and writes two figures comparing
## four multiple comparisons corrections (Bonferroni, Dunn-Sidak, Holm,
## single-step max-|T|) applied to the full model (no step()
## simplification), for comparison with Fig. 1a:
##   output/fig1_corrections.png          -- experimentwise type I error
##                                           rate (at least one term significant)
##   output/fig1_corrections_perparam.png -- per-parameter type I error rate
##                                           E[V]/k (V = number of significant terms)
##   output/fig1_corrections_2panel.png   -- experimentwise, Bonferroni and
##                                           max-|T| only (Holm is identical to
##                                           Bonferroni here, Dunn-Sidak nearly so)

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

source(here::here("forstmeier_sim", "R", "graphics_utils.R"))

results <- readRDS(here::here("forstmeier_sim", "output", "corrections_results.rds"))

method_labels <- c(
  bonferroni = "(a) Bonferroni",
  dunnsidak = "(b) Dunn–Šidák",
  holm = "(c) Holm",
  maxT = "(d) Single-step max|T|"
)

## corrections_results.rds screens three model-simplification scenarios;
## these plots cover only the unselected (full-model) one -- see
## fig1_corrections_by_scenario.png for the selected cases.

## One row per replicate and method: the experimentwise indicator (sig) and
## the per-parameter rate (perparam = number significant / k).
long_method <- function(data, method) {
  data |>
    transmute(N, m, interactions, k, method = method_labels[[method]],
              sig = .data[[paste0("sig_", method, "_unselected")]],
              perparam = .data[[paste0("n_sig_", method, "_unselected")]] / k)
}

## Gaussian-approximation 95% CIs: p +/- 1.96*se. For the per-parameter
## rate (a mean of per-replicate rates in [0, 1], not of 0/1 indicators)
## the se is computed from the replicate values' own sd.
summ <- bind_rows(lapply(names(method_labels), long_method, data = results)) |>
  group_by(N, m, interactions, k, method) |>
  summarise(
    n = n(),
    prop_sig = mean(sig),
    se_sig = sqrt(prop_sig * (1 - prop_sig) / n),
    prop_perparam = mean(perparam),
    se_perparam = sd(perparam) / sqrt(n),
    .groups = "drop"
  ) |>
  mutate(
    method = factor(method, levels = method_labels),
    interactions_label = ifelse(interactions, "With 2-way interactions", "Without interactions")
  )

with_ci <- function(data, type) {
  data |>
    mutate(
      prop = .data[[paste0("prop_", type)]],
      se = .data[[paste0("se_", type)]],
      ci_lo = pmax(0, prop - 1.96 * se),
      ci_hi = pmin(1, prop + 1.96 * se)
    )
}

theme_set(theme_bw())

## Shared layers: lines, points and CI ribbons coloured by interactions,
## linetype by N, one facet per method.
corrections_plot <- function(data) {
  ggplot(data, aes(m, prop, colour = interactions_label, linetype = factor(N))) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = interactions_label,
                    group = interaction(interactions_label, N)),
                colour = NA, alpha = 0.3, show.legend = FALSE) +
    geom_line() +
    geom_point() +
    scale_colour_manual(values = okabe_ito, name = NULL) +
    scale_fill_manual(values = okabe_ito, guide = "none") +
    scale_linetype_discrete(name = "N") +
    m_k_scale_x()
}

summ_sig <- with_ci(summ, "sig")

experimentwise_plot <- function(data) {
  corrections_plot(data) +
    ## ymin = -Inf rather than 0: 0 is not representable on a logit scale
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.05,
             fill = "grey80", alpha = 0.5) +
    geom_hline(yintercept = 0.05, colour = "grey40", linetype = "dashed", linewidth = 1) +
    facet_wrap(~method) +
    ## limits fitted to the CI range, with breaks finer than logit_breaks()'s
    ## decades/halves, which would stretch the axis to the next 0.01/0.1
    scale_y_continuous(trans = "logit", limits = range(c(data$ci_lo, data$ci_hi, 0.05)),
                       breaks = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.1)) +
    ylab("Proportion with type I errors (experimentwise)")
}

ggsave(here::here("forstmeier_sim", "output", "fig1_corrections.png"),
       experimentwise_plot(summ_sig), width = 10, height = 8, dpi = 150)

two_panel_labels <- c("(a) Bonferroni", "(b) Single-step max|T|")
summ_sig_2panel <- summ_sig |>
  filter(method %in% method_labels[c("bonferroni", "maxT")]) |>
  mutate(method = factor(method, levels = method_labels[c("bonferroni", "maxT")],
                         labels = two_panel_labels))

ggsave(here::here("forstmeier_sim", "output", "fig1_corrections_2panel.png"),
       experimentwise_plot(summ_sig_2panel), width = 10, height = 4.5, dpi = 150)

## Per-parameter version. Reference curve alpha/k: under the global null
## each Bonferroni test rejects with probability exactly alpha/k, so this
## is Bonferroni's expected per-parameter rate. Holm and Dunn-Sidak reject
## every term Bonferroni does, so their rates can only be at or above it.
summ_perparam <- with_ci(summ, "perparam")
alpha_over_k <- summ_perparam |>
  distinct(m, interactions, k) |>
  mutate(ref = 0.05 / k)

fig1_corrections_perparam <- corrections_plot(summ_perparam) +
  geom_line(data = alpha_over_k, aes(x = m, y = ref, group = interactions),
            inherit.aes = FALSE, colour = "grey50", linewidth = 2, alpha = 0.5) +
  facet_wrap(~method) +
  scale_y_logit(c(summ_perparam$ci_lo, summ_perparam$ci_hi)) +
  ylab("Per-parameter type I error rate")

ggsave(here::here("forstmeier_sim", "output", "fig1_corrections_perparam.png"),
       fig1_corrections_perparam, width = 10, height = 8, dpi = 150)
