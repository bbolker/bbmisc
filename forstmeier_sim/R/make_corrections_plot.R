## Reads output/corrections_results.rds and writes
## output/fig1_corrections.png: type I error rate under four significance
## criteria (uncorrected, Bonferroni-type/Sidak, Holm, single-step max-|T|)
## applied to the full model (no step() simplification), for comparison
## with Fig. 1a.

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

source(here::here("forstmeier_sim", "R", "graphics_utils.R"))

results <- readRDS(here::here("forstmeier_sim", "output", "corrections_results.rds"))

method_labels <- c(
  sig_raw  = "(a) Uncorrected",
  sig_bonf = "(b) Dunn–Šidák",
  sig_holm = "(c) Holm",
  sig_maxT = "(d) Single-step max|T|"
)

## corrections_results.rds screens three model-simplification scenarios;
## this plot covers only the unselected (full-model) one, for comparison
## with Fig. 1a -- see fig1_corrections_by_scenario.png for the selected
## cases.
summ_raw <- results |>
  group_by(N, m, interactions, k) |>
  summarise(
    n = n(),
    sig_raw  = mean(sig_raw_unselected),
    sig_bonf = mean(sig_bonf_unselected),
    sig_holm = mean(sig_holm_unselected),
    sig_maxT = mean(sig_maxT_unselected),
    .groups = "drop"
  )

## Gaussian-approximation 95% CI for a binomial proportion: p +/- 1.96*se
long_method <- function(data, prop_col, method_name) {
  data |>
    transmute(N, m, interactions, k, method = method_name, prop = .data[[prop_col]], n) |>
    mutate(
      se = sqrt(prop * (1 - prop) / n),
      ci_lo = pmax(0, prop - 1.96 * se),
      ci_hi = pmin(1, prop + 1.96 * se)
    )
}

summ <- bind_rows(
  long_method(summ_raw, "sig_raw", method_labels[["sig_raw"]]),
  long_method(summ_raw, "sig_bonf", method_labels[["sig_bonf"]]),
  long_method(summ_raw, "sig_holm", method_labels[["sig_holm"]]),
  long_method(summ_raw, "sig_maxT", method_labels[["sig_maxT"]])
) |>
  mutate(
    method = factor(method, levels = method_labels),
    interactions_label = ifelse(interactions, "With 2-way interactions", "Without interactions")
  )

## theoretical independence expectation alpha' = 1-(1-alpha)^k is only a
## meaningful reference for the uncorrected panel -- the three corrections
## are designed to sit near the flat alpha=0.05 zone regardless of k
grey_lines <- summ |>
  filter(method == method_labels[["sig_raw"]]) |>
  distinct(m, interactions, k, method) |>
  mutate(alpha_prime = 1 - (1 - 0.05)^k)

## flat alpha=0.05 reference, shown only on the corrected panels (the
## uncorrected panel already has the theoretical alpha' curve for reference)
nominal_ref <- summ |>
  filter(method != method_labels[["sig_raw"]]) |>
  distinct(method)

theme_set(theme_bw())

fig1_corrections <- ggplot(summ, aes(m, prop, colour = interactions_label, linetype = factor(N))) +
  ## ymin = -Inf rather than 0: 0 is not representable on a logit scale
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.05,
           fill = "grey80", alpha = 0.5) +
  geom_line(data = grey_lines, aes(x = m, y = alpha_prime, group = interactions),
            inherit.aes = FALSE, colour = "grey50", linewidth = 2, alpha = 0.5) +
  geom_hline(data = nominal_ref, aes(yintercept = 0.05), inherit.aes = FALSE,
             colour = "grey40", linetype = "dashed", linewidth = 1) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = interactions_label, group = interaction(interactions_label, N)),
              colour = NA, alpha = 0.3, show.legend = FALSE) +
  geom_line() +
  geom_point() +
  facet_wrap(~method, scales = "free_y") +
  scale_colour_manual(values = okabe_ito, name = NULL) +
  scale_fill_manual(values = okabe_ito, guide = "none") +
  scale_linetype_discrete(name = "N") +
  scale_x_continuous(name = "Number of explanatory variables", breaks = 1:6,
                     sec.axis = k_sec_axis()) +
  scale_y_continuous(trans = "logit", breaks = logit_breaks()) +
  ylab("Proportion of models with type I errors")

ggsave(here::here("forstmeier_sim", "output", "fig1_corrections.png"), fig1_corrections,
       width = 10, height = 8, dpi = 150)
