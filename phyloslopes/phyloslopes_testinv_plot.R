library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)

## Okabe-Ito palette minus black and yellow
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

base_breaks <- function(n = 10) {
  function(x) axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
}

recode_type <- function(x) {
  recode_values(x,
                "elapsed" ~ "elapsed time (sec)",
                "size_bytes" ~ "object size (Mb)",
                "sparsity" ~ "1-sparsity",
                default = x)
}

dd <- readRDS("phyloslopes_testinv.rds") |>
  bind_rows() |>
  select(n_tips, rep, step, elapsed, size_bytes, sparsity) |>
  mutate(across(size_bytes, \(x) x / 1024^2)) |>
  pivot_longer(-c(n_tips, rep, step), names_to = "type") |>
  mutate(across(type, recode_type))

dd$type <- factor(dd$type,
                   levels = c("elapsed time (sec)", "object size (Mb)", "1-sparsity"))

## put vcv and solve first, keep the rest in their existing order
step_levels <- unique(dd$step)
step_levels <- c("vcv", "solve", setdiff(step_levels, c("vcv", "solve")))
dd$step <- factor(dd$step, levels = step_levels)

p <- ggplot(dd, aes(n_tips, value, colour = step, fill = step)) +
  geom_point() +
  ## only 6 distinct n_tips values, so the gam default k = 10 knots
  ## fails for every group; k = 5 stays below that
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5)) +
  scale_colour_manual(values = okabe_ito) +
  scale_fill_manual(values = okabe_ito) +
  facet_wrap(~type, scales = "free") +
  scale_y_log10(breaks = base_breaks()) +
  scale_x_log10(breaks = base_breaks()) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave("phyloslopes_testinv_plot.png", p, width = 10, height = 4)

