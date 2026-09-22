library(ggplot2); theme_set(theme_bw())
library(tidyverse)
dd <- read.csv("package_counts_by_year.csv") |>
  mutate(across(package, ~reorder(factor(.), X = mixed_count,
                                  FUN = \(x) tail(x,1)))) |>
  pivot_longer(ends_with("count"), names_to = "type")
ggplot(dd, aes(year, value, colour = package)) +
  facet_wrap(~type, scale = "free") +
  geom_line() +
  scale_x_continuous(breaks = seq(2015, 2025, by = 5))
