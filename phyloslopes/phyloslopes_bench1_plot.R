## Plot of the phyloslopes_bench1.R timing results
## (phylo_bench.rds / phylo_bench_slopes.rds).

library(ggplot2); theme_set(theme_bw())

bench <- readRDS("phylo_bench.rds")
bench_slopes <- readRDS("phylo_bench_slopes.rds")

## combined plot: hand-built (not autoplot()) violin plot, times on x,
## methods on y, faceted top/bottom by model class with free x scales
## (the two classes differ by ~1-2 orders of magnitude in fitting time).
## phyr_compare is excluded here (its timings are so much faster than
## everything else that it squashes the other violins) but stays in
## `bench`/phylo_bench.rds -- see its median/5%/95% printed by phyloslopes_bench1.R
bench_df_intercept <- as.data.frame(bench)
bench_df_intercept <- droplevels(bench_df_intercept[bench_df_intercept$expr != "phyr_compare", ])
bench_df <- rbind(
  transform(bench_df_intercept, model_type = "random-intercept only"),
  transform(as.data.frame(bench_slopes), model_type = "random-slopes (+ intercept)")
)
bench_df$time_ms <- bench_df$time / 1e6
bench_df$model_type <- factor(bench_df$model_type,
                              levels = c("random-slopes (+ intercept)",
                                         "random-intercept only"))

bench_plot <- ggplot(bench_df, aes(x = time_ms, y = expr)) +
  geom_violin(fill = "gray") +
  scale_x_log10() +
  facet_wrap(~ model_type, ncol = 1, scales = "free") +
  labs(x = "time (ms)", y = NULL,
      title = "phyloslopes: fitting time by method")
ggsave("phyloslopes_bench.png", bench_plot, width = 8, height = 7)

## intercept-only model, standalone (no faceting); phyr_compare already
## excluded from bench_df_intercept above
bench_plot_intercept <- ggplot(bench_df_intercept, aes(x = time / 1e6, y = expr)) +
  geom_violin(fill = "gray") +
  scale_x_log10() +
  labs(x = "time (ms)", y = NULL,
      title = "phyloslopes: fitting time by method (random-intercept only)")
ggsave("phyloslopes_bench1.png", bench_plot_intercept, width = 8, height = 4)
