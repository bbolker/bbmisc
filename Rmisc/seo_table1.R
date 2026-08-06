library(tidyverse); theme_set(theme_bw())
library(patchwork)
zmargin <- theme(panel.spacing = grid::unit(0, "lines"))

dd <- read_csv("seo_table1.csv") |>
    ## drop constant columns
    select(-c(k, Delta, `Lambda^2`)) |>
    pivot_longer(-matches("n[12]{2}"),
                 names_to = "variable", values_to = "value") |>
    separate(variable, into = c("type", "correct", "method"), sep = "\\.") |>
    mutate(vars = paste(n11, n21, n12, n22, sep = ", "), ntot = n11 + n21 + n12 + n22) |>
    ## order variable sets by total sample size
    mutate(vars = factor(vars, levels = unique(vars)[order(unique(ntot))]),
           method = factor(method, levels = c("GS", "emp", "plug", "N3", "N5", "N7")))


## 'b', 'c' are biased vs corrected; 'x' applies when neither of those is relevant (e.g. GS, plug-in)
## use open points for x, closed points with different fills for b and c?

## baseline plot (use as basis for downstream plots)
theme_set(theme_bw(base_size = 20))
gg0 <- ggplot(dd, aes(x = method, y = value, colour = correct, shape = correct)) +
    zmargin +
    scale_shape_manual(values = c(16,17, 1)) + scale_colour_manual(values = c("red", "blue", "black"))

gg0 + geom_point() + facet_grid(type~vars, scales = "free")


## other things to play with:
## for a variety of reasons it would be better to create separate plots for coverage and
##   variance, put them together with `patchwork` or `cowplot::plot_grid`
##     * coverage might be plotted on logit scale
##     * for coverage, add reference line at 0.95, binomial CIs on reference line *or* CIs on individual points
##     * for variance, add reference line at GS for each panel?
## * better colour scheme (Okabe-Ito?)

dds <- split(dd, dd$type)

nsim <- 2000
cval <- 0.95
qq <- qnorm(c(0.025, 0.975))
## manually set breaks for logit scale
brkvec <- c(0.25, 0.5, 0.75, 0.9, 0.95)
## Gaussian-approx binom CIs around nominal coverage
binom_sd <- sqrt(cval*(1-cval)/nsim)

## Gaussian-approx binom CIs around observed coverage values
dds$cover2 <- (dds$cover
    |> mutate(se = sqrt(value*(1-value)/nsim),
              lwr = value + qq[1]*se,
              upr = value + qq[2]*se))

plot_cover <- gg0 %+% dds$cover2 +
    geom_pointrange(aes(ymin = lwr, ymax = upr),
                    position = position_dodge(width=0.25)) +
    facet_wrap(~vars, nrow = 1) +
    labs(title = "Coverage", y  = "") + 
    scale_y_continuous(trans = "logit", breaks = brkvec) +
    geom_hline(yintercept = cval, lty = 2) +
    annotate("rect", alpha = 0.25, color = NA, fill = "black",
             xmin = -Inf, xmax = Inf, ymin = cval-1.96*binom_sd, ymax = cval + 1.96*binom_sd)

GS <- dds$var |> filter(method == "GS") |> select(vars, value)
plot_var <- gg0 %+% filter(dds$var, method != "GS") +
    geom_point() +
    facet_wrap(~vars, nrow = 1) +
    labs(title = "Variance", y = "") +
    geom_hline(data = GS, aes(yintercept = value), color = "black", lty = 2) +
    expand_limits(y=0)

## needs patchwork package loaded!
plot_cover / plot_var

## FIXME: add a filled-white foreground for "x" points (pch = 22, pt.bg = "white") ?
