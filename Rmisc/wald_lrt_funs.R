## Computational and plotting functions shared by wald_lrt_examples.qmd and
## wald_lrt_examples_lowp.R. Needs glmmTMB, car, tidyverse and futurize
## (and patchwork, for the plots) to be loaded; the plot functions also
## use zmargin, defined by the calling script. The '## ---- label' lines
## mark chunks for knitr::read_chunk() (the qmd pulls each section into
## the empty chunk with the same label); they are ordinary comments when
## this file is source()d.

## ---- nb-funs
## extract drop1 (LRT) and Anova (Wald) p-values
get_p_values <- function(model, w = 1) {
  d <- drop1(model, test="Chisq")
  p_drop <- d[["Pr(>Chi)"]][[w+1]]
  a <- Anova(model)
  p_anova <- a[["Pr(>Chisq)"]][[w]]
  tibble(p_drop, p_anova)
}

## simplified null-hypothesis simulation (no zi);
## species effect (fixed effects parameters 2-7) is set to zero
simfun <- function(n = nrow(Salamanders)) {
    dd <- transform(Salamanders,
            count = simulate_new(~ spp + mined + (1|site),
            newdata = Salamanders,
            newparams = list(beta = c(-0.5, rep(0,6), 1.4),
            betadisp = 0.5,
            theta = -1),
            family = nbinom2)[[1]])
    ## randomly subsample (sorting isn't strictly necessary)
    dd <- dd[sort(sample(nrow(Salamanders), size = n)), ]
    return(dd)
}

fitfun <- function(data = simfun()) {
  model <- glmmTMB(count ~ spp + mined + (1|site),
                   family=nbinom2, data=data)
  return(model)
}

sumfun <- function(model = fitfun()) {
  get_p_values(model)
}

full_sim <- function(n = nrow(Salamanders), nsim = 1000, seed = 101) {
  set.seed(seed)
  tt <- system.time(
    res <- replicate(nsim, sumfun(fitfun(simfun(n = n))), simplify = FALSE) |>
      futurize() |>
      bind_rows()
  )
  attr(res, "time") <- tt
  res
}

rpt_time <- function(x) { attr(x, "time")[["elapsed"]] |> round() }

## ---- binom-exact-funs
## log priors on each group's logit: Jeffreys (Firth) and its
## glmmTMB-compatible t approximation
lp_firth <- function(eta) -log(cosh(eta/2))
t_scale <- round(pi * sqrt(5/7), 2)
lp_t <- function(eta) dt(eta/t_scale, df = 7, log = TRUE) - log(t_scale)
lp_none <- function(eta) 0

## penalized log-likelihood for two groups (y0, y1 successes out of n)
## given the baseline log-odds alpha, log-odds ratio psi and log prior lp
ll_2g <- function(alpha, psi, y0, y1, n, lp = lp_none) {
  dbinom(y0, n, plogis(alpha), log = TRUE) + lp(alpha) +
    dbinom(y1, n, plogis(alpha + psi), log = TRUE) + lp(alpha + psi)
}

## (penalized) maximum for a single group's logit; the model is saturated,
## so the overall maximum is the sum of the per-group maxima
gmax_2g <- function(y, n, lp = lp_none, alpha_max = 30) {
  optimize(\(eta) dbinom(y, n, plogis(eta), log = TRUE) + lp(eta),
           c(-alpha_max, alpha_max), maximum = TRUE)
}

## penalized Wald CI: MAP estimate +/- z * SE from the penalized
## information, n p q - lp''(eta) for each group (lp'' by finite differences)
pen_wald_ci_2g <- function(y0, y1, n, lp, h = 1e-4) {
  info <- function(m) {
    eta <- m$maximum
    lp2 <- (lp(eta + h) - 2 * lp(eta) + lp(eta - h)) / h^2
    n * plogis(eta) * plogis(-eta) - lp2
  }
  m0 <- gmax_2g(y0, n, lp)
  m1 <- gmax_2g(y1, n, lp)
  se <- sqrt(1/info(m0) + 1/info(m1))
  (m1$maximum - m0$maximum) + c(-1, 1) * qnorm(0.975) * se
}

## (penalized) profile CI for one outcome; the limits are searched for
## over [-bmax, bmax], as in the simulations
prof_ci_2g <- function(y0, y1, n, lp = lp_none, bmax = 20, alpha_max = 30) {
  crit <- qchisq(0.95, 1)
  m0 <- gmax_2g(y0, n, lp, alpha_max)
  m1 <- gmax_2g(y1, n, lp, alpha_max)
  ll_max <- m0$objective + m1$objective
  prof_dev <- function(psi) {
    ll <- optimize(ll_2g, c(-alpha_max, alpha_max), psi = psi,
                   y0 = y0, y1 = y1, n = n, lp = lp, maximum = TRUE)$objective
    2 * (ll_max - ll) - crit
  }
  ## without a prior the estimate is infinite under separation: start the
  ## search just inside the boundary on the side the data point to
  ctr <- pmin(pmax(m1$maximum - m0$maximum, -0.99 * bmax), 0.99 * bmax)
  ## no sign change between the estimate and the limit means the
  ## profile interval is unbounded on that side
  root <- function(lim) {
    if (prof_dev(lim) < 0) return(sign(lim) * Inf)
    uniroot(prof_dev, sort(c(ctr, lim)))$root
  }
  c(root(-bmax), root(bmax))
}

## unpenalized and penalized estimates, Wald and profile CIs for one outcome
ci_2g <- function(y0, y1, n = 50) {
  sep <- y0 %in% c(0, n) || y1 %in% c(0, n)
  ## infinite under separation (NaN if both groups are all 0s or all 1s)
  est <- qlogis(y1/n) - qlogis(y0/n)
  pen_est <- function(lp) gmax_2g(y1, n, lp)$maximum - gmax_2g(y0, n, lp)$maximum
  se <- sqrt(1/y0 + 1/(n - y0) + 1/y1 + 1/(n - y1))
  wald <- if (sep) c(-Inf, Inf) else est + c(-1, 1) * qnorm(0.975) * se
  cis <- rbind(wald,
               prof_ci_2g(y0, y1, n),
               pen_wald_ci_2g(y0, y1, n, lp = lp_firth),
               prof_ci_2g(y0, y1, n, lp = lp_firth),
               pen_wald_ci_2g(y0, y1, n, lp = lp_t),
               prof_ci_2g(y0, y1, n, lp = lp_t))
  # tibble::tibble() used explicitly to avoid 'future' warnings
  tibble::tibble(y0, y1,
                 type = rep(c("Wald", "profile"), 3),
                 penalty = rep(c("none", "Firth", "t prior"), each = 2),
                 est = rep(c(est, pen_est(lp_firth), pen_est(lp_t)), each = 2),
                 lwr = cis[, 1],
                 upr = cis[, 2])
}

## exact one-sided coverage at true log-odds ratio beta
exact_cov <- function(beta, cis, n = 50, p0 = 0.2) {
  cis |>
    mutate(w = dbinom(y0, n, p0) * dbinom(y1, n, plogis(qlogis(p0) + beta))) |>
    group_by(type, penalty) |>
    summarise(lower = sum(w * (lwr <= beta)),
              upper = sum(w * (upr >= beta)),
              .groups = "drop") |>
    mutate(beta = beta)
}

## estimates and CIs for all possible outcomes with group size n
## (these don't depend on the baseline probability p0, which only
## enters through the weights in exact_cov() and exact_medians())
exact_cis_n <- function(n) {
  expand_grid(y0 = 0:n, y1 = 0:n, n = n) |>
    pmap(ci_2g) |>
    futurize() |>
    bind_rows() |>
    mutate(n = n)
}

## exact coverage over a grid of true log-odds ratios (cis for one n)
exact_cov_n <- function(cis, betavec = seq(0, 6, by = 0.02), p0 = 0.2) {
  n <- cis$n[[1]]
  map_dfr(betavec, exact_cov, cis = cis, n = n, p0 = p0) |>
    mutate(n = n)
}

## weighted median (Inf allowed; NaN, which has negligible probability
## for p0 = 0.2, is dropped)
wmedian <- function(x, w) {
  ok <- !is.na(x)
  x <- x[ok]
  w <- w[ok]
  o <- order(x)
  x[o][which(cumsum(w[o]) >= sum(w)/2)[1]]
}

## exact medians of the estimates and CI limits at true log-odds ratio beta
exact_medians <- function(beta, cis, p0 = 0.2) {
  cis |>
    mutate(w = dbinom(y0, n, p0) * dbinom(y1, n, plogis(qlogis(p0) + beta))) |>
    group_by(n, type, penalty) |>
    summarise(across(c(est, lwr, upr), \(x) wmedian(x, w)),
              .groups = "drop") |>
    mutate(beta = beta)
}

## ---- binom-plot-funs
## label facets by group size, with horizontal row strips
n_labeller <- labeller(n = \(x) paste0("n = ", x, "\nper group"))
hstrip <- theme(strip.text.y = element_text(angle = 0))
pen_levs <- c("none", "Firth", "t prior")

## median estimates and CI limits (output of exact_medians()) vs true beta
plot_medians <- function(med_2g, dodge_w = 0.8) {
  ## a single grouping factor fixes the dodge order (penalty, then CI type)
  meth_levs <- c(outer(c("profile", "Wald"), pen_levs, \(t, p) paste(p, t)))
  med_2g <- med_2g |>
    mutate(penalty = factor(penalty, levels = pen_levs),
           type = factor(type, levels = c("profile", "Wald"), labels = c("LRT", "Wald")),
           method = factor(paste(penalty, if_else(type == "LRT", "profile", "Wald")),
                           levels = meth_levs),
           ## x positions matching position_dodge(dodge_w), for the arrows
           xd = beta + (as.integer(method) - (length(meth_levs) + 1)/2) *
             dodge_w / length(meth_levs))
  ylims <- range(unlist(select(med_2g, est, lwr, upr)) |> keep(is.finite)) +
    c(-0.5, 0.5)
  ## infinite limits are clipped to the panel edge and get arrowheads;
  ## infinite estimates are dropped
  med_inf <- bind_rows(
    filter(med_2g, upr == Inf) |> mutate(y = ylims[2] - 0.3, yend = ylims[2]),
    filter(med_2g, lwr == -Inf) |> mutate(y = ylims[1] + 0.3, yend = ylims[1])
  )
  med_fin <- med_2g |>
    mutate(lwr = pmax(lwr, ylims[1]), upr = pmin(upr, ylims[2]),
           est = if_else(is.finite(est), est, NA_real_))
  ## geom_linerange() + geom_point() rather than geom_pointrange(), which
  ## would drop the whole interval when the estimate is missing; every
  ## layer keeps all six methods at each beta so the dodging lines up
  ggplot(med_fin, aes(beta, est, colour = penalty, shape = type, group = method)) +
    geom_segment(data = distinct(med_2g, n, beta),
                 aes(x = beta - 0.45, xend = beta + 0.45, y = beta, yend = beta),
                 inherit.aes = FALSE, colour = "black") +
    geom_linerange(aes(ymin = lwr, ymax = upr),
                   position = position_dodge(width = dodge_w)) +
    geom_point(position = position_dodge(width = dodge_w), size = 2, na.rm = TRUE) +
    geom_segment(data = med_inf, aes(x = xd, xend = xd, y = y, yend = yend),
                 arrow = arrow(length = grid::unit(4, "pt")), show.legend = FALSE) +
    scale_shape_manual(values = c(LRT = 16, Wald = 17), name = NULL) +
    facet_wrap(~ n, labeller = n_labeller) + zmargin +
    coord_cartesian(ylim = ylims, expand = FALSE) +
    labs(x = "true log-odds ratio",
         y = "log-odds ratio (median estimate and CI limits)")
}

## one-sided (left) and two-sided (right) exact coverage vs true beta
## (output of exact_cov_n()); pen_labs gives the positions of the direct
## penalty labels (columns beta, coverage, n, in the lower-tail panels)
plot_coverage <- function(cov_2g, pen_labs = NULL, nsim_ref = 2000) {
  tail_labs <- c(lower = "lower limit below true value",
                 upper = "upper limit above true value")
  cov_band <- 0.975 + c(-1, 1) * qnorm(0.975) * sqrt(0.975 * 0.025 / nsim_ref)
  cov_2g_long <- cov_2g |>
    pivot_longer(c(lower, upper), names_to = "tail", values_to = "coverage") |>
    mutate(tail = factor(tail_labs[tail], levels = tail_labs),
           penalty = factor(penalty, levels = pen_levs),
           type = factor(type, levels = c("profile", "Wald")),
           ## coverage reaches 1 (to rounding), which would stretch or break a
           ## logit axis, so values above 0.995 are drawn at 0.995
           coverage = pmin(coverage, 0.995))
  p_tails <- ggplot(cov_2g_long, aes(beta, coverage, colour = penalty)) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = cov_band[1], ymax = cov_band[2],
             fill = "gray", alpha = 0.4) +
    geom_hline(yintercept = 0.975, colour = "gray40") +
    geom_line(aes(linetype = type)) +
    scale_y_continuous(transform = "logit", limits = c(NA, 0.995),
                       breaks = c(0.92, 0.95, 0.975, 0.99, 0.995)) +
    scale_linetype_manual(values = c(Wald = "dashed", profile = "solid"),
                          breaks = c("Wald", "profile"),
                          labels = c("Wald", "LRT"), name = NULL) +
    facet_grid(n ~ tail, labeller = n_labeller) + zmargin + hstrip +
    ## inset linetype legend in the lower left corner of the top left panel
    ## (inside positions are relative to the whole panel area; the rows of
    ## panels have equal heights and no spacing)
    theme(legend.position = "inside",
          legend.position.inside = c(0.005, 1 - 1/n_distinct(cov_2g$n) + 0.005),
          legend.justification = c(0, 0),
          legend.background = element_rect(fill = "white", colour = NA)) +
    labs(x = "true log-odds ratio", y = "one-sided coverage (logit scale)")
  if (!is.null(pen_labs)) {
    ## direct labels take the place of the colour legend
    pen_labs <- pen_labs |>
      mutate(penalty = factor(penalty, levels = pen_levs),
             tail = factor(tail_labs[["lower"]], levels = tail_labs))
    p_tails <- p_tails + geom_text(data = pen_labs, aes(label = label)) +
      guides(colour = "none")
  }

  ## two-sided coverage: missing below and missing above are disjoint events
  cov_band2 <- 0.95 + c(-1, 1) * qnorm(0.975) * sqrt(0.95 * 0.05 / nsim_ref)
  cov_2g_two <- cov_2g |>
    mutate(coverage = lower + upper - 1,
           ## single-level column factor, for a top strip matching p_tails
           tail = "true value within limits",
           penalty = factor(penalty, levels = pen_levs),
           type = factor(type, levels = c("profile", "Wald")))
  p_two <- ggplot(cov_2g_two, aes(beta, coverage, colour = penalty)) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = cov_band2[1], ymax = cov_band2[2],
             fill = "gray", alpha = 0.4) +
    geom_hline(yintercept = 0.95, colour = "gray40") +
    geom_line(aes(linetype = type)) +
    scale_linetype_manual(values = c(Wald = "dashed", profile = "solid")) +
    ## breaks kept away from the panel edges, so labels of adjacent rows
    ## don't collide
    scale_y_continuous(breaks = seq(0.88, 0.98, by = 0.02)) +
    facet_grid(n ~ tail) + zmargin +
    ## row strips would repeat those of the one-sided plot
    theme(strip.text.y = element_blank(), strip.background.y = element_blank(),
          legend.position = "none") +
    labs(x = "true log-odds ratio", y = "two-sided coverage (linear scale)")

  p_tails + p_two + plot_layout(widths = c(2, 1))
}
