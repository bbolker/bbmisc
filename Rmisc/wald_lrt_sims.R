## Simulation-based power and coverage for the two-group logistic
## regression example in wald_lrt_examples.qmd (not currently used there;
## the document uses exact calculations instead). Runs glmmTMB fits for
## nsim simulated data sets at each true log-odds ratio, for n = 50 and
## n = 100 per group, and plots power.

library(glmmTMB)
library(tidyverse)
library(futurize)

## avoid overcommitting threads
library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)

## adjust this to something that makes sense on your machine
plan(multicore, workers = min(28, parallel::detectCores() - 1))

theme_set(theme_bw())
zmargin <- theme(panel.spacing = grid::unit(0, "pt"))
hstrip <- theme(strip.text.y = element_text(angle = 0))
n_labeller <- labeller(n = \(x) paste0("n = ", x, "\nper group"))
## Okabe-Ito minus black and yellow
oi_cols <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
scale_colour_discrete <- function(...) scale_colour_manual(..., values = oi_cols)

rpt_time <- function(x) { attr(x, "time")[["elapsed"]] |> round() }

## simulate two groups of size n; baseline prob p0, log-odds ratio beta
simfun_b <- function(beta, n = 50, p0 = 0.2) {
  g <- factor(rep(0:1, each = n))
  eta <- qlogis(p0) + beta * (g == "1")
  data.frame(g, y = rbinom(2 * n, size = 1, prob = plogis(eta)))
}

## Wald and LRT p-values; Wald and profile CIs for the log-odds ratio.
## The profile CI uses root-finding (TMB::tmbroot) over the absolute
## range [-bmax, bmax]: the default search range is +/- 7 Wald SEs,
## which is useless under (quasi-)separation, where the SE is huge.
## tmbroot's parm.range is a displacement *from the estimate* (unlike
## tmbprofile's, which is absolute), hence the subtraction. A bound with
## no root in the range (NA) means the profile interval is unbounded on
## that side, so it is set to +/- Inf
sumfun_b <- function(data, beta, bmax = 20) {
  m1 <- suppressWarnings(glmmTMB(y ~ g, family = binomial, data = data))
  m0 <- glmmTMB(y ~ 1, family = binomial, data = data)
  cc <- summary(m1)$coefficients$cond["g1", ]
  ## at extreme separation the Hessian can be non-positive-definite, so
  ## the Wald SE is NaN: treat this as the limiting (infinite-SE) case
  bad_se <- !is.finite(cc[["Std. Error"]])
  if (bad_se) cc[c("Std. Error", "Pr(>|z|)")] <- c(Inf, 1)
  wald_ci <- cc[["Estimate"]] + c(-1, 1) * qnorm(0.975) * cc[["Std. Error"]]
  prof_ci <- suppressWarnings(
    confint(m1, parm = "g1", method = "uniroot", estimate = FALSE,
            parm.range = c(-bmax, bmax) - cc[["Estimate"]]))[1, ]
  prof_ci <- ifelse(is.na(prof_ci), c(-Inf, Inf), prof_ci)
  tibble(method = c("Wald", "LRT/profile"),
         p = c(cc[["Pr(>|z|)"]],
               pchisq(2 * (logLik(m1) - logLik(m0)), df = 1, lower.tail = FALSE)),
         lwr = c(wald_ci[1], prof_ci[1]),
         upr = c(wald_ci[2], prof_ci[2]),
         miss_below = beta < lwr,  ## CI entirely above the true value
         miss_above = beta > upr,  ## CI entirely below the true value
         wald_bad_se = c(bad_se, NA))
}

binom_sim <- function(beta, nsim = 2000, n = 50, p0 = 0.2) {
  replicate(nsim, sumfun_b(simfun_b(beta, n = n, p0 = p0), beta),
            simplify = FALSE) |>
    futurize() |>
    bind_rows() |>
    mutate(beta = beta)
}

## simulations across a range of log-odds ratios for group size n
binom_sims <- function(n = 50, betavec = 0:6, seed = 20260927, ...) {
  set.seed(seed)
  tt <- system.time(
    res <- map_dfr(betavec, binom_sim, n = n, ...) |>
      mutate(n = n)
  )
  attr(res, "time") <- tt
  res
}

## binom_sims() sets the seed (default 20260927) before each run;
## futurize() gives the parallel replicates independent RNG streams
res_b50 <- binom_sims(n = 50)
res_b100 <- binom_sims(n = 100)
cat("simulation times (s): n = 50:", rpt_time(res_b50),
    "; n = 100:", rpt_time(res_b100), "\n")

res_b <- bind_rows(res_b50, res_b100)
res_b_sum <- res_b |>
  group_by(n, beta, method) |>
  summarise(n_na = sum(is.na(p)),
            n_bad_se = sum(wald_bad_se),
            reject = mean(p < 0.05, na.rm = TRUE),
            coverage = mean(!(miss_below | miss_above), na.rm = TRUE),
            miss_below = mean(miss_below, na.rm = TRUE),
            miss_above = mean(miss_above, na.rm = TRUE),
            .groups = "drop")

## power (at beta = 0, this is the type I error rate)
p_power <- ggplot(res_b_sum, aes(beta, reject, colour = method)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0.05, linetype = 2) +
  facet_grid(n ~ ., labeller = n_labeller) + zmargin + hstrip +
  labs(x = "true log-odds ratio", y = "proportion rejected (p < 0.05)")
print(p_power)
