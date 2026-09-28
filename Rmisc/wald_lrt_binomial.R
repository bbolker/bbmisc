## setup
library(glmmTMB)
library(car)
suppressPackageStartupMessages(library(tidyverse))
library(futurize)

## avoid overcommitting threads
library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)

## original example from Christoph Scherber
m2 <- glmmTMB(count ~ spp + mined + (1|site),
   zi=~spp + mined,
   family=nbinom2, data=Salamanders)
drop1(m2,test="Chisq")
Anova(m2)

## extract drop1 (LRT) and Anova (Wald) p-values
get_p_values <- function(model, w = 1) {
  d <- drop1(model, test="Chisq")
  p_drop <- d[["Pr(>Chi)"]][[w+1]]
  a <- Anova(model)
  p_anova <- a[["Pr(>Chisq)"]][[w]]
  tibble(p_drop, p_anova)
}
## test
get_p_values(m2)

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

set.seed(101)
sumfun()

## adjust this to something that makes sense on your machine
plan(multicore, workers = min(28, parallel::detectCores() -1 ))

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
  
## 136 seconds
res <- full_sim()
colMeans(res < 0.05)
attr(res, "time")

## sumfun() found very little difference in p-values for this
## particular case, so there may be little difference in type-1 error
## here; try to find a more extreme case for comparison?

## consider coverage (for a non-null case)
## instead of/in addition to type 1 error ?

res_small <- full_sim(n=100)

## ---- binomial example: Wald vs LRT/profile away from the null ----
##
## Two-group logistic regression (baseline probability 0.2, 50 per
## group); the log-odds ratio (beta) ranges from 0 (null) to large
## values, where the second group's probability approaches 1 and the
## log-likelihood becomes strongly asymmetric (Hauck-Donner territory).
## With a sparser baseline (e.g. 0.05) the LRT itself becomes
## anticonservative under the null, muddying the comparison

theme_set(theme_bw())
zmargin <- theme(panel.spacing = grid::unit(0, "pt"))
## Okabe-Ito minus black and yellow
oi_cols <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
scale_colour_discrete <- function(...) scale_colour_manual(..., values = oi_cols)

## Hauck-Donner: fix the data in group 0 (10/50 successes) and increase
## the number of successes in group 1. The LRT statistic keeps
## increasing, while the Wald statistic peaks and then falls off
## (to ~0 at complete separation, k = 50)
hd_stats <- function(k, y0 = 10, n = 50) {
  dd <- data.frame(g = factor(0:1), y = c(y0, k), n = n)
  m1 <- glm(cbind(y, n - y) ~ g, family = binomial, data = dd)
  m0 <- update(m1, . ~ 1)
  tibble(k,
         Wald = coef(summary(m1))["g1", "z value"]^2,
         LRT = m0$deviance - m1$deviance)
}

hd_res <- map_dfr(11:50, hd_stats)
hd_long <- pivot_longer(hd_res, -k, names_to = "test", values_to = "statistic")
gg_hd <- ggplot(hd_long, aes(k, statistic, colour = test)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = qchisq(0.95, 1), linetype = 2) +
  labs(x = "successes in group 1 (out of 50; group 0 = 10/50)",
       y = expression(chi^2 ~ "statistic"))
print(gg_hd)

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

set.seed(20260927)
betavec <- 0:6
tt_b <- system.time(
  res_b <- map_dfr(betavec, binom_sim)
)
tt_b

res_b_sum <- res_b |>
  group_by(beta, method) |>
  summarise(n_na = sum(is.na(p)),  ## failed fits
            n_bad_se = sum(wald_bad_se),  ## NaN Wald SEs (counted as p = 1)
            reject = mean(p < 0.05, na.rm = TRUE),
            coverage = mean(!(miss_below | miss_above), na.rm = TRUE),
            miss_below = mean(miss_below, na.rm = TRUE),
            miss_above = mean(miss_above, na.rm = TRUE),
            .groups = "drop")
print(res_b_sum, n = Inf)

## power (at beta = 0, this is the type I error rate)
gg_power <- ggplot(res_b_sum, aes(beta, reject, colour = method)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0.05, linetype = 2) +
  labs(x = "true log-odds ratio", y = "proportion rejected (p < 0.05)")
print(gg_power)

## non-coverage split by tail: a well-calibrated 95% CI should miss
## ~2.5% of the time on each side
res_b_tails <- res_b_sum |>
  select(beta, method, miss_below, miss_above) |>
  pivot_longer(starts_with("miss"), names_to = "tail", values_to = "prop") |>
  mutate(tail = factor(tail, levels = c("miss_below", "miss_above"),
                       labels = c("CI above true value", "CI below true value")))
gg_tails <- ggplot(res_b_tails, aes(beta, prop, colour = method)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0.025, linetype = 2) +
  facet_wrap(~ tail) + zmargin +
  labs(x = "true log-odds ratio", y = "proportion of CIs missing the true value")
print(gg_tails)
