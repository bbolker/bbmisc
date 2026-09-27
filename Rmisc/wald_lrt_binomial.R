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
simfun <- function() {
    dd <- transform(Salamanders,
            count = simulate_new(~ spp + mined + (1|site),
            newdata = Salamanders,
            newparams = list(beta = c(-0.5, rep(0,6), 1.4),
            betadisp = 0.5,
            theta = -1),
            family = nbinom2)[[1]])
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

## 136 seconds
set.seed(101)
system.time(
  res <- replicate(1000, sumfun(), simplify = FALSE) |>
    futurize() |>
    bind_rows()
)

colMeans(res < 0.05)

## sumfun() found very little difference in p-values for this
## particular case, so there may be little difference in type-1 error
## here; try to find a more extreme case for comparison?

## consider coverage (for a non-null case)
## instead of/in addition to type 1 error ?
