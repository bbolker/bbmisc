## ----pkgs, message = FALSE----------------------------------------------------
library(lme4)
library(brms)
## options(mc.cores = min(4, parallel::detectCores()-1))
options(mc.cores = 1)
options(brms.backend = "cmdstanr")
if (require("unix")) {
    rlimit_as(14e9)  ## 1e9 ~ 1 GB
}
library(broom.mixed)
library(tidybayes)

prior1 <- prior(normal(0,10), class = b) +
    prior(cauchy(0,2), class = sd)
fit1 <- brm(count ~ zAge + zBase * Trt + (1|patient),
            data = epilepsy,
            family = poisson(), prior = prior1)
     
## generate a summary of the results
summary(fit1)

data("sleepstudy", package = "lme4")
form1 <- Reaction ~ Days + (Days|Subject)

## ----prior5, cache = TRUE, warning = FALSE, dependson = "prior1"--------------
b_prior5 <- c(set_prior("normal(200, 10)", "Intercept"),
              set_prior("normal(0, 8)", "b"),
              set_prior("student_t(10, 0, 3)", "sd"),
              set_prior("student_t(10, 0, 3)", "sigma")
             )

## ----lmer_fit-----------------------------------------------------------------
m_lmer <- lmer(form1, sleepstudy)


## ----brms_fit_reg, cache = TRUE-----------------------------------------------
b_reg <- 
    brm(form1, sleepstudy, prior = b_prior5,
        seed = 101,               ## reproducibility
        ## silent = 2, refresh = 0,  ## be vewy vewy quiet ...
        ## handle divergence
        ## (currently makes memory use balloon ...)
        control = list(adapt_delta = 0.95)
        )


## ----brms_fit_default, cache = TRUE-------------------------------------------
b_default <-
    brm(form1, sleepstudy,
        seed = 101,              ## reproducibility
        ## silent = 2, refresh = 0  ## be vewy vewy quiet ...
        )

save(list=c("b_reg", "b_default", "m_lmer"),
     file = "examples1.rda")

