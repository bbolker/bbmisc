library(glmmTMB)
library(lme4)
library(nlme)

## simulate data corresponding to a two-group t-test with unequal variances
dd <- data.frame(
  f = factor(rep(c("a","b"), each = 200)),
  g = factor(rep(1:20, each = 10))
)

true_params <- list(beta = c(0, 4),
                    theta = 1,
                    betadisp = c(0, 1))
dd$y <- simulate_new(~f + (1|g),
                  dispformula = ~f,
                  newdata = dd,
                  family = gaussian,
                  newparams = true_params)[[1]]

## in glmmTMB with dispformula
m1 <- glmmTMB(y ~ f + (1|g),
              dispformula = ~f,
              family = gaussian,
              data = dd)

## in lmer: use a dummy variable (0/1) based on the
## group, with an observation-level random effects
## this only works if we use the group with the *smallest*
## residual variance as the baseline level
## Have to override some of lmer's sanity checks
dd$obs <- seq(nrow(dd))
m2 <- lmer(y ~ f + (1|g) + (0 + dummy(f, "b")|obs), data = dd,
           control = lmerControl(check.nobs.vs.nlev = "ignore",
                                 check.nobs.vs.nRE = "ignore"),
           REML = FALSE)

## in lme: use varIdent()
m3 <- lme(y ~ f,
          random = ~1|g,
          data = dd,
          weights = varIdent(form = ~1|f),
          method = "ML")
          
## compare effects ...
ae <- function(...) {
  all.equal(..., tolerance = 1e-6, check.attributes = FALSE)
}


ae(sigma(m2), exp(fixef(m1)$disp)[1])
ae(sqrt(sigma(m2)^2 + drop(VarCorr(m2)$obs)), exp(sum(fixef(m1)$disp)))

ae(sigma(m2), sigma(m3))
ae(coef(m3$modelStruct$varStruct), fixef(m1)$disp[2])

## in the non-mixed case, we can compare gls() with a weights= argument
## to t.test(..., var.equal = FALSE)

## see also https://rpubs.com/bbolker/factorvar (maybe -- it's ugly
## and tangential)
