## goal: experiment with ordinal models in RTMB, in particular to
## fit complex-intercept models (sensu Scandola and Tidoni 2024) with
## desired form

library(ordinal)
library(broom)
library(broom.mixed)
library(RTMB)
library(reformulas)
library(Matrix)
library(purrr)
library(dplyr)

dd <- expand.grid(A=factor(0:1), B = factor(0:1), id = 1:300)
## FIXME: simulate from a real model, rather than randomly?
dd$y <- ordered(sample(1:4, size = nrow(dd), replace = TRUE))

## 7 seconds
system.time(
    fit_clmm_RE <- clmm(y ~ A*B + (A*B|id), data = dd)
)


## these all fail
try(clmm(y ~ A*B + (1|id/(A*B)), data = dd,
         control = clmm.control(trace = 1)))

try(clmm(y ~ A*B + (1|id/(A:B)), data = dd,
         control = clmm.control(trace = 1)))

dd$AB <- with(dd, paste0(A,B))
try(clmm(y ~ A*B + (1|id/(A:B)), data = dd,
         control = clmm.control(trace = 1)))

dd$id_AB <- with(dd, paste0(id,A,B))

## can do it if we spell it out
system.time(
    fit2 <- clmm(y ~ A*B + (1|id) + (1|id_AB), data = dd)
)

ord_nll0 <- function(par) {
    getAll(par, tmbdata)
    nobs <- nrow(X)
    theta <- cumsum(c(beta0[1], exp(beta0[-1])))  ## length (J-1)
    eta <- drop(X %*% beta)      ## length n
    ## plogis() gets in trouble because of RTMB magic
    ## sort out outer()?
    gamma0 <- 1/(1+exp(-(-1*outer(eta, theta, FUN = "-")))) ## n x J
    logprob <- log(apply(cbind(0, gamma0, 1), 1, diff))
    -sum(logprob[cbind(y, 1:nobs)])
}

## * would loops be better for TMB? this may be overvectorizing
## read stuff in ordinal vignette about parameterization of intercept
## * if implementing this in glmmTMB, how are we going to organize eta vs theta?

form <- y ~ A*B + ((A*B)|id)
resp <- dd$y
mf <- model.frame(subbars(form), data = dd)
X <- model.matrix(~A*B, data = dd)[,-1] ## drop intercept
rt <- mkReTrms(findbars(form), mf)
Z <- t(rt$Zt)
image(rt$Lambdat[1:20,1:20])
par0 <- list(
    ## len(intercept) = J-1
    beta0 = rep(0, length(levels(dd$y))-1)
    , beta = rep(0, ncol(X))
    ## ,b = rep(0, ncol(Z))
)
tmbdata <- tibble::lst(X, y = dd$y)
ord_nll0(par0)

ff <- MakeADFun(ord_nll0, par0, silent = TRUE)
class(ff) <- "TMB"
ff$fn()

fit_polr <-MASS::polr(y~A*B, data = dd)
fit_clm <- ordinal::clm(y~A*B, data = dd)
fit <- with(ff, nlminb(par, fn, gr))

glance <- function(x, ...) {
    broom::glance(x, ...) |>
        mutate(across(logLik, c))
}

cmp_models <- function(...) {
    tibble::lst(...) |>
        purrr::map_dfr(glance, .id = "model") |>
        select(model, AIC) |>
        mutate(across(AIC, ~ . - min(.))) |>
        tidyr::pivot_wider(names_from = model, values_from = AIC)
}

cmp_models(fit_polr, fit_clm, fit_RTMB = ff)

blksize <- length(rt$cnms$id)
tmbdata_RE <- c(tmbdata, list(Z = Z, blksize = blksize))
par0_RE <- c(par0, list(vcpars = rep(0, blksize*(blksize+1)/2), b = rep(0, ncol(Z))))
             
## same as above, but with random effects
ord_nll_re <- function(par) {
    getAll(par, tmbdata_RE)
    nobs <- nrow(X)
    theta <- cumsum(c(beta0[1], exp(beta0[-1])))  ## length (J-1)
    eta <- drop(X %*% beta + Z %*% b)      ## length n
    ## plogis() gets in trouble because of RTMB magic
    ## sort out outer()?
    gamma0 <- 1/(1+exp(-(-1*outer(eta, theta, FUN = "-")))) ## n x J
    logprob <- log(apply(cbind(0, gamma0, 1), 1, diff))
    nll <- -sum(logprob[cbind(y, 1:nobs)])
    ## random effects
    us <- unstructured(blksize)
    cc <- us$corr(vcpars[-(1:blksize)])
    sdvec <- exp(vcpars[1:blksize])
    ## note misleading error message in RTMB:::dscale
    ## if (length(scale) != nc) stop("Vector 'scale' must be compatible with *rows* of 'x'")
    nllpen <- -sum(dmvnorm(t(matrix(b, nrow = blksize)), Sigma = cc, scale = sdvec, log = TRUE))
    nll + nllpen
}
ord_nll_re(par0_RE)
ff_RE <- MakeADFun(ord_nll_re, par0_RE, silent = TRUE, random = "b")
class(ff_RE) <- "TMB"
ff_RE$fn()
fit_RE_nlminb <- with(ff_RE, nlminb(par, fn, gr,
                                    control = list(eval.max = 1e4, iter.max = 1e4)))
## annoying that pars are stored internally ...
ff_RE2 <- MakeADFun(ord_nll_re, par0_RE, silent = TRUE, random = "b")
class(ff_RE2) <- "TMB"
fit_RE_optim <- with(ff_RE2, optim(par, fn, gr, method = "BFGS"))

## slight difference in AIC, but not too bad?
cmp_models(fit_clmm_RE, fit_RTMB_RE = ff_RE, fit_RTMB_RE_optim = ff_RE2)

tt1 <- tidy(fit_clmm_RE)
tt2 <- tidy(ff_RE) |> filter(grepl("^beta", term))

## now try compound-symmetric model ...

## same as above, but with random effects (only one term)
ord_nll_re_cs <- function(par) {
    getAll(par, tmbdata_RE)
    nobs <- nrow(X)
    theta <- cumsum(c(beta0[1], exp(beta0[-1])))  ## length (J-1)
    ## random effects:
    ##  do these *first* so $simulate() can work
    ## may need to replace plogis() for simulation?
    sdval <- exp(vcpars[1])
    corrval <- plogis(vcpars[2]) ## unneeded
    cc <- matrix(plogis(logit_corr), nrow = blksize, ncol = blksize)
    diag(cc) <- 1
    nllpen <- -sum(dmvnorm(t(matrix(b, nrow = blksize)), Sigma = cc, scale = sdval, log = TRUE))
    eta <- drop(X %*% beta + Z %*% b)      ## length n
    gamma0 <- 1/(1+exp(-(-1*outer(eta, theta, FUN = "-")))) ## n x J
    logprob <- log(apply(cbind(0, gamma0, 1), 1, diff))
    nll <- -sum(logprob[cbind(y, 1:nobs)])

    nll + nllpen
}
par0_RE_cs <- c(par0, list(vcpars = rep(0,2), b = rep(0, ncol(Z))))

form2 <- y ~ A*B + ((0+A:B)|id)
form3 <- y ~ A*B + (1|id/(A:B))

resp <- dd$y
mf <- model.frame(subbars(form), data = dd)
X <- model.matrix(~A*B, data = dd)[,-1] ## drop intercept
rt <- mkReTrms(findbars(form), mf)
Z <- t(rt$Zt)
