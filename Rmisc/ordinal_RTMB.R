library(ordinal)
dd <- expand.grid(A=factor(0:1), B = factor(0:1), id = 1:300)
dd$y <- ordered(sample(1:4, size = nrow(dd), replace = TRUE))

fit <- clmm(y ~ A*B + (A*B|id), data = dd,
            control = clmm.control(trace = 1))

try(clmm(y ~ A*B + (1|id/(A*B)), data = dd,
         control = clmm.control(trace = 1)))

try(clmm(y ~ A*B + (1|id/(A:B)), data = dd,
         control = clmm.control(trace = 1)))

dd$AB <- with(dd, paste0(A,B))
try(clmm(y ~ A*B + (1|id/(A:B)), data = dd,
         control = clmm.control(trace = 1)))

dd$id_AB <- with(dd, paste0(id,A,B))

fit2 <- clmm(y ~ A*B + (1|id) + (1|id_AB), data = dd,
             control = clmm.control(trace = 1))


library(RTMB)
library(reformulas)
library(Matrix)
ord_nll <- function(par) {
    getAll(par, tmbdata)
    nobs <- nrow(X)
    theta <- cumsum(exp(beta0))  ## length (J-1)
    eta <- drop(X %*% beta)      ## length n
    ## plogis() gets in trouble because of RTMB magic
    ## sort out outer()?
    gamma0 <- 1/(1+exp(-(-1*outer(eta, theta, FUN = "-")))) ## n x J
    logprob <- log(apply(cbind(0, gamma0, 1), 1, diff))
    -sum(logprob[cbind(y, 1:nobs)])
}

## would loops be better for TMB?
## read stuff in ordinal vignette about parameterization of intercept

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
tmbdata <- tibble::lst(X, Z, y = dd$y)
ord_nll(par0)

ff <- MakeADFun(ord_nll, par0, silent = TRUE)
ff$fn()

MASS::polr(y~A*B, data = dd)
fit <- with(ff, nlminb(par, fn, gr))
