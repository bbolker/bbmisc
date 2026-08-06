library(nnet)
library(RTMB)

## illustrating constructing/fitting a multinomial model with RTMB

## https://en.wikipedia.org/wiki/Compositional_data
## softmax function is the inverse of the center log ratio transform
## (log(x(i)/GM(x)) where GM(x) is the geometric mean of the x(i))
softmax <- function(x) {
    ex <- exp(x)
    ex/sum(ex)
}
softmax0 <- function(x) {
    softmax(c(0,x))
}

## simulate data
##  B is a matrix of coefficients: if we have p parameters (including intercept,
##  i.e. ncol(X)) and k categories, then B is a p x (k-1) matrix (we only need
##  k-1 parameters to define k probabilities
set.seed(101)
k <- 3
ng <- 10
ntot <- 100
dd <- data.frame(x = rnorm(ntot), y = rnorm(ntot), g = factor(rep(1:ng, length.out = ntot)))
X <- model.matrix(~ x + y, data = dd)
np <- ncol(X)
B <- matrix(rnorm(np*(k-1)), nrow = np)
## could use reformulas for a more complicated (e.g. random-slopes) model
Z <- Matrix::sparse.model.matrix(~ g - 1, data = dd)
b <- matrix(rnorm(ng*(k-1)), nrow = ng)
probs <- t(apply(X %*% B + Z %*% b, 1, softmax0))
stopifnot(all(abs(rowSums(probs)-1.0) < 1e-6))

Y <- t(apply(probs, 1, \(p) rmultinom(1, size = 10, prob = p)))
size <- rowSums(Y)

tmb_data <- list(X = X, Y = Y, Z  = Z, size = size)
pars <- list(B = B, b = b, logsd = 0)

fn0 <- function(pars) {
    getAll(pars, tmb_data)
    probs <- t(apply(X %*% B + Z %*% b, 1, softmax))
    liks <- sapply(1:nrow(Z),
                   \(i) lgamma(size[i] + 1) +
                        sum(Z[i,] * log(probs[i,]) - lgamma(Z[i,] + 1)))
    b_lik <- sum(dnorm(b, 0, exp(logsd), log = TRUE))
    -1*(b_lik + sum(liks))
}

fn0(pars)
adfun0 <- MakeADFun(fn0, pars, silent = TRUE)
adfun0$fn()

adfun <- MakeADFun(fn0, pars, random = "b", silent = TRUE)
adfun$fn()

fit <- with(adfun, nlminb(par, fn, gr))
pp <- adfun$env$parList()
Bhat <- matrix(pp$B, nrow = np)
Bhat - B

## hmm, more work to do
