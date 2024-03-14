## quick and dirty stuff for n-component (1-D) Gaussian mixtures:
## density, CDF, random deviate generator
dmix <- function(x, p, m, s) {
    m <- Map(\(p, m, s) p*dnorm(x, m, s), p, m, s)
    rowSums(do.call(cbind, m))
}

pmix <- function(x, p, m, s) {
    m <- Map(\(p, m, s) p*pnorm(x, m, s), p, m, s)
    rowSums(do.call(cbind, m))
}

rmix <- function(n, p, m, s) {
    n <- drop(rmultinom(1, size =n, prob = p))
    r <- Map(\(n, m, s) rnorm(n, m, s), n, m, s)
    unlist(r)
}

## test/draw some curves
xvec <- seq(-1, 3, length.out = 101)
d1 <- dmix(xvec,
        p = c(0.2, 0.8),
        m = c(1, 2),
        s = c(0.3, 0.3))
plot(xvec, d1, type = "l")

true_pars <- list(p = c(0.2, 0.8),
                  m = c(1, 2),
                  s = c(0.3, 0.3))
## sample
set.seed(101)
r1 <- do.call("rmix",
              c(list(n=1000),
                true_pars))
## get quantiles
qvec <- seq(0, 1, by = 0.1)
qq <- quantile(r1, qvec)

## fit a monotonic spline to the CDF (quantiles are the x-variate)
s1f <- splinefun(qq, qvec, method = "hyman")

s2f <- splinefun(qq, qvec, method = "monoH.FC")

## spline() doesn't work, not sure why we have to jump through this hoop
## if we want values (not used any more though)
## s2 <- list(x = s1$x, y = s2f(s1$x))

## what we have now is a smooth representation of the cumulative
## distribution function.
## we can evaluate the derivatives directly ... (see below)

## fit 2-component Gaussian mixture model

## JD's order statistics stuff would make sense if we wanted to account
##  for error in a known-small sample
## source("https://raw.githubusercontent.com/dushoff/notebook/master/orderStats.R")
## given a large/unknown sample size, least-squares fitting to CDF is
## probably good enough

## fit predicted CDF to quantile values (easier than sorting out
## quantile function for Gaussian mixture)
skel <- list(prob = NA, mean = rep(NA, 2), sd = rep(NA,2))
## relist(1:5, skel)
lsqfun <- function(par) {
    p <- relist(par, skel)
    with(p, sum((pmix(qq, prob, mean, sd) - qvec)^2))
}
lsqfit <- optim(c(0.4, 1,2, 0.4, 0.5), lsqfun)

## results

## cvec <- c("#000000", colorblindr::palette_OkabeIto)

cvec <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
          "#D55E00", "#CC79A7", "#999999")
lfun <- function(pos = "topright") {
    legend(pos, col = cvec,
       bg = "white",
       lty = 1,
       lwd = 2,
       legend = c("true", "mono-spline (hyman)", "mono-spline (monoH.FC)",
                  "mixture model"))
}


par(las=1, bty = "l")
## compare cumulative dist predictions
plot(qq, qvec)
with(true_pars, curve(pmix(x, p, m, s),
      lwd = 2, add = TRUE, col = cvec[1]))
curve(s1f, add = TRUE, col = cvec[2], lwd = 2)
curve(s2f, add = TRUE, col = cvec[3], lwd = 2)
with(relist(lsqfit$par, skel),
     curve(pmix(x, prob, mean, sd), add = TRUE, col= cvec[4], lwd = 2))
lines(ecdf(r1))
abline(v=qq, lty = 2)
lfun(pos = "bottomright")

par(las = 1, bty = "l", yaxs="i")
hist(r1, freq = FALSE, breaks = 40, main = "", ylim = c(0, 1.3))
with(true_pars, curve(dmix(x, p, m, s),
      lwd = 2, add = TRUE, col = cvec[1])) ## true distribution
curve(s1f(x, 1L), add = TRUE, n = 1001, col = cvec[2], lwd =2)
curve(s2f(x, 1L), add = TRUE, n = 1001, col = cvec[3], lwd = 2)
with(relist(lsqfit$par, skel),
     curve(dmix(x, prob, mean, sd), add = TRUE, col= cvec[4], lwd = 2))
abline(v=qq, lty = 2)
lfun("topleft")

## Conclusions:
## *  mixture model is best for this case (naturally, since it matches
## the true model!); however, quantile spacing is such that we pretty
## much miss the lower mode
## * monoH.FC has some odd cusps (i.e. second derivative of the spline
##   is not continuous), but seems better controlled than hyman spline



