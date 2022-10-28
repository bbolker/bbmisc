## https://cran.r-project.org/web/packages/gnorm/vignettes/gnormUse.html
#' @param lwr lower (loose) bound
#' @param upr upper (loose) bound
#' @param tail_prob prob outside {lwr, upr}
#' @param ctr_prob probability of being in the middle 50\%  of {lwr, upr}
#' @examples
#' get_gnorm()  ## default; gaussian with sd=1
#' get_gnorm(tail_prob=0.01, ctr_prob=0.6)
#' ## set prior approx between 1 and 1000 (log scale)
#' (p <- get_gnorm(lwr=log(1), upr=log(1000), tail_prob=0.01, ctr_prob=0.55))
#' fx <- function(x) do.call("dgnorm", c(list(exp(x)), as.list(p)))
#' curve(fx, from=log(0.5), to=log(1200), n=501)
#' abline(v=c(0,log(1000)), lty=2)
#' brks <- outer(c(1,2,5),c(1,10,100))
#' axis(side=3, at=log(brks), labels=brks)
get_gnorm <- function(lwr=-1, upr=1, tail_prob=2*pnorm(lwr),
                      ctr_prob=abs(diff(pnorm(c(-1,1)*lwr/2)))) {
  require("gnorm")
  ## default tail_prob/ctr_prob assume lwr/upr symmetric around 0 ...
  ## start from Gaussian
  ## desired alpha
  sd <- abs(upr-lwr)/(-2*qnorm(tail_prob/2))
  ## convert to sd (?pgnorm)
  ## conversion factor: sqrt(1/(gamma(3/2)/(gamma(1/2))))
  alpha <- sd*sqrt(2)
  mu <- (lwr+upr)/2 ## symmetric, we don't have to estimate this
  start <- c(alpha=alpha, beta=2)
  tfun <- function(x) {
    ## compute probability within range
    pfun <- function(r) abs(diff(vapply(r,
         function(z) do.call("pgnorm",c(list(z, mu=mu), as.list(x))),
         FUN.VALUE=numeric(1))))
    tail_obs <- 1-pfun(c(upr,lwr))
    ctr_range <- c((mu+lwr)/2, (mu+upr)/2)
    ctr_obs <- pfun(ctr_range)
    return((tail_prob-tail_obs)^2 + (ctr_prob-ctr_obs)^2)
  }
  return(c(mu=mu,optim(par=start,fn=tfun)$par))
}

plot_gnorm <- function(..., add = FALSE, xlim = NULL, ylim = NULL,
                       lcol = 1, fill = NULL) {
    p <- get_gnorm(...)
    L <- list(...)
    fx <- function(x) do.call("dgnorm", c(list(x, as.list(p))))
    if (!add) {
        if (is.null(xlim)) xlim <- do.call("qgnorm",
                                           c(list(c(L$tail_prob/2,
                                           (1-L$tail_prob/2))),
                                           as.list(p)))
        cc <-curve(fx, from = xlim[1], to = xlim[2])
    } else {
        cc <- curve(fx, add = TRUE)
    }
}
              
#' curve(fx, from=log(0.5), to=log(1200), n=501)
#' abline(v=c(0,log(1000)), lty=2)

}

if (FALSE) {
  png("dgn.png")
  curve(do.call("dgnorm", c(list(x), as.list(res))), from=-2, to=2, ylab="")
  dev.off()
}

