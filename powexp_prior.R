## https://cran.r-project.org/web/packages/gnorm/vignettes/gnormUse.html
library(gnorm)
#' @param lwr lower (loose) bound
#' @param upr upper (loose) bound
#' @param tail_prob prob outside {lwr, upr}
#' @param ctr_prob probability of being in the middle 50\%  of {lwr, upr}
get_gnorm <- function(lwr=-1, upr=1, tail_prob=2*pnorm(lwr),
                      ctr_prob=abs(diff(pnorm(c(-1,1)*lwr/2)))) {
  ## default tail_prob/ctr_prob assume lwr/upr symmetric around 0 ...
  ## start from Gaussian
  ## desired alpha
  sd <- abs(upr-lwr)/(-2*qnorm(tail_prob/2))
  ## convert to sd (?pgnorm)
  alpha <- sd/sqrt(gamma(3/2)/(gamma(1/2)))
  mu <- (lwr+upr)/2 ## symmetric, we don't have to estimate this
  start <- c(alpha=alpha, beta=2)
  tfun <- function(x) {
    ## compute probability within range
    pfun <- function(r) abs(diff(vapply(r,
         function(z) do.call("pgnorm",c(list(z, mu=mu), as.list(x))),
         FUN.VALUE=numeric(1))))
    tail_obs <- 1-pfun(c(upr,lwr))
    ctr_range <- mu+c(-1,1)/2*lwr
    ctr_obs <- pfun(ctr_range)
    return((tail_prob-tail_obs)^2 + (ctr_prob-ctr_obs)^2)
  }
  return(optim(par=start,fn=tfun)$par)
}


(res <- get_gnorm(tail_prob=0.01, ctr_prob=0.6))
png("dgn.png")
curve(do.call("dgnorm", c(list(x), as.list(res))), from=-2, to=2, ylab="")
dev.off()

