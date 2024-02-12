## mostly trying to avoid tidyverse but ...
library(ggplot2); theme_set(theme_bw())
                 


## simulate a t-test-like comparison: what is the distribution of outcomes/
## probability of Dushoff outcomes
## https://dushoff.github.io/ResearchSandbox/clarStrength.Rout.pdf

## start by determing required n to get 80% power for a SESOI of 1 with SD 1
pp <- power.t.test(delta = 1, sd = 1, power = 0.8)
n <- ceiling(pp$n)  ## 17
tt <- power.t.test(delta = 1, sd = 1, n = n)
tt$power ## 0.807

## Q: what parameters are identifiable/jointly confounded in
## classic power analysis? e.g. delta and sd must be exchangeable
## since one is a scale parameter?

## two-group sim for equal-var t-test power/outcome calculation
simfun <- function(n, delta=1, sd=1, conf.level = 0.95, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    x <- replicate(2, rnorm(n, mean = delta, sd = sd), simplify = FALSE)
    tt <- t.test(x[[1]], x[[2]], conf.level = conf.level)
    with(tt, c(est = unname(-1*diff(estimate)),
                        lwr = conf.int[1], upr = conf.int[2]))
}    

levs <- c("large/positive",
          "large?/positive",
          "small/positive",
          "small/positive?",
          "not (large&negative)",
          "unclear")

catfun <- function(x, delta=1) {
    m <- mean(x[2:3])
    ci <- x[2:3]
    ## adjust for symmetry
    ## FIXME: eventually, do this some other way
    ci <- ci*sign(m)
    m <- abs(m)
    lwr <- ci[1]
    upr <- ci[2]
    ## case-when?
    if (lwr>delta) return(levs[1])
    if (m > delta & lwr>0) return(levs[2])
    if (lwr>0 & upr < delta) return(levs[3])
    if (lwr<0 & upr < delta) return(levs[4])
    if (lwr>(-delta) & upr>delta) return(levs[5])
    if (lwr<(-delta) & upr>delta) return(levs[6])
}


set.seed(101)
system.time({
    ## using auto-simplify of replicate ...
    dd1 <- as.data.frame(t(replicate(10000, simfun(n=17))))
    dd1$cat <- apply(dd1, 1, catfun) |> factor(levels = levs)
})

## flip sign
dd1 <- transform(dd1, lwr = sign(est)*lwr, upr = sign(est)*upr,
                 est = abs(est))
## ordering for caterpillar plot
dd1$n <- reorder(factor(seq(nrow(dd1))), dd1$est)

print(ggplot(dd1, aes(n, est))
      + geom_pointrange(aes(colour = cat, ymin = lwr, ymax = upr),
                        alpha = 0.5)
      + scale_x_discrete(labels = NULL, breaks = NULL))

props <- table(dd1$cat) |> prop.table()
print(power <- sum(props[1:3]))  ## why only 0.5, not 0.8?
## (something about two-sidedness?)

## note: this is the body of the t-test power calculation
## where tsample is 1 for one-sample, 2 for two-sample
## tside is 1 for one-sided, 2 for two-sided
## (n is the number of samples *per group*)

## can we generalize this?

## nu <- pmax(1e-07, n - 1) * tsample   ## df
## qu <- qt(sig.level/tside, nu, lower.tail = FALSE)
## pt(qu, nu, ncp = sqrt(n/tsample) * delta/sd, lower.tail = FALSE) + 
##     pt(-qu, nu, ncp = sqrt(n/tsample) * delta/sd, 
##        lower.tail = TRUE)

## our problem has *two* scales (we should be able to manipulate
##  the effect size and the cutoff separately!
## making effect size = SESOI = cutoff is (over?)simplifying the problem
