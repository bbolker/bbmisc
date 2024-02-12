## simulate a t-test-like comparison: what is the distribution of outcomes/
## probability of Dushoff
## https://dushoff.github.io/ResearchSandbox/clarStrength.Rout.pdf

## start by determing what 80% power for a SESOI of 1 looks like:
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
    t.test(x[[1]], x[[2]], conf.level = conf.level)$conf.int
}    

catfun <- function(ci, delta=1) {
    m <- mean(ci) ## central point
    ## adjust for symmetry
    ci <- ci*sign(m)
    m <- abs(m)
    lwr <- ci[1]
    upr <- ci[2]
    if (lwr>delta) return("large/positive")
    if (m > delta & lwr>0) return("large?/positive")
    if (lwr>0 & upr < delta) return("small/positive")
    if (lwr<0 & upr < delta) return("small/positive?")
    if (lwr>(-delta) & upr>delta) return("not large&negative")
    if (lwr<(-delta) & upr>delta) return("unclear")
}

levs <- c("large/positive",
          "large?/positive",
          "small/positive",
          "small/positive?",
          "not large&negative",
          "unclear")

set.seed(101)
system.time(test1 <- replicate(10000, catfun(simfun(n=17))))
test1 <- factor(test1, levels = levs)
props <- table(test1) |> prop.table()
power <- sum(props[1:3])  ## why only 0.5, not 0.8?
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
