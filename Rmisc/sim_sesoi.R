## mostly trying to avoid tidyverse but ...
library(ggplot2); theme_set(theme_bw() + theme(panel.spacing=grid::unit(0, "lines")))
                 
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
    x <- rnorm(2*n, mean = rep(c(0,delta), each =n), sd = sd)
    tt <- t.test(x[1:n], x[-(1:n)], conf.level = conf.level, var.equal = TRUE)
    with(tt, c(est = unname(-1*diff(estimate)),
                        lwr = conf.int[1], upr = conf.int[2]))
}    

levs <- c("large/clear sign",
          "unclear magnitude/clear sign",
          "small/clear sign",
          "small/unclear sign",
          "NOT (large and opposite est)",
          "unclear")

catfun <- function(x, s=1) {
    lwr <- x[2]
    upr <- x[3]
    ## adjust for symmetry? ci <- ci*sign(m); m <- abs(m)
    ## case-when?
    if (lwr>s || upr<(-s)) return(levs[1])
    if ((upr>s && lwr>0 && lwr<s) || (lwr<(-s) && upr<0 && upr>(-s))) return(levs[2])
    if ((lwr>0 && upr<s) || (upr<0 && lwr>(-s))) return(levs[3])
    if ((lwr>(-s) && lwr<0 && upr>0 && upr<s) ||
        (upr<s   && upr>0 && lwr>0 && upr>(-s))) return(levs[4])
    if ((lwr<0 && lwr>(-s) && upr>s) || (upr>0 && upr<s && lwr<(-s))) return(levs[5])
    if (lwr<(-s) && upr>s) return(levs[6])
}


set.seed(101)
system.time({
    ## using auto-simplify of replicate ...
    dd1 <- as.data.frame(t(replicate(10000, simfun(n=17))))
    dd1$cat <- apply(dd1, 1, catfun) |> factor(levels = levs)
})

## flip sign
## dd1 <- transform(dd1, lwr = sign(est)*lwr, upr = sign(est)*upr,
## est = abs(est))
## ordering for caterpillar plot
## alternately could arrange() and then fct_inorder
dd1$n <- reorder(factor(seq(nrow(dd1))), dd1$est)


print(gg0 <- ggplot(dd1, aes(n, est))
      + geom_pointrange(aes(colour = cat, ymin = lwr, ymax = upr),
                        alpha = 0.5)
      + geom_hline(yintercept = 0)
      + geom_hline(yintercept = c(-1, 1), linetype = 2)
      + scale_x_discrete(labels = NULL, breaks = NULL))

print(gg0 + facet_wrap(~cat))
## plot looks plausible

print(props <- table(dd1$cat) |> prop.table())
print(power <- sum(props[1:3]))

## what if the true effect size is twice the SESOI/powered effect size?
dd1$cat2 <- apply(dd1, 1, catfun, s=0.5) |> factor(levels = levs)
print(props2 <- table(dd1$cat2) |> prop.table())

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

## non-central t based on two different centers??

## our problem has *two* scales (we should be able to manipulate
##  the effect size and the cutoff separately!
## making effect size = SESOI = cutoff is (over?)simplifying the problem

