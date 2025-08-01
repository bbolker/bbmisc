<## mostly trying to avoid tidyverse but ...
library(ggplot2); theme_set(theme_bw() + theme(panel.spacing=grid::unit(0, "lines")))
source("sim_sesoi_funs.R")

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


proptest <- function(x, s = 1) {
    lwr <- x[,2]
    upr <- x[,3]
    c(lwr_gt_0 = mean(lwr>0),
      lwr_gt_s = mean(lwr>s),
      lwr_gt_negs = mean(lwr>(-s)),
      upr_gt_0 = mean(upr>0),
      upr_gt_s = mean(upr>s),
      upr_gt_negs = mean(upr>(-s)))
}

set.seed(101)
system.time({
    ## using auto-simplify of replicate ...
    dd1 <- as.data.frame(t(replicate(50000, simfun(n=17))))
    dd1$cat <- apply(dd1, 1, catfun) |> factor(levels = levs)
})


tabfun <- function(..., nsim = 1000) {
  res <- lapply(seq.int(nsim), function(i) simfun(...)) |> do.call(what=rbind)
  dd1 <- as.data.frame(res)
  dd1$cat <- apply(dd1, 1, catfun) |> factor(levels = levs)
  table(dd1$cat) |> prop.table()
}
set.seed(101)
tabfun(n=17)
tabfun(n=100)

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

## example plot: 80% power
## prob (lower CI > 0)
n <- 17; delta <- 1; sd <- 1; s <- 0.5; nu <- 2*n-2; conf.level <- 0.95
curve(dt(x, df = nu, ncp = sqrt(n/2) * delta/sd), from = -2, 10)
lims <- qt((1+conf.level)/2, df = nu)
abline(v = c(-1, 1)*lims, lty = 2)
pt(lims, df = nu, ncp = sqrt(n/2) * delta/sd, lower.tail = FALSE) ## 0.807
## (the other tail is negligible
pt(-lims, df = nu, ncp = sqrt(n/2) * delta/sd, lower.tail = TRUE) ## ~ 1e-6
## prob (lower CI > delta [SESOI])

curve(dt(x, df = 5, ncp = 0), from = -5, to = 10)
abline(v = qt(c(0.025, 0.975), df = 5), lty = 2)

calc_catpower <- function(n, delta=1, sd=1, s=1, conf.level = 0.95, debug = FALSE,
                          retval = c("catprops", "indprops")) {
    retval <- match.arg(retval)
    nu <- 2*n-2
    ## three comparisons: ncp of (delta, delta-SESOI, delta+SESOI)
    lims <- qt((1+conf.level)/2, df = nu)

    prob_lwr_gt_0 <- pt(lims, df = nu, ncp = sqrt(n/2) * delta/sd, lower.tail = FALSE)
    ## equal to value for -lims, -delta, lower tail (upr < 0 for negative delta):
    ##  pt(-lims, df = nu, ncp = sqrt(n/2) * -delta/sd, lower.tail = TRUE)
    prob_upr_lt_0 <- pt(-lims, df = nu, ncp = sqrt(n/2) * delta/sd, lower.tail = TRUE)
    prob_lwr_gt_s <- pt(lims, df = nu, ncp = sqrt(n/2) * (delta-s)/sd, lower.tail = FALSE)
    prob_lwr_gt_s <- pt(lims, df = nu, ncp = sqrt(n/2) * (delta-s)/sd, lower.tail = FALSE)
    prob_lwr_gt_negs <- pt(lims, df = nu, ncp = sqrt(n/2) * (delta+s)/sd, lower.tail = FALSE)
    prob_upr_gt_0 <- pt(-lims, df = nu, ncp = sqrt(n/2) * delta/sd, lower.tail = FALSE)
    prob_upr_gt_s <- pt(-lims, df = nu, ncp = sqrt(n/2) * (delta-s)/sd, lower.tail = FALSE)
    ## suppress full-precision warning? (1-3.430589e-14)
    prob_upr_gt_negs <- pt(-lims, df = nu, ncp = sqrt(n/2) * (delta+s)/sd, lower.tail = FALSE)
    if (retval == "indprops") {
        return(c(lwr_gt_0 = prob_lwr_gt_0, lwr_gt_s = prob_lwr_gt_s, lwr_gt_negs = prob_lwr_gt_negs,
                 upr_gt_0 = prob_upr_gt_0, upr_gt_s = prob_upr_gt_s, upr_gt_negs = prob_upr_gt_negs))
    }
    res <- c(prob_lwr_gt_s,
             prob_upr_gt_s * prob_lwr_gt_0 * (1-prob_lwr_gt_s),
             prob_lwr_gt_0 * (1-prob_upr_gt_s),  ## small/clear sign
             prob_lwr_gt_negs * (1-prob_lwr_gt_0) * prob_upr_gt_0 * (1-prob_upr_gt_s),
             (1-prob_lwr_gt_0) * prob_lwr_gt_negs * prob_upr_gt_s,
             (1-prob_lwr_gt_negs) * prob_upr_gt_s) |> setNames(levs)
    return(res)
}

(ii <- calc_catpower(n = 17, retval = "indprops"))
## flip because effect is actually negative in simulations
## is flipping sign always correct? I'm not sure.
proptest(transform(dd1, lwr = -upr, upr = -lwr))
## these probabilities look correct, so the discrepancies must (???) be in the way they are combined ...

cc <- calc_catpower(n=17)
cbind(cc, props)
diffs <- cc-props
diffs[abs(diffs)>1e-3]
sum(cc)  ## too large?? by 1.0054
## errors are of order 0.02 (opposite tails?)
## analytical answer is too large for "small/clear sign" (element 3), "NOT (large and opposite)" (element 5)
##  (can't be a two-tailed problem?)

## small/clear sign:
## ((lwr>0 && upr<s) || (upr<0 && lwr>(-s)))
with(dd1, mean(upr<0))
ii[["lwr_gt_0"]]
with(dd1, mean(lwr>(-s)))  ## NOT INDEPENDENT!! Ugh.
chisq.test(with(dd1, table(upr<0, lwr>(-s))))
1-ii[["upr_gt_s"]]  ## 0.025 vs 0.0028?
## 0.975 vs 0.99972
ii[["upr_gt_s"]]
with(dd1, mean(lwr<(-s)))
with(dd1, mean((upr<0 & lwr>(-s)) | (lwr>0 & upr<s)))

library(mvtnorm)
conf.level <- .95; nu <- 32; sd <- 1; s <- 1
lims <- qt((1+conf.level)/2, df = nu)
## making corr nearly perfect -- is this the right way? what 'type' do I want?
pmvt(lower=c(lims,-Inf),
     upper=c(Inf,lims),
     df = nu,
     delta = c(sqrt(n/2)*delta/sd, sqrt(n/2)*(delta-s)/sd), corr = matrix(c(1, 0.999, 0.999, 1), 2))

pmvt(lower=lims, upper = Inf,
     df = nu,
     delta = sqrt(n/2)*delta/sd)

## the key point is that we want the same value of X in both cases (I think). This paper seems geared toward what we want:
## 
## Owen, D. B. 1965. “A Special Case of a Bivariate Non-Central $t$-Distribution.” Biometrika 52 (3/4): 437–46. https://doi.org/10.2307/2333696.

## I think it's easier to look at
## Julious, Steven A. 2004. “Sample Sizes for Clinical Trials with Normal Data.” Statistics in Medicine 23 (12): 1921–86. https://doi.org/10.1002/sim.1783.
## which includes power calculations for equivalence tests, especially eq 23:
## 1-beta = Probt(-t_{1-\alpha},n_A(r+1)-2, n_A(r+1)-2, \tau_2)
##        - Probt(-t_{1-\alpha},n_A(r+1)-2, n_A(r+1)-2, \tau_1)
## where tau_1,2 = ((mu_A-mu_B) ± d)\sqrt{r n_A}/(sqrt(r+1) sigma^2)
## n_A is the sample size of one group and n_B = r*n_A
## (r=1 is the balanced case)
