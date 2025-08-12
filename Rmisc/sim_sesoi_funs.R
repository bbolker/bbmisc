#' two-group sim for equal-var t-test power/outcome calculation
#' @param n number of simulations
#' @param delta difference between means
#' @param standard deviation of observations
#' @param conf.level confidence level
#' @param seed random-number seed
simfun <- function(n, delta=1, sd=1, conf.level = 0.95, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    x <- rnorm(2*n, mean = rep(c(0,delta), each =n), sd = sd)
    tt <- t.test(x[1:n], x[-(1:n)], conf.level = conf.level, var.equal = TRUE)
    with(tt, c(est = unname(-1*diff(estimate)),
                        lwr = conf.int[1], upr = conf.int[2]))
}    

## should be able to do this much faster if we're sticking to equal-sample size, etc. etc. etc.?

## how many cases should we distinguish?
## (1) show the effect is small or large
##    * care less about the sign if it's small?

levs <- c("large/positive",
          "unclear/positive",
          "small/positive",
          "small/unclear",
          "NOT (large and negative)",
          "What")

#' categorize outcomes
#' @param x a 3-element vector with 'lower' and 'upper' as the second and third elements
#' @param s sesoi (critical value distinguishing small/large effect sizes)
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

#' compute fractions in each category
#' @param x a matrix-like object with columns 2 and 3 equal to lower/upper CIs
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

tabfun <- function(..., nsim = 1000) {
  res <- lapply(seq.int(nsim), function(i) simfun(...)) |> do.call(what=rbind)
  dd1 <- as.data.frame(res)
  dd1$cat <- apply(dd1, 1, catfun) |> factor(levels = levs)
  table(dd1$cat) |> prop.table()
}
