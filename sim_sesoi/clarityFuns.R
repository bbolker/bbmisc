library(shellpipes)
library(readr)

catTable <- read_csv(file="
	lwr,upr,Short,Mag/Sign
	1,1,LN,large/negative
	1,2,N,unclear/negative
	1,3,xLP,unclear/unclear
	1,4,Z,unclear/unclear
	2,2,SN,small/negative
	2,3,S,small/unclear
	2,4,xLN,unclear/unclear
	3,3,SP,small/positive
	3,4,P,unclear/positive
	4,4,LP,large/positive
")

simfun <- function(n, delta=1, sd=1, conf.level = 0.95, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    x <- rnorm(2*n, mean = rep(c(0,delta), each =n), sd = sd)
    tt <- t.test(x[1:n], x[-(1:n)], conf.level = conf.level, var.equal = TRUE)
    res <- with(tt, c(est = unname(-1*diff(estimate)),
                      lwr = conf.int[1], upr = conf.int[2]))
    return(res)
}    

f_lengthen <- function(x) {
  (x 
    |> as.data.frame()
    |> dplyr::mutate(n = nvec, .before = 1)
    |> tidyr::pivot_longer(col = -n)
    |> dplyr::mutate(name = factor(name, levels = levs))
  )
}

## should be able to do this much faster if we're sticking to equal-sample size, etc. etc. etc.?

## how many cases should we distinguish?
## (1) show the effect is small or large
##    * care less about the sign if it's small?

levs <- c("large/clear sign",
          "unclear magnitude/clear sign",
          "small/clear sign",
          "small/unclear sign",
          "NOT (large and opposite est)",
          "unclear")

#' categorize outcomes
#' @param x a 3-element vector with 'lower' and 'upper' as the second and third elements
#' @param s sesoi (critical value distinguishing small/large effect sizes)
catfun <- function(x, s=1, refs=c(-1, 0, 1), levNames=NULL) {
    lwr <- x[2]
    upr <- x[3]
	 marks <- s*refs
	 return(
	 	c(1, 4, 6, 7)[1+sum(lwr>marks)] + sum(upr>marks)
	)
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

tabfun <- function(nsim, levNames=NULL, ...) {
  res <- lapply(seq.int(nsim), function(i) simfun(...)) |> do.call(what=rbind)
  dd1 <- as.data.frame(res)
  dd1$cat <- apply(dd1, 1, catfun)
  res <- table(dd1$cat) |> prop.table()
  return(res)
}

saveEnvironment()
