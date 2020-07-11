library(tidyverse)

## all standardized
n <- 26
set.seed(101)
ngauss <- 4             ## number of Gaussian reps
dfvec <- c(2,5,10,20)   ## t-distribution df
shapevec <- c(1,2,5,10) ## Gamma shape parameters
res0 <- c(replicate(ngauss, rnorm(n), simplify=FALSE),
         lapply(shapevec, rgamma, n=n, scale=1),
         lapply(dfvec, rt, n=n))

## Shapiro-Wilk tests (attach to results below)
sw_fun <- function(x) unlist(shapiro.test(x)[c("statistic","p.value")])
sw_pvals <- do.call(rbind,lapply(res0,sw_fun))

## names
nm0 <- c(sprintf("gauss(rep=%d)",seq(ngauss)),
         sprintf("Gamma(shape=%d)",shapevec),sprintf("t(df=%d)",dfvec))

## scramble!
oo <- sample(length(res0))

## scrambled versions
nm <- nm0[oo]
res <- res0[oo]

## generate plots (coded A-L in scrambled order)
ff <- sprintf("qq_%s.png",LETTERS[seq_along(res)])
for (i in seq_along(res)) {
    a <- LETTERS[i]
    png(ff[[i]])
    par(las=1,bty="l")
    qqnorm(res[[i]],main=a)
    qqline(res[[i]])
    dev.off()
}

## construct combined data frame for ggplot
ord_df <- tibble(model=nm0,
                 code=LETTERS[order(oo)],  ## scrambled order
                 num=seq_along(nm0))
names(res0) <- nm0
dd <- (purrr:::map_dfr(res0,~tibble(x=.),.id="model")
    %>% full_join(ord_df,by="model")
    %>% arrange(num)
    %>% mutate_at("model",~forcats::fct_inorder(factor(.)))
)
saveRDS(dd, file="qq_data.rds")

zip("qq.zip",files=ff)

## compute type-1 error and power for t-tests and F-tests
## run many simulations with the each of the argument sets specified above
## (running for all 4 Gauss reps is a waste of time, but keeps things parallel)

## type-1: partition an iid sample into two halves and test them against each other
## power : partition sample, displace one of them far enough to provide 80% power in
##         the Gaussian case

## setup for estimating appropriate shift/ratio to get 80% power for t- and F-tests
pdel0 <- power.t.test(n=13,sd=1,power=0.8)$delta
prat0 <- 2.3  ## F-test power ~ 0.8, by trial and error
sumfun <- function(rfun=rnorm,args=list(), n=26, nsim=10000, pdel=pdel0, prat=prat0,
                   verbose=FALSE) {
    if (verbose) cat(".")
    x <- replicate(nsim,scale(do.call(rfun,c(args,list(n=n)))),simplify=FALSE)
    t_fun <- function(z,effect=0) t.test(z[1:(n/2)],effect+z[(n/2+1):n], var.equal=TRUE)$p.value
    t_1 <- mean(vapply(x,t_fun,effect=0, FUN.VALUE=numeric(1))<0.05)
    if (verbose) cat(".")
    t_pow <- mean(vapply(x,t_fun,effect=pdel,FUN.VALUE=numeric(1))<0.05)
    if (verbose) cat(".")
    F_fun <- function(z, ratio=1) var.test(z[1:(n/2)],ratio*z[(n/2+1):n])$p.value
    F_1 <- mean(vapply(x,F_fun,ratio=1, FUN.VALUE=numeric(1))<0.05)
    if (verbose) cat(".")
    F_pow <- mean(vapply(x,F_fun,ratio=prat, FUN.VALUE=numeric(1))<0.05)
    if (verbose) cat(".")
    sw_vals <- lapply(x,sw_fun)
    if (verbose) cat(".")
    W_avg <- mean(vapply(sw_vals,"[",1,FUN.VALUE=numeric(1)))
    p_avg <- mean(vapply(sw_vals,"[",2,FUN.VALUE=numeric(1)))
    if (verbose) cat(".\n")
    return(c(t_1,t_pow,F_1,F_pow,W_avg,p_avg))
}

## in original order
set.seed(101)
answers <- c(replicate(ngauss,sumfun(verbose=TRUE),simplify=FALSE),
  lapply(dfvec,function(df) sumfun(rfun=rt, args=list(df=df),verbose=TRUE)),
  lapply(shapevec,function(shape) sumfun(rfun=rgamma, args=list(scale=1,shape=shape),
                                         verbose=TRUE))
  )

ans0 <- as.data.frame(do.call(rbind,answers))
dimnames(ans0) <- list(nm0, c("type1_t","power_t","type1_F","power_F",
                              "avg_SW_W","avg_SW_pval"))
ans0$code <- LETTERS[seq_along(res)][order(oo)]  ## codes (scrambled!)
## include Shapiro-Wilk stat and p-value for individual samples
ans0 <- data.frame(ans0,sw_pvals,check.names=FALSE)
saveRDS(ans0,file="qq_ans.rds")


### JUNK
## investigating type-1 error from t-statistic (obsolete)
if (FALSE) {
    simfun <- function(n=50,dist="rt", dist.args=list(df=10)) {
        dd <- data.frame(f=factor(rep(1:2,n)),
                         x=do.call(dist,c(list(2*n), dist.args)))
        t.test(x~f, data=dd)$p.value
    }

    set.seed(101)
    dfvec <- 2:20
    type1_t <- sapply(dfvec,
                      function(df) mean(replicate(1000,simfun(n=10, dist.args=list(df=df)))<0.05))
    plot(dfvec,type1_t)


    set.seed(101)
    shapevec <- 2:20
    type1_gamma <- sapply(shapevec,
                          function(shape) mean(replicate(1000,
                                                         simfun(n=10, dist="rgamma", list(scale=5, shape=shape)))<0.05))
    plot(shapevec,type1_gamma)
}

