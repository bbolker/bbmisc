library(shellpipes)

## It would work to use a .R dependency and say `sourceFiles()`
## But using .rda has always seemed cleaner to me
## See also Makefile for how I did this without editing the funs file

loadEnvironments()
set.seed(101)
nvec <- c(5:10, (2:9)*10, 100, 200)
sim1 <- lapply(nvec, tabfun, delta=0.5, nsim = 10000) |> do.call(what = rbind)

saveEnvironment()
