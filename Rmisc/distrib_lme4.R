## start by doing everything in serial, without headaches
library(lme4)

nc <- 4
dd <- split(sleepstudy, rep(1:nc, length.out = nrow(sleepstudy)))

mkdevfun <- function(data, form, ...) {
    lmod <- lFormula(form, data, ...)
    do.call(mkLmerDevfun, lmod)
}

form <- Reaction ~ Days + (Days | Subject)
clusterExport(cl, c("mkdevfun", "dd", "form"))
invisible(clusterEvalQ(cl, library(lme4)))

## https://stackoverflow.com/questions/61383366/parallelisation-within-a-function-in-r-avoiding-automatic-export-of-objects-t



## doing this in parallel is sort of pointless
## again, ideally would want to run mkdevfun() separately on
## each worker, applied to the worker's own data shard
devfuns <- parLapply(cl, dd, mkdevfun, form = form)

comb_devfun <- function(theta) {
    force(theta)
    clusterExport(cl, "theta", envir = environment())
    parLapply(cl, devfuns, function(f) f(theta))
}

try(comb_devfun(c(1, 1, 1, 1)))


## sudo apt install openmpi-common openmpi-bin libopenmpi-dev

## mpirun  -oversubscribe -np 1 R --vanilla
library(Rmpi)
library(snow)
mpi.spawn.Rslaves()

## mpi.send.Robj to send data subset to slaves
## https://scicomp.ethz.ch/wiki/Distributed_computing_in_R_with_Rmpi

## https://stackoverflow.com/questions/62538349/rmpi-mpi-remote-exec-cant-access-user-defined-function
library(lme4)
library(parallel)

nc <- 4
cl <- makeCluster(nc)
dd <- split(sleepstudy, rep(1:nc, length.out = nrow(sleepstudy)))
## with parallel pkg, can't easily export particular
## shards to particular workers;
## have to export all of the pieces (which defeats the purpose of
## distributed memory computing ...)
## see ?lme4::modular

mkdevfun <- function(data, form, ...) {
    lmod <- lFormula(form, data, ...)
    do.call(mkLmerDevfun, lmod)
}

form <- Reaction ~ Days + (Days | Subject)
clusterExport(cl, c("mkdevfun", "dd", "form"))
invisible(clusterEvalQ(cl, library(lme4)))

## https://stackoverflow.com/questions/61383366/parallelisation-within-a-function-in-r-avoiding-automatic-export-of-objects-t



## doing this in parallel is sort of pointless
## again, ideally would want to run mkdevfun() separately on
## each worker, applied to the worker's own data shard
devfuns <- parLapply(cl, dd, mkdevfun, form = form)

comb_devfun <- function(theta) {
    force(theta)
    clusterExport(cl, "theta", envir = environment())
    parLapply(cl, devfuns, function(f) f(theta))
}

try(comb_devfun(c(1, 1, 1, 1)))

## this doesn't work, so try it in serial (as proof of concept)

comb_devfun_serial <- function(theta) {
    sum(vapply(devfuns, function(f) f(theta), FUN.VALUE = numeric(1)))
}

comb_devfun_serial(c(1,1,1,1))
