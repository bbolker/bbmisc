## source("sim_sesoi/sim_sesoi_funs.R")
library(shellpipes)

loadEnvironments()


set.seed(101)
nsim <- 1e5
system.time(res <- lapply(seq.int(nsim), function(i) simfun(delta=0, sd = 1, n = 2)) |> do.call(what=rbind))
colMeans(res)

system.time(res2 <- simfun2(delta=0, sd = 1, n = 2, nsim = nsim))
colMeans(res2)
colMeans(res2)/colMeans(res)

## with n = 2, delta =0,  sd = 1
## res2:
##          est          lwr          upr 
##  0.001758344 -4.331556609  4.335073297 
## vs resL
#       est          lwr          upr 
##  0.001735941 -3.803351493  3.806823374
## about 13% too wide with n=2 (per group)

