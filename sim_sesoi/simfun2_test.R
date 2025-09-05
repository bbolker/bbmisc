## source("sim_sesoi/sim_sesoi_funs.R")
library(shellpipes)
library(dplyr)
library(tidyr)

loadEnvironments()


set.seed(101)
nsim <- 1e5
system.time(res <- lapply(seq.int(nsim), function(i) simfun(delta=0, sd = 1, n = 2)) |> do.call(what=rbind))
colMeans(res)

system.time(res2 <- simfun_fast(delta=0, sd = 1, n = 2, nsim = nsim))
colMeans(res2)
colMeans(res2)/colMeans(res)

tt0 <- tabfun(n=17, nsim =  10000, fast = FALSE)
tt1 <- tabfun(n=17, nsim =  10000, fast = TRUE)
cbind(tt0, tt1)

res0 <- lapply(seq.int(10000), function(i) simfun(n=17)) |> do.call(what=rbind)
res1 <- simfun_fast(nsim = 10000, n=17)

comb <- list(slow = res0, fast = res1) |>
  purrr::map(as.data.frame) |>
  purrr::map_dfr(pivot_longer, cols = everything(), .id = "method")

ggplot(comb, aes(x=value, fill=method)) +
  geom_density(alpha=0.5) +
  facet_wrap(~name, scales = "free")

## with n = 2, delta =0,  sd = 1
## res2:
##          est          lwr          upr 
##  0.001758344 -4.331556609  4.335073297 
## vs resL
#       est          lwr          upr 
##  0.001735941 -3.803351493  3.806823374
## about 13% too wide with n=2 (per group)

