library(shellpipes); manageConflicts()
library(dplyr)

loadEnvironments()

set.seed(101)
t0 <- system.time(
  tt0 <- tabfun(n=17, nsim =  1e2, levNames="Mag/Sign")
)
print(tt0)
