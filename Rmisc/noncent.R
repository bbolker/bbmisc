
seed <- 32
n <- 1e5
mu <- 2
nu <- 10

set.seed(seed)

z <- rnorm(n, mu)
V <- rchisq(n, nu)

nc <- z/sqrt(V/nu)

hist(nc)
