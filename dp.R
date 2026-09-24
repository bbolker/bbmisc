library(poibin)
library(ggplot2); theme_set(theme_bw())
library(purrr)
library(dplyr)

## stuff related to the Dirichlet process prior
## machinery from Claude (needs poibin, Poisson-binomial)
## https://en.wikipedia.org/wiki/Poisson_binomial_distribution

alphavec <- c(1,2,4)
nvec <- c(20, 200, 2000)
dd <- expand.grid(alpha = alphavec, n = nvec) |>
  mutate(s = 1:n(), .before = 1)

ddirproc <- function(alpha, n) {
  ## needs numerical accuracy: alpha + seq_len(n) -1 in denominator can fail
  p0 <- alpha / (alpha + 0:(n-1))
  poibin::dpoibin(kk = 0:n, pp = p0)
}

dirprocmean <- function(alpha, n) {
  alpha*(digamma(alpha+n) - digamma(alpha))
}


## from Claude
dirprocmode <- function(alpha, n) {
  mu <- dirprocmean(alpha, n)
  k <- floor(mu)
  if (k >= n) return(n)
  ## Darroch (1964): mode is floor(mu) or ceiling(mu); these bounds settle it
  if (mu < k + 1/(k + 2)) return(k)
  if (mu > k + 1 - 1/(n - k + 1)) return(k + 1)
  ## ambiguous zone: compare the two candidate probabilities directly
  p <- ddirproc(alpha, n)[k + 1:2]
  k + (p[2] > p[1])
}

dd2 <- purrr::map_dfr(
                1:nrow(dd), 
                \(i) with(dd[i,],
                          tibble(s=i, x = 0:n, prob = ddirproc(alpha, n)))) |>
  full_join(x=dd, by = "s") |>
  mutate(pars = factor(s, labels = sprintf("alpha=%d, n = %d", dd$alpha, dd$n)))

## distributions
ggplot(dd2, aes(x=x, xend = x, y = 0, yend = prob, colour = pars)) +
  geom_segment() +
  geom_point() +
  geom_line(lty = 2, aes (y = prob)) +
  scale_x_log10(limit = c(NA, 50))

## mean
dd$m <- sapply(1:nrow(dd),
               \(i) with(dd[i,], dirprocmean(alpha, n)))

## what alpha do we need to get a mean of 2 with n=200?
solve_alpha <- function(target, n = 200) {
  uniroot(\(a) dirprocmean(a, n) - target,
          lower = 1e-5, upper = 10)$root
}

solve_alpha(2)
solve_alpha(1.1)

## apparently the mode is always either trunc(mean) or trunc(mean)+1

dirprocmode(0.017, 200)
alpha_vec <- 10^seq(-2, 1, length.out = 101)
mode_vec <- sapply(alpha_vec, dirprocmode, n = 200)
mean_vec <- sapply(alpha_vec, dirprocmean, n = 200)

pdf("dp.pdf")
par(las=1, bty = "l")
plot(alpha_vec, mode_vec, type = "s", log = "xy",
     xlab = expression(alpha),
     ylab = "mode and mean",
     main = "Mean and mode of Dirichlet process for n=200")
## grid()
lines(alpha_vec, mean_vec, col = 2)
text(0.2, 5, "mean", col = 2)
text(1, 2, "mode")
dev.off()


