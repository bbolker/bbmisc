library(invgamma)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size=15))

## Whatever
qspan <- seq(0, 1, by=0.01)
q <- (qspan[-1] + qspan[-length(qspan)])/2

##
shape <- 2
mean <- 2

ratemax <- 10
timemax <- 5
steps <- 500
svec <- seq(1/2, steps-1/2)/steps
print(svec)

rden <- tibble(rate = ratemax*svec
	, density=dgamma(rate, shape=shape, scale=mean/shape)
)

tden <- tibble(time = timemax*svec
	, density=dinvgamma(time, shape=shape, scale=mean/shape)
)

print(ggplot(rden)
	+ aes(rate, density)
	+ geom_line()
)

print(ggplot(tden)
	+ aes(time, density)
	+ geom_line()
)

######################################################################

## Joint approach (above may be deprecated)

ran <- 100
steps <- 500

comb <- tibble(lrate = log(ran)*seq(-steps, steps)/steps
	, rate = exp(lrate)
	, time = 1/rate
	, rden = dgamma(rate, shape=shape, scale=mean/shape)
	, tden = dinvgamma(time, shape=shape, scale=mean/shape)
	, drat = rate^2*rden/tden
)

print(comb, n=Inf)
