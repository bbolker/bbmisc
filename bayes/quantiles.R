library(invgamma)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size=15))

shape <- 2
mean <- 2
ratemax <- 10
timemax <- 5

ran <- 100
steps <- 500

comb <- tibble(lrate = log(ran)*seq(-steps, steps)/steps
	, rate = exp(lrate)
	, time = 1/rate
	, rden = dgamma(rate, shape=shape, scale=mean/shape)
	, tden = dinvgamma(time, shape=shape, scale=mean/shape)
	, drat = rate^2*rden/tden
	, pp = pgamma(rate, shape=shape, scale=mean/shape)
	, tp = pinvgamma(time, shape=shape, scale=mean/shape)
)

qcomb <- comb |> filter(pp>0.025 & pp < 0.975)
tqcomb <- comb |> filter(tp>0.025 & tp < 0.975)

rdplot <- (ggplot(comb)
	+ aes(rate, rden)
	+ geom_line()
	+ xlim(c(0, ratemax))
)

print(ggplot(comb)
	+ aes(time, tden)
	+ geom_line()
	+ xlim(c(0, timemax))
)

print(rdplot)
print(rdplot
	+ geom_segment(data=qcomb, aes(x=rate, y=rden, xend=rate, yend=0))
)
print(rdplot
	+ geom_segment(data=tqcomb, aes(x=rate, y=rden, xend=rate, yend=0))
)
