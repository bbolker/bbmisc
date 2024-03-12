library(invgamma)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size=15))

shape <- 2
mean <- 2
ratemax <- 10
timemax <- 5

ran <- 100
steps <- 500

## drat checks the transformation logic, should be 1
comb <- tibble(lrate = log(ran)*seq(-steps, steps)/steps
	, rate = exp(lrate)
	, time = 1/rate
	, rden = dgamma(rate, shape=shape, scale=mean/shape)
	, tden = dinvgamma(time, shape=shape, scale=mean/shape)
	, drat = rate^2*rden/tden
	, pp = pgamma(rate, shape=shape, scale=mean/shape)
	, tp = pinvgamma(time, shape=shape, scale=mean/shape)
)
print(comb, n=Inf)

## Quantile intervals (these are the same)
qcomb <- comb |> filter(pp>0.025 & pp < 0.975)
tqcomb <- comb |> filter(tp>0.025 & tp < 0.975)

## HPD wrong logic! Need to find a density cutoff and then calculate area.

svec <- seq(1/2, steps-1/2)/steps
rdensity=dgamma(ratemax*svec, shape=shape, scale=mean/shape)
tdensity=dinvgamma(timemax*svec, shape=shape, scale=mean/shape)

rpdcut <- quantile(rdensity, probs=0.05)
tpdcut <- quantile(tdensity, probs=0.05)

print(quantile(svec, probs=c(0.05, 0.95)))
## print((comb$tden<tpdcut))

print(c(rpdcut=rpdcut, tpdcut=tpdcut))

hist(rdensity)

rpdcomb <- comb |> filter(rden <= rpdcut)
tpdcomb <- comb |> filter(tden <= tpdcut)

print(tpdcomb, n=Inf)

rdplot <- (ggplot(comb)
	+ aes(rate, rden)
	+ geom_line()
	+ xlim(c(0, ratemax))
)

tdplot <- (ggplot(comb)
	+ aes(time, tden)
	+ geom_line()
	+ xlim(c(0, timemax))
)


print(rdplot)

## rates <- data.frame(rate=svec*ratemax, rden=rdensity)
## print(rdplot + geom_point(data=rates))

print(rdplot
	+ geom_segment(data=qcomb, aes(x=rate, y=rden, xend=rate, yend=0))
)
print(rdplot
	+ geom_segment(data=tqcomb, aes(x=rate, y=rden, xend=rate, yend=0))
)
print(rdplot
	+ geom_segment(data=rpdcomb, aes(x=rate, y=rden, xend=rate, yend=0))
)
print(rdplot
	+ geom_segment(data=tpdcomb, aes(x=rate, y=rden, xend=rate, yend=0))
)

print(tdplot)
print(tdplot
	+ geom_segment(data=qcomb, aes(x=time, y=tden, xend=time, yend=0))
)
print(tdplot
	+ geom_segment(data=tqcomb, aes(x=time, y=tden, xend=time, yend=0))
)
print(tdplot
	+ geom_segment(data=rpdcomb, aes(x=time, y=tden, xend=time, yend=0))
)
print(tdplot
	+ geom_segment(data=tpdcomb, aes(x=time, y=tden, xend=time, yend=0))
)
