library(invgamma)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size=15))

## Parameters
shape <- 10
mean <- 2
ratemax <- 10
timemax <- 5

ran <- 100
steps <- 500

scale=mean/shape

## Make a combined rate/time-based density profile

## drat checks the transformation logic, should be 1
comb <- tibble(lrate = log(ran)*seq(-steps, steps)/steps
	, rate = exp(lrate)
	, time = 1/rate
	, rden = dgamma(rate, shape=shape, scale=scale)
	, tden = dinvgamma(time, shape=shape, scale=scale)
	, drat = rate^2*rden/tden
	, pp = pgamma(rate, shape=shape, scale=scale)
	, tp = pinvgamma(time, shape=shape, scale=scale)
)
print(comb, n=Inf)

## Quantile intervals (these are the same)
qcomb <- comb |> filter(pp>0.025 & pp < 0.975)
tqcomb <- comb |> filter(tp>0.025 & tp < 0.975)

######################################################################

## HPD functions

## Density at a quantile
dq <- function(x, qfun=qfun, dfun=dfun, ...){
	return(dfun(qfun(x, ...), ...))
}

## Density difference between tails
deltaDensity <- function(prop, qfun, dfun, alpha, ...){
	left <- dq(alpha*prop, qfun, dfun, ...)
	right <- dq(1-alpha*(1-prop), qfun, dfun, ...)
	return(right-left)
}

## Balance tails and find an HPD cutoff
dcut <- function(qfun, dfun, alpha=0.05, eps=1e-3, ...){
	u <- uniroot(deltaDensity
		, lower=eps, upper=1-eps, qfun=qfun, dfun=dfun, alpha=alpha, ...
	)
	prop <- u$root
	stopifnot(is.numeric(prop))
	return(dq(alpha*prop, qfun, dfun, ...))
	
	## Not reached
	return(c(
		prop
		, dq(alpha*prop, qfun, dfun, ...)
		, dq(1-alpha*(1-prop), qfun, dfun, ...)
	))
}

## Works OK, but less accurate than expected
## print(dcut(qnorm, dnorm))
## print(dcut(qgamma, dgamma, shape=2, scale=1))

## HPD calc
rpdcut <- dcut(qgamma, dgamma, shape=shape, scale=scale)
tpdcut <- dcut(qinvgamma, dinvgamma, shape=shape, scale=scale)

rpdcomb <- comb |> filter(rden >= rpdcut)
tpdcomb <- comb |> filter(tden >= tpdcut)

rdplot <- (ggplot(comb)
	+ aes(rate, rden)
	+ geom_line()
	+ xlim(c(0, ratemax))
	+ xlab("rate (per day)")
	+ ylab("density (day)")
)

tdplot <- (ggplot(comb)
	+ aes(time, tden)
	+ geom_line()
	+ xlim(c(0, timemax))
	+ xlab("time (day)")
	+ ylab("density (per day)")
)

qrdplot <- (rdplot
	+ geom_segment(data=qcomb, aes(x=rate, y=rden, xend=rate, yend=0))
	+ ggtitle("Quantile-based credible interval")
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

qtdplot <- (tdplot
	+ geom_segment(data=qcomb, aes(x=time, y=tden, xend=time, yend=0))
	+ ggtitle("Quantile-based credible interval")
)
cqtdplot <- (tdplot
	+ geom_segment(data=tqcomb, aes(x=time, y=tden, xend=time, yend=0))
	+ ggtitle("Quantile-based credible interval from the rate scale")
)

print(cqtdplot)

quit()
print(tdplot
	+ geom_segment(data=rpdcomb, aes(x=time, y=tden, xend=time, yend=0))
)
print(tdplot
	+ geom_segment(data=tpdcomb, aes(x=time, y=tden, xend=time, yend=0))
)
