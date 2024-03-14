library(invgamma)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size=15))

## Parameters
shape <- 4
mean <- 2

ratemax <- 6
timemax <- 3
yspace <- 0.0

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

RatePlot <- (ggplot(comb)
	+ aes(rate, rden)
	+ geom_line()
	+ xlim(c(0, ratemax))
	+ xlab("rate (per day)")
	+ ylab("density (day)")
	+ ggtitle("Our posterior")
)

TimePlot <- (ggplot(comb)
	+ aes(time, tden)
	+ geom_line()
	+ xlim(c(0, timemax))
	+ xlab("time (day)")
	+ ylab("density (per day)")
	+ ggtitle("The same posterior (after non-linear transformation)")
)

qRatePlot <- (RatePlot
	+ geom_ribbon(data=qcomb, aes(x=rate, ymax=rden), ymin=-yspace
		, alpha=0.3
	)
	+ ggtitle("Quantile-based credible interval")
)
cqRatePlot <- (RatePlot
	+ geom_ribbon(data=tqcomb, aes(x=rate, ymax=rden), ymin=-yspace
		, alpha=0.3
	)
	+ ggtitle("Quantile-based credible interval from the time scale")
)
hdRatePlot <- (RatePlot
	+ geom_ribbon(data=rpdcomb, aes(x=rate, ymax=rden), ymin=-yspace
		, alpha=0.3
	)
	+ ggtitle("Highest density credible interval")
)
chdRatePlot <- (RatePlot
	+ geom_ribbon(data=tpdcomb, aes(x=rate, ymax=rden), ymin=-yspace
		, alpha=0.3
	)
	+ ggtitle("Highest density credible interval from the time scale")
)

qTimePlot <- (TimePlot
	+ geom_ribbon(data=qcomb, aes(x=time, ymax=tden), ymin=-yspace
		, alpha=0.3
	)
	+ ggtitle("Quantile-based credible interval")
)
cqTimePlot <- (TimePlot
	+ geom_ribbon(data=tqcomb, aes(x=time, ymax=tden), ymin=-yspace
		, alpha=0.3
	)
	+ ggtitle("Quantile-based credible interval from the rate scale")
)

hdTimePlot <- (TimePlot
	+ geom_ribbon(data=tpdcomb, aes(x=time, ymax=tden), ymin=-yspace
		, alpha=0.3
	)
	+ ggtitle("Highest density credible interval")
)
chdTimePlot <- (TimePlot
	+ geom_ribbon(data=rpdcomb, aes(x=time, ymax=tden), ymin=-yspace
		, alpha=0.3
	)
	+ ggtitle("Highest density credible interval from the rate scale")
)

rpdcomb |> pull(time) |> max() |> print()
rpdcomb |> pull(rate) |> min() |> print()

print(RatePlot)
print(TimePlot)

print(qRatePlot)
print(qTimePlot)
print(cqTimePlot)

print(hdRatePlot)
print(hdTimePlot)
print(chdTimePlot)
