library(shellpipes); manageConflicts()
library(dplyr)

library(ggplot2); theme_set(theme_bw() + theme(panel.spacing=grid::unit(0, "lines")))
library(directlabels)

loadEnvironments()

## Okabe-Ito minus black and yellow
oi3 <- palette.colors(9)[-c(1, 5)]
out_scale <- scale_colour_manual(name = "outcome category",  values = oi3)

## ----dushoff_effects----------------------------------------------------------
largeEffect=1.3; smallEffect=0.5; tinyEffect = 0.1; smallVar = 0.2; medVar = 0.4; largeVar = 0.7; hugeVar = 1.2
span = 2
vget = Vectorize(get)
vr <- (
	read.table(header=TRUE, strip.white=TRUE, sep=":", text="
		pic : val : unc : atext : ntext
		PL : largeEffect : smallVar : Clearly large|and positive : different
		PU : largeEffect : largeVar : Clearly positive,|maybe large : different
		PS : smallEffect : smallVar : Clearly positive|and not large : different
		US : tinyEffect : medVar : Maybe positive,|clearly small : different
		UU : smallEffect : largeVar : Not both (large|and negative) : different
		nopower : tinyEffect : hugeVar : Should have|done a power|analysis first : different
	")
)
vf <- (vr
	|> mutate(NULL
		, pic = factor(pic, levels=rev(pic))
		, val = vget(val)
		, unc = vget(unc)
		, lwr = val-unc
		, upr = val+unc
		, atext = gsub('\\|', '\n', atext)
	)
)
print(ggplot(vf)
	+ aes(val, pic)
	+ geom_pointrange(aes(xmin=lwr, xmax=upr))
    + theme(panel.grid.minor.x = element_blank())  ## delete vertical grid lines
	+ geom_vline(xintercept=c(0, -1, 1), lty = c(1, 2, 2))
	+ scale_x_continuous(limits=c(-span, span)
		, breaks = -1:1
		, labels = c("cutoff\n(SESOI)", 0 , "cutoff\n(SESOI)")
	)
	+ scale_y_discrete(labels=rev(vf$atext))
  + labs(x = "", y = "")
)

## ----pow1---------------------------------------------------------------------
pp <- power.t.test(delta = 1, sd = 1, power = 0.8)
n <- ceiling(pp$n)  ## 17
tt <- power.t.test(delta = 1, sd = 1, n = n)
pow1 <- round(tt$power,3) ## 0.807

## ----simfun0------------------------------------------------------------------
set.seed(101)
t0 <- system.time(
  tt0 <- tabfun(n=17, nsim =  10000)
)
pow1_sim <- sum(tt0[1:3])
stopifnot(all.equal(pow1_sim, tt$power, tolerance = 2e-3))
print(tt0)

## ----sim1, cache=FALSE--------------------------------------------------------
set.seed(101)
nvec <- c(5:10, (2:9)*10, 100, 200)
mm <- lapply(nvec, tabfun, delta=0.5, nsim = 10000) |> do.call(what = rbind)

## ----sim1_sum-----------------------------------------------------------------
mmw <- f_widen(mm)

## ----sim1_plot, fig.width = 8-------------------------------------------------
gg1 <- ggplot(mmw, aes(n, value, colour = name)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  labs(y = "proportion", x = "sample size per group") +
  out_scale
## see https://tdhock.github.io/directlabels/docs/index.html
##  for direct labeling choices
print(direct.label(gg1, "top.bumptwice"))

## ----mm_tail------------------------------------------------------------------
drop(tail(mm,1))

## ----sim2, cache=FALSE--------------------------------------------------------
mm2 <- lapply(nvec, tabfun, nsim = 10000, delta = 1.5) |> do.call(what = rbind)

## ----sim2_plot, fig.width = 8-------------------------------------------------
mmw2 <- f_widen(mm2)
print(direct.label(gg1 %+% mmw2, "top.bumptwice"))

## ----si-----------------------------------------------------------------------
sessionInfo()

