## problem from Jake Szamosi
## pre-processing/anonymization in RatioMods_pre.R
library(glmmTMB)
library(ggplot2); theme_set(theme_bw())
zmargin = theme(panel.spacing=grid::unit(0,"lines")) ## squash together
library(emmeans)    
library(broom.mixed)  ## tidying models
library(dotwhisker)   ## coefficient plots
library(DHARMa)       ## diagnostics
library(plotrix)      ## for corner.label()
library(robustbase)   ## for lmrob()
library(GGally)       ## ggpairs()
library(bbmle)        ## AICctab()
library(tidyverse)

cd_df <- readRDS("cd_df.rds")
## z: predictor
## D: counts of "denominator" species
## C: counts of "numerator" species
## ratio: what it says

named_list = lme4:::namedList

## bare-bones: no conf.int, conf.level, etc.
tidy.lmrob <- function(x) {
    res <- (x
        %>% summary()
        %>% coef()
        %>% as.data.frame()
        %>% rownames_to_column("term")
        %>% as_tibble()
        %>% rename(estimate="Estimate",
                   std.error="Std. Error",
                   statistic="t value",
                   p.value="Pr(>|t|)")
    )
    return(res)
}


## GLMs
fam_list = named_list(poisson,nbinom1,nbinom2,gamma=Gamma(link="log"))
mod_list = map(fam_list,
                ~glmmTMB(C ~ z + offset(log(D)),
                         data = cd_df, family = .))
## add log-linear model
mod_list = c(mod_list,
             list(
                 lm_log=lm(log(C) ~ z + offset(log(D)), data = cd_df),
                 lm_log2=lm(log(ratio) ~ z , data = cd_df)
             )
             )

## turns out most of the fancy downstream stuff I want to use doesn't
## work with robustbase::lmrob models. Ugh.
## more trouble than it's worth to include in mod_list
lmrob_log = lmrob(log(C) ~ z + offset(log(D)), data = cd_df)

## NB1 and Poisson give nearly identical results. They assume the
## same mean-variance relationship, up to scaling: Poisson, V=mu,
## NB1, V=phi*mu (for some estimated phi>=1).
## Therefore, points will be weighted equally in the
## Poisson and NB1 regressions.

## NB2 and Gamma also give nearly identical results.
## var-mean relationships are **not** identical: NB2, var=mu*(1+mu/k),
## Gamma, var=mu^2/shape (I think).  However: for mu>>k, NB2 var is
## *approx* mu^2/k, so it looks like a Gamma.
##
## It's mildly surprising that log-linear model doesn't give similar
## results to Gamma (it also has var proportional to mu^2), but it
## doesn't ... 

ss = suppressMessages(map(mod_list, DHARMa::simulateResiduals))
plot_fun = function(x,n) {
    plotResiduals(x)
    plotrix::corner.label(x=1,y=1,n,cex=2)
    invisible(plotQQunif(x))
}
## parameters to squash everything together
## need a tall, skinny aspect ratio for this to be at all readable
op = par(mfrow=c(length(ss),2),mar=c(3,1,1,0),mgp=c(1,0.5,0), no.readonly=TRUE)
map2(ss,names(ss),plot_fun)
par(op) ## restore parameters

## I don't guarantee the DHARMa plots -- sometimes
##  I get weird results  that I don't understand -- but
##  it seems pretty clear that the lm results suck less than
##  the rest

em = (mod_list
    %>% map(emmeans,
            specs='z',
            cov.reduce=FALSE,
            offset=0)
    %>% map(as.data.frame)
    %>% bind_rows(.id="Model")
)

cd_est = ggplot(em, aes(z, emmean)) +
	geom_point(data = cd_df, aes(y = log(ratio), size=D), alpha = 0.5) +
	geom_line(aes(colour = Model))+
	geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = Model),
				alpha = 0.3) +
	scale_colour_brewer(palette = 'Dark2') +
    scale_fill_brewer(palette = 'Dark2')+
    labs(y="log(ratio)")

library(directlabels)
cd_est_d <- direct.label(cd_est,method="first.bumpup")+expand_limits(x=-3)
ggsave(cd_est_d,file="nbfig1.png",width=6,height=6,dpi=150)

## It looks like Pois is missing here. It's not, it's just entirely overlapping
## NB1, and also has extremely narrow confidence intervals.
## Similarly, Gamma is almost identical to NB2

cd_est + facet_wrap(~Model) + zmargin

## I think this works:
pred_res = map(mod_list,~augment(.,type.residual="pearson")) %>%
    bind_rows(.id="Model") %>%
    group_by(Model) %>%
    mutate(sample=seq(n())) %>%
    ungroup()

## Poisson is crazy
## BMB: this is because we are using Pearson residuals,
## i.e.  (raw resid)/sqrt(V_est)
## because the Poisson estimates ridiculously small variances,
## the Pearson residuals are ridiculously inflated

resid_plt = ## ggplot(pred_res, aes(Pred, Resid)) +
    ggplot(pred_res, aes(.fitted, .resid)) +
    geom_point() +
    geom_smooth() +
    geom_hline(yintercept = 0) +
    facet_grid(.~Model, scales = 'free_x') +
    zmargin

resid_plt
resid_plt %+% filter(pred_res, Model!="poisson")
## all models have severe bias problems except for lm_log
## nbinom2 and gamma are again nearly identical

pred_all = (pred_res
			%>% select(sample,Model,.fitted)
			%>% spread(Model, .fitted))

head(pred_all)

## There's also a suspiciously clean linear relation between all the pairs of
## estimates
## BMB: why is that suspicious? or am I getting different answers?

## order columns to group models with similar resultso
ggpairs(pred_all,columns=c(2,5,4,6,3))+ zmargin


ggplot(pred_all, aes(nbinom1, poisson)) +
	geom_point() +
	geom_abline(slope = 1, intercept = 0) +
	geom_smooth(method = 'lm')


## AICc table
## we're transforming the response variable from F to log(F/B):
## Jacobian = d(log(F/B))/dF = 1/(F/B)*(1/B)*dF/dF = 1/F
## or alternatively
## d(log(F/B)) = d(log(F)-log(B))/dF = d(log(F))/dF-d(log(B))/dF = 1/F
## need to add -2 * sum(log(J)) = -2*sum(log(1/F)) = 2*sum(log(F))

## Akaike, Hirotugu. “On the Likelihood of a Time Series Model.” Journal of the Royal Statistical Society. Series D (The Statistician) 27, no. 3/4 (1978): 217–35. https://doi.org/10.2307/2988185.

## compute AICc values and hack lm_log value appropriately
aa = bbmle::AICctab(mod_list,base=TRUE)
nm = attr(aa,"row.names")
cc <- 2*sum(log(cd_df$C))
aa$AICc[nm=="lm_log"] = aa$AICc[nm=="lm_log"] + cc
aa$dAICc <- aa$AICc-min(aa$AICc) ## recompute dAICc
aa

## lm_log < (nbinom2, gamma) << nbinom1 <<< poisson

### So, to summarize, the Poisson model:
## - is not giving any kind of error or warning
## - is producing nearly identical point estimates to NB1
## - is producing residuals several orders of magnitude bigger than the other models
## - is producing SEs near zero (with p-values to match)
## - has a comically large AIC compared to the two NB models.
##    BMB: it thinks the fit is terrible because it is underestimating
##     the residual variance


## BMB shouldn't look at this until we're all done ...
## boring since there is only one parameter
mod_tab <- (mod_list
    %>% map(tidy)
    ## add lmrob results
    %>% c(list(lmrob_log=tidy(lmrob_log)))
    %>% bind_rows(.id="model")
)
dotwhisker::dwplot(mod_tab) + coord_fixed(ratio=0.5) +
    geom_vline(xintercept=0,lty=2)+
    guides(colour = guide_legend(reverse = TRUE))
## miraculously doesn't need reordering except for legend reversal

