library("glmmTMB")
library("car")
m2 <- glmmTMB(count ~ spp + mined + (1|site),
   zi=~spp + mined,
   family=nbinom2, data=Salamanders)
drop1(m2,test="Chisq")
Anova(m2)

## simplify (no zi)
simfun <- function() {
    dd <- transform(Salamanders,
            count = simulate_new(~ spp + mined + (1|site),
            newdata = Salamanders,
            newparams = list(beta = c(-0.5, rep(0,6), 1.4),
            betadisp = 0.5,
            theta = -1),
            family = nbinom2)[[1]])
    return(dd)
}

