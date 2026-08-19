---
title: "Scale-dependent parameter estimation"
author: "Ben Bolker"
date: "20 Aug 2021"
bibliography: scale_dep.bib
---

Scale-dependent parameter estimation is a wonderful idea that seems to be super-useful and yet super-underused in ecology. The basic idea is that we have a response variable (e.g., a population density) measured at some points in space ($y(z)$) and one or more predictor variables measured at the same, *or other*, points in space ($x(z)$). If we assume that the values of $x$ determine the values of $y$ *at exactly the same point in space*, then this becomes a standard regression problem. If the values of $x$ affect *nearby* values of $y$, however, this doesn't work. The deterministic (linear) relationship underlying the relation, rather than $E[y(z)] = \beta_0 + \beta_1 x(z)$, is $E[y(z)] = \beta_0 + \beta_1 \int K(|z'-z|) x(z') \, dz$. In other words, $y$ depends on a *spatial convolution* of $x$ with a *kernel* $K$ (sometimes written as $(K * x)(z)$. 

If we happen to know $K$, we can compute the convolution directly and put the result into some regression framework. But what if we don't know $K$?

One standard approach in spatial ecology (REFS) is to assume a shape of the kernel $K$ (typically uniform or "top-hat", i.e. all locations at a distance less than some radius $\rho$ affect the response at $y$ equally, locations farther away have no effect). Then for each scale or range $\rho$, compute the convolution (i.e. the average of the environmental at all points within a distance $\rho$ of each response), fit the model, and compute some goodness-of-fit statistic ($R^2$, AIC, etc.). Pick the $\rho$ with the best goodness-of-fit and use it in the final model. In my opinion this is fine as far as it goes, but (1) it's expensive (you need to compute a convolution for every possible scale and fit a model for each one); (2) it assumes a particular and unrealistic shape for $K$; (3) it's hard to account for uncertainty in the scale, or do inference on the scale itself (although these problems could certainly be overcome, e.g. by doing model averaging across models with different scales).

I have previously spent a lot of time and effort until about 2012 or so trying to come up with correlation-based methods to extract information about scale and shape of convolution kernels, using the idea of *deconvolution* and Fourier transform methods [@bolker_spatial_2009;@bolker_inferring_2011;@seabloom_spatial_2005]. I spent a long time on this stuff, but mostly didn't get it published. It is theoretically elegant, and is in principle useful for combining the effects of different spatial processes (dispersal, competition, environmental filtering), but there are lots of technical hurdles that kept it from being practical. 

Since then I have learned about *scale-dependent parameter estimation*, which gets at the idea of estimating $K$ in a different way. It's related to, or could be considered a variant of *varying-coefficient* or *geographic* regression [@gelfand_spatial_2003;@finley_bayesian_2020]: in varying-coefficient regression, we assume that $E[y(z_1)] = \beta_0(z_1) + \beta_1(z_1) x(z_1)$, where $\beta_0$ and $\beta_1$ are smoothly varying functions of spatial position. In other words, the regression coefficients are spatial variables, e.g. the effects of $x$ on $y$ could change (smoothly) from place to place. In scale-dependent parameter estimation, in contrast, we use the same basic tools to estimate $E[y(z_1)] = \beta_0 + \beta_1 (K*x)(z_1)$ (the same equation that we introduced at the beginning).

The basic trick is to compute the average values of $x$ within *annuli* around each response point (i.e. rings where $\rho < |z_1-z_2| < \rho + \Delta \rho$), include all of these columns in the set of predictors, and use the machinery of GAMs to penalize the fit so that parameters of $\beta(\rho)$ vary smoothly with $\rho$.

I got this idea from presentations of Thomas Cornulier and Alexandre Villers [@cornulier_modelling_2015;@cornulier_modelling_2016] (also not published!). They cite an earlier (less technically elegant) version that uses a similar idea to estimate the *temporal* scale of environmental influences [@sims_identifying_2007].

I asked eco-stats tweeps about this in a [twitter thread](https://twitter.com/bolkerb/status/1320452421492051968) earlier this year but didn't get much response.

I think this is a fantastic method that deserves to be better developed. Please let me know if you're interested in working on it! (Or if you find a good reference showing that it has already been well developed and explicated by someone, somewhere ...)

## References

