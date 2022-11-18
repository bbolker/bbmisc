# Brain dump on MCMC for epidemics

We would like to be able to estimate parameter values for dynamic models with *all* of the following characteristics:

* discrete states: we usually want to count at the level of individuals. Especially for beginnings/ends of epidemics, and outbreaks in small populations, finite-size effects (increased sampling noise at low prevalence and fadeout/extinction processes) are important
* continuous time: epidemics 'really' run into continuous time; even though time scales of epidemic processes are usually longer than a day, some processes can be close to this time scale, and discreteness can cause annoying dynamical instabilities [Ref Mollison and Ud Din?]
* both observation and process error
   * note that 'process error' can occur at two weakly separable scales, i.e. 'sampling-level' (demographic noise, either 'simple' [Poisson noise/Poisson-process branching events] or overdispersed [Hooke processes, Gamma-white noise processes [Ionides and King], negative binomial/beta-binomial epidemic sampling]) or stochastic time-varying parameter values, especially transmission rates
* 'plug-and-play' analysis of complex epidemic models
* convenient inference, especially Bayesian
* regularization/priors