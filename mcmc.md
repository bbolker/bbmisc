---
title: "MCMC notes"
author: "Ben Bolker"
date: "6 October 2021"
bibliography: mcmc.bib
---

This is a **very** broad outline/brain dump.

## Basic ideas

[@gilks_markov_1995; @van_ravenzwaaij_simple_2018; @bolker_ecological_2008]

- need to construct a Markov chain on the parameter space
- *stationarity*, *ergodicity*, *irreducibility*
- *detailed balance*/reversibility: $\pi_i P_{ij} = \pi_j P_{ji}$

## Sampling strategies

- Metropolis-Hastings
- Gibbs (conditional conjugate sampling)
- slice
- 'generalized Gibbs' (x-within-Gibbs)
- hybrid/Hamiltonian MCMC (uses gradient info)
- reversible-jump MCMC (variable dimensionality)
- sequential MC (time series etc.)
- Metropolis-coupled MCMC (MC^3) (tempering) [@altekar_parallel_2004]

## Platforms

- BUGS/JAGS
- Stan (+ `rethinking::ulam`)
- Nimble
- PyMC3
- TMB + `tmbstan`
- [greta](https://greta-stats.org)
- special-case front ends (`MCMCglmm`, `brms`, `rstanarm`, ...)

\(1) Which samplers are available? (2) Simplicity vs flexiblity (2) *Procedural* or *graphical* model definition?

## Diagnostics and visualization

- trace plots
- R-hat, improved R-hat, effective sample size, ... [@vehtari_rank-normalization_2019; @lambert_r_2020]
- simulation-based calibration [@talts_validating_2020]

## Troubleshooting

- tuning/adaptation [@rosenthal_optimal_2011]
   - shape of candidate distribution should match shape of posterior
   - acceptance probability 0.1 to 0.6, approx 0.24 in high dim
- change samplers
- reparameterize
- make priors stronger
- run longer!

See @gelman_bayesian_2020

## References
