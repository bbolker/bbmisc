## clark_2024_ex.R -- reproduces the worked example from @clarkPhylogeneticSmoothing2024
## ("Phylogenetic smoothing using mgcv", https://ecogambler.netlify.app/blog/phylogenetic-smooths-mgcv/,
## see phyloslopes.bib), verbatim except for wrapping in a script and timing
## the fit. Run for a size/timing comparison against the (much bigger, 49
## species x tp(10)) te()-vs-nllfun_tensor comparison in phyloslopes_combo.R:
## this example is 12 species (omega is 12x12) x an 8-level mrf time margin
## (96 total tensor-product coefficients, N = 450 observations, mrf x mrf),
## vs. our 49 species (Sphylo is 49x49) x a 10-df tp() margin (490 total
## coefficients, N = 49 observations, mrf x tp).
##
## Run with: Rscript clark_2024_ex.R

suppressMessages({
  library(ape)
  library(mgcv)
  library(mvnfast)
  library(ggplot2)
  library(dplyr)
  library(MRFtools)
})

sim_gp = function(N, alpha, rho){
  Sigma <- alpha ^ 2 *
    exp(-0.5 * ((outer(1:N, 1:N, "-") / rho) ^ 2)) +
    diag(1e-9, N)
  mvnfast::rmvn(1,
                mu = rep(0, N),
                sigma = Sigma)[1,]
}

set.seed(101)

N_species <- 12
tree <- rcoal(N_species, tip.label = paste0('sp_', 1:N_species))
species_names <- tree$tip.label

N <- 50
shared <- sim_gp(N, alpha = 1, rho = 8) + 10

warp1 <- sim_gp(N, alpha = 2, rho = 20) + 10
warp2 <- sim_gp(N, alpha = 2, rho = 20) + 10
weights1 <- as.vector(scale(rTraitCont(tree)))
weights2 <- as.vector(scale(rTraitCont(tree)))

dat <- do.call(rbind,
               lapply(seq_len(N_species),
                      function(i){
                        sp_trend <- warp1 * weights1[i] +
                          warp2 * weights2[i] + shared
                        obs <- rnorm(N,
                                     mean = as.vector(scale(sp_trend)),
                                     sd = 0.35)
                        if(i %in% c(3, 7)){
                          weight <- 0
                          obs <- NA
                        } else {
                          weight <- 1
                        }
                        data.frame(species = species_names[i],
                                   weight = weight,
                                   time = 1:N,
                                   truth = as.vector(scale(sp_trend)),
                                   y = obs)
                      }))
dat$species <- factor(dat$species, levels = species_names)

dat %>%
  dplyr::mutate(y = dplyr::case_when(
    time <= N-5 ~ y,
    time > N-5 ~ NA,
    TRUE ~ y
  )) -> dat

cat("nrow(dat):", nrow(dat), "\n")

omega <- solve(vcv(tree))
cat("dim(omega):", dim(omega), "\n")

rw_penalty <- mrf_penalty(object = 1:max(dat$time),
                          type = 'linear')
dat$time_factor <- factor(1:max(dat$time))

cat("total te() coefficients: k[1]*k[2] =", 8 * N_species, "\n\n")

cat("=== timing gam() fit ===\n")
t_gam <- system.time({
  mod <- gam(y ~ s(time, k = 10) +
                 te(time_factor, species,
                    bs = c("mrf", "mrf"),
                    k = c(8, N_species),
                    xt = list(list(penalty = rw_penalty),
                              list(penalty = omega))),
             data = dat,
             drop.unused.levels = FALSE,
             method = "REML")
})
print(t_gam)
cat("\nlogLik:", as.numeric(logLik(mod)), "\n")
