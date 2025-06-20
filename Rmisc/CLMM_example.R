#CLMM

## BMB: don't do this
## e.g. see https://forum.posit.co/t/project-oriented-workflow-setwd-rm-list-ls-and-computer-fires/3549
## setwd("C:/Users/anmol/OneDrive/Documents/suicide stats 25/r SUICIDE")
library(readxl)
library(ordinal)
library(bbmle) ## for AICtab()
library(broom.mixed)
library(parallel)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(emmeans)
## library(dotwhisker) ## not used [raw ggplot instead]
## library(marginaleffects) ## doesn't work yet ...

dataset <- NULL
infn <- "SUICIDE STAT CLEANED FOR THESIS.xlsx"
if (file.exists(infn)) {
  dataset <- read_excel(infn)
}

if (!is.null(dataset)) {
dataset <- transform(dataset,
  RESP = factor(RESP, 
                levels = 1:4,
                labels = c("Strongly Disagree", "Disagree", "Agree", "Strongly Agree"),
                ordered = TRUE),
  TB = factor(TB),
  PB = factor(PB),
  LV = factor(LV),
  ## construct interactions by hand, clmm() can't handle these
  PB_LV_TB = interaction(PB, LV, TB),
  RESPID_PB_LV_TB = interaction(RESPID, PB, LV, TB),
  DE_PB_LV_TB = interaction(DE, PB, LV, TB)
)

##SUBSET DATA FOR EACH ANALYSIS
## See below, I've done this on-the-fly within the model statement
dataQ1 <- subset(dataset, QUES == 1)
dataQ2 <- subset(dataset, QUES == 2)
dataQ3 <- subset(dataset, QUES == 3)
}

fn  <- "clmm_fits.rds"
if (file.exists(fn)) {
  fits <- readRDS(fn)
} else {
  
  ## we only have 4 jobs to run (so far), but each one is multi-threaded, so
  ##   need more cores in cluster ... ???
  ##
  ## I happen to have > 16 cores on my computer.
  cl <- makeCluster(min(16, detectCores()-1))

  ## all possible (so far) choices for random effect terms
  re_forms <- c(
    ## fully crossed (maximal model)
    full = "(PB*LV*TB|RESPID) + (PB*LV*TB|DE)",
    ## two-way interactions among factors
    ##  (only reduces from 8x8 to 7x7 covariance matrix)
    twoway = "((PB+LV+TB)^2|RESPID) + ((PB+LV+TB)^2|DE)",
    ## main effects only
    maineff = "(PB+LV+TB|RESPID) + (PB+LV+TB|DE)",
    ## complex random intercepts à la Scandola and Tidoni
    ## (need to spell these out, clmm can't expand in-place)
    CRI = "(1|RESPID) + (1|RESPID_PB_LV_TB) + (1|DE) + (1|DE_PB_LV_TB)",
    CRI_minus = "(1|RESPID) + (1|DE) + (1|DE_PB_LV_TB)"
  )

  fitfun <- function(re_form, question = 1) {
    ## add fixed effects (full three-way interaction) to specified RE
    form <- reformulate(c("PB*LV*TB", re_form), response = "RESP")
    dd <- subset(dataset, QUES == question)
    t1 <- system.time(
      fit <- clmm(form,
                  data = dd,
                  ## 'trace' gives noisy output (but disappears if run in parallel)
                  control = clmm.control(trace = 1,
                                         iter.max = 1000, eval.max = 1000)
                  )
    )
    attr(fit, "time") <- t1
    fit
  }

  ## fit everything in parallel - have to load packages and export data
  invisible(clusterEvalQ(cl, library(ordinal)))
  clusterExport(cl, c("dataset", "re_forms"))
  fits <- parLapply(cl, re_forms, fitfun)

  stopCluster(cl)
  saveRDS(fits, file = fn)
}

if (exists("dataQ1")) {
  ## checking distribution/experimental design
  with(dataQ1, table(table(PB, LV, TB, RESPID)))  ## 2 observations for every combination
  ## for example ...
  subset(dataQ1, LV ==1 & TB == 0 & PB == 1 & RESPID ==84)
}

## both RE terms from the full model
##  are singular (determinant approx 0)
sapply(VarCorr(fits[["full"]]), det)

## check if any of the RE terms is singular
isSingular <- function(x, tol = 1e-8, collapse = TRUE) {
  dets <- vapply(VarCorr(x), det, FUN.VALUE = numeric(1))
  sing <- dets < tol
  if (!collapse) sing else any(sing)
}

get_time <- function(x, w = "elapsed") {
  attr(x, "time")[[w]]
}

AICtab(fits)
## choosing REs on the basis of minimum AIC is reasonable ...
## e.g. Matuschek et al 2017 J Mem Language,
##      Moritz et al 2023 Ecology Letters
##
sapply(fits, isSingular) ## although *all* models are singular
lapply(fits, isSingular, collapse = FALSE)

sapply(fits, get_time)  ## elapsed time in seconds (72 mins for full model)
## CPU time (¿not quite sure this is the right element to look at?)
sapply(fits, get_time, "sys.self")  

## proceed with simplest (lowest-AIC) model
fCRI <- fits[["CRI_minus"]]
## RE standard devs -- still v small for one term
get_sds <- function(fit) {
  sapply(VarCorr(fit), attr, "stddev")
}
print(get_sds(fCRI))
## ?? not sure why last two variances are identical ... ??
## maybe the model is over-parameterized after all?

## would normally worry about scaling by 2SD, but all predictors are
##  categorical/binary
tt <- tidy(fCRI, effects = "fixed", conf.int = TRUE)  |>
  filter(coef.type == "location") |>
  rename(lwr = "conf.low", upr = "conf.high") |>
  ## order effects by magnitude
  mutate(across(term, ~ reorder(factor(.), estimate)))

## first two fits give non-pos-def Hessians, have to work harder
## (move this fix upstream to the `tidy` method?)
tfun <- function(x) {
  npd <- inherits(try(vcov(x), silent = TRUE), "try-error")
  if (npd) {
    cc <- coef(x)
    cc <- cc[!grepl("|", names(cc), fixed = TRUE)]
    ret <- tibble(term = names(cc), estimate = cc,
                  lwr = NA_real_, upr = NA_real_)
    return(ret)
  }
  tidy(x, effects = "fixed", conf.int = TRUE)  |>
                      filter(coef.type == "location") |>
                      rename(lwr = "conf.low", upr = "conf.high")
}
tt_all <- purrr::map_dfr(fits, tfun, .id = "model") |>
  mutate(across(term, ~ reorder(factor(.), estimate)),
         across(model, forcats::fct_inorder))


ggplot(tt_all, aes(estimate, term, color = model)) +
  geom_pointrange(aes(xmin=lwr, xmax = upr),
                  position = position_dodge(width = 0.3)) +
  geom_vline(xintercept = 0, lty = 2)

## these are in units of log-odds

## final model only
## could use dotwhisker::dwplot() but I like to tweak things a little bit more
ggplot(tt, aes(estimate, term)) +
  geom_pointrange(aes(xmin=lwr, xmax = upr)) +
  geom_vline(xintercept = 0, lty = 2)
## note for interpretation purposes
##   that these effects are using *treatment contrasts*

## can look at expected marginal means, for every level
plot(emmeans(fCRI, ~PB*LV*TB))
## or for one particular main effect
(em_PB <- emmeans(fCRI, ~PB))
pairs(em_PB)



##  TODO:
##   * re-run for other questions, e.g.
##      fits_Q2 <- parLapply(cl, re_forms, fitfun, question = 2)
##   * try out RTMB?
##   * `marginaleffects` package: doesn't know about `ordinal` yet ...
