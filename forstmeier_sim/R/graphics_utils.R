
## extend = TRUE (the default) brackets the data range with one candidate
## break just below lo and just above hi (when available), rather than
## only ever showing breaks strictly inside [lo, hi]. Without it, a data
## range whose endpoint falls between two candidates (e.g. 0.02, between
## the 0.01 and 0.05 candidates) shows no break near that endpoint at all,
## leaving a visible gap between the nearest labelled gridline and the
## actual edge of the (expanded) panel. Pass extend = FALSE to restore the
## strictly-inside-range behaviour.
logit_breaks <- function(n = 6, extend = TRUE) {
  function(x) {
    rng <- range(x[is.finite(x) & x > 0 & x < 1])
    lo <- rng[1]; hi <- rng[2]
    k <- 1:10
    decades <- 10^-k              # .1 .01 .001 ...
    halves  <- 5 * 10^-k          # .5 .05 .005 ...
    lower <- sort(unique(c(decades, halves)))
    lower <- lower[lower <= 0.5]
    full  <- sort(unique(c(lower, 1 - lower)))
    keep <- full[full >= lo & full <= hi]
    if (extend) {
      ## Denser 1-9 x 10^-k grid, used only to find the closest bracketing
      ## candidate just outside [lo, hi]. The sparser decades/halves grid
      ## used for in-range breaks can leave a big gap between candidates
      ## (e.g. nothing between 0.1 and 0.5) that the panel's default
      ## expansion isn't wide enough to actually reach, silently clipping
      ## the bracket break rather than showing it.
      dense <- sort(unique(as.vector(outer(1:9, 10^-k))))
      dense <- dense[dense <= 0.5]
      dense <- sort(unique(c(dense, 1 - dense)))
      below <- dense[dense < lo]
      above <- dense[dense > hi]
      if (length(below) > 0) keep <- c(max(below), keep)
      if (length(above) > 0) keep <- c(keep, min(above))
    }
    if (length(keep) > n) {
      decade_pts <- sort(unique(c(decades, 1 - decades)))
      thinner <- keep[keep %in% decade_pts | keep == 0.5]
      if (length(thinner) >= 3) keep <- thinner
    }
    if (length(keep) < 3) {          # range too tight around 0.5
      keep <- scales::extended_breaks(n = n)(c(lo, hi))
      keep <- keep[keep > 0 & keep < 1]
    }
    keep
  }
}

## A logit-scaled y scale for the range of `values` (e.g. a data frame's
## ci_lo/ci_hi columns, or just its proportion column), with breaks from
## logit_breaks() and -- critically -- the scale's own `limits` explicitly
## widened to include them. Just passing `breaks = logit_breaks()` to
## scale_y_continuous() is not enough for extend = TRUE's bracket
## candidates to actually render: ggplot2 censors any break lying outside
## a continuous scale's data-derived limits to NA regardless of what the
## breaks function returns, no matter how much panel expansion padding
## surrounds it. Computing the breaks eagerly here, from the actual data,
## lets the limits be widened to match before the scale is even built.
scale_y_logit <- function(values, n = 6, extend = TRUE, ...) {
  rng <- range(values[is.finite(values) & values > 0 & values < 1])
  brks <- logit_breaks(n = n, extend = extend)(rng)
  ggplot2::scale_y_continuous(trans = "logit", breaks = brks,
                               limits = range(c(rng, brks)), ...)
}

## Okabe-Ito palette, excluding black and yellow
okabe_ito <- c("#E69F00","#56B4E9","#009E73","#0072B2","#D55E00","#CC79A7")

scale_colour_okabeito <- function(...) ggplot2::scale_colour_manual(values = okabe_ito, ...)
scale_fill_okabeito   <- function(...) ggplot2::scale_fill_manual(values = okabe_ito, ...)
scale_color_okabeito  <- scale_colour_okabeito

## Log10 scale for the raw variable x (e.g. R0), positioning points on
## log10(x - offset) so values near `offset` (e.g. R0 near 1) are spread out,
## with breaks chosen by the base-R log-tick algorithm on the shifted values
## and labelled in the original (unshifted) units, e.g. 1.01, 1.05, 1.1, ...
trans_log10_shifted <- function(offset = 1) {
  scales::trans_new(
    paste0("log10 (x - ", offset, ")"),
    transform = function(x) log10(x - offset),
    inverse   = function(x) offset + 10^x)
}

scale_x_log10_shifted <- function(offset = 1, n = 10, ...) {
  ggplot2::scale_x_continuous(
    trans = trans_log10_shifted(offset),
    breaks = function(lims) offset + axisTicks(log10(range(lims) - offset), log = TRUE, n = n),
    labels = function(b) format(b, trim = TRUE, drop0trailing = TRUE),
    ...)
}

scale_y_log10_shifted <- function(offset = 1, n = 10, ...) {
  ggplot2::scale_y_continuous(
    trans = trans_log10_shifted(offset),
    breaks = function(lims) offset + axisTicks(log10(range(lims) - offset), log = TRUE, n = n),
    labels = function(b) format(b, trim = TRUE, drop0trailing = TRUE),
    ...)
}

## Primary x-axis folding the number of model parameters k = m + choose(m, 2)
## into each number-of-predictors (m) tick label as "m\n(k)", rather than a
## secondary axis on top (an earlier design, since replaced everywhere) --
## exact for models that include all two-way interactions; for
## main-effects-only models k = m trivially. Specific to forstmeier_sim,
## not part of the shared plague_virulence original this file was copied
## from.
m_k_scale_x <- function(breaks = 1:6) {
  ggplot2::scale_x_continuous(
    name = "Number of explanatory (predictor) variables",
    breaks = breaks,
    labels = sprintf("%d\n(%d)", breaks, choose(breaks, 2) + breaks)
  )
}

