
logit_breaks <- function(n = 6) {
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

