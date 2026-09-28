
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

## Session-wide default colour/fill palettes, e.g.
##   old <- set_ggplot_palettes(); ...plots...; options(old)
## `discrete` is a vector of colours (ggplot2 falls back to its default hue
## palette if a variable has more levels than colours); `continuous` is the
## name of a scale_colour_*()/scale_fill_*() family, e.g. "viridis"
## (https://stackoverflow.com/questions/53750310)
set_ggplot_palettes <- function(discrete = okabe_ito, continuous = "viridis") {
  invisible(options(ggplot2.discrete.colour = discrete,
                    ggplot2.discrete.fill = discrete,
                    ggplot2.continuous.colour = continuous,
                    ggplot2.continuous.fill = continuous))
}

## Log-axis breaks chosen by base R's rule (axisTicks(..., log = TRUE)):
## powers of 10 for wide ranges, 1-2-5 (or finer) sequences for narrow ones
## (https://stackoverflow.com/questions/14255533)
base_breaks <- function(n = 10) {
  function(x) axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
}

## log10 axes with base_breaks() and plain-number labels (prettyNum()
## formats each label separately, avoiding a shared scientific format)
scale_x_log10_base <- function(n = 10, ...) {
  ggplot2::scale_x_log10(breaks = base_breaks(n), labels = prettyNum, ...)
}
scale_y_log10_base <- function(n = 10, ...) {
  ggplot2::scale_y_log10(breaks = base_breaks(n), labels = prettyNum, ...)
}

## log10 axes with powers-of-10 breaks labelled as 10^k (plotmath superscripts)
## (https://stackoverflow.com/questions/15178081)
scale_x_log10_pow <- function(...) {
  ggplot2::scale_x_log10(breaks = scales::breaks_log(),
                         labels = scales::label_log(), ...)
}
scale_y_log10_pow <- function(...) {
  ggplot2::scale_y_log10(breaks = scales::breaks_log(),
                         labels = scales::label_log(), ...)
}

## Remove minor grid lines, which often sit awkwardly between log-scale breaks
no_minor_grid <- function() {
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

## Direct labels with every label shifted by (dx, dy) cm, e.g.
##   direct_label_shift(p, "first.qp", dx = -0.1)
## dl.trans() captures unevaluated expressions, so the shifts are inlined
## with bquote() rather than passed as variables
## (https://stackoverflow.com/questions/17839153)
dl_shift <- function(dx = 0, dy = 0) {
  eval(bquote(directlabels::dl.trans(x = x + .(dx), y = y + .(dy))))
}
direct_label_shift <- function(p, method = "last.qp", dx = 0, dy = 0) {
  directlabels::direct.label(p, list(method, dl_shift(dx, dy)))
}

## Make tiles/rasters fill each panel to its edges (no axis expansion);
## for tiles that fill every facet when facets cover different ranges,
## also use scales = "free" in facet_wrap()/facet_grid()
## (https://stackoverflow.com/questions/7001710)
fill_panels <- function(...) {
  ggplot2::coord_cartesian(expand = FALSE, ...)
}

