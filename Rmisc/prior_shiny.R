## https://shiny.posit.co/r/reference/shiny/0.13.1/plotoutput.html
## https://mastering-shiny.org/action-graphics.html


#' draw prior samples from a normal distribution
#' @param nm (character) name of the parameter
#' @param vals (numeric) vector of values (ci range)
#' @param level (numeric) confidence level for the normal distribution
#' @param n (numeric) number of samples to draw
norm_samples <- function(nm, vals, level = 0.95, n = 1000) {
  breadth <- 2* qnorm((1+level)/2)
  sd <- diff(vals)/breadth
  m <- mean(vals)
  res <- rnorm(n, m, sd)
  if (grepl("^log_", nm)) res <- exp(res)
  res
}

#' @inheritParams norm_samples
#' @param vals vector of (min, max)
unif_samples <- function(nm, vals, n = 1000) {
  res <- runif(n, min = vals[1], max = vals[2])
  if (grepl("^log_", nm)) res <- exp(res)
  res
}

#' @param prior_ranges named list of range vectors
#' @param seed random-number seed
#' @param sampfun function for sampling (e.g. norm_samples, unif_samples)
draw_prior_samples <- function(prior_ranges, seed = NULL, sampfun = norm_samples, ...) {
  if (!is.null(seed)) set.seed(seed)
  Map(\(nm, vals) sampfun(nm, vals, ...),
      names(prior_ranges), prior_ranges) |>
    do.call(what = cbind)
}

#' @param ncurve Number of curves to draw
#' @param xvec (numeric) vector of x-values
#' @param prior_vals (data.frame) 
run_shiny <- function(ranges,
                      curve_fun,
                      xvec = seq(15, 30, length = 51),
                      ncurve = 100,
                      seed = NULL) {
  require(shiny)
  require(ggplot2); theme_set(theme_bw()) ## FIXME, don't hard-code?

  if (!is.null(seed)) { set.seed(seed) }
  curve_ids <- seq(ncurve)

  prior_vals <- draw_prior_samples(ranges, n = ncurve)
  prior_curves <- purrr::map_dfr(seq(ncurve),
                               ~ dplyr::tibble(xvec,
                                               val = curve_fun(x = xvec, prior_vals[.x, ])),
                                               .id = "curve_id")

  
  shinyApp(
    ui = fluidPage(
      ## textInput( 
      ##   "formula", 
      ##   "Function", 
      ##   placeholder = "Enter text..."
      ## ), 
      plotOutput("plot", brush = "plot_brush"),
      actionButton("reset_selection", "Reset selection"),
      actionButton("save_selection", "Save selection")
    ),
  
  server = function(input, output, session) {
    selected <- reactiveVal(rep(FALSE, nrow(prior_vals)))

    observeEvent(input$plot_brush, {
      brushed <- brushedPoints(prior_curves, input$plot_brush)
      any_brushed <- unique(brushed$curve_id)
      selected(curve_ids  %in% any_brushed | selected())
    })
    
    observeEvent(input$reset_selection, {
      selected(rep(FALSE, ncurve))
    })

    observeEvent(input$save_selection, {
      write.csv(prior_vals[selected(),], file = "output.csv")
    })

    output$plot <- renderPlot({
    prior_curves$selected <- prior_curves$curve_id %in% which(selected())
    ggplot(prior_curves, aes(xvec, val)) + 
      geom_line(aes(group=curve_id, colour = selected, alpha = selected)) +
      scale_colour_manual(breaks = c("TRUE", "FALSE"), values = c("black", "red")) +
      scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c(0.05, 1))
    }, res = 96)
  }
  )
}

## FIXME: generalize all of this
## source("tpc_funs.R") ## prior_ranges, tpc_prob
## prior_ranges <- c(prior_ranges, list(C=c(1,1)))
## run_shiny(ranges = prior_ranges, curve_fun = \(x, p) tpc_prob(pars = p, Tvec = x))

## Example from McElreath 2d ed. Figure 4.5

mcel_prior_ranges <- list(a = c(138, 218), b = c(-20, 20))
mcel_curve_fun <- \(x, p) p[1] + p[2]*(x-mean(x))
run_shiny(mcel_prior_ranges, mcel_curve_fun, xvec = 30:70, seed = 101)

## FIXME:
## user-specified formula (?? eval ??? [safety??])
## user-specified parameter ranges? (normal, uniform?)
## user-specified output file name?
## zoom in on selected range (coord_cartesian) [optional based on radio/action button?]
## print MVN summaries (or other) of selected values?
## x, y-axis labels?

## utilities to plot selected values (pairs plot), find mean/covariance matrix (or other distributions?)


