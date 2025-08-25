library(shiny)

## not sure how to source() in a shiny env?

#' two-group sim for equal-var t-test power/outcome calculation
#' @param n number of simulations
#' @param delta difference between means
#' @param standard deviation of observations
#' @param conf.level confidence level
#' @param seed random-number seed
simfun <- function(n, delta=1, sd=1, conf.level = 0.95, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    x <- rnorm(2*n, mean = rep(c(0,delta), each =n), sd = sd)
    tt <- t.test(x[1:n], x[-(1:n)], conf.level = conf.level, var.equal = TRUE)
    with(tt, c(est = unname(-1*diff(estimate)),
                        lwr = conf.int[1], upr = conf.int[2]))
}    

## how many cases should we distinguish?
## (1) show the effect is small or large
##    * care less about the sign if it's small?

## true effect is positive, small/large
## 
levs <- c("large/clear sign",
          "unclear magnitude/clear sign",
          "small/clear sign",
          "small/unclear sign",
          "NOT (large and opposite est)",
          "unclear")

#' categorize outcomes
catfun <- function(x, s=1) {
    lwr <- x[2]
    upr <- x[3]
    ## adjust for symmetry? ci <- ci*sign(m); m <- abs(m)
    ## case-when?
    if (lwr>s || upr<(-s)) return(levs[1])
    if ((upr>s && lwr>0 && lwr<s) || (lwr<(-s) && upr<0 && upr>(-s))) return(levs[2])
    if ((lwr>0 && upr<s) || (upr<0 && lwr>(-s))) return(levs[3])
    if ((lwr>(-s) && lwr<0 && upr>0 && upr<s) ||
        (upr<s   && upr>0 && lwr>0 && upr>(-s))) return(levs[4])
    if ((lwr<0 && lwr>(-s) && upr>s) || (upr>0 && upr<s && lwr<(-s))) return(levs[5])
    if (lwr<(-s) && upr>s) return(levs[6])
}

proptest <- function(x, s = 1) {
    lwr <- x[,2]
    upr <- x[,3]
    c(lwr_gt_0 = mean(lwr>0),
      lwr_gt_s = mean(lwr>s),
      lwr_gt_negs = mean(lwr>(-s)),
      upr_gt_0 = mean(upr>0),
      upr_gt_s = mean(upr>s),
      upr_gt_negs = mean(upr>(-s)))
}


## TO DO
## - indicate locations of mean, median, mode, match colours between
## table and points:
## http://stackoverflow.com/questions/22684162/r-shiny-color-dataframe
# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
  titlePanel("Extended power simulation"),

  ## Sidebar with slider inputs
  sidebarPanel(
    sliderInput("n",
                "number of samples",
                min = 5,
                max = 100,
                value = 20,
                step = 1),
    sliderInput("delta",
                "difference between means",
                min = 0,
                max = 5,
                value = 1,
                step = 0.25),
    sliderInput("sd",
                "standard deviation",
                min = 0.1,
                max = 2,
                value = 1,
                step = 0.05),
    sliderInput("nsim",
                "number of simulations",
                min = 100,
                max = 50000,
                value = 1000,
                step = 100),
    actionButton("redo", "Rerun sim")
  ),
  mainPanel(
    tableOutput('table')
  )
)

server <- function(input, output) {

  res <- reactive({
    input$redo
    dd1 <- as.data.frame(t(replicate(input$nsim, simfun(n=input$n, delta = input$delta, sd = input$sd))))
    dd1$cat <- apply(dd1, 1, catfun) |> factor(levels = levs)
    list(tab = table(dd1$cat) |> prop.table())
  })
  output$table <- renderTable(res()$tab)
}

shinyApp(ui = ui, server = server)
