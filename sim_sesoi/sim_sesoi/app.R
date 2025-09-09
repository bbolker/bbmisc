library(shiny)
## intended to be run from head of bbmisc repo
source("sim_sesoi/sim_sesoi_funs.R")

# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application title
  titlePanel("Disaggregated power simulation"),

  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      actionButton("redo","Re-run sim"),
      textInput("delta",
                "delta",
                value = "1"),
      textInput("s",
                "SESOI",
                value = "1"),
      textInput("nmin",
                "min n",
                value = "10"),
      textInput("nmax",
                "max n",
                value = "200"),
      textInput("nsteps",
                "n steps",
                value = "20"),
      textInput("nsim",
                "simulations per n",
                value = "10000")
    ) ## sidebarPanel
   ,
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("powerPlot")
    ) ## mainPanel
  ) ## sidebarLayout
) # fluidPage


# Define server logic required to draw a histogram
## server <- function(input, output) {
##   simplot <- NULL
##   vals <- reactiveValues()
##   observe({
##     vals <- lapply(input, as.numeric)
##     nvec <- round(seq(vals$nmin, vals$nmax, length.out = vals$nsteps))
##     sims <- lapply(nvec, tabfun, delta= vals$delta, s = vals$s, nsim = vals$nsim)
##     simdf <- do.call(rbind, sims)
##     attr(sims, "pars") <- with(vals, c(delta = delta, s = s))
##     simplot <- plotfun(f_lengthen(sims))
##   })
##   observeEvent(input$redo, {
##     output$powerPlot <- renderPlot(simplot)
##   })
## }

## version from Claude
server <- function(input, output) {
  library(ggplot2)
  theme_set(theme_bw())
  # Initialize reactive values properly
  vals <- reactiveValues()
  
  # Create a reactive expression for the simulation data
  sim_data <- reactive({
    # Validate inputs exist before proceeding
    req(input$nmin, input$nmax, input$nsteps, input$delta, input$s, input$nsim)
    
    # Convert inputs to numeric
    vals$nmin <- as.numeric(input$nmin)
    vals$nmax <- as.numeric(input$nmax)
    vals$nsteps <- as.numeric(input$nsteps)
    vals$delta <- as.numeric(input$delta)
    vals$s <- as.numeric(input$s)
    vals$nsim <- as.numeric(input$nsim)
    
    # Create sequence and run simulations
    nvec <- round(seq(vals$nmin, vals$nmax, length.out = vals$nsteps))
    sims <- lapply(nvec, tabfun, delta = vals$delta, s = vals$s, nsim = vals$nsim)
    simdf <- do.call(rbind, sims) |> f_lengthen(nvec)
    
    return(simdf)
  })

  # Render the plot initially and when redo button is pressed
  output$powerPlot <- renderPlot({
    plotfun(sim_data())
  })
  
  # Optional: Update plot when redo button is pressed
  observeEvent(input$redo, {
    # Force re-evaluation of reactive expressions
    # This will automatically update the plot if inputs have changed
    output$powerPlot <- renderPlot({
      plotfun(sim_data())
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
