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
      actionButton("redo","Run sim"),
      textInput("delta",
                "delta;",
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
server <- function(input, output) {
  simplot <- NULL
  vals <- reactiveValues()
  observe({
    vals <- lapply(input, as.numeric)
    nvec <- round(seq(vals$nmin, vals$nmax, length.out = vals$nsteps))
    sims <- lapply(nvec, tabfun, delta= vals$delta, s = vals$s, nsim = vals$nsim)
    simdf <- do.call(rbind, sims)
    attr(sims, "pars") <- with(vals, c(delta = delta, s = s))
    simplot <- plotfun(f_lengthen(sims))
  })
  observeEvent(input$redo, {
    output$powerPlot <- renderPlot(simplot)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
