## http://shiny.rstudio.com/tutorial/lesson1/
library("shiny")

# Define UI for application that draws a histogram
shinyUI(fluidPage(

  # Application title
  titlePanel("Sickle-cell equilibrium frequency"),

  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("fit.het",
                  "heterozygote (AS) advantage",
                  min = 1,
                  max = 2,
                  value = 1.15),
      sliderInput("fit.hom",
                  "homozygote (SS) fitness",
                  min = 0,
                  max = 1,
                  value = 0)
    ),


    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("main_plot",width="400px"),
        p(paste("The top plot shows the fitness of the wild-type",
                "(A, blue) and sickle-cell alleles (S, orange)",
                "as a function of sickle-cell allele frequency in",
                "the population; S fitness is highest when S is rare,",
                "while A fitness is highest when S is common")),
        p(paste("The bottom plot shows the stable frequency (i.e,",
                "the value of S frequency where the lines cross in the upper plot)",
                "as a function of heterozygote advantage (the fitness of",
                "a heterozygote relative to the wild-type); the dashed",
                "lines show the stable frequency for the current value",
                "of heterozygote advantage."))
)
      
)))
