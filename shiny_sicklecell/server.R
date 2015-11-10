## adapted from:
## https://github.com/rstudio/shiny-examples/blob/master/083-front-page/server.R


if (FALSE) {
    ## local test
    library("shiny")
    runApp("shiny_sickleCell")

    ## deploy (from MacOS side/with R 3.2.2)
    shinyapps::deployApp("shiny_sickleCell")
}

shinyServer(function(input, output) {
  output$main_plot <- renderPlot({
      fit.het <- input$fit.het
      fit.hom <- input$fit.hom
      ## layout parameters
      par(mfrow=c(2,1),mgp=c(2.75,0.75,0),mar=c(3.5,4,0.5,0),oma=c(3,0,0,0))
      par(las=1,bty="l",yaxs="i")
      colvec <- c("orange","blue")
      ## linear decline: from fit.het to fit.hom
      ## S fitness: fit.het at 0% S; fit.hom at 100% S
      Sfit <- function(x) fit.het*(1-x)+fit.hom*x
      Afit <- function(x) fit.het*x+(1-x)*1
      ## stable frequency:
      ## fhet*(1-x)+fhom*x == fhet*x + (1-x)
      ## (-fhet+fhom)*x + fhet == (fhet-1)*x + 1
      ## (-2*fhet+1+fhom)*x = 1-fhet
      ## (fhet-1)/(2*fhet-1-fhom)
      stabfun <- function(x) (x-1)/(2*x-1-fit.hom)
      par(xpd=NA)
      curve(Sfit(x),from=0,to=1,ylim=c(0,2),
            xaxs="i",yaxs="i",
            xlab="frequency of S",
            ylab="average fitness",col=colvec[1])
      curve(Afit(x),add=TRUE,col=colvec[2],from=0)
      par(xpd=FALSE)
      abline(h=fit.het,lty=2)
      abline(v=stabfun(fit.het),lty=2)
      text(c(0.7,0.7),c(Sfit(0.7),Afit(0.7)),col=colvec,
           c("S fitness","A fitness"),pos=4)
      par(xpd=NA)
      curve(stabfun(x),from=1,to=1.5,
            xlab="heterozygote advantage",
            ylab="stable frequency")
      par(xpd=FALSE)
      segments(fit.het,0,fit.het,stabfun(fit.het),lty=2)
      segments(0,stabfun(fit.het),fit.het,stabfun(fit.het),lty=2)
  })
})

## shinyapps::deployApp('path/to/your/app')
