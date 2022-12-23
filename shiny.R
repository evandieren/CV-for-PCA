library(shiny)
library(MASS)
library(locpol)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(mvtnorm)

source("SimulationStudy.R")

# Define UI
ui <- fluidPage(
  
  # Application title
  titlePanel("Cross-Validation on PCA"),
  
  # Sidebars
  sidebarLayout(
    sidebarPanel(
      sliderInput("n",
                  "Observations",
                  min = 20,
                  max = 100,
                  value = 50,
                  step = 5),
      sliderInput("p",
                  "Variables",
                  min = 3,
                  max = 15,
                  value = 8,
                  step = 1),
      sliderInput("r",
                  "Dimension of truncated data set",
                  min = 3,
                  max = 15,
                  value = 3,
                  step = 1),
      sliderInput("K",
                  "CV Folds",
                  min = 5,
                  max = 15,
                  value = 5,
                  step = 5),
      sliderInput("sim",
                  "Simulation runs",
                  min = 3,
                  max = 20,
                  value = 10,
                  step = 1),
      sliderInput("noise",
                  "Noise",
                  min = 0.001,
                  max = 0.1,
                  value = 0.01,
                  step = 0.001),

    ),
    
    # Resulting Plot
    mainPanel(
      plotOutput("distPlot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  output$distPlot <- renderPlot({
    n <- as.numeric(input$n)
    p <- as.numeric(input$p)
    K <- as.numeric(input$k)
    r <- as.numeric(input$r)
    sim <- as.numeric(input$sim)
    noise <- as.numeric(input$noise)
    chosen <- SimulationStudy(WrongPCAImproved, "Wrong PCA Improved", n, p, K, r, sim, noise, meth=T, eigen=F)

    #par(mfrow=c(3,3))
    #for (i in 1:3) {plot(1:p, colMeans(chosen[[4]][[i]]), xlab="Rank r", ylab="Error", type = "b", pch = 19, lty = 1, col = 1)}
    #for (i in 1:3) {plot(1:p, colMeans(chosen[[5]][[i]]), xlab="Rank r", ylab="Error", type = "b", pch = 19, lty = 1, col = 1)}
    #for (i in 1:3) {plot(1:p, colMeans(chosen[[6]][[i]]), xlab="Rank r", ylab="Error", type = "b", pch = 19, lty = 1, col = 1)}
    plot(1:2, 1:2)

  })#,height = 600, width = 600 )
}

# Run the application 
debug(shinyApp)
shinyApp(ui = ui, server = server)


