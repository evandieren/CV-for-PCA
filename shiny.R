library(shiny)
library(MASS)
library(locpol)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(mvtnorm)

source("WrongPCA.R")
source("WrongPCAImproved.R")
source("Missing_data.R")
source("KDEApproach.R")
source("MatrixCompletion.R")
source("SimulationStudy2.R")

# Define UI
ui <- fluidPage(
  
  # Application title
  titlePanel("CV on PCA"),
  
  # Sidebars
  sidebarLayout(
    sidebarPanel(
      sliderInput("n",
                  "Sample size:",
                  min = 100,
                  max = 1000,
                  value = 100,
                  step = 100),
      sliderInput("r",
                  "Real r:",
                  min = 2,
                  max = 10,
                  value = 3,
                  step = 1),
      sliderInput("K",
                  "folds:",
                  min = 2,
                  max = 20,
                  value = 5,
                  step = 1),
      sliderInput("sim",
                  "number of simulations:",
                  min = 3,
                  max = 20,
                  value = 10,
                  step = 1),
      sliderInput("p",
                  "possible p:",
                  min = 3,
                  max = 20,
                  value = 7,
                  step = 1),
      sliderInput("noi",
                  "noise:",
                  min = 0.01,
                  max = 1,
                  value = 0.05,
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
    r=input$r
    n=input$n
    K=input$k
    sim=input$sim
    p=input$p
    noi=input$noi
    Simul=SimulationStudy(WrongPCAImproved, n, p, K, r, sim, noi )
    
  })
}
# Run the application 
shinyApp(ui = ui, server = server)


