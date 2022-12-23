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
                  min = 50,
                  max = 500,
                  value = 100,
                  step = 10),
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
                  max = 10,
                  value = 5,
                  step = 5),
      sliderInput("sim",
                  "Simulation runs",
                  min = 1,
                  max = 20,
                  value = 5,
                  step = 1),
      # sliderInput("noise",
      #             "Noise",
      #             min = 0.001,
      #             max = 0.1,
      #             value = 0.01,
      #             step = 0.001),
      selectInput("method", 
                  "Method of CV: ",
                  choices = names(list("Wrong PCA"=WrongPCA, "Wrong PCA Improved"=WrongPCAImproved, 
                                       "MissingData"=MissingData, "Matrix Completion"=MatrixCompletion, 
                                       "KDE Approach"=KDEApproach)),
                  selected = 'Wrong PCA'),

    ),
    
    # Resulting Plot
    mainPanel(
      plotOutput("distPlot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  ls <- list("Wrong PCA"=WrongPCA, "Wrong PCA Improved"=WrongPCAImproved, "MissingData"=MissingData, "Matrix Completion"=MatrixCompletion, "KDE Approach"=KDEApproach)
  
  output$distPlot <- renderPlot({
    n <- input$n
    p <- input$p
    K <- input$K
    r <- input$r
    sim <- input$sim
    method <- input$method
    
    print(ls[which(names(ls)==method)])
    print(method)
    print(str(method))
    # noise <- as.numeric(input$noise)
    noise <- c(0.02,0.001,0.005)
    chosen <- SimulationStudy(ls[[method]], method, n, p, K, r, sim, noise, meth=T, eigen=F)
    
    plot_wrt_covariance <- function(method_str,noise){
      
      #load(paste0("datasets_plots/",method_str,".Rdata"))
      
      errors <- list(chosen[[4]],chosen[[5]],chosen[[6]]) # D_0,D_1,D_2
      p <- dim(errors[[1]][[1]])[2]
      
      color_p <- c("orange", "darkblue","darkgray")
      # compute ylim
      ymin <- 0
      ymax <- -Inf
      for (j in 1:3){
        candid <- log(colMeans(errors[[j]][[noise]]))
        if (ymin > min(candid)){
          ymin <- min(candid)*1.05
        }
        if (ymax < max(candid)){
          ymax <- max(candid)*1.05 
        }
      }
      
      plot(1:p, log(colMeans(errors[[1]][[noise]])), xlab="Rank r", ylab="Error in logscale", main=TeX(paste0(method_str," error with Noise n°", noise), bold=T),type="b",cex.main=1.3, cex.lab=1.5,cex.axis=1.4,lwd = 2,ylim = c(ymin,ymax),col=color_p[1])
      for (j in 2:3){
        lines(1:p, log(colMeans(errors[[j]][[noise]])), xlab="Rank r", ylab="Error",type="b",lwd = 2,col=color_p[j])
      }
      
      legend("topright", legend=c(TeX("$D_0$"),TeX("$D_1$"),TeX("$D_2$")),
             col=color_p,lty = 1, cex=0.8)
      
    }
    
    par(mfrow=c(1,3))
    plot_wrt_covariance(method,1)
    plot_wrt_covariance(method,2)
    plot_wrt_covariance(method,3)
    # title(paste0("Error of ", method), line = - .9, outer = TRUE)

    
  },height = 400, width = 750)
}

# Run the application 
shinyApp(ui = ui, server = server)


