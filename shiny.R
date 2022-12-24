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
  p("This Shiny app is intended to give an impression on how the different 
            paramters in a PCA analysis come into play. Keep in mind that simulations take a bit of time
            and that the truncated rank of the data set has to be smaller then the amount of
             variables. The noise in % indicates the amount of noise added with respect to the last untruncated
    singular value of the initial base data set."),
  # Sidebars
  sidebarLayout(
    sidebarPanel(
      selectInput("method", 
                  "Method of CV",
                  choices = names(list("Wrong PCA"=WrongPCA, "Wrong PCA Improved"=WrongPCAImproved, 
                                       "Missing Data"=MissingData, "Matrix Completion"=MatrixCompletion, 
                                       "KDE Approach"=KDEApproach)),
                  selected = 'Wrong PCA'),
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
      sliderInput("hnoise",
                  "High noise in %",
                  min = 0.01,
                  max = 0.2,
                  value = 0.02,
                  step = 0.01),
      sliderInput("lnoise",
                  "Low noise in %",
                  min = 0.0001,
                  max = 0.0101,
                  value = 0.001,
                  step = 0.001),
      sliderInput("dnoise",
                  "Differing noise in %",
                  min = 0.001,
                  max = 0.0201,
                  value = 0.005,
                  step = 0.002),


    ),
    
    # Resulting Plot
    mainPanel(
      plotOutput("distPlot"),
    ),
  ),
)

# Define server logic
server <- function(input, output) {
  ls <- list("Wrong PCA"=WrongPCA, "Wrong PCA Improved"=WrongPCAImproved, "Missing Data"=MissingData, "Matrix Completion"=MatrixCompletion, "KDE Approach"=KDEApproach)
  
  output$distPlot <- renderPlot({
    method <- input$method
    n <- input$n
    p <- input$p
    K <- input$K
    r <- input$r
    sim <- input$sim
    hn <- input$hnoise
    ln <- input$lnoise
    dn <- input$dnoise
    
    noise <- c(hn, ln, dn)
    set.seed(1312)
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
             col=color_p,lty = 1, cex=1.2)
    }
    
    plot_three <- function(method_str,dataset_number){
      load(paste0("datasets_plots/",method_str,".Rdata"))
      errors <- chosen[[4+dataset_number]]
      p <- dim(errors[[1]])[2]
      color_p <- c("orange", "darkblue","darkgray")
      for (j in 1:3){
        plot(1:p, colMeans(errors[[j]]), xlab="Rank r", ylab="Error", main=TeX(paste0(method_str," error on $D^",j,"_",dataset_number,"$"), bold=T),type="b",lwd = 2,col=color_p[j])
      }
    }
    
    if (method == "KDE Approach") {
      par(mfrow=c(3,3))
      plot_three("KDE Approach",0)
      plot_three("KDE Approach",1)
      plot_three("KDE Approach",2)
      
    }
    else{
      par(mfrow=c(3,1))
      plot_wrt_covariance(method,1)
      plot_wrt_covariance(method,2)
      plot_wrt_covariance(method,3)
      # title(paste0("Error of ", method), line = - .9, outer = TRUE)
    }
    
  },height = 650, width = 550)
}

# Run the application 
shinyApp(ui = ui, server = server)


