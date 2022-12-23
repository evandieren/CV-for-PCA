library(shiny)
library(MASS)
library(locpol)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(mvtnorm)

source("Methods/WrongPCA.R")
source("Methods/WrongPCAImproved.R")
source("Methods/Missing_data.R")
source("Methods/KDEApproach.R")
source("Methods/MatrixCompletion.R")

dataset <- function(df,r,noi,n,p){
  df_svd <- svd(df)
  df_svd$d[-(1:r)] <- 0
  last_sv <- df_svd$d[r]
  df <- df_svd$u %*% diag(df_svd$d) %*% t(df_svd$v)
  df_uni_noise <- df + rmvnorm(n = n, mean = rep(0, p), sigma = noi*last_sv*diag(p))
  return(list("Noise"=df_uni_noise))
}

# Define UI
ui <- fluidPage(
  
  # Application title
  titlePanel("CV on PCA"),
  
  # Sidebars
  sidebarLayout(
    sidebarPanel(
      sliderInput("n",
                  "Sample size (n):",
                  min = 100,
                  max = 1000,
                  value = 100,
                  step = 100),
      sliderInput("r",
                  "Rank of the truncated dataset (r):",
                  min = 2,
                  max = 10,
                  value = 3,
                  step = 1),
      sliderInput("K",
                  "Number of folds (k):",
                  min = 2,
                  max = 20,
                  value = 5,
                  step = 1),
      sliderInput("sim",
                  "Number of simulations:",
                  min = 3,
                  max = 20,
                  value = 10,
                  step = 1),
      sliderInput("p",
                  "Number of variables (p):",
                  min = 3,
                  max = 20,
                  value = 7,
                  step = 1),
      sliderInput("noi",
                  "Amount of noise :",
                  min = 0.01,
                  max = 1,
                  value = 0.05,
                  step = 0.001),
      selectInput("method", 
                  "Method of CV: ",
                  choices = c('WrongPCA', 'WrongPCAImproved', 'MissingData', 'MatrixCompletion', 'KDEApproach'),
                  selected = 'WrongPCA'),

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
    r=as.numeric(input$r)
    n=as.numeric(input$n)
    K=as.numeric(input$K)
    sim=as.numeric(input$sim)
    p=as.numeric(input$p)
    noi=as.numeric(input$noi)
    method=input$method
  
    lsEigen0 <- replicate(1, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
    
    lsmethod0 <- replicate(1, matrix(rep(0, sim*p), nrow=sim, ncol = p), simplify=F)
    
    choose0 <- replicate(1, rep(0, p), simplify=F)
    
    
    for (i in 1:sim) {
      samples <- matrix(sample(1:n), ncol = K)
      df <- rmvnorm(n = n, mean = rep(0, p), sigma = diag(p)) 
      data0 <- dataset(df,r,noi,n,p)
      
      for (da in 1:1) {
        lsEigen0[[da]][i,] <- svd(data0[[da]])$d
        
        # WrongPCA, WrongPCAImproved, MissingData, MatrixCompletion, KDEApproach
        if (method == 'KDEApproach'){
          lsmethod0[[da]][i,] <- KDEApproach(data0[[da]], samples)
        }
        if (method == 'WrongPCA'){
          lsmethod0[[da]][i,] <- WrongPCA(data0[[da]], samples)
        }
        if (method == 'WrongPCAImproved'){
          lsmethod0[[da]][i,] <- WrongPCAImproved(data0[[da]], samples)
        }
        if (method == 'MissingData'){
          lsmethod0[[da]][i,] <- MissingData(data0[[da]], samples)
        }
        if (method == 'MatrixCompletion'){
          lsmethod0[[da]][i,] <- MatrixCompletion(data0[[da]], samples)
        }
        
        min_err0 <- which.min(lsmethod0[[da]][i,])
        choose0[[da]][min_err0] <- choose0[[da]][min_err0] +1
        
      }
    }
    
    
    par(mfrow=c(1,2))
    
    plot1 = plot(1:p, colMeans(lsEigen0[[1]]), xlab="kth Eigenvalue", ylab="Value", main="Scree plot data set ", type = "b", pch = 19, lty = 1, col = 1)
    
    plot2 = plot(1:p, colMeans(lsmethod0[[1]]), xlab="Rank r", ylab="Value", main= "Error data set ", type = "b", pch = 19, lty = 1, col = 1)
    
    })
}
# Run the application 
shinyApp(ui = ui, server = server)
