#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(MASS)
library(locpol)
library(ggplot2)
library(ggpubr)

# Define UI
ui1 <- fluidPage(
  
  # Application title
  titlePanel("h_opt(x) exploration by sample size, location and beta parameters"),
  
  # Sidebars
  sidebarLayout(
    sidebarPanel(
      sliderInput("n",
                  "Sample size:",
                  min = 10,
                  max = 10000,
                  value = 5000),
      sliderInput("location",
                  "Location:",
                  min = 0,
                  max = 1,
                  value = 0.5),
      sliderInput("shape1",
                  "Shape 1 beta parameter:",
                  min = 0,
                  max = 10,
                  value = 2),
      sliderInput("shape2",
                  "Shape 2 beta parameter:",
                  min = 0,
                  max = 10,
                  value = 2)
    ),
    
    # Resulting Plot
    mainPanel(
      plotOutput("distPlot"),
    )
  )
)

# Define server logic
server1 <- function(input, output) {
  
  output$distPlot <- renderPlot({
    n <- input$n
    shape1 <- input$shape1
    shape2 <- input$shape2
    x <- input$location
    X <- rbeta(n,shape1,shape2)
    Y <- sin(1/(X/3+0.1))+rnorm(n)
    sig_sq = 1
    # Chosen kernel : Gaussian
    int_z_sqrd_kernel <- sqrt(2*pi)
    int_sqrd_kernel <- sqrt(pi)
    sec_der_m <- function(x) {
      N <- (6*x+1.8)*cos(3/(x+0.3))-9*sin(3/(x+0.3))
      D <- (x+0.3)**4
      return(N/D)
    }
    
    h_opt <- function(x) {
      N <- sig_sq**2*int_sqrd_kernel
      D <- (sec_der_m(x)*int_z_sqrd_kernel)**2*dbeta(x, shape1, shape2)
      return( n**(-1/5)*(N/D)**(1/5))
    }
    
    x_plot <- seq(from=0,to=1,by=0.01)
    h_val_x <- h_opt(x)
    par(mfrow=c(3,1))
    plot(x_plot,h_opt(x_plot),type="l",xlab="Location (x)",ylab="Optimal bandwith")
    abline(h=h_val_x,v=x,col="red",lty=2)
    legend('topleft',legend=parse(text = as.character(h_val_x)),bty='n')
    plot(x_plot,sin(1/(x_plot/3+0.1)),type="l",xlab="x",ylab="m(x)")
    plot(x_plot,dbeta(x_plot,shape1,shape2),type="l",xlab="Location (x)",ylab="Beta Density")
    abline(v=x,col="red",lty=2)
  },height = 600, width = 600 )
}

# Run the application 
shinyApp(ui = ui1, server = server1)
