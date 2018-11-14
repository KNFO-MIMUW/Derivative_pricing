library(shiny)

MC_discrete <- function(n) {
  dW_t <- rnorm(n[6])
  W_t <- 1:n[6]
  W_t[1] = n[1]*exp((n[2]-(n[3]^2)/2)*n[4]/n[6]+ n[3]*sqrt(n[4]/n[6])*dW_t[1])
  for(i in 2:n[6]){
    W_t[i] = W_t[i-1] * (exp((n[2] - ((n[3]^2)/2))*n[4]/n[6] + n[3]*sqrt(n[4]/n[6])*dW_t[i]))
  } 
  X = exp(-n[2]*n[4])*(max(mean(W_t)-n[5],0)) - n[9]*(exp(-n[2]*n[4])*(max(funct2(W_t)-n[5],0)) - n[8])
  return(X)
}

funct1 <- function(n) {
  dW_t <- rnorm(n[6])
  W_t <- 1:n[6]
  W_t[1] = n[1]*exp((n[2]-(n[3]^2)/2)*n[4]/n[6]+ n[3]*sqrt(n[4]/n[6])*dW_t[1])
  for(i in 2:n[6]){
    W_t[i] = W_t[i-1] * (exp((n[2] - ((n[3]^2)/2))*n[4]/n[6] + n[3]*sqrt(n[4]/n[6])*dW_t[i]))
  }
  return(W_t)
}


funct2 = function(y, na.rm = TRUE) {
  exp(sum(log(y[y > 0]), na.rm = na.rm) / length(y))
}


asian_call <- function(S_0, r, sigma, dT, K, n, M) {
  v0 <- c(S_0, r, sigma, dT, K, n, M)
  St_1 <- matrix(rep(v0, M), 7, 10000)
  W_t1 <- apply(St_1, 2, funct1)
  X_1 <- exp(-r * dT) * pmax(colMeans(W_t1) - K, 0)
  X_2 <- exp(-r * dT) * pmax(apply(W_t1, 2, funct2) - K, 0)
  theta <- cov(X_1, X_2) / var(X_2)
  
  price1 <- exp(-r*dT+(r-(sigma^2)/2)*dT*(n+1)/(2*n)+(sigma^2)*dT*(n+1)*(2*n+1)/(12*(n^2)))*S_0*pnorm((log(S_0/K)+(r-(sigma^2)/2)*dT*(n+1)/(2*n)+(sigma^2)*dT*(n+1)*(2*n+1)/(6*(n^2)))/(sqrt(dT*(n+1)*(2*n+1)/6)*sigma/n))-exp(-r*dT)*K*pnorm((log(S_0/K)+(r-(sigma^2)/2)*dT*(n+1)/(2*n))/(sqrt(dT*(n+1)*(2*n+1)/6)*sigma/n))
  v <- c(S_0, r, sigma, dT, K, n, M, price1, theta)
  S_t <- matrix(rep(v, M), 9, M)
  Z1 <- apply(S_t, 2 , MC_discrete)
  
  price2 <- exp(-(r+(sigma^2)/6)*dT/2)*S_0*pnorm((log(S_0/K)+(r+(sigma^2)/6)*dT/2)/(sigma*sqrt(dT/3))) - exp(-r*dT)*K*pnorm( (log(S_0/K)+(r+(sigma^2)/6)*dT/2)/(sigma*sqrt(dT/3)) - sigma*sqrt(dT/3))
  v <- c(S_0, r, sigma, dT, K, n, M, price2, theta)
  S_t <- matrix(rep(v, M), 9, M)
  Z2 <- apply(S_t, 2 , MC_discrete)
  avg <- mean(Z1)
  Beta1 <- sum(1.96 * sqrt((avg - Z1) ^ 2) / (sqrt((M - 1) * M)))
  # print("Srednia dyskretna. Estymacja C0:")
  # print(avg)
  # print("Przedzial ufnosci (95%):")
  # print(avg - Beta1)
  # print(avg + Beta1)
  avg2 <- mean(Z2)
  Beta2 <- sum(1.96 * sqrt((avg2 - Z2) ^ 2) / (sqrt((M - 1) * M)))
  # print("Srednia ciagla. Estymacja C0:")
  # print(avg2)
  # print("Przedzial ufnosci (95%):")
  # print(avg2 - Beta2)
  # print(avg2 + Beta2)
  wynik <- c(avg, avg - Beta1,avg + Beta1,avg2,avg2 - Beta2, avg2 - Beta2)
  names(wynik) <-c("Srednia dyskretna", "CI: Dolna granica ", "CI: Gorna granica","Srednia ciagla", "CI: Dolna granica", "CI: Gorna granica")
  return(wynik)
}


MonteCarlo <- function(S_0, r, sigma, dT, K, M){
  
  vector <- 1:M
  for(i in 1:M){
    vector1 <- runif(1)
    vector2 <- runif(1)
    vector[i] <- sqrt(-2*log(vector1))*cos(2*pi*vector2)
  }
  V_t = pmax(exp(-r*dT)*(K-S_0*exp((r-(sigma^2)/2)*dT+sigma*sqrt(dT)*vector)),0)
  aver <-mean(V_t)
  value <- exp(-r*dT)*pnorm(-(log(S_0/K)+r*dT-1/2*dT*sigma^2)/(sigma*sqrt(dT)))*K - S_0*pnorm(-(log(S_0/K)+r*dT+1/2*dT*sigma^2)/(sigma*sqrt(dT)))
  error <- (abs(value-mean(V_t)))/value
  beta <-sum(1.96*(aver-V_t)^2/((M-1)*sqrt(M)))
  # print("V_t:")
  # print(aver)
  # print("Przedzial ufnosci (95%):")
  # print(aver - beta)
  # print(aver + beta)
  # print("Blad wzgledny:")
  # print(error)
  wynik <- c(aver, error, aver - beta, aver + beta)
  names(wynik) <- c("MC AVG", "MC Error", "MC CI Dolna granica", "MC CI Gorna granica")
  return(wynik)
}


halton <- function(k){
  
  f<-1
  r1<-0
  m<-k
  while(k>0){
    f<-f/2
    r1<- r1+f*(k%%2)
    k<-k%/%2
  }
  g<-1
  r2<-0
  
  while(m>0){
    g<-g/3
    r2<- r2+g*(m%%3)
    m<-m%/%3
  }
  result <- 1:2
  result[1] <- r1
  result[2] <- r2
  return(result)
}

QuasiMonteCarlo <- function(S_0, r, sigma, dT, K, M){
  
  vector <- 1:M
  for(i in 1:M){
    randomv <- halton(i)
    vector[i] <- sqrt(-2*log(randomv[1]))*cos(2*pi*randomv[2]) 
  }
  V_t = pmax(exp(-r*dT)*(K-S_0*exp((r-(sigma^2)/2)*dT+sigma*sqrt(dT)*vector)),0)
  
  value <- exp(-r*dT)*pnorm(-(log(S_0/K)+r*dT-1/2*dT*sigma^2)/(sigma*sqrt(dT)))*K - S_0*pnorm(-(log(S_0/K)+r*dT+1/2*dT*sigma^2)/(sigma*sqrt(dT)))
  aver <-mean(V_t)
  beta <-sum(1.96*(aver-V_t)^2/((M-1)*sqrt(M))) 
  error <- abs((value-mean(V_t)))/value
  # print("V_t:")
  # print(aver)
  # print("Blad wzgledny:")
  # print(error)
  wynik <- c(aver, error)
  names(wynik) <- c("QM AVG", "QM Relative Error")
  return(wynik)
}




heston <- function(S0, V0, r, sigma, a, b, rho, dT, K, n, M) {
  vect <- c(S0, V0, r, sigma, a, b, rho, dT, K, n, M)
  St <- matrix(rep(vect, M), 11, M)
  Zt <- apply(St, 2 , milstein)
  C0 = pmax(exp(-r * dT) * (Zt - K), 0)
  srednia <- mean(C0)
  betainterval <- sum(1.96 * sqrt((srednia - C0) ^ 2) / (sqrt((M - 1) *
                                                                M)))
  #print("Estymacja C0:")
  #print(srednia)
  #print("95% Przedzial ufnosci:")
  #print(srednia - betainterval)
  #print(srednia + betainterval)
  wynik <- c(srednia,srednia - betainterval,srednia + betainterval)
  names(wynik) <-c("Srednia", "Dolna granica", "Gorna granica")
  print(wynik)
  return(wynik)
}


milstein <- function(v) {
  Z1 <- rnorm(v[10])
  Z2 <- rnorm(v[10])
  dWt = Z1
  dVt = Z1 * v[7] + Z2 * sqrt(1 - v[7] ^ 2)
  Vt <- 1:v[10]
  Wt <- 1:v[10]
  Vt[1] = v[2] + ((v[5] * (v[6] - v[2]) - (v[4] ^ 2) / 4) * v[8] / v[10]) + v[4] *
    sqrt(v[2]) * sqrt(v[8] / v[10]) * dVt[1] + (dVt[1] ^ 2) * (v[8] / v[10]) *
    (v[4] ^ 2) / 4
  Wt[1] = v[1] * exp((v[3] - max(v[2], 0) / 2) * v[8] / v[10] + sqrt(max(v[2], 0)) *
                       sqrt(v[8] / v[10]) * dWt[1])
  
  for (i in 2:v[10]) {
    Vt[i] = Vt[i - 1] + ((v[5] * (v[6] - max(Vt[i - 1], 0)) - (v[4] ^ 2) / 4) *
                           v[8] / v[10]) + v[4] * sqrt(max(Vt[i - 1], 0)) * sqrt(v[8] / v[10]) * dVt[i] +
      (dVt[i] ^ 2) * (v[8] / v[10]) * (v[4] ^ 2) / 4
    Wt[i] = Wt[i - 1] * (exp((v[3] - (max(
      Vt[i - 1], 0
    ) / 2)) * v[8] / v[10] + sqrt(max(Vt[i - 1], 0) * v[8] / v[10]) * dWt[i]))
  }
  return(Wt[v[10]])
}






# Define UI for dataset viewer app ----
ui <- fluidPage(
  # App title ----
  titlePanel("Option Price Calc"),
  
  # Sidebar layout with input and output definitions ----
  fluidRow(
    
    # Sidebar panel for inputs ----
    
  
      # Input: Select a dataset ----
    column(3,  
    selectInput("dataset", "Choose a model:",
                  choices = c("Heston", "CIR", "BS")))),
  conditionalPanel( condition = "input.dataset == 'Heston'",
                    
  fluidRow(
      column(3,            
      # Input: Specify the S0 ----
      numericInput("S0", "S0:", 100)),

      # Input: Specify the V0 ----
      column(3,     
      numericInput("V0", "V0:", 0.09)),
      
      column(3,     
      # Input: Specify the r ----
      numericInput("r", "r:", 0.05))),
  fluidRow(
    column(3,
      # Input: Specify the sigma ----
      numericInput("sigma", "sigma:", 1)),
      
    column(3,
      # Input: Specify the a ----
      numericInput("a", "a:", 2)),
      
    column(3,
      # Input: Specify the b ----
      numericInput("b", "b:", 0.09))),
  
  fluidRow(    
    column(3,
      # Input: Specify the rho ----
      numericInput("rho", "rho:", -0.3)),
      
    column(3,
      # Input: Specify the dT ----
      numericInput("dT", "Time:", 5)),
    
    column(3,
      # Input: Specify the K ----
      numericInput("K", "K:", 100))),
    
  fluidRow(
    column(3,
      # Input: Specify the n ----
      numericInput("n", "n:", 500)),
      
    column(3,
      # Input: Specify the M ----
      numericInput("M", "M:", 1000)))
  ),
  conditionalPanel( condition = "input.dataset == 'CIR'",
                    
                    fluidRow(
                      column(3,            
                             # Input: Specify the S0 ----
                             numericInput("cirS0", "S0:", 50)),
                      
                      # Input: Specify the cirr ----
                      column(3,     
                             numericInput("cirr", "r:", 0.05)),
                      
                      column(3,     
                             # Input: Specify the cirsigma ----
                             numericInput("cirsigma", "sigma:", 0.3))),
                    fluidRow(
                      column(3,
                             # Input: Specify the cirdT ----
                             numericInput("cirdT", "dT:", 2)),
                      
                      column(3,
                             # Input: Specify the cirK ----
                             numericInput("cirK", "K:", 50)),
                      
                      column(3,
                             # Input: Specify the cirn ----
                             numericInput("cirn", "n:", 200))),
                    
                    fluidRow(    
                      column(3,
                             # Input: Specify the cirM ----
                             numericInput("cirM", "M:", 10000)))
                      

                    ),
    conditionalPanel( condition = "input.dataset == 'BS'",
                      
                      fluidRow(
                        column(3,            
                               # Input: Specify the BSS0 ----
                               numericInput("BSS0", "S0:", 50)),
                        
                        # Input: Specify the BSr ----
                        column(3,     
                               numericInput("BSr", "r:", 0.05)),
                        
                        column(3,     
                               # Input: Specify the BSsigma ----
                               numericInput("BSsigma", "sigma:", 0.3))),
                      fluidRow(
                        column(3,
                               # Input: Specify the BSdT ----
                               numericInput("BSdT", "Time:", 0.5)),
                        
                        column(3,
                               # Input: Specify the BSK ----
                               numericInput("BSK", "K:", 50)),
                        
                        column(3,
                               # Input: Specify the BSM ----
                               numericInput("BSM", "M:", 1000)))                      
      
                      ),
      
      # Include clarifying text ----
      #helpText("Note: while the data view will show only the specified",
      #        "number of observations, the summary will still be based",
      #         "on the full dataset."),
      
      # Input: actionButton() to defer the rendering of output ----
      # until the user explicitly clicks the button (rather than
      # doing it immediately when inputs change). This is useful if
      # the computations required to render output are inordinately
      # time-consuming.
  fluidRow(
    column(3,
           actionButton("update", "Calculate"))
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
  

      # Output: Header + table of distribution ----
      h4("Option Price"),
      tableOutput("PriceAnalysis")

    )
    
  
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output, session) {
  
  # Return the requested dataset ----
  # Note that we use eventReactive() here, which depends on
  # input$update (the action button), so that the output is only
  # updated when the user clicks the button
  datasetInput <- eventReactive(input$update, {
    switch(input$dataset,
           "rock" = rock,
           "pressure" = pressure,
           "cars" = cars)
  }, ignoreNULL = FALSE)
  
  argumentsUpdate <- eventReactive(input$update, {
    if(input$dataset == 'Heston'){
    
    
    h <- heston(input$S0, input$V0, input$r, input$sigma, input$a, input$b,
                input$rho, input$dT, input$K, isolate(input$n), input$M)
    }
    if(input$dataset == 'CIR'){
      
      
      h <- asian_call(input$cirS0, input$cirr, input$cirsigma, input$cirdT, input$cirK, input$cirn, input$cirM)
      
    }
    if(input$dataset == 'BS'){
      
      
      h <- c(MonteCarlo(input$BSS0,input$BSr,input$BSsigma,input$BSdT,input$BSK,input$BSM),
             QuasiMonteCarlo(input$BSS0,input$BSr,input$BSsigma,input$BSdT,input$BSK,input$BSM))
      
      
    }
    return(h)
  }, ignoreNULL = FALSE)
  # Show the first "n" observations ----
  # The use of isolate() is necessary because we don't want the table
  # to update whenever input$obs changes (only when the user clicks
  # the action button)
  output$PriceAnalysis <- renderTable({
   head(t(argumentsUpdate()))
   })
  
}

# Create Shiny app ----
shinyApp(ui, server)
