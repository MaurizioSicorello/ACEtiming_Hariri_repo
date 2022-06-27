library(psych)


#####################################################################################
# This script contains three important functions 

# The function ACE.struct() takes information on the data frame as arguments and outputs a dataframe of the variable "age" and the
# correlations between trauma during that age and a dependent variable

# The function ACE.eval() takes the output of ACE.struct() and calculates the best fitting trend-line up to a third-order polynomial

# ACE.plot() takes the output of ACE.eval() and creates a plot


#####################################################################################

# function 1 should output age and effect size for a given variable
# arguments: 
#DV = name of dependent variable 
#first_pred = name of first predictors in quotation marks (columns of dataframe have to be sorted!)
#start_age = first age to evaluate
#end_age = last age to evaluate

#möglicherweise könnte ich noch eine Sortierfunktion einbauen, damit die Jahre auf jeden Fall richtig sortiert sind


ACE.struct <- function(DV, first_pred, start_age = 1, end_age = 17, data = df){
  
  df <- data 
  age <- c(start_age : end_age)
  start_col <- which(names(df) == first_pred)
  
  ACE_effect <- NULL
  
  for(i in 1:(end_age - start_age + 1)){
    
    ACE_effect[i] <- cor(df[, start_col + i - 1], data[, DV])
  }
  
  data.frame("age" = age, "effects" = ACE_effect)
}




# function 2: helper function to calculate corrected AIC (called AICc) from an LM
AICc <- function(lm.model){
  
  k <- length(lm.model$coefficients)
  n <- length(lm.model$residuals)
  AIC <- AIC(lm.model)
  
  AICc <- AIC + (2*k^2 + 2*k)/(n - k - 1)
  AICc
}



# function 3: chooses the best fitting model between zero order and 3rd order polynomial
ACE.eval <- function(ACE.data){
  
  age1 <- ACE.data$age
  age2 <- (age1-mean(age1))^2
  age3 <- (age1-mean(age1))^3
  
  effects <- ACE.data$effects
  
  pol0 <- lm(effects ~ 1)
  pol1 <- lm(effects ~ age1)
  pol2 <- lm(effects ~ age1 + age2)
  pol3 <- lm(effects ~ age1 + age2 + age3)
  
  AICc_vector <- c(AICc(pol0), AICc(pol1), AICc(pol2), AICc(pol3))
  
  best.polynomial <- switch(which(AICc_vector == min(AICc_vector)), "Zero Order Polynomial", "First Order Polynomial", "Second Order Polynomial", "Third Order Polynomial")
  best.model <- switch(which(AICc_vector == min(AICc_vector)), pol0, pol1, pol2, pol3)
  best.AICc <- switch(which(AICc_vector == min(AICc_vector)), AICc(pol0), AICc(pol1), AICc(pol2), AICc(pol3))
  
  output <- list("Best.fit" = best.polynomial, "Model.Parameters" = best.model, "AICc" = best.AICc, "All.AICc" = AICc_vector, "data" = ACE.data)
  output
}





# function 4

ACE.plot <- function(ACE.model, N = 66, y_axis = 0.4, zero_line = TRUE, CI_margins = TRUE, trend_line = TRUE, title = NULL){
  
  #confidence interval of correlation r = 0f for margins
  SE <- 1/sqrt(N - 2)
  UL <- fisherz2r(1.96*SE)
  LL <- fisherz2r(-1.96*SE)
  
  #main plot
  plot(effects~age, type = "b", ylim = c(-y_axis, y_axis), data = ACE.model$data, main = title)
  
  if(zero_line == TRUE){
  #zero line
  abline(a = 0, b = 0, col = "red", lwd = 1.5)
  }
  
  if(CI_margins == TRUE){
  #confidence margins
  abline(a = UL, b = 0, col = "black", lty = 2)
  abline(a = LL, b = 0, col = "black", lty = 2)
  }
  
  if(trend_line == TRUE){
  #best fitted trend
  pred <- predict(ACE.model$Model.Parameters, newdata = ACE.model$data)
  with(ACE.model$data, lines(x = age, y = pred, col = "blue"))
  }
}










