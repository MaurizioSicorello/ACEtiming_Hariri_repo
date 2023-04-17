N <- 100


############################################
# function to simulate and evaluate regression data
simReg <- function(N, effSize1, effSize2, noise){

  x1 <- rnorm(N)
  x2 <- rnorm(N)
  
  x_mean <- (x1+x2)/2
  
  y <- effSize1*x1 + effSize2*x2 + rnorm(N, 0, noise)
  
  r1 <- lm(y ~ x_mean)
  r2 <- lm(y ~ x1+x2)
  
  an <- anova(r1, r2)
  
  return(an$`Pr(>F)`[2])
}  



############################################
# run model comparisons

iterations <- 1000
p_out <- numeric(iterations)

for(i in 1:iterations){
  
  p_out[i] <- simReg(1000, 0.27, 0.3, 0.3)
  
}

sum(p_out <= 0.05)/length(p_out)


