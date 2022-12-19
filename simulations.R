
######################
# to do
# the two comparison models are not equivalent


# ALTERNATIVE SOLUTION:
# https://cborchers.com/2021/04/22/permutation-test-for-f-score-differences-in-python/
# https://arxiv.org/pdf/cs/0008005.pdf



######################
# packages
library("MASS")   
library("caret")
library("here")
library("party")

source(here("functions", "randomForest_functions.R"))

######################
# create functions


###########
# function to create data uncorrelated data and predictable outcome

createData <- function(nPred = 160, nSubj = 141, permute = FALSE){
  
  data <- mvrnorm(n = nSubj,        # Create random data
                  mu = rep(0, nPred),
                  Sigma = matrix(diag(nPred), nrow = nPred)
  )
  data <- as.data.frame(data)
  names(data) <- paste0(rep("x", nPred), 1:nPred)
  
  #b_coef <- c(0.1, 0.2, 0.3, -0.2, -0.3, rep(0, nPred-5))
  b_coef <- c(0.3, rep(0, nPred-1))
  
  y <- as.matrix(data) %*% b_coef + rnorm(nSubj)
  
  if(permute == TRUE){y <- sample(y)}
  
  data <- data.frame(y, data)
  
  return(data)
  
}


# dataframe to store results
iterations = 2
dfOut <- data.frame(matrix(nrow = iterations, ncol = 4))



for(i in 1:iterations){
  
  # create data
  df <- createData(nPred = 160)
  
  # define predictor sets
  set_small <- names(df)[2:11]
  set_large <- names(df)[7:length(names(df))]
  predSets <- list(set_small, set_large)
  
  RFmodel_sq <- RFmain("y",
                    predictorSets = predSets,
                    mtryArg = "sqroot")
  
  RFmodel_ba <- RFmain("y",
                       predictorSets = predSets,
                       mtryArg = "bagging")
  
  dfOut[i, ] <- c(RFmodel_sq, RFmodel_ba)
  
  
}









