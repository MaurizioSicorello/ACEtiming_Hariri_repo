##########################
# main random forest function

RFmain <- function(DV, 
                             data = df,
                             predictorSets = NULL,
                             numtree = 1000,
                             mtryArg = "sqroot",
                             permute_DV = FALSE, 
                             include.imp = FALSE,
                             saveModels = FALSE){
  
  multiMod <- ifelse(is.list(predictorSets) & length(predictorSets) > 1, TRUE, FALSE)
  
  if(multiMod){
    numSets <- length(predictorSets)
  }else{
    numSets <- 1
  }
  
  RFresults <- numeric(numSets)
  RFmodels <- list()
  
  if(include.imp == T & multiMod){
    warning("variable importance only supported if one single model is given")
  }
  
  if(permute_DV == TRUE){
    data[, DV] <- sample(data[, DV])
  }
  
  for(i in 1:numSets){
    
    if(multiMod){
      RF.data <- data.frame(data[, DV], data[, predictorSets[[i]]])
    }else{
      RF.data <- data.frame(data[, DV], data[, predictorSets])
    }
    
    missings <- sum(!complete.cases(RF.data))
    if(missings > 0){
      warning("Partially incomplete data in model ", i, "\n Dropping ", missings, " observations")
      
      RF.data <- RF.data[complete.cases(RF.data), ]
    }
    
    names(RF.data)[1] <- DV
    
    mtryPar <- switch(mtryArg,
                      default5 = 5,
                      sqroot = round(sqrt(ncol(RF.data)-1), 0),
                      bagging = ncol(RF.data) -1)

    f <- as.formula(paste(DV, ".", sep = " ~ "))
    cfModel <- cforest(f, data = RF.data, controls = cforest_unbiased(ntree = numtree, mtry = mtryPar))
    
    R_squared <- 1 - sum((RF.data[,DV]-predict(cfModel, OOB = T))^2)/sum((RF.data[,DV]-mean(RF.data[,DV]))^2)
    RFresults[i] <- R_squared
    
    if(include.imp == T & !(is.list(predictorSets) & length(predictorSets) > 1)){
      imp <- varimp(cfModel, conditional = TRUE)
    }
    
    RFmodels[[i]] <- cfModel
    
  }
  
  ## DEPRECATED: only if differences between models should be exported 
  # diffScores <- unmatrix(
  #     as.matrix(
  #       dist(RFresults)
  #       )
  #     )
  
  if(exists("imp")){
    results <- list("RF_results" = RFresults, "varImp" = imp)
  }else if(saveModels == TRUE){
    results <- list("Accuracies" = RFresults, "RFmodels" = RFmodels)
  }else{
    results <- RFresults
  }
  
  return(results)
  
}


#############
# repeated random forest for single models, including variable importance

RFmain_repeat_singleModel <- function(DV, 
                     data = df,
                     predictorSets = NULL,
                     numtree = 1000,
                     mtryArg = "sqroot",
                     permute_DV = FALSE, 
                     include.imp = TRUE, 
                     repeats = 5,
                     seed = 1000){
  
  accVec <- numeric(repeats)
  impVec <- matrix(nrow = repeats, ncol = length(predictorSets))
  
  for(i in 1:repeats){
    
    set.seed(i+1)
    
    RFmodel <- RFmain(DV = DV,
                     data = data,
                     predictorSets = predictorSets,
                     numtree = numtree,
                     mtryArg = mtryArg,
                     permute_DV = permute_DV,
                     include.imp = include.imp
                     )
    
    accVec[i] <- RFmodel$RF_results
    impVec[i,] <- RFmodel$varImp
    
    if(i == 1){
      colnames(impVec) <- names(RFmodel$varImp)
    }
    
  }
  
  results <- list("Accuracy" = mean(accVec), "Importance" = colMeans(impVec))
  
  return(results)
  
}




#############
# repeated random forest for multiple models (without variable importance)

RFmain_repeat_multiModel <- function(DV, 
                                      data = df,
                                      predictorSets = NULL,
                                      numtree = 1000,
                                      mtryArg = "sqroot",
                                      permute_DV = FALSE, 
                                      include.imp = FALSE, 
                                      repeats = 5,
                                      seed = 1000){
  
  accOut <- matrix(nrow = repeats, ncol = length(predictorSets))
  
  for(i in 1:repeats){
    
    set.seed(i+1)
    
    RFmodel <- RFmain(DV = DV,
                      data = data,
                      predictorSets = predictorSets,
                      numtree = numtree,
                      mtryArg = mtryArg,
                      permute_DV = permute_DV,
                      include.imp = include.imp)
    
    accOut[i,] <- RFmodel
    
  }
  
  results <- colMeans(accOut)
  
  return(results)
  
}



#############
# Permutation test for single models

RFperm <- function(DV, 
                   data = df,
                   predictorSets = NULL,
                   numtree = 1000,
                   mtryArg = "sqroot",
                   nPerm = 1000
                   ){
  
  permOut <- data.frame(
    matrix(
      nrow = nPerm,
      ncol = length(predictorSets)+1,
      dimnames = list(NULL, c("Accuracy", predictorSets))
    ))
  
  for(i in 1:nPerm){
    
    cat("\r permutation ", i, " of ", nPerm)
    
    RFmodel <- RFmain(DV = DV,
                      data = data,
                      predictorSets = predictorSets,
                      numtree = numtree,
                      mtryArg = mtryArg,
                      permute_DV = TRUE,
                      include.imp = TRUE)
    
    permOut[i,] <- c(RFmodel$RF_results, RFmodel$varImp)
    
  }

  return(permOut)
  
}


#############
# R squared helper function

R_squared_cust <- function(response, prediction){
  
  R2 <- 1-sum((response-prediction)^2)/sum((response-mean(response))^2)
  return(R2)
}



#############
# Permutation test for model comparisons
# (only two models permitted)
# remember to check for missing values! all models need to contain the same observations in the same order

RFcompare <- function(Diff, responseVar, model1, model2, nPerms = 1000){
  
  predModel1 <- predict(model1, OOB = TRUE)
  predModel2 <- predict(model2, OOB = TRUE)
  
  predictionModels <- cbind(predModel1, predModel2)
  
  if(length(predModel1) != length(predModel2)){
    stop("the two models contain a different number of predictions")
  }
  
  diffOut <- numeric(nPerms)
  
  for(i in 1:nPerms){
    
    randInd <- sample(rep(c(1:2), length.out = nrow(predictionModels)))
    predShuffled1 <- predictionModels[cbind(seq_len(nrow(predictionModels)), randInd)]
    
    randInd_inverse <- ifelse(randInd == 1, 2, 1)
    predShuffled2 <- predictionModels[cbind(seq_len(nrow(predictionModels)), randInd_inverse)]
    
    diffOut[i] <- R_squared_cust(responseVar, predShuffled1) - R_squared_cust(responseVar, predShuffled2)
    
  }
  
  p_value <- sum(diffOut >= abs(Diff))/length(diffOut)
  return(p_value)
  
}


#############
# return p-avlues

returnP <- function(values, permutations){
  
  pOut <- numeric(length(values))
  
  for(i in 1:length(values)){
    
    pOut[i] <- sum(permutations[,i] >= values[i])/nrow(permutations)
    
  }
  return(pOut)
  
}


