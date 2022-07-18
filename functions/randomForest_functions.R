##########################
# Function to conduct random forest regression
ACErandomForest <- function(DV,
                            data = df,
                            first.age = "KERF_Sum_3", 
                            last.age = "KERF_Sum_17", 
                            covariates = c("KERF_Sum", "KERF_Multi", "KERF_Duration", "Group", "Age", "Sex"), 
                            numtree = 1000,
                            include.imp = T,
                            plot.imp = T){
  
  # Make smaller df with relevant variables
  data.DVcov <- data[, c(DV, covariates)]
  data.timing <- data[, which(names(data) == first.age):which(names(data) == last.age)]
  RF.data <- data.frame(data.timing, data.DVcov)
  
  # Do random forest regression
  f <- as.formula(paste(DV, ".", sep = " ~ "))
  cfModel <- cforest(f, data = RF.data, controls = cforest_unbiased(ntree = numtree))
  
  # Performance measures
  R_squared <- 1 - sum((data[,DV]-predict(cfModel, OOB = T))^2)/sum((data[,DV]-mean(data[,DV]))^2)
  rmse <- sqrt(mse(predict(cfModel, OOB = T), data[,DV]))
  rmse_baseline <- rmse(data[,DV], mean(data[,DV]))
  if(include.imp == T){imp <- varimp(cfModel, conditional = TRUE)}else{
    plot.imp <- FALSE
    imp <- NA
  }
  
  # Output model and performance in list
  output <- list(cfModel, R_squared, rmse, rmse_baseline, imp)
  names(output) <- c("Model", "Variance_explained", "RMSE", "RMSE_baseline", "variableImportance")
  
  # Plot variable importance
  if(plot.imp == T){
    limits <- c(5*min(imp), 2*max(imp))
    barplot(imp, xlab = "Features", ylab = "Importance", ylim = limits, las = 2)
    abline(h = abs(min(imp)), lty = 2)
  }
  
  output
}




##########################
# Function for permutation test of variance explained and variable importance
RF.permutation <- function(DV,
                           data = df,
                           first.age = "KERF_Sum_3", 
                           last.age = "KERF_Sum_17", 
                           covariates = c("KERF_Sum", "KERF_Multi", "KERF_Duration", "Group", "Age", "Sex"), 
                           numtree = 1000,
                           numperm = 1000, 
                           directory = "results"){
  
  Results.df <- data.frame()
  
  #Start loop
  for(i in 1:numperm){
    
    # Progress report
    cat("Permutation: ", "(", i, "/", numperm, ")\n")
    
    # resample DV
    data.resampled <- data
    data.resampled[, DV] <- sample(data[ ,DV])
    
    # run random forest
    RF.model.perm <- ACErandomForest(DV, 
                                     data = data.resampled, 
                                     first.age = first.age,
                                     last.age = last.age,
                                     covariates = covariates,
                                     numtree = numtree,
                                     plot.imp = F)
    
    # extract R_squared and variable importance
    results.row <- c(RF.model.perm$Variance_explained, RF.model.perm$variableImportance)
    names(results.row)[1] <- "R_squared"
    
    # save to dataframe
    if(i == 1){
      Results.df <- results.row
    }else{
      Results.df <- rbind(Results.df, results.row)
    }
    
  }
  
  # save dataframe to results folder
  Results.df <- as.data.frame(Results.df)
  row.names(Results.df) <- NULL
  
  #if(dir.exists(directory) == F){
  #  directory = dirname(rstudioapi::getActiveDocumentContext()$path)
  #  cat("\n Ouput folder cannot be located. Saving results to script location")
  #}
  
  write.csv(Results.df, file = here("results", paste0("RFperm_", DV, ".csv")), row.names = F)
}




##########################
# Function for p-values from empirical results against permutation-based distributions


# helper function to get the frequency of values larger (or smaller for RMSE) than an empirical value
pValuePerm <- function(xEmp, xPerm, larger = T){
  
  xEmp <- as.numeric(xEmp)
  
  if(larger == T){
    p <- sum(xPerm >= xEmp)/length(xPerm)
  }else{
    p <- sum(xPerm <= xEmp)/length(xPerm)
  }
  p
}

# Main function
Rf.permResults <- function(ROI, RF.model, directory = paste0(wd,"/results")){
  
  # read csv
  df.perm <- read.csv(paste0(directory, "/RFperm_", ROI, ".csv"))
  
  # p-value for variance explained
  R_squared_P <- pValuePerm(RF.model$Variance_explained, df.perm[, 1])
  
  # Loop which calculates p-values for variable importances
  varImp <- RF.model$variableImportance
  varImp_P <- numeric(length = length(varImp))
  for(i in 1:length(varImp)){
    varImp_P[i] <- pValuePerm(varImp[i], df.perm[, i+1])
  }
  names(varImp_P) <- names(varImp)
  
  output <- list(R_squared_P, varImp_P)
  names(output) <- c("R_squared_p", "variableImportance_p")
  
  output
}