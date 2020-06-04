###### Mean and covariance matrix estimation (our graphical method - in Appendix)
# Name: Mean_covariances_estimations_graphical
# Description: #it estimates the mean and the covariance matrix of a data matrix containing MNAR variables using our graphical method. 
               #it returns a list containing the estimated mean and covariance matrix.
# Arguments: 
# YNA: data matrix containing missing values of size n*p.
# indMissVar: indexes of the missing variables. 
# indRegVar: indexes of the pivot variables.
# opt: "agg" (by default) if we aggregate many possible combinations of variables to perform the regressions, "noagg" if we randomly choose one combination.
###### 

Mean_covariances_estimations_graphical <- function(YNA,indMissVar,indRegVar,opt){
  
  indMissVar <- sort(indMissVar)
  indRegVar <- sort(indRegVar)
  
  if (length(indRegVar)<2){
    
    stop("Error: the number of possible regression variables must be greater than 2.")
    
  }else{
    combinations <- permutations(length(indRegVar),2,indRegVar)
    if (opt!="agg"){
      combinations <- matrix(1000,nrow=length(indRegVar),ncol=2)
      for (k in 1:length(indRegVar)){
        combinations[k,2] <- indRegVar[k]
        combinations[k,1] <- sample(setdiff(indRegVar,indRegVar[k]),1)
      }
    }
    ## Table of the results for the regressions. 
    reg_results <- NULL
    reg_results_variance <- NULL
    for (k in indMissVar){
      indvalNA <- which(is.na(YNA[,k]))
      YNAcc <- YNA[-indvalNA,]
      reg_results[[as.character(k)]] <- matrix(1000,nrow=nrow(combinations),ncol=3)
      if (opt=="agg"){
        indRegVar_forvariance <- indRegVar
      }else{
        indRegVar_forvariance <- combinations[1,2]
      }
      reg_results_variance[[as.character(k)]] <- cbind(matrix(1000,nrow=length(indRegVar_forvariance),ncol=2),indRegVar_forvariance)
      for (l in 1:nrow(combinations)){
        reg_results[[as.character(k)]][l,] <- coefficients(lm(YNAcc[, combinations[l,1]] ~  YNAcc[, c(k,combinations[l,2])])) 
      }
      if (opt=="agg"){
        combi_choosed <- NULL
        combinations_mean <- combinations 
      }else{
        combi_choosed <- sample(nrow(combinations),1)
        combinations_mean <- combinations[combi_choosed,]
      }
      for (l in 1:length(indRegVar_forvariance)){
        reg_results_variance[[as.character(k)]][l,1:2] <- coefficients(lm(YNAcc[, indRegVar_forvariance[l]] ~  YNAcc[, k])) 
      }
    }
  }
  
  mean_results <- apply(YNA,2,mean)
  mean_results[is.na(mean_results)] <- 0
  cov_results <- var(YNA)
  cov_results[is.na(cov_results)] <- 0
  
  for (jmiss in indMissVar){
    indvalNA <- which(is.na(YNA[,jmiss]))
    YNAcc <- YNA[-indvalNA,]
    
    
    mean_results <- Mean_graphical(jmiss,YNA,combinations_mean,combi_choosed,reg_results,mean_results,opt)
    cov_results <- Variance_graphical(jmiss,YNA,combinations,reg_results,reg_results_variance,cov_results)
    cov_results <- Covariance_graphical(jmiss,YNA,indRegVar,combinations,reg_results,cov_results)
    
    for (jmiss2 in indMissVar){
      if(jmiss2<jmiss){
        cov_results <- Covariance_2missing_graphical(jmiss,jmiss2,YNA,indRegVar,reg_results,cov_results)
      }
    }
  }
  
  
  return(list(mean=mean_results,cov=cov_results))
  
  
}


