###### Mean estimation 
# Name: Mean_graphical
# Description: #it estimates the mean of one missing variable considering MNAR mechanism. 
               #it returns the upgraded vector of the estimated means. 
# Arguments: 
# jmiss: index of the missing variable.
# YNA: data matrix containing missing values of size n*p.
# combinations: possible combinations of variables to perform the regressions for the mean estimation.
# combi_choosed: NULL if opt="agg", the choosed combination if opt="noagg". 
# reg_results: results of the regressions. 
# mean_results: vector of the estimated means of size p. 
# opt: "agg" (by default) if we aggregate many possible combinations of variables to perform the regressions, "noagg" if we randomly choose one combination.
###### 

Mean_graphical <- function(jmiss,YNA,combinations,combi_choosed,reg_results,mean_results,opt){
  
  if (opt=="agg"){
  Mean <- c()
  for (i in 1:nrow(combinations)){
    Mean <- c(Mean,mean(YNA[, combinations[i,1]]) - reg_results[[as.character(jmiss)]][i,1] - mean(YNA[,combinations[i,2]])*reg_results[[as.character(jmiss)]][i,3] / reg_results[[as.character(jmiss)]][i,2])
  }
  mean_results[jmiss] <- median(Mean)
  }else{
    mean_results[jmiss] = mean(YNA[, combinations[1]]) - reg_results[[as.character(jmiss)]][combi_choosed,1] - mean(YNA[,combinations[2]])*reg_results[[as.character(jmiss)]][combi_choosed,3] / reg_results[[as.character(jmiss)]][combi_choosed,2]
  }
  
  return(mean_results)
}


###### Variance estimation
# Name: Variance_graphical
# Description: #it estimates the variance f one missing variable considering MNAR mechanism. 
               #it returns the upgraded covariance matrix.
# Arguments: 
# jmiss: index of the missing variable.
# YNA: data matrix containing missing values of size n*p.
# combinations: possible combinations of variables to perform the regressions for the variance estimation.
# reg_results: results of the regressions. 
# reg_results_variance: results of the specific regressions for the variance estimatoin. 
# cov_results: matrix of the estimated covariances of size p*p (the estimated covariances for the different combinations of variables have been aggregated).
###### 

Variance_graphical <- function(jmiss,YNA,combinations,reg_results,reg_results_variance,cov_results){
  
  Var <- c()
  for (i in 1:nrow(combinations)){
    comb <- combinations[i,]
    row_reg_results_variance <- which(reg_results_variance[[as.character(jmiss)]][,3]==comb[2])
    coeff_div <- reg_results_variance[[as.character(jmiss)]][row_reg_results_variance,2]
    coeff <- (1 / reg_results[[as.character(jmiss)]][i,2]) * ((cov(YNA[, comb[1]], YNA[, comb[2]]) / var(YNA[, comb[2]])) -  reg_results[[as.character(jmiss)]][i,3])
    Var <- c(Var,(coeff * var(YNA[, comb[2]])) / coeff_div)
  }
  
  cov_results[jmiss,jmiss] <- median(Var)
  
  return(cov_results) 
}


###### Covariances (with pivot variables) estimation
# Name: Covariance_graphical
# Description: #it estimates the variance and covariances (with pivot variables) of one missing variable considering MNAR mechanism. 
              #it returns the upgraded covariance matrix.
# Arguments: 
# jmiss: index of the missing variable.
# YNA: data matrix containing missing values of size n*p.
# indRegVar: indexes of the pivot variables.
# combinations: possible combinations of variables to perform the regressions for the covariances estimation.
# reg_results: results of the regressions. 
# cov_results: matrix of the estimated covariances of size p*p (the estimated covariances for the different combinations of variables have been aggregated).
###### 

Covariance_graphical <- function(jmiss,YNA,indRegVar,combinations,reg_results,cov_results){
  
  for (k in indRegVar){
    reg_results[[as.character(jmiss)]]
    Cov <- c()
    row_reg_results <- which(combinations[,2]==k)
    for (l in row_reg_results){
      b <- reg_results[[as.character(jmiss)]][l,2]
      c <- reg_results[[as.character(jmiss)]][l,3]
      coeff <- (1 / reg_results[[as.character(jmiss)]][l,2]) * ((cov(YNA[, combinations[l,1]], YNA[, k]) / var(YNA[, k])) - reg_results[[as.character(jmiss)]][l,3])
      Cov <- c(Cov,coeff*var(YNA[,k]))
    }
    cov_results[jmiss,k] <- median(Cov)
    cov_results[k,jmiss] <- cov_results[jmiss,k]
  }
  
  return(cov_results) 
  
}


###### Covariances (with not pivot variables) estimation
# Name: Covariance_2missing_graphical
# Description: #it estimates the covariances (with not pivot variables) of one missing variable considering MNAR mechanism. 
               #it returns the upgraded covariance matrix.
# Arguments: 
# jmiss: index of the missing variable.
# jmiss2: 
# YNA: data matrix containing missing values of size n*p.
# indRegVar: indexes of the pivot variables.
# reg_results: results of the regressions. 
# cov_results: matrix of the estimated covariances of size p*p (the estimated covariances for the different combinations of variables have been aggregated).
###### 

Covariance_2missing_graphical <- function(jmiss,jmiss2,YNA,indRegVar,reg_results,cov_results){
  
  Cov <- c()
  variables_reg_minusy <- permutations(2,2,c(jmiss,jmiss2))
  for (k in indRegVar){
    YNAcc_allvar <- YNA
    for (ind in c(jmiss,jmiss2,k)){
      indvalNA <- which(is.na(YNAcc_allvar[,ind]))
      if(length(indvalNA)>0){
        YNAcc_allvar <- YNAcc_allvar[-indvalNA,] 
      }
    }
    for (l in 1:2){
      coeff_reg <- coefficients(lm(YNAcc_allvar[,k] ~ YNAcc_allvar[,variables_reg_minusy[l,]]))
      coeff <- (1 / coeff_reg[2]) * (cov_results[k,variables_reg_minusy[l,2]] / cov_results[variables_reg_minusy[l,2],variables_reg_minusy[l,2]] - coeff_reg[3])
      Cov <- c(Cov,coeff*cov_results[variables_reg_minusy[l,2],variables_reg_minusy[l,2]])
    }
    cov_results[jmiss,jmiss2] <- median(Cov)
    cov_results[jmiss2,jmiss] <- cov_results[jmiss,jmiss2]
  }

  return(cov_results) 
  
}
