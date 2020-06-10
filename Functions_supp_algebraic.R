library(MASS)


###### Mean estimation 
# Name: Mean_algebraic
# Description: #it estimates the mean of one missing variable considering MNAR mechanism. 
               #it returns the upgraded vector of the estimated means. 
# Arguments: 
# jmiss: index of the missing variable.
# r: number of latent variables.
# YNA: data matrix containing missing values of size n*p.
# combinations: possible combinations of variables to perform the regressions for the mean estimation.
# combination_choosed: NULL if opt="agg", the choosed combination if opt="noagg". 
# reg_results: results of the regressions. 
# mean_results: vector of the estimated means of size p. 
# opt: "agg" (by default) if we aggregate many possible combinations of variables to perform the regressions, "noagg" if we randomly choose one combination.
###### 

Mean_algebraic <- function(jmiss,r,YNA,combinations,combination_choosed,reg_results,mean_results,opt){
  
  if (opt=="agg"){
  Mean <- c()
  
  for (i in 1:nrow(combinations)){
    Mean <- c(Mean,mean_results[combinations[i,1]] - reg_results[[as.character(jmiss)]][i,1] - sum(mean_results[combinations[i,2:r]]*reg_results[[as.character(jmiss)]][i,3:(r+1)]) / reg_results[[as.character(jmiss)]][i,2])
  }
  
  mean_results[jmiss] <- median(Mean)
  
  }else{
    mean_results[jmiss] <- mean_results[combinations[1]] - reg_results[[as.character(jmiss)]][combination_choosed,1] - sum(mean_results[combinations[2:r]]*reg_results[[as.character(jmiss)]][combination_choosed,3:(r+1)]) / reg_results[[as.character(jmiss)]][combination_choosed,2]
  }
  
  return(mean_results)
}



###### Variance and covariances (with pivot variables) estimation
# Name: Variance_covariance_algebraic
# Description: #it estimates the variance and covariances (with pivot variables) of one missing variable considering MNAR mechanism. 
               #it returns the upgraded list of the estimated covariances.
# Arguments: 
# jmiss: index of the missing variable.
# r: number of latent variables.
# indMissVar: indexes of the missing variables. 
# YNA: data matrix containing missing values of size n*p.
# combinations: possible combinations of variables to perform the regressions for the variance estimation and the covariances between a MNAR variable and a pivot one. 
# combinations_sortrows: rows in the matrix called combinations corresponding to the ordered combinations given one nonordered combination. The first line gives the rows in the matrix called combinations correponding to the possible ordered combinations from the nonordered combination in the first line of combinations_uniq. 
# combinations_uniq: possible nonordered combinations of variables (nonordered: for instance 1 2 3 and 2 1 3 are the same nonordered combination).
# reg_results: results of the regressions. 
# mean_results: vector of the estimated means of size p. 
# cov_results: matrix of the estimated covariances of size p*p (the estimated covariances for the different combinations of variables have been aggregated).
# cov_stock: list of the estimated covariances (for the different combinations of variables).
# opt: "agg" (by default) if we aggregate many possible combinations of variables to perform the regressions, "noagg" if we randomly choose one combination.
###### 

Variance_covariance_algebraic <- function(jmiss,r,indMissVar,YNA,combinations,combinations_sortrows,combinations_uniq,reg_results,mean_results,cov_results,cov_stock,opt){

  
  if (opt == "agg"){
    
    for (combi in 1:dim(combinations_uniq)[1]){
      combi_ord <- combinations_sortrows[combi,]
      for (k in combi_ord){
        YNAcc_allvar <- YNA
        for (u in c(jmiss,combinations[k,])){
          indvalNA <- which(is.na(YNAcc_allvar[,u]))
          if (length(indvalNA)>0){
          YNAcc_allvar <- YNAcc_allvar[-indvalNA,]
          }
        } 
        c_sortrows_minus_chosen <- c(k,setdiff(combinations_sortrows[combi,],k))
        M1 <- matrix(0,nrow=r+1,ncol=r+1)
        M2 <- numeric(r+1)
        M1[1,1] <- reg_results[[as.character(jmiss)]][k,2]^2 
        M1[1,3:(r+1)] <- 2*reg_results[[as.character(jmiss)]][k,3:(r+1)]*reg_results[[as.character(jmiss)]][k,2]
        for (ind in 2:r){
          M1[ind,] <- -c(reg_results[[as.character(jmiss)]][c_sortrows_minus_chosen[ind-1],2:ind],-1,reg_results[[as.character(jmiss)]][c_sortrows_minus_chosen[ind-1],(ind+1):(r+1)])
        }
        M1[r+1,] <- -c(reg_results[[as.character(jmiss)]][c_sortrows_minus_chosen[r],2:(r+1)],-1)
        variables_reg <- combinations[k,]
        cov_cc <- cov(YNAcc_allvar[,variables_reg[1]],YNAcc_allvar[,c(jmiss,variables_reg[2:r])])
        Q_cc <- var(YNAcc_allvar[,variables_reg[1]])-cov_cc%*%solve(var(YNAcc_allvar[,c(jmiss,variables_reg[2:r])]))%*%t(cov_cc)
        if (r==2){
          M2[1] <- cov_results[variables_reg[1],variables_reg[1]] - Q_cc - t(reg_results[[as.character(jmiss)]][k,3:(r+1)])*cov_results[variables_reg[2:r],variables_reg[2:r]]*reg_results[[as.character(jmiss)]][k,3:(r+1)]
          for (ind in 2:(r+1)){
            M2[ind] <- (reg_results[[as.character(jmiss)]][c_sortrows_minus_chosen[ind-1],]%*%c(1,mean_results[jmiss],mean_results[combinations[c_sortrows_minus_chosen[ind-1],][2:r]]) - mean_results[variables_reg[ind-1]]) * mean_results[jmiss] 
          }    
        }else{
          M2[1] <- cov_results[variables_reg[1],variables_reg[1]] - Q_cc - t(reg_results[[as.character(jmiss)]][k,3:(r+1)])%*%cov_results[variables_reg[2:r],variables_reg[2:r]]%*%reg_results[[as.character(jmiss)]][k,3:(r+1)] 
          for (ind in 2:(r+1)){
            M2[ind] <- (reg_results[[as.character(jmiss)]][c_sortrows_minus_chosen[ind-1],]%*%c(1,mean_results[jmiss],mean_results[combinations[c_sortrows_minus_chosen[ind-1],][2:r]]) - mean_results[variables_reg[ind-1]]) * mean_results[jmiss] 
          }
        }

        Res <- ginv(M1,tol = 10^-20)%*%M2

        cov_stock[[paste(as.character(jmiss),as.character(jmiss))]] <- c(cov_stock[[paste(as.character(jmiss),as.character(jmiss))]],Res[1])
        for (ind in 1:r){
          cov_stock[[paste(as.character(jmiss),as.character(combinations[k,][ind]))]]<-c(cov_stock[[paste(as.character(jmiss),as.character(combinations[k,][ind]))]],Res[ind+1])
        }
      }
    }
  
  }else{
    YNAcc_allvar <- YNA
    for (u in c(jmiss,combinations)){
      indvalNA <- which(is.na(YNAcc_allvar[,u]))
      if (length(indvalNA)>0){
        YNAcc_allvar <- YNAcc_allvar[-indvalNA,]
      }
    }
    M1 <- matrix(0,nrow=r+1,ncol=r+1)
    M2 <- numeric(r+1)
    combi_ord <- combinations_sortrows
    M1[1,1] <- reg_results[[as.character(jmiss)]][combi_ord[1],2]^2 
    M1[1,3:(r+1)] <- 2*reg_results[[as.character(jmiss)]][combi_ord[1],3:(r+1)]*reg_results[[as.character(jmiss)]][combi_ord[1],2]
    for (ind in 2:r){
      M1[ind,] <- -c(reg_results[[as.character(jmiss)]][combi_ord[ind-1],2:ind],-1,reg_results[[as.character(jmiss)]][combi_ord[ind-1],(ind+1):(r+1)])
    }
    M1[r+1,] <- -c(reg_results[[as.character(jmiss)]][combi_ord[r],2:(r+1)],-1)
    variables_reg <- combinations
    cov_cc <- cov(YNAcc_allvar[,variables_reg[1]],YNAcc_allvar[,c(jmiss,variables_reg[2:r])])
    Q_cc <- var(YNAcc_allvar[,variables_reg[1]])-cov_cc%*%solve(var(YNAcc_allvar[,c(jmiss,variables_reg[2:r])]))%*%t(cov_cc)
    if (r==2){
      M2[1] <- cov_results[variables_reg[1],variables_reg[1]] - Q_cc - t(reg_results[[as.character(jmiss)]][combi_ord[1],3:(r+1)])*cov_results[variables_reg[2:r],variables_reg[2:r]]*reg_results[[as.character(jmiss)]][combi_ord[1],3:(r+1)] #+ sum(sigma^2*reg_results[[as.character(jmiss)]][combi_ord[1],2:(r+1)]^2)
      for (ind in 2:(r+1)){
        M2[ind] <- (reg_results[[as.character(jmiss)]][combi_ord[ind-1],]%*%c(1,mean_results[jmiss],mean_results[variables_reg[2:r]]) - mean_results[variables_reg[ind-1]]) * mean_results[jmiss] #- sigma^2*reg_results[[as.character(jmiss)]][combi_ord[ind-1],2]
      }    
    }else{
      M2[1] <- cov_results[variables_reg[1],variables_reg[1]] - Q_cc - t(reg_results[[as.character(jmiss)]][combi_ord[1],3:(r+1)])%*%cov_results[variables_reg[2:r],variables_reg[2:r]]%*%reg_results[[as.character(jmiss)]][combi_ord[1],3:(r+1)] #+ sum(sigma^2*reg_results[[as.character(jmiss)]][combi_ord[1],2:(r+1)]^2)
      for (ind in 2:(r+1)){
        M2[ind] <- (reg_results[[as.character(jmiss)]][combi_ord[ind-1],]%*%c(1,mean_results[jmiss],mean_results[variables_reg[2:r]]) - mean_results[variables_reg[ind-1]]) * mean_results[jmiss] #- sigma^2*reg_results[[as.character(jmiss)]][combi_ord[ind-1],2]
      }
    }

    Res <- ginv(M1,tol = 10^-20)%*%M2
    
    cov_stock[[paste(as.character(jmiss),as.character(jmiss))]] <- Res[1]
    for (ind in 1:r){
      cov_stock[[paste(as.character(jmiss),as.character(combinations[combi_ord][ind]))]]<-Res[ind+1]
    }
  }
  
  return(cov_stock)
  
}


###### Covariances (with not pivot variables) estimation
# Name: Covariance_2missing_algebraic 
# Description: #it estimates the covariances (with not pivot variables) of one missing variable considering MNAR mechanism. 
               #it returns the upgraded list of the estimated covariances.
# Arguments: 
# jmiss: index of the missing variable.
# jmiss2: index of the second missing variable (for which the covariance is computed).
# r: number of latent variables.
# indMissVar: indexes of the missing variables. 
# combinations_cov: possible combinations of variables to perform the regressions for the covariance estimation between two not pivot variables. 
# cov_results: matrix of the estimated covariances of size p*p (the estimated covariances for the different combinations of variables have been aggregated).
# YNA: data matrix containing missing values of size n*p.
# opt: "agg" (by default) if we aggregate many possible combinations of variables to perform the regressions, "noagg" if we randomly choose one combination.
###### 

Covariance_2missing_algebraic <- function(jmiss,jmiss2,r,indMissVar,combinations_cov,cov_results,YNA,opt){

  
  if (opt=="agg"){
    Cov <- c()
    for (cb in 1:dim(combinations_cov)[1]){
      combi <- combinations_cov[cb,]
      YNAcc_allvar <- YNA
      for (k in c(jmiss,jmiss2,combi)){
        indvalNA <- which(is.na(YNAcc_allvar[,k]))
        if(length(indvalNA)>0){
          YNAcc_allvar <- YNAcc_allvar[-indvalNA,] 
        }
      }
      if (r>2){
        variables_reg_minusy <- c(jmiss,jmiss2,combi[2:(r-1)])
        variables_reg <- combi[2:(r-1)]
        coeff_reg <- coefficients(lm(YNAcc_allvar[,combi[1]] ~ as.matrix(YNAcc_allvar[,variables_reg_minusy])))
        variances_vec <- diag(cov_results[variables_reg_minusy,variables_reg_minusy])
        term <- 0
        for (k in 1:length(variables_reg_minusy)){
          for (l in 1:length(combi[2:(r-1)])){
            if (k!=l+2){
              term <- term + 2*coeff_reg[k+1]*coeff_reg[l+3]*cov_results[variables_reg_minusy[k],combi[2:(r-1)][l]]
            }          
          }
        }
      }else{
        variables_reg_minusy <- c(jmiss,jmiss2)  
        variables_reg <- c()
        coeff_reg <- coefficients(lm(YNAcc_allvar[,combi[1]] ~ as.matrix(YNAcc_allvar[,variables_reg_minusy])))
        variances_vec <- diag(cov_results[variables_reg_minusy,variables_reg_minusy])
        term <- 0
      }
      cov_cc <- cov(YNAcc_allvar[,combi[1]],YNAcc_allvar[,c(jmiss,jmiss2,variables_reg)])
      Q_cc <- var(YNAcc_allvar[,combi[1]])-cov_cc%*%solve(var(YNAcc_allvar[,c(jmiss,jmiss2,variables_reg)]))%*%t(cov_cc)
      K <- 2*coeff_reg[2]*coeff_reg[3]
      bruit <- 0
      Cov <- c(Cov,(cov_results[combi[1],combi[1]]-Q_cc-t(coeff_reg[2:(r+1)]^2)%*%variances_vec-term-bruit)/K)
    }
    
    cov_results[jmiss,jmiss2] <- median(Cov)
    cov_results[jmiss2,jmiss] <- cov_results[jmiss,jmiss2]
  }else{
    combi <- combinations_cov
    YNAcc_allvar <- YNA
    for (k in c(jmiss,jmiss2,combi)){
      indvalNA <- which(is.na(YNAcc_allvar[,k]))
      if(length(indvalNA)>0){
        YNAcc_allvar <- YNAcc_allvar[-indvalNA,] 
      }
    }
    if (r>2){
      variables_reg_minusy <- c(jmiss,jmiss2,combi[2:(r-1)])
      variables_reg <- combi[2:(r-1)]
      coeff_reg <- coefficients(lm(YNAcc_allvar[,combi[1]] ~ as.matrix(YNAcc_allvar[,variables_reg_minusy])))
      variances_vec <- diag(cov_results[variables_reg_minusy,variables_reg_minusy])
      term <- 0
      for (k in 1:length(variables_reg_minusy)){
        for (l in 1:length(combi[2:(r-1)])){
          if (k!=l+2){
            term <- term + 2*coeff_reg[k+1]*coeff_reg[l+3]*cov_results[variables_reg_minusy[k],combi[2:(r-1)][l]]
          }  
        }
      }
    }else{
      variables_reg_minusy <- c(jmiss,jmiss2)  
      variables_reg <- c()
      coeff_reg <- coefficients(lm(YNAcc_allvar[,combi[1]] ~ as.matrix(YNAcc_allvar[,variables_reg_minusy])))
      variances_vec <- diag(cov_results[variables_reg_minusy,variables_reg_minusy])
      term <- 0
    }
    cov_cc <- cov(YNAcc_allvar[,combi[1]],YNAcc_allvar[,c(jmiss,jmiss2,variables_reg)])
    Q_cc <- var(YNAcc_allvar[,combi[1]])-cov_cc%*%solve(var(YNAcc_allvar[,c(jmiss,jmiss2,variables_reg)]))%*%t(cov_cc)
    K <- 2*coeff_reg[2]*coeff_reg[3]
    bruit <- 0 
    Cov <- (cov_results[combi[1],combi[1]]-Q_cc-t(coeff_reg[2:(r+1)]^2)%*%variances_vec-term-bruit)/K

    cov_results[jmiss,jmiss2] <- Cov
    cov_results[jmiss2,jmiss] <- cov_results[jmiss,jmiss2]
  }
  return(cov_results)
}
