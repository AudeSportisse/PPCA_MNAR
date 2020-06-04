###### Mean and covariance matrix estimation (our method - Algorithm 1)
# Name: Mean_covariances_estimations_algebraic
# Description: #it estimates the mean and the covariance matrix of a data matrix containing MNAR variables using our method. 
               #it returns a list containing the estimated mean and covariance matrix.
# Arguments: 
# YNA: data matrix containing missing values of size n*p.
# indMissVar: indexes of the missing variables. 
# indRegVar: indexes of the pivot variables.
# r: number of latent variables.
# YNAimp: complete data matrix of size n*p imputed in advance if opt_data="real".
#         NULL if opt_data="simu" or "MCAR".
# opt: "agg" (by default) if we aggregate many possible combinations of variables to perform the regressions, "noagg" if we randomly choose one combination.
# opt_data: "simu" (by default), "real" or "MCAR": 
#           "simu": for the pivot variables (they are assumed to be observed), the estimated mean and covariances are computed using the empirical quantities. 
#           "real": for the pivot variables, the estimated mean and covariances are computed using the empirical quantities for Yimp.
#           "MCAR": for the pivot variables (they are assumed to be MCAR), the estimated mean and covariances are computed using the empirical quantities for the complete case. 
###### 


Mean_covariances_estimations_algebraic <- function(YNA,indMissVar,indRegVar,r,YNAimp=NULL,opt="agg",opt_data="simu"){
  
  indMissVar <- sort(indMissVar)
  indRegVar <- sort(indRegVar)
  
  if (length(indRegVar)<r){
    
    stop("Error: the number of possible regression variables must be greater than the rank.")
    
  }else{
    
    if (opt=="agg"){
      
    ##Possible regressions
    nrow_pervar <- choose(length(indRegVar)-1,r-1)
    combinations <- matrix(1000,nrow=nrow_pervar*length(indRegVar),ncol=r)
    for (i in 1:length(indRegVar)){
      combinations[((i-1)*nrow_pervar+1):((i-1)*nrow_pervar+nrow_pervar),1] <- indRegVar[i]
      if (length(indRegVar)>r){
        combinations[((i-1)*nrow_pervar+1):((i-1)*nrow_pervar+nrow_pervar),2:r] <- t(combn(setdiff(indRegVar,indRegVar[i]),r-1))
      }else{
        combinations[((i-1)*nrow_pervar+1):((i-1)*nrow_pervar+nrow_pervar),2:r] <- setdiff(indRegVar,indRegVar[i])
      }
    }
    combinations_mean <- combinations
    combinations_ord <- t(apply(combinations,1,sort))
    
    ## Table of the results for the regressions. 
    reg_results <- NULL
    
      for (k in indMissVar){
        reg_results[[as.character(k)]] <- matrix(1000,nrow=nrow(combinations),ncol=r+1)
        for (l in 1:nrow(combinations)){
          YNAcc <- YNA
          for (u in c(k,combinations[l,])){
            indvalNA <- which(is.na(YNAcc[,u]))
            if (length(indvalNA)>0){
              YNAcc <- YNAcc[-indvalNA,]
            }
          }
          YNAcc=as.matrix(YNAcc)
          reg_results[[as.character(k)]][l,] <-coefficients(lm(YNAcc[, combinations[l,1]] ~  YNAcc[, c(k,combinations[l,2:r])])) 
        }
      }
      combination_choosed <- NULL
      
      ## For the variance and covariances with pivot variables
      combinations_uniq <- unique(combinations_ord)
      combinations_sortrows <- matrix(1000,nrow=dim(combinations_uniq)[1],ncol=r)
      for (i in 1:dim(combinations_uniq)[1]){
        combinations_sortrows[i,] <- which(apply(combinations_ord,1,identical,y=combinations_uniq[i,]))
      }
      combinations_variance <- combinations
      
      ## For the covariances between 2 MNAR variables
      if (r>2){
        nrow_pervar_cov <- choose(length(indRegVar)-1,r-2)
        combinations_cov <- matrix(1000,nrow=nrow_pervar_cov*length(indRegVar),ncol=r-1)
        for (i in 1:length(indRegVar)){
          combinations_cov[((i-1)*nrow_pervar_cov+1):((i-1)*nrow_pervar_cov+nrow_pervar_cov),1] <- indRegVar[i]
          combinations_cov[((i-1)*nrow_pervar_cov+1):((i-1)*nrow_pervar_cov+nrow_pervar_cov),2:(r-1)] <- t(combn(setdiff(indRegVar,indRegVar[i]),r-2))
        }
      }else{
        combinations_cov <- matrix(indRegVar,nrow=length(indRegVar),ncol=1)
      }
    }else{
      
      ##Choose one regression
      
      var_used <- sample(indRegVar,r,replace=FALSE)
      combinations <- matrix(1000,nrow=length(var_used),ncol=r)
      for (i in 1:length(var_used)){
        combinations[i,1] <- var_used[i]
        combinations[i,2:r] <- setdiff(var_used,var_used[i])
      }
      
      ##For mean:
      combination_choosed <- sample(nrow(combinations),1)
      combinations_mean <- combinations[combination_choosed,]
      combinations_variance <- combinations[combination_choosed,]
      
      ## Table of the results for the regressions. 
      reg_results <- NULL
      
      for (k in indMissVar){
        indvalNA <- which(is.na(YNA[,k]))
        YNAcc <- YNA[-indvalNA,]
        reg_results[[as.character(k)]] <- matrix(1000,nrow=nrow(combinations),ncol=r+1)
        for (l in 1:nrow(combinations)){
          reg_results[[as.character(k)]][l,] <-coefficients(lm(YNAcc[, combinations[l,1]] ~  YNAcc[, c(k,combinations[l,2:r])])) 
        }
      }
      combinations_sortrows <- c(combination_choosed,setdiff(1:nrow(combinations),combination_choosed))
      combinations_uniq <- NULL
      
      combinations_cov <- sample(combinations_mean,r-1,replace=FALSE)
    }
  }

  ## MCAR and for real data
  if (opt_data == "real"){
    mean_results <- apply(YNAimp,2,mean,na.rm = TRUE)
    mean_results[indMissVar] <- 1000
  }else if(opt_data == "MCAR"){
    mean_results <- apply(YNA,2,mean,na.rm = TRUE)
  }else{
    mean_results <- apply(YNA,2,mean)
    mean_results[is.na(mean_results)] <- 1000
  }
  
  ## for real data
  if (opt_data == "real"){
    cov_results <- var(YNAimp)
    for (l in indMissVar){
      cov_results[l,] <- 1000
      cov_results[,l] <- 1000
    }
  }else if(opt_data == "MCAR"){
    cov_results <- matrix(1000,ncol=ncol(YNA),nrow=ncol(YNA))
    for (l in 1:ncol(YNA)){
      cov_results[l,] <- apply(YNA,2,cov,YNA[,l],use="complete.obs")
    }
  }else{
    cov_results <- var(YNA)
    cov_results[is.na(cov_results)] <- 1000
  }
  cov_stock <- NULL
  
  for (jmiss in indMissVar){
    
    cov_stock[[paste(as.character(jmiss),as.character(jmiss))]]<-numeric()
    for (j in indRegVar){
      cov_stock[[paste(as.character(jmiss),as.character(j))]]<-numeric()
    }
    
    mean_results <- Mean_algebraic(jmiss,r,YNA,combinations_mean,combination_choosed,reg_results,mean_results,opt)
    
    cov_stock <- Variance_covariance_algebraic(jmiss,r,indMissVar,YNA,combinations_variance,combinations_sortrows,combinations_uniq,reg_results,mean_results,cov_results,cov_stock,opt)
    
    cov_results[jmiss,jmiss] <- median(cov_stock[[paste(as.character(jmiss),as.character(jmiss))]])
    for (j in indRegVar){
      cov_results[jmiss,j] <- median(cov_stock[[paste(as.character(jmiss),as.character(j))]])
      cov_results[j,jmiss] <- cov_results[jmiss,j]
    }
    
    if (opt!="agg"){
      for (j in setdiff(indRegVar,combinations_mean)){
        cov_results <- Covariance_2missing_algebraic(jmiss,j,r,indMissVar,combinations_cov,cov_results,YNA,opt)
      }
    }else if( (length(indRegVar)+length(indMissVar))!=ncol(YNA)){
      for (j in setdiff(1:ncol(YNA),c(indMissVar,indRegVar))){
        cov_results <- Covariance_2missing_algebraic(jmiss,j,r,indMissVar,combinations_cov,cov_results,YNA,opt)
      }
    }
    
    for (jmiss2 in indMissVar){
      if(jmiss2<jmiss){
        cov_results <- Covariance_2missing_algebraic(jmiss,jmiss2,r,indMissVar,combinations_cov,cov_results,YNA,opt)
      }
    }
  }
  
  
  return(list(mean=mean_results,cov=cov_results))
  
  
}


