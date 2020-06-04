###### PPCA estimation and imputation of missing values
# Name: Results_imputation_fixedeffect
# Description: #it imputes the missing values of an incomplete matrix, assuming that the data matrix contains missing variables and has been generated under a PPCA model with a known number of latent variable r and a known noise level sigma.
               #It uses an estimation of the covariance matrix and the mean.        
               #It returns a list containing the imputation error. 
# Arguments:
# Covar_estimated: estimated covariance matrix of size p*p.
# Mean_estimated: estimated mean vector of size p.
# YNA: data matrix containing missing values of size n*p.
# Y: complete data matrix of size n*p.
# indMissVar: indexes of the missing variables. 
# M: missing-data indicator, binary matrix of size n*p. 
# r: number of latent variables. 
# sigma: noise level (standart deviation).
###### 

Results_imputation_fixedeffect <- function(Covar_estimated,Mean_estimated,YNA,Y,indMissVar,M,r,sigma){ 
  n <- nrow (YNA)
  p <- ncol(YNA)
  
  ressvd <- svd(Covar_estimated - sigma ^ 2 * diag(1, ncol = p, nrow = p))
  BMNAR <- sqrt(diag(ressvd$d)[1:r, 1:p]) %*% t(ressvd$u)
  
  Cov <- t(BMNAR)%*%BMNAR + sigma ^ 2 * diag(1, ncol = p, nrow = p)
  for (j in 1:length(indMissVar)){
    enscol <- 1:p
    for (jbis in indMissVar){
      enscol <- enscol[enscol!=jbis]
    }
    remp1 <- c()
    remp3bis <- c()
    for (jter in enscol){
      remp1 <- c(remp1,Cov[jter,indMissVar[j]])
      remp3bis <- c(remp3bis,mean(YNA[,jter]))
    }
    remp2 <- solve(Cov[enscol, enscol]) 
    remp3 <- t(YNA[, enscol] - matrix(remp3bis, nrow = n, ncol = (p-length(indMissVar)), byrow = TRUE))
    remp <- Mean_estimated[indMissVar[j]] + remp1 %*% remp2 %*% remp3
    missing <- which(is.na(YNA[,indMissVar[j]]))
    YNA[missing, indMissVar[j]] <- remp[missing]
  }
  mse <- MSE(YNA*(1-M),Y*(1-M)) 
  return(list(MSE=mse))
}


###### PPCA estimation
# Name: Results_fixedeffect
# Description: #it gives the imputation error , assuming that the data matrix contains missing variables and has been generated under a PPCA model with a known number of latent variable r and a known noise level sigma.
               #It uses an imputed matrix (implying also an estimation of the covariance matrix and an imputed matrix).   
              #It returns a list containing the imputation error. 
# Arguments:
# Covar_estimated: estimated covariance matrix of size p*p.
# Yimp: imputed data matrix of size n*p.
# Y: complete data matrix of size n*p.
# M: missing-data indicator, binary matrix of size n*p. 
# r: number of latent variables. 
# sigma: noise level (standart deviation).
###### 

Results_fixedeffect <- function(Covar_estimated,Yimp,Y,M,r,sigma){
  
  ressvd <- svd(Covar_estimated - sigma ^ 2 * diag(1, ncol = p, nrow = p))
  Bestim <- sqrt(diag(ressvd$d)[1:r, 1:p]) %*% t(ressvd$u)
  mse <- MSE(Yimp*(1-M),Y*(1-M))
  
  return(list(MSE=mse))
}





