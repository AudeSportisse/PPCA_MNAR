###### Measuring the performance
# Name: MSE
# Description: it returns the normalized MSE to measure the quality of the imputation. 
# Arguments: 
# X1: imputed matrix.
# X2: complete matrix (of same size of X1).
###### 

MSE <- function(X1, X2){ return(norm(as.matrix(X1 - X2),type="F")/norm(as.matrix(X2),type="F"))}


###### Imputation by the mean
# Name: ImputeMean0
# Description: #it imputes an incomplete matrix by the mean
               #It returns an imputed matrix by the mean. 
# Arguments:
# tab: incomplete matrix. 
###### 

ImputeMean0 <- function(tab){
  m <- apply(tab, 2, mean, na.rm = TRUE)
  tab <- sapply(1:ncol(tab), function(x) ifelse(is.na(tab[,x]), m[x], tab[,x]))
  tab <- as.data.frame(tab)
  return(tab)
}

###### 
# Name: ImputeMean
# Description: #it imputes an incomplete matrix by the mean but all values in the column are replaced by the average of the column. 
               #It returns the imputed matrix.
# Arguments:
# tab: incomplete matrix. 
###### 

ImputeMean <- function(tab){
  m <- apply(tab, 2, mean, na.rm = TRUE)
  tab <- sapply(1:ncol(tab), function(x) replicate(nrow(tab),m[x]))
  tab <- as.data.frame(tab)
  return(tab)
}


###### PPCA estimation and imputation of missing values
# Name: Results_PPCA_imputation
# Description: #it estimates the loading matrix and imputes the missing values of an incomplete matrix, assuming that the data matrix contains missing variables and has been generated under a PPCA model with a known number of latent variable r and a known noise level sigma.
               #It uses an estimation of the covariance matrix and the mean.        
               #It returns a list containing the correlation coefficient (RV coefficient) between the estimated loading matrix and the true one, and also the imputation error. 
# Arguments:
# Covar_estimated: estimated covariance matrix of size p*p.
# Mean_estimated: estimated mean vector of size p.
# YNA: data matrix containing missing values of size n*p.
# Y: complete data matrix of size n*p.
# B: true loading matrix of size r*p.
# indMissVar: indexes of the missing variables. 
# M: missing-data indicator, binary matrix of size n*p. 
# r: number of latent variables. 
# sigma: noise level (standart deviation).
###### 

Results_PPCA_imputation <- function(Covar_estimated,Mean_estimated,YNA,Y,B,indMissVar,M,r,sigma){ 
  n <- nrow (YNA)
  p <- ncol(YNA)
  
  ressvd <- svd(Covar_estimated - sigma ^ 2 * diag(1, ncol = p, nrow = p))
  BMNAR <- sqrt(diag(ressvd$d)[1:r, 1:p]) %*% t(ressvd$u)
  Corr <- coeffRV(t(BMNAR),t(B))$rv 
  
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
  return(list(Corr=Corr,MSE=mse))
}


###### PPCA estimation
# Name: Results_PPCA
# Description: #it estimates the loading matrix, assuming that the data matrix contains missing variables and has been generated under a PPCA model with a known number of latent variable r and a known noise level sigma.
               #It uses an imputed matrix (implying also an estimation of the covariance matrix and an imputed matrix).   
               #It returns a list containing the correlation coefficient (RV coefficient) between the estimated loading matrix and the true one, and also the imputation error. 
# Arguments:
# Covar_estimated: estimated covariance matrix of size p*p.
# Yimp: imputed data matrix of size n*p.
# Y: complete data matrix of size n*p.
# B: true loading matrix of size r*p.
# M: missing-data indicator, binary matrix of size n*p. 
# r: number of latent variables. 
# sigma: noise level (standart deviation).
###### 

Results_PPCA <- function(Covar_estimated,Yimp,Y,B,M,r,sigma){
  
  ressvd <- svd(Covar_estimated - sigma ^ 2 * diag(1, ncol = p, nrow = p))
  Bestim <- sqrt(diag(ressvd$d)[1:r, 1:p]) %*% t(ressvd$u)
  Correlation <- coeffRV(t(Bestim),t(B))$rv
  mse <- MSE(Yimp*(1-M),Y*(1-M))
  
  return(list(Corr=Correlation,MSE=mse))
}



###### PPCA with MAR variables, our method
# Name: Mean_covariances_estimations_MAR
# Description: #it estimates the mean and the covariance matrix of a data matrix containing MAR variables using our method. 
               #it returns a list containing the estimated mean and covariance matrix.
# Arguments: 
# YNA: data matrix containing missing values of size n*p.
# indMissVar: indexes of the missing variables. 
# indRegVar: indexes of the pivot variables. 
# r: number of latent variables. 
######

Mean_covariances_estimations_MAR <- function(YNA,indMissVar,indRegVar,r){
  
  MeanMAR <- apply(YNA,2,mean)
  CovMAR <- var(YNA)
  
  for (jmiss in indMissVar){
    
    indvalNA <- which(is.na(YNA[,jmiss]))
    YNAcc <- YNA[-indvalNA,]
    
    indcolMAR <- sample(indRegVar,r,replace=FALSE)
    indselec <- sample(indcolMAR,1)
    indcolMAR <- indcolMAR[indcolMAR!=indselec]
    indcolMAR_all <- c(indcolMAR,indselec)
    regcc <- lm(YNAcc[, jmiss] ~  YNAcc[,indcolMAR_all])
    sumMean <- c()
    for (ind in 1:length(indcolMAR_all)){ sumMean <- c(sumMean,coefficients(regcc)[ind+1] * mean(YNA[,indcolMAR_all[ind]]))}
    MeanMAR[jmiss] <- coefficients(regcc)[1] + sum(sumMean)
    
    #Variance estimation 
    if(length(indcolMAR)==1){indselec2=indcolMAR}else{indselec2 <- sample(indcolMAR,1)}
    coeff_reg <- coefficients(lm(YNAcc[, jmiss] ~  YNAcc[, c(indselec,indselec2)]))
    M21 <- c(cov(YNAcc[, indselec], YNAcc[, jmiss]), cov(YNAcc[, indselec2], YNAcc[, jmiss])) 
    MMAR <-solve(var(YNAcc[, c(indselec, indselec2)])) 
    add1 <- coeff_reg[2] ^ 2 * var(YNA[, indselec]) + coeff_reg[3] ^ 2 * var(YNA[, indselec2])
    add2 <- 2 * coeff_reg[2] * coeff_reg[3] * cov(YNA[, indselec], YNA[, indselec2])
    CovMAR[jmiss,jmiss] <- var(YNAcc[, jmiss]) - t(M21) %*% MMAR %*% M21 + add1 + add2  
    
    
    for (j in indRegVar){
      indRegVar_minusj <- indRegVar[indRegVar!=j]
      if (length(indRegVar_minusj)==1){indselec=indRegVar_minusj}else{indselec <- sample(indRegVar_minusj,1,replace=FALSE)}
      coeff_reg <- coefficients(lm(YNAcc[, jmiss] ~  YNAcc[, c(j,indselec)]))
      ResCov <- coeff_reg[2]*(var(YNA[,j])+mean(YNA[,j])^2)+coeff_reg[1]*mean(YNA[,j])+coeff_reg[3]*(cov(YNA[,j],YNA[,indselec])+mean(YNA[,j])*mean(YNA[,indselec]))-mean(YNA[,j])*MeanMAR[jmiss] 
      CovMAR[j,jmiss] <- ResCov
      CovMAR[jmiss,j] <- CovMAR[j,jmiss]
    } 
    
    for (jmiss2 in indMissVar){
      if(jmiss2 < jmiss){
        if (length(indRegVar)==1){indselec=indRegVar}else{indselec <- sample(indRegVar,1,replace=FALSE)}
        coeff_reg <- coefficients(lm(YNAcc[, jmiss] ~  YNAcc[, c(jmiss2,indselec)]))
        ResCov <- coeff_reg[2]*(CovMAR[jmiss2,jmiss2]+MeanMAR[jmiss2]^2)+coeff_reg[1]*MeanMAR[jmiss2]+coeff_reg[3]*(CovMAR[indselec,jmiss2]+MeanMAR[jmiss]*mean(YNA[,indselec]))-MeanMAR[jmiss2]*MeanMAR[jmiss]
        CovMAR[jmiss,jmiss2] <- ResCov
        CovMAR[jmiss2,jmiss] <- CovMAR[jmiss,jmiss2]
      } 
    } 
    
  }
  
  return(list(mean=MeanMAR,cov=CovMAR))
  
}


###### SoftImpute
# Name: Mean_covariances_estimations_imputation_softImpute
# Description: #it imputes missing values using the algorithm softImpute. 
               #it returns a list containing the estimated mean and covariance matrix (deduced from the imputation) and the imputed matrix. 
# Arguments: 
# YNA: data matrix containing missing values of size n*p.
# Y: complete data matrix of size n*p.
# M: missing-data indicator, binary matrix of size n*p. 
# indMissVar: indexes of the missing variables. 
######

Mean_covariances_estimations_imputation_softImpute <- function(YNA,Y,M,indMissVar){
  
  RES <- NULL
  gridlambda1 <- seq(0, lambda0(YNA)*0.1, length = 100)
  for (k in 1:length(gridlambda1)){
    fit1 <- softImpute(as.matrix(YNA),rank=min(dim(YNA))-1,lambda=gridlambda1[k],maxit = 10000,type="svd")
    if (fit1$d[1]==0){
      X1.soft <- as.matrix(ImputeMean(YNA))
    }else if(length(fit1$d)==1){
      X1.soft <- (fit1$u*fit1$d)%*%t(fit1$v)
    }else{
      X1.soft <- (fit1$u%*%diag(fit1$d))%*%t(fit1$v)
    }
    RES[k] <- MSE(X1.soft*(1-M),Y*(1-M))
  }
  fit1 <- softImpute(as.matrix(YNA),rank=min(dim(YNA))-1,lambda=gridlambda1[which.min(RES)],maxit = 10000,type="svd")
  if (fit1$d[1]==0){
    X1.soft <- as.matrix(ImputeMean(YNA))
  }else if(length(fit1$d)==1){
    X1.soft <- (fit1$u*fit1$d)%*%t(fit1$v)
  }else{
    X1.soft <- (fit1$u%*%diag(fit1$d))%*%t(fit1$v)
  }
  
  YSoft <- YNA
  MeanSoft <- apply(YNA,2,mean)
  for (jmiss in indMissVar){
    missing <- which(is.na(YNA[,jmiss]))
    YSoft[missing,jmiss] <- X1.soft[missing,jmiss]
    MeanSoft[jmiss] <- mean(YSoft[, jmiss])
  }
  CovSoft <- var(YSoft)
  
  return(list(mean=MeanSoft,cov=CovSoft,Yimp=YSoft))
  
}


###### PPCA with MAR data, EM algorithm
# Name: Estimations_PPCA_MAREM
# Description: #it estimates the loading matrix and imputes missing values, assuming that the data matrix contains MAR variables and has been generated under a PPCA model with a known number of latent variable r and a known noise level sigma.
               #it returns a list containing the estimated mean and covariance matrix (deduced from the loading matrix estimation and from the imputation), the correlation coefficient (RV coefficient) between the estimated loading matrix and the true one, and also the imputation error. 
# Arguments: 
# YNA: data matrix containing missing values of size n*p.
# Y: complete data matrix of size n*p.
# M: missing-data indicator, binary matrix of size n*p. 
# B: true loading matrix of size r*p.
# indMissVar: indexes of the missing variables. 
# r: number of latent variables. 
# sigma: noise level (standart deviation).
######

Estimations_PPCA_MAREM <- function(YNA,Y,M,B,indMissVar,r,sigma){
  
  p <- ncol(YNA)
  coeff_init=0
  ccompt0 = 0
  while(coeff_init<0.4 & ccompt0<50){
    B_estimNew <-  matrix(rnorm(p*r), nrow = r, ncol = p)
    coeff_init=coeffRV(t(B_estimNew),t(B))$rv 
    ccompt0 = ccompt0 + 1
  }
  Mean_estimated <- rnorm(p,sd=1,mean=0)
  seuilEM <- 10^-4
  diff <- 100
  ccompt <- 0
  
  # Algo EM
  while (diff > seuilEM & ccompt < 300){ #ccompt<300
    resPPCA_MAREM_it <- PPCA_MAREM_it(B_estimNew,B,Mean_estimated,YNA,sigma,r)
    B_estimNew <- resPPCA_MAREM_it$B_estimNew
    Mean_estimated <- resPPCA_MAREM_it$Mean_estimated
    diff <- resPPCA_MAREM_it$diff
    ccompt <- ccompt + 1
  }
  
  Y_MAREM <- YNA
  Cov <- t(B_estimNew)%*%B_estimNew + sigma ^ 2 * diag(1, ncol = p, nrow = p)
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
    remp3 <- t(YNA[, enscol] - remp3bis)
    remp <- Mean_estimated[indMissVar[j]] + remp1 %*% remp2 %*% remp3
    missing <- which(is.na(YNA[,indMissVar[j]]))
    Y_MAREM[missing, indMissVar[j]] <- remp[missing]
  }
  
  MeanMAREM <- apply(Y_MAREM,2,mean)
  for (j in indMissVar){
    MeanMAREM[j] <- mean(Y_MAREM[,j])
  }
  CovMAREM <- var(Y_MAREM)
  CorrelationMAREM <- coeffRV(t(B_estimNew),t(B))$rv 
  MSEMAREM <- MSE(Y_MAREM*(1-M),Y*(1-M))
  
  return(list(mean=MeanMAREM,cov=CovMAREM,Corr=CorrelationMAREM,MSE=MSEMAREM))
}

######
# Name: PPCA_MAREM_it
# Description: it computes one iterate of the EM algorithm designed for the PPCA model in presence of MAR data. 
# Arguments: 
#B_estimOld: estimated coefficient matrix (at the iteration) of size r*p
#B: true coefficient matrix of size r*p.
#Mean_estimated: estimated mean vector of size p.
#YNA: the data matrix containing missing values of size n*p.
#sigma: noise level (standart deviation). 
#r: number of latent variables. 
######

PPCA_MAREM_it <- function(B_estimOld,B,Mean_estimated,YNA,sigma,r){
  
  n <- nrow(YNA)
  p <- ncol(YNA)
  r <- nrow(B_estimOld)
  B_estimNew <- B_estimOld
  
  #Expectation
  Mean_W <- matrix(0,ncol=n,nrow=r)
  Sigma_W <- list()
  for (k in 1:n){
    colObs <- which(is.na(YNA[k,])==0)
    sumB <- 0
    for (j in colObs){
      sumB <- sumB + B_estimOld[,j]%*%t(B_estimOld[,j])
    }
    Sigma_W[[k]] <- sigma^2*ginv(sigma^2*diag(1, ncol = r, nrow = r) + sumB)
    sumMean <- 0
    for (j in colObs){
      sumMean <- sumMean + B_estimOld[,j]*(YNA[k,j]-Mean_estimated[j])
    }
    Mean_W[,k] <- 1/(sigma^2)*Sigma_W[[k]]%*%sumMean
  }
  
  #Maximization
  for (j in 1:p){
    sum_Meanj <- 0
    rowObs <- which(is.na(YNA[,j])==0)
    for (k in rowObs){
      sum_Meanj <- sum_Meanj + (YNA[k,j]-t(B_estimOld[,j])%*%Mean_W[,k])
    }
    Mean_estimated[j] <- 1/(length(rowObs))*sum_Meanj
  }
  for (j in 1:p){
    MatB <- 0
    sumBMax <- 0
    rowObs <- which(is.na(YNA[,j])==0)
    for (k in rowObs){
      MatB <- MatB + Mean_W[,k]%*%t(Mean_W[,k]) + Sigma_W[[k]]
      sumBMax <- sumBMax + Mean_W[,k]*(YNA[k,j]-Mean_estimated[j])
    }
    B_estimNew[,j] <- ginv(MatB)%*%sumBMax
  }
  diff <- norm(B_estimNew - B_estimOld, type = "F") / (norm(B_estimOld, type = "F") + 10 ^ -3)
  
  cor <- coeffRV(t(B_estimNew),t(B))$rv 
  
  return(list(B_estimNew=B_estimNew,Mean_estimated=Mean_estimated,diff=diff,cor=cor))
}


