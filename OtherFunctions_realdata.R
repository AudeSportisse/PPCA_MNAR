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
# Name: Results_PPCA_imputation_realdata
# Description: #it estimates the loading matrix and imputes the missing values of an incomplete matrix, assuming that the data matrix contains missing variables and has been generated under a PPCA model with a known number of latent variable r and a known noise level sigma.
               #It uses an estimation of the covariance matrix and the mean.        
               #It returns the imputed matrix.
# Arguments:
# Covar_estimated: estimated covariance matrix of size p*p.
# Mean_estimated: estimated mean vector of size p.
# YNA: data matrix containing missing values of size n*p.
# indMissVar: indexes of the missing variables. 
# r: number of latent variables. 
# sigma: noise level (standart deviation).
###### 

Results_PPCA_imputation_realdata <- function(Covar_estimated,Mean_estimated,YNA,indMissVar,r,sigma){ 
  n <- nrow (YNA)
  p <- ncol(YNA)
  
  ressvd <- svd(Covar_estimated - sigma ^ 2 * diag(1, ncol = p, nrow = p))
  BMNAR <- sqrt(diag(ressvd$d)[1:r, 1:p]) %*% t(ressvd$u)
  
  YNA_imputed <- YNA
  Cov <- t(BMNAR)%*%BMNAR + sigma ^ 2 * diag(1, ncol = p, nrow = p)
  YNA_imputed <- YNA
  for (j in 1:length(indMissVar)){
    for (i in 1:nrow(YNA)){
      if (is.na(YNA[i,indMissVar[j]])){
        enscol <- 1:p
        for (jbis in indMissVar){
          enscol <- enscol[enscol!=jbis]
        }
        for (jbis in enscol){
          if(is.na(YNA[i,jbis])){
            enscol <- enscol[enscol!=jbis]
          }
        }
        remp1 <- c()
        remp3bis <- c()
        for (jter in enscol){
          remp1 <- c(remp1,Cov[jter,indMissVar[j]])
          remp3bis <- c(remp3bis,mean(YNA[,jter],na.rm=TRUE)) #mean(YNAimp[,jter])
        }
        remp2 <- solve(Cov[enscol, enscol]) 
        remp3 <- t(YNA[i, enscol] - remp3bis)
        remp <- Mean_estimated[indMissVar[j]] + t(as.matrix(remp1)) %*% remp2 %*% as.matrix(remp3)
        YNA_imputed[i, indMissVar[j]] <- remp
      }
    }
  }
  
  return(YNA_imputed)
}




###### SoftImpute
# Name: Mean_covariances_estimations_imputation_softImpute_realdata
# Description: #it imputes missing values using the algorithm softImpute. 
               #it returns a list containing the imputed matrix. 
# Arguments: 
# YNA: data matrix containing missing values of size n*p.
# Y: complete data matrix of size n*p.
# M: missing-data indicator, binary matrix of size n*p. 
# indMissVar: indexes of the missing variables. 
######

Mean_covariances_estimations_imputation_softImpute_realdata <- function(YNA,M,indMissVar){
  
  fit1 <- softImpute(as.matrix(YNA),rank=min(dim(YNA))-1,lambda=lambda0(YNA)*0.08,maxit = 10000,type="svd")
  if (fit1$d[1]==0){
    X1.soft <- as.matrix(ImputeMean(YNA))
  }else if(length(fit1$d)==1){
    X1.soft <- (fit1$u*fit1$d)%*%t(fit1$v)
  }else{
    X1.soft <- (fit1$u%*%diag(fit1$d))%*%t(fit1$v)
  }
  
  YSoft <- YNA
  for (jmiss in indMissVar){
    missing <- which(is.na(YNA[,jmiss]))
    YSoft[missing,jmiss] <- X1.soft[missing,jmiss]
  }
  
  return(list(Yimp=YSoft))
  
}

###### PPCA with MAR data, EM algorithm
# Name: Estimations_PPCA_MAREM_realdata
# Description: #it estimates the loading matrix and imputes missing values, assuming that the data matrix contains MAR variables and has been generated under a PPCA model with a known number of latent variable r and a known noise level sigma.
               #it returns a list containingthe imputed matrix. 
# Arguments: 
# YNA: data matrix containing missing values of size n*p.
# indMissVar: indexes of the missing variables. 
# r: number of latent variables. 
# sigma: noise level (standart deviation).
######

Estimations_PPCA_MAREM_realdata <- function(YNA,indMissVar,r,sigma){
  
  p <- ncol(YNA)
  ressvd <- svd(var(YNA,na.rm=TRUE)  - sigma ^ 2 * diag(1, ncol = p, nrow = p))
  B_estimNew  <- sqrt(diag(ressvd$d)[1:r, 1:p]) %*% t(ressvd$u)
  Mean_estimated <- rnorm(p,sd=1,mean=0)
  seuilEM <- 10^-5
  diff <- 100
  ccompt <- 0
  
  # Algo EM
  while (diff > seuilEM & ccompt < 200){
    resPPCA_MAREM_it <- PPCA_MAREM_it_realdata(B_estimNew,Mean_estimated,YNA,sigma,r)
    B_estimNew <- resPPCA_MAREM_it$B_estimNew
    Mean_estimated <- resPPCA_MAREM_it$Mean_estimated
    diff <- resPPCA_MAREM_it$diff
    ccompt <- ccompt + 1
  }
  
  Y_MAREM <- YNA
  Cov <- t(B_estimNew)%*%B_estimNew + sigma ^ 2 * diag(1, ncol = p, nrow = p)
  for (j in 1:length(indMissVar)){
    for (i in 1:nrow(YNA)){
      if (is.na(YNA[i,indMissVar[j]])){
        enscol <- 1:p
        for (jbis in indMissVar){
          enscol <- enscol[enscol!=jbis]
        }
        for (jbis in enscol){
          if(is.na(YNA[i,jbis])){
            enscol <- enscol[enscol!=jbis]
          }
        }
        remp1 <- c()
        remp3bis <- c()
        for (jter in enscol){
          remp1 <- c(remp1,Cov[jter,indMissVar[j]])
          remp3bis <- c(remp3bis,mean(YNA[,jter],na.rm=TRUE))
        }
        remp2 <- solve(Cov[enscol, enscol]) 
        remp3 <- t(YNA[i, enscol] - remp3bis)
        remp <- Mean_estimated[indMissVar[j]] + t(as.matrix(remp1)) %*% remp2 %*% as.matrix(remp3)
        Y_MAREM[i, indMissVar[j]] <- remp
      }
    }
  }
  
  
  return(list(Y_estim=Y_MAREM))
}



######
# Name: PPCA_MAREM_it_realdata
# Description: it computes one iterate of the EM algorithm designed for the PPCA model in presence of MAR data. 
# Arguments: 
#B_estimOld: estimated coefficient matrix (at the iteration) of size r*p
#Mean_estimated: estimated mean vector of size p.
#YNA: the data matrix containing missing values of size n*p.
#sigma: noise level (standart deviation). 
#r: number of latent variables. 
######

PPCA_MAREM_it_realdata <- function(B_estimOld,Mean_estimated,YNA,sigma,r){
  
  n <- nrow(YNA)
  p <- ncol(YNA)
  r <- nrow(B_estimOld)
  B_estimNew <- B_estimOld
  
  #Expectation
  Moy_W <- matrix(0,ncol=n,nrow=r)
  Sigma_W <- list()
  for (k in 1:n){
    colObs <- which(is.na(YNA[k,])==0)
    sumB <- 0
    for (j in colObs){
      sumB <- sumB + B_estimOld[,j]%*%t(B_estimOld[,j])
    }
    Sigma_W[[k]] <- sigma^2*ginv(sigma^2*diag(1, ncol = r, nrow = r) + sumB)
    sumMoy <- 0
    for (j in colObs){
      sumMoy <- sumMoy + B_estimOld[,j]*(YNA[k,j]-Mean_estimated[j])
    }
    Moy_W[,k] <- 1/(sigma^2)*Sigma_W[[k]]%*%sumMoy
  }
  
  #Maximization
  for (j in 1:p){
    sum_moyj <- 0
    rowObs <- which(is.na(YNA[,j])==0)
    for (k in rowObs){
      sum_moyj <- sum_moyj + (YNA[k,j]-t(B_estimOld[,j])%*%Moy_W[,k])
    }
    Mean_estimated[j] <- 1/(length(rowObs))*sum_moyj
  }
  for (j in 1:p){
    MatB <- 0
    sumBMax <- 0
    rowObs <- which(is.na(YNA[,j])==0)
    for (k in rowObs){
      MatB <- MatB + Moy_W[,k]%*%t(Moy_W[,k]) + Sigma_W[[k]]
      sumBMax <- sumBMax + Moy_W[,k]*(YNA[k,j]-Mean_estimated[j])
    }
    B_estimNew[,j] <- ginv(MatB)%*%sumBMax
  }
  diff <- norm(B_estimNew - B_estimOld, type = "F") / (norm(B_estimOld, type = "F") + 10 ^ -3)
  
  
  return(list(B_estimNew=B_estimNew,Mean_estimated=Mean_estimated,diff=diff))
}

