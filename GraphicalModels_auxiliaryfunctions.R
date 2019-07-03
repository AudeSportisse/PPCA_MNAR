######
#Preliminary functions
######


#MSE (normalized)
MSE <- function(X1, X2){ return(norm(as.matrix(X1 - X2),type="F")/norm(as.matrix(X2),type="F"))}


#Imputation by the mean
ImputeMean0 <- function(tab){
  m <- apply(tab, 2, mean, na.rm = TRUE)
  tab <- sapply(1:ncol(tab), function(x) ifelse(is.na(tab[,x]), m[x], tab[,x]))
  tab <- as.data.frame(tab)
  return(tab)
}


# Imputation and estimation by the mean
ImputeMean <- function(tab){
  m <- apply(tab, 2, mean, na.rm = TRUE)
  tab <- sapply(1:ncol(tab), function(x) replicate(nrow(tab),m[x]))
  tab <- as.data.frame(tab)
  return(tab)
}


######
# Name: Results_graph
# Date: 25/06/2019
# Description: it estimates the coefficient matrix B using SVD and imputes values in the data matrix Y using the conditionnal expectation given the observed variables. It returns the coefficient RV for the estimation of B and the normalized MSE.
# Arguments: 
#Covar: estimated covariance matrix.
#Moy: estimated means of the missing variables. 
#YNA: matrix containing missing values.
#Y: matrix. 
#B: coefficient matrix. 
#####


Results_graph <- function(Covar,Moy,YNA,Y,B,indcolNA,M,r,p){ 
  ressvd <- svd(Covar - sigma ^ 2 * diag(1, ncol = p, nrow = p))
  BMNAR <- sqrt(diag(ressvd$d)[1:r, 1:p]) %*% t(ressvd$u)
  Corr <- coeffRV(t(BMNAR),t(B))$rv
  
  Cov <- t(BMNAR)%*%BMNAR + sigma ^ 2 * diag(1, ncol = p, nrow = p)
  for (j in 1:length(indcolNA)){
    enscol <- 1:p
    for (jbis in indcolNA){
      enscol <- enscol[enscol!=jbis]
    }
    remp1 <- c()
    remp3bis <- c()
    for (jter in enscol){
      remp1 <- c(remp1,Cov[jter,indcolNA[j]])
      remp3bis <- c(remp3bis,mean(YNA[,jter]))
    }
    remp2 <- solve(Cov[enscol, enscol]) 
    remp3 <- t(YNA[, enscol] - matrix(remp3bis, nrow = n, ncol = (p-length(indcolNA)), byrow = TRUE))
    remp <- Moy[j] + remp1 %*% remp2 %*% remp3
    missing <- which(is.na(YNA[,indcolNA[j]]))
    YNA[missing, indcolNA[j]] <- remp[missing]
  }
  mse <- MSE(YNA*(1-M),Y*(1-M)) 
  return(list(Corr,mse))
}


######
# Name: Mean_graphical_algebraic
# Date: 25/06/2019
# Description: it estimates the mean of the missing variable Yj using the graphical method.
# Arguments: 
#indcolNA: indexes of the missing variables.
#indcolReg:  indexes of the observed variables which are used for the regressions. 
#j: index of the missing variable.
#YNA: the data matrix containing missing values. 
#YNAcc: the reduced data matrix containing the observations such that the variable j is always observed
#r: number of latent variables. 
#####


Mean_graphical_algebraic <- function(indcolNA,indcolReg,j,YNA,YNAcc,r){
  
  Mean <- c()
  coeff_mean <- c()
  
  for (u0 in 1:length(indcolReg)){
    
    indselec1 <- indcolReg[u0]
    indcolReg_mselec1 <- indcolReg[indcolReg!=indselec1]
    if (r-1-length(indcolReg_mselec1)==0){
      regcc <- lm(YNAcc[, indselec1] ~  YNAcc[,c(indcolReg_mselec1,indcolNA[j])]) 
      sumMean <- c()
      for (ind in 1:length(indcolReg_it)){ sumMean <- c(sumMean,coefficients(regcc)[ind+1] * mean(YNA[,indcolReg_mselec1[ind]]))}
      Mean <- c(Mean,(mean(YNA[, indselec1]) - coefficients(regcc)[1] - sum(sumMean)) / coefficients(regcc)[length(indcolReg_mselec1)+2])
      coeff_mean <- c(coeff_mean,coefficients(regcc)[length(indcolReg_mselec1)+2])
    }else{
      Choice_indcolReg <- t(combn(indcolReg_mselec1,r-1))
      for (u1 in 1:length(indcolReg_mselec1)){
        
        indcolReg_it <- Choice_indcolReg[u1,]
        
        regcc <- lm(YNAcc[, indselec1] ~  YNAcc[,c(indcolReg_it,indcolNA[j])]) 
        sumMean <- c()
        for (ind in 1:length(indcolReg_it)){ sumMean <- c(sumMean,coefficients(regcc)[ind+1] * mean(YNA[,indcolReg_it[ind]]))}
        Mean <- c(Mean,(mean(YNA[, indselec1]) - coefficients(regcc)[1] - sum(sumMean)) / coefficients(regcc)[length(indcolReg_it)+2])
        coeff_mean <- c(coeff_mean,coefficients(regcc)[length(indcolReg_it)+2])
        
      }
    }
  }
  
  return(list(Mean=Mean,coeff_mean=coeff_mean))
  
}


######
# Name: Variance_graphical
# Date: 25/06/2019
# Description: it estimates the variance of the missing variable Yj using the graphical method.
# Arguments: 
#indcolNA: indexes of the missing variables.
#indcolReg:  indexes of the observed variables which are used for the regressions. 
#j: index of the missing variable.
#YNA: the data matrix containing missing values. 
#YNAcc: the reduced data matrix containing the observations such that the variable j is always observed
#r: number of latent variables. 
#####


Variance_graphical <- function(indcolNA,indcolReg,j,YNA,YNAcc,r){
  
  Var <- c()
  coeff_var <- c()
  
  for (u0 in 1:length(indcolReg)){
    
    indselec1 <- indcolReg[u0]
    indcolReg_mselec1 <- indcolReg[indcolReg!=indselec1]
    
    for (u in 1:length(indcolReg_mselec1)){
      
      indselec2 <- indcolReg_mselec1[u]
      
      regcc <- lm(YNAcc[, indselec1] ~  YNAcc[, c(indcolNA[j],indselec2)])
      regcc2 <- lm(YNAcc[,indselec2] ~ YNAcc[, indcolNA[j]])
      b <- coefficients(regcc)[2]
      c <- coefficients(regcc)[3]
      abar <- coefficients(regcc2)[2]
      a <- (1 / b) * ((cov(YNA[, indselec1], YNA[, indselec2]) / var(YNA[, indselec2])) - c)
      Var <- c(Var,(a * var(YNA[, indselec2])) / abar)
      coeff_var <- c(coeff_var,min(abs(abar),abs(b)))
      
    }
  }
  
  return(list(Var=Var,coeff_var=coeff_var))
  
  
}


######
# Name: Covariance_obsmiss_graphical
# Date: 25/06/2019
# Description: it estimates the covariance between the missing variable Yj and the observed variable Yjother using the graphical method. 
# Arguments: 
#indcolNA: indexes of the missing variables.
#indcolReg:  indexes of the observed variables which are used for the regressions. 
#j: index of the first missing variable.
#jother: index of the second missing variable.
#YNA: the data matrix containing missing values. 
#YNAcc: the reduced data matrix containing the observations such that the variable j is always observed.
#r: number of latent variables. 
#####


Covariance_obsmiss_graphical <- function(indcolNA,indcolReg,j,jother,YNA,YNAcc,r){
  
  cov <- c()
  coeff_cov <- c()
  indcolReg_it <- indcolReg[indcolReg!=jother]
  for (indselec in indcolReg_it){
    regcc <- lm(YNAcc[, indselec] ~  YNAcc[, c(indcolNA[j],jother)])
    b <- coefficients(regcc)[2]
    c <- coefficients(regcc)[3]
    a <- (1 / b) * ((cov(YNA[, indselec], YNA[, jother]) / var(YNA[, jother])) - c)
    cov <- c(cov,a*var(YNA[,jother]))
    coeff_cov <- c(coeff_cov,b)
  }
  
  return(list(cov=cov,coeff_cov=coeff_cov))
  
}


######
# Name: Covariance_2missing_graphical
# Date: 25/06/2019
# Description: it estimates the covariance between two missing variables Yj and Yjother using the graphical method. 
# Arguments: 
#Variance: computed variances of the missing variables. 
#CovarianceMat: computed covariances associated to the missing variables. 
#meth: chosen method (best, agg, noagg).
#indcolNA: indexes of the missing variables.
#indcolReg:  indexes of the observed variables which are used for the regressions. 
#j: index of the first missing variable.
#jother: index of the second missing variable.
#YNA: the data matrix containing missing values. 
#YNAcc: the reduced data matrix containing the observations such that the variable j is always observed.
#r: number of latent variables. 
#####


Covariance_2missing_graphical <- function(Variance,CovarianceMat,meth,indcolNA,indcolReg,j,jother,YNA,YNAcc,r){
  
  cov <- c()
  coeff_cov <- c()
  for (indselec in indcolReg){
    indvalNA2 <- which(is.na(YNAcc[,indcolNA[jother]]))
    YNAccbis <- YNAcc[-indvalNA2,]
    
    regcc <- lm(YNAccbis[, indselec] ~  YNAccbis[, c(indcolNA[j],indcolNA[jother])])
    b <- coefficients(regcc)[2]
    c <- coefficients(regcc)[3]
    a <- (1 / b) * (CovarianceMat[indselec,indcolNA[jother]] / Variance[jother] - c)
    cov <- c(cov,a*Variance[jother])
    coeff_cov <- c(coeff_cov,b)
  }
  
  if (meth=="best"){
    res  <- cov[which.max(abs(coeff_cov))]
  }  
  if(meth=="agg"){
    res  <- median(cov)
  }
  if(meth=="noagg"){
    res  <- sample(cov,1)
  }
  
  return(res)

}


######
# Name: Variance_algebraic
# Date: 25/06/2019
# Description: it estimates the variance of the missing variable Yj using the algebraic method.
# Arguments: 
#Mean: computed mean of the missing variable. 
#meth: chosen method (best, agg, noagg).
#shrink: 1 to return the results for the algebraic MNAR method by using a shrinkage for the matrix inversions, 0 (default) otherwise.
#indcolNA: indexes of the missing variables.
#indcolReg:  indexes of the observed variables which are used for the regressions. 
#j: index of the missing variable.
#YNA: the data matrix containing missing values. 
#YNAcc: the reduced data matrix containing the observations such that the variable j is always observed
#CovTheo: theoretical covariance matrix.
#r: number of latent variables.
#####


Variance_algebraic <- function(Mean,meth,shrink,indcolNA,indcolReg,j,YNA,YNAcc,CovTheo,r){
  
  Vars <- c()
  Var <- c()
  Choice_indcolReg <- t(combn(indcolReg,r))
  
  for (ind2 in 1:dim(Choice_indcolReg)[1]){
    
    Choice2_indcolReg <- Choice_indcolReg[ind2,]
    
    for (ind3 in 1:length(Choice2_indcolReg)){
      
      indcol_it <- c(Choice2_indcolReg[ind3],Choice2_indcolReg[Choice2_indcolReg!=Choice2_indcolReg[ind3]])
      
      Mres <- matrix(0,nrow=r+1,ncol=r+1)
      vec <- numeric(r+1)
      for (ind in 1:r){
        indcol_it2 <- indcol_it[indcol_it!=indcol_it[ind]]
        regcc <-  lm(YNAcc[, indcol_it[ind]] ~ YNAcc[, indcolNA[j]] + YNAcc[, indcol_it2])
        coeff_regcc <- coefficients(regcc)
        if (ind==1){
          Mres[1,1] <- coeff_regcc[2]^2
          for (ind2 in 3:(r+1)){
            Mres[1,ind2] <- 2*coeff_regcc[ind2]*coeff_regcc[2]
          }
          Q <- var(YNAcc[,indcol_it[ind]])-t(c(cov(YNAcc[, indcol_it[ind]], YNAcc[, indcolNA[j]]), cov(YNAcc[, indcol_it[ind]], YNAcc[, indcol_it2])))%*%solve(var(YNAcc[,c(indcolNA[j],indcol_it2)]))%*%c(cov(YNAcc[, indcol_it[ind]], YNAcc[, indcolNA[j]]), cov(YNAcc[, indcol_it[ind]], YNAcc[, indcol_it2]))
          if (r==2){
            #vec[1] <- var(YNA[, indcol_it[ind]]) - Q - t(coeff_regcc[3:(r+1)])*var(YNA[,indcol_it2])*coeff_regcc[3:(r+1)]
            vec[1] <- var(YNA[, indcol_it[ind]]) - Q - t(coeff_regcc[3:(r+1)])*var(YNA[,indcol_it2])*coeff_regcc[3:(r+1)] + sum(sigma^2*coeff_regcc[2:(r+1)]^2)
          }else{
            #vec[1] <- var(YNA[, indcol_it[ind]]) - Q - t(coeff_regcc[3:(r+1)])%*%var(YNA[,indcol_it2])%*%coeff_regcc[3:(r+1)]
            vec[1] <- var(YNA[, indcol_it[ind]]) - Q - t(coeff_regcc[3:(r+1)])%*%var(YNA[,indcol_it2])%*%coeff_regcc[3:(r+1)] + sum(sigma^2*coeff_regcc[2:(r+1)]^2)
          }
        }
        if ((ind+2)<=r+1){
          Mres[ind+1,] <- -c(coeff_regcc[2:(ind+1)],-1,coeff_regcc[(ind+2):(r+1)]) 
        }else{
          Mres[ind+1,] <- -c(coeff_regcc[2:(ind+1)],-1)
        }
        if (r>2){
          #vec[ind+1] <- (coeff_regcc%*%c(1,Mean[j],apply(YNA[,indcol_it2],2,mean)) - mean(YNA[,indcol_it[ind]])) * Mean[j]
          vec[ind+1] <- (coeff_regcc%*%c(1,Mean[j],apply(YNA[,indcol_it2],2,mean)) - mean(YNA[,indcol_it[ind]])) * Mean[j] - sigma^2*coeff_regcc[2]
        }else{
          #vec[ind+1] <- (coeff_regcc%*%c(1,Mean[j],mean(YNA[,indcol_it2])) - mean(YNA[,indcol_it[ind]])) * Mean[j] 
          vec[ind+1] <- (coeff_regcc%*%c(1,Mean[j],mean(YNA[,indcol_it2])) - mean(YNA[,indcol_it[ind]])) * Mean[j] - sigma^2*coeff_regcc[2]
        }        
      }
      #Shrinkage
      if (shrink==1){
        listLambda=seq(0,0.01,0.0001)
        ChoiceLambda=c()
        for (lambda in listLambda){
          Mres_shrink=Mres+lambda*diag(1,nrow=r+1,ncol=r+1)
          Res <- ginv(Mres_shrink,tol = 10^-20)%*%vec
          ChoiceLambda <- c(ChoiceLambda,abs(Res[1]-CovTheo[indcolNA[j],indcolNA[j]]))
        }
        Mres_shrink=Mres+listLambda[which.min(ChoiceLambda)]*diag(1,nrow=r+1,ncol=r+1)
        Res <- ginv(Mres_shrink,tol = 10^-20)%*%vec
        Vars <- c(Vars,Res[1])
      }
      #No shrinkage
      Res <- ginv(Mres,tol = 10^-20)%*%vec
      Var <- c(Var,Res[1])
    }
  }
  
  if (shrink==1){
    if (meth=="best" | meth=="agg"){
      res <- c(median(Var),median(Vars))
    }else{
      res <- c(sample(Var,1),median(Vars,1))
    } 
  }else{
    if (meth=="best" | meth=="agg"){
      res <-  median(Var)
    }else{
      res <-  sample(Var,1)
    } 
  }
  
  return(res)
  
}


######
# Name: Covariance_obsmiss_algebraic
# Date: 25/06/2019
# Description: it estimates the covariance between the missing variable Yj and the observed variable Yindselec using the algebraic method.
# Arguments: 
#Mean: computed mean of the missing variable. 
#meth: chosen method (best, agg, noagg).
#shrink: 1 to return the results for the algebraic MNAR method by using a shrinkage for the matrix inversions, 0 (default) otherwise.
#indselec: index of the missing variable. 
#indcolNA: indexes of the missing variables.
#indcolReg:  indexes of the observed variables which are used for the regressions. 
#j: index of the missing variable.
#YNA: the data matrix containing missing values. 
#YNAcc: the reduced data matrix containing the observations such that the variable j is always observed
#CovTheo: theoretical covariance matrix.
#r: number of latent variables.
#####


Covariance_obsmiss_algebraic <- function(Mean,meth,shrink,indselec,indcolNA,indcolReg,j,YNA,YNAcc,CovTheo,r){
  
  cov <- c()
  covs <- c()
  
  indcolReg_mselec <- indcolReg[indcolReg!=indselec]
  if (length(indcolReg_mselec)==r-1){
    Choice_indcolReg <- indcolReg_mselec
    indcol_it <- c(indselec,Choice_indcolReg)
    Mres <- matrix(0,nrow=r+1,ncol=r+1)
    vec <- numeric(r+1)
    for (ind in 1:r){
      indcol_it2 <- indcol_it[indcol_it!=indcol_it[ind]]
      regcc <-  lm(YNAcc[, indcol_it[ind]] ~ YNAcc[, indcolNA[j]] + YNAcc[, indcol_it2])
      coeff_regcc <- coefficients(regcc)
      if (ind==1){
        Mres[1,1] <- coeff_regcc[2]^2
        for (ind2 in 3:(r+1)){
          Mres[1,ind2] <- 2*coeff_regcc[ind2]*coeff_regcc[2]
        }
        Q <- var(YNAcc[,indcol_it[ind]])-t(c(cov(YNAcc[, indcol_it[ind]], YNAcc[, indcolNA[j]]), cov(YNAcc[, indcol_it[ind]], YNAcc[, indcol_it2])))%*%solve(var(YNAcc[,c(indcolNA[j],indcol_it2)]))%*%c(cov(YNAcc[, indcol_it[ind]], YNAcc[, indcolNA[j]]), cov(YNAcc[, indcol_it[ind]], YNAcc[, indcol_it2]))
        if (r==2){
          #vec[1] <- var(YNA[, indcol_it[ind]]) - Q - t(coeff_regcc[3:(r+1)])*var(YNA[,indcol_it2])*coeff_regcc[3:(r+1)]
          vec[1] <- var(YNA[, indcol_it[ind]]) - Q - t(coeff_regcc[3:(r+1)])*var(YNA[,indcol_it2])*coeff_regcc[3:(r+1)] + sum(sigma^2*coeff_regcc[2:(r+1)]^2)
        }else{
          #vec[1] <- var(YNA[, indcol_it[ind]]) - Q - t(coeff_regcc[3:(r+1)])%*%var(YNA[,indcol_it2])%*%coeff_regcc[3:(r+1)]
          vec[1] <- var(YNA[, indcol_it[ind]]) - Q - t(coeff_regcc[3:(r+1)])%*%var(YNA[,indcol_it2])%*%coeff_regcc[3:(r+1)] + sum(sigma^2*coeff_regcc[2:(r+1)]^2)
        }
      }
      if ((ind+2)<=r+1){
        Mres[ind+1,] <- -c(coeff_regcc[2:(ind+1)],-1,coeff_regcc[(ind+2):(r+1)]) 
      }else{
        Mres[ind+1,] <- -c(coeff_regcc[2:(ind+1)],-1)
      }
      if(r>2){
        #vec[ind+1] <- (coeff_regcc%*%c(1,Mean[j],apply(YNA[,indcol_it2],2,mean)) - mean(YNA[,indcol_it[ind]])) * Mean[j]
        vec[ind+1] <- (coeff_regcc%*%c(1,Mean[j],apply(YNA[,indcol_it2],2,mean)) - mean(YNA[,indcol_it[ind]])) * Mean[j] - sigma^2*coeff_regcc[2]
      }else{
        #vec[ind+1] <- (coeff_regcc%*%c(1,Mean[j],mean(YNA[,indcol_it2])) - mean(YNA[,indcol_it[ind]])) * Mean[j]
        vec[ind+1] <- (coeff_regcc%*%c(1,Mean[j],mean(YNA[,indcol_it2])) - mean(YNA[,indcol_it[ind]])) * Mean[j] - sigma^2*coeff_regcc[2]
      }          
    }
    #Shrinkage
    if (shrink==1){
      listLambda=seq(0,0.01,0.0001)
      ChoiceLambda=c()
      for (lambda in listLambda){
        Mres_shrink=Mres+lambda*diag(1,nrow=r+1,ncol=r+1)
        Res <- ginv(Mres_shrink,tol = 10^-20)%*%vec
        ChoiceLambda <- c(ChoiceLambda,abs(Res[2]-CovTheo[indselec,indcolNA[j]]))
      }
      Mres_shrink=Mres+listLambda[which.min(ChoiceLambda)]*diag(1,nrow=r+1,ncol=r+1)
      Res <- ginv(Mres_shrink,tol = 10^-20)%*%vec
      covs <- c(covs,Res[2])
    }
    #No shrinkage
    Res <- ginv(Mres,tol = 10^-20)%*%vec
    cov <- c(cov,Res[2])
  }else{
    Choice_indcolReg <- t(combn(indcolReg_mselec,r-1))
  
    for (indbis in 1:dim(Choice_indcolReg)[1]){
      
      if (length(indcolReg_mselec)==1){
        indcol_it <- c(indselec,Choice_indcolReg[indbis])
      }else{
        indcol_it <- c(indselec,Choice_indcolReg[indbis,])
      }
      Mres <- matrix(0,nrow=r+1,ncol=r+1)
      vec <- numeric(r+1)
      for (ind in 1:r){
        indcol_it2 <- indcol_it[indcol_it!=indcol_it[ind]]
        regcc <-  lm(YNAcc[, indcol_it[ind]] ~ YNAcc[, indcolNA[j]] + YNAcc[, indcol_it2])
        coeff_regcc <- coefficients(regcc)
        if (ind==1){
          Mres[1,1] <- coeff_regcc[2]^2
          for (ind2 in 3:(r+1)){
            Mres[1,ind2] <- 2*coeff_regcc[ind2]*coeff_regcc[2]
          }
          Q <- var(YNAcc[,indcol_it[ind]])-t(c(cov(YNAcc[, indcol_it[ind]], YNAcc[, indcolNA[j]]), cov(YNAcc[, indcol_it[ind]], YNAcc[, indcol_it2])))%*%solve(var(YNAcc[,c(indcolNA[j],indcol_it2)]))%*%c(cov(YNAcc[, indcol_it[ind]], YNAcc[, indcolNA[j]]), cov(YNAcc[, indcol_it[ind]], YNAcc[, indcol_it2]))
          if (r==2){
            #vec[1] <- var(YNA[, indcol_it[ind]]) - Q - t(coeff_regcc[3:(r+1)])*var(YNA[,indcol_it2])*coeff_regcc[3:(r+1)]
            vec[1] <- var(YNA[, indcol_it[ind]]) - Q - t(coeff_regcc[3:(r+1)])*var(YNA[,indcol_it2])*coeff_regcc[3:(r+1)] + sum(sigma^2*coeff_regcc[2:(r+1)]^2)
          }else{
            #vec[1] <- var(YNA[, indcol_it[ind]]) - Q - t(coeff_regcc[3:(r+1)])%*%var(YNA[,indcol_it2])%*%coeff_regcc[3:(r+1)]
            vec[1] <- var(YNA[, indcol_it[ind]]) - Q - t(coeff_regcc[3:(r+1)])%*%var(YNA[,indcol_it2])%*%coeff_regcc[3:(r+1)] + sum(sigma^2*coeff_regcc[2:(r+1)]^2)
          }
        }
        if ((ind+2)<=r+1){
          Mres[ind+1,] <- -c(coeff_regcc[2:(ind+1)],-1,coeff_regcc[(ind+2):(r+1)]) 
        }else{
          Mres[ind+1,] <- -c(coeff_regcc[2:(ind+1)],-1)
        }
        if(r>2){
          #vec[ind+1] <- (coeff_regcc%*%c(1,Mean[j],apply(YNA[,indcol_it2],2,mean)) - mean(YNA[,indcol_it[ind]])) * Mean[j]
          vec[ind+1] <- (coeff_regcc%*%c(1,Mean[j],apply(YNA[,indcol_it2],2,mean)) - mean(YNA[,indcol_it[ind]])) * Mean[j] - sigma^2*coeff_regcc[2]
        }else{
          #vec[ind+1] <- (coeff_regcc%*%c(1,Mean[j],mean(YNA[,indcol_it2])) - mean(YNA[,indcol_it[ind]])) * Mean[j]
          vec[ind+1] <- (coeff_regcc%*%c(1,Mean[j],mean(YNA[,indcol_it2])) - mean(YNA[,indcol_it[ind]])) * Mean[j] - sigma^2*coeff_regcc[2]
        }          
      }
      #Shrinkage
      if (shrink==1){
        listLambda=seq(0,0.01,0.0001)
        ChoiceLambda=c()
        for (lambda in listLambda){
          Mres_shrink=Mres+lambda*diag(1,nrow=r+1,ncol=r+1)
          Res <- ginv(Mres_shrink,tol = 10^-20)%*%vec
          ChoiceLambda <- c(ChoiceLambda,abs(Res[2]-CovTheo[indselec,indcolNA[j]]))
        }
        Mres_shrink=Mres+listLambda[which.min(ChoiceLambda)]*diag(1,nrow=r+1,ncol=r+1)
        Res <- ginv(Mres_shrink,tol = 10^-20)%*%vec
        covs <- c(covs,Res[2])
      }
      #No shrinkage
      Res <- ginv(Mres,tol = 10^-20)%*%vec
      cov <- c(cov,Res[2])
    }
  }
  
  if (shrink==1){
    if (meth=="best" | meth=="agg"){
      res <- c(median(cov),median(covs))
    }else{
      res <- c(sample(cov,1),median(covs,1))
    } 
  }else{
    if (meth=="best" | meth=="agg"){
      res <-  median(cov)
    }else{
      res <-  sample(cov,1)
    } 
  }
  
  return(res)
  
}


######
# Name: Covariance_2missing_algebraic
# Date: 25/06/2019
# Description: it estimates the covariance between two missing variables Yj and Yjother using the algebraic method.
# Arguments: 
#Mean: computed mean of the missing variable. 
#meth: chosen method (best, agg, noagg).
#shrink: 1 to return the results for the algebraic MNAR method by using a shrinkage for the matrix inversions, 0 (default) otherwise.
#indcolNA: indexes of the missing variables.
#indcolReg:  indexes of the observed variables which are used for the regressions. 
#j: index of the first missing variable.
#jother: index of the second missing variable.
#YNA: the data matrix containing missing values. 
#YNAcc: the reduced data matrix containing the observations such that the variable j is always observed
#CovTheo: theoretical covariance matrix.
#r: number of latent variables.
#####


Covariance_2missing_algebraic <- function(Variance,CovarianceMat,meth,shrink,indcolNA,indcolReg,j,jother,YNA,YNAcc,CovTheo,r){
  
  cov <- c()
  coeff_cov <- c()
  for (u0 in 1:length(indcolReg)){
    indselec <- indcolReg[u0]
    indcolReg_init <- indcolReg[indcolReg!=indselec]
    Choice_indcolReg_init <- t(combn(indcolReg_init,r-2))
    for (u1 in 1:dim(Choice_indcolReg_init)[1]){
      indcol_it <- Choice_indcolReg_init[u1,]
      
      indvalNA2 <- which(is.na(YNAcc[,indcolNA[jother]]))
      YNAccbis <- YNAcc[-indvalNA2,]
      
      regcc <- lm(YNAccbis[, indselec] ~  YNAccbis[, c(indcolNA[j],indcolNA[jother],indcol_it)])
      b <- coefficients(regcc)[2]
      c <- coefficients(regcc)[3]
      Q <- var(YNAccbis[,indselec])-t(c(cov(YNAccbis[, indselec], YNAccbis[, indcolNA[j]]), cov(YNAccbis[, indselec], YNAccbis[, indcolNA[jother]]), cov(YNAccbis[, indselec], YNAccbis[, indcol_it])))%*%solve(var(YNAccbis[,c(indcolNA[j],indcolNA[jother],indcol_it)]))%*%c(cov(YNAccbis[, indselec], YNAccbis[, indcolNA[j]]), cov(YNAccbis[, indselec], YNAccbis[, indcolNA[jother]]),cov(YNAccbis[, indselec], YNAccbis[, indcol_it]))
      t1 <- 0
      t2 <- 0
      if(length(indcol_it)>=1){
        for (ubis in 1:length(indcol_it)){
          t1 <- t1-2*b*coefficients(regcc)[ubis+3]%*%CovarianceMat[indcolNA[j],indcol_it[ubis]]
          t2 <- t2-2*c*coefficients(regcc)[ubis+3]%*%CovarianceMat[indcolNA[jother],indcol_it[ubis]]
        }
      }
      if (length(indcol_it)==1){
        #Res <- (1/(2*b*c))*(-b^2*Variance[j]+var(YNA[,indselec])-Q-c^2*Variance[jother]-t(coefficients(regcc)[4:length(coefficients(regcc))])*var(YNA[,indcol_it])*coefficients(regcc)[4:length(coefficients(regcc))]+t1+t2)
        Res <- (1/(2*b*c))*(-b^2*Variance[j]+var(YNA[,indselec])-Q-c^2*Variance[jother]-t(coefficients(regcc)[4:length(coefficients(regcc))])*var(YNA[,indcol_it])*coefficients(regcc)[4:length(coefficients(regcc))]+t1+t2+sum(sigma^2*coefficients(regcc)[2:(r+1)]^2)) 
      }else if(length(indcol_it)>1){
        #Res <- (1/(2*b*c))*(-b^2*Variance[j]+var(YNA[,indselec])-Q-c^2*Variance[jother]-t(coefficients(regcc)[4:length(coefficients(regcc))])%*%var(YNA[,indcol_it])%*%coefficients(regcc)[4:length(coefficients(regcc))]+t1+t2)
        Res <- (1/(2*b*c))*(-b^2*Variance[j]+var(YNA[,indselec])-Q-c^2*Variance[jother]-t(coefficients(regcc)[4:length(coefficients(regcc))])%*%var(YNA[,indcol_it])%*%coefficients(regcc)[4:length(coefficients(regcc))]+t1+t2+sum(sigma^2*coefficients(regcc)[2:(r+1)]^2))
      }else{
        #Res <- (1/(2*b*c))*(-b^2*Variance[j]+var(YNA[,indselec])-Q-c^2*Variance[jother])
        Res <- (1/(2*b*c))*(-b^2*Variance[j]+var(YNA[,indselec])-Q-c^2*Variance[jother]+sum(sigma^2*coefficients(regcc)[2:(r+1)]^2))
      }
      cov <- c(cov,Res)
      coeff_cov <- c(coeff_cov,2*b*c)
    }
  }

  if (meth=="best"){
    res <-  cov[which.min(coeff_cov)]
  }else if(meth=="agg"){
    res <-  median(cov)
  }else{
    res <-  sample(cov,1)
  }
  
  return(res)
  
}


######
# Name: PPCA_MAREM_it
# Date: 25/06/2019
# Description: it computes one iterate of the EM algorithm designed for the PPCA model in presence of MAR data. 
# Arguments: 
#B_estimOld: estimated coefficient matrix. 
#B: true coefficient matrix. 
#moy_estim: estimated mean vector of the data matrix. 
#YNA: the data matrix containing missing values. 
#sigma: noise level (standart deviation). 
#####


PPCA_MAREM_it <- function(B_estimOld,B,moy_estim,YNA,sigma){
  
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
      sumMoy <- sumMoy + B_estimOld[,j]*(YNA[k,j]-moy_estim[j])
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
    moy_estim[j] <- 1/(length(rowObs))*sum_moyj
  }
  for (j in 1:p){
    MatB <- 0
    sumBMax <- 0
    rowObs <- which(is.na(YNA[,j])==0)
    for (k in rowObs){
      MatB <- MatB + Moy_W[,k]%*%t(Moy_W[,k]) + Sigma_W[[k]]
      sumBMax <- sumBMax + Moy_W[,k]*(YNA[k,j]-moy_estim[j])
    }
    B_estimNew[,j] <- ginv(MatB)%*%sumBMax
  }
  diff <- norm(B_estimNew - B_estimOld, type = "F") / (norm(B_estimOld, type = "F") + 10 ^ -3)
  
  cor <- coeffRV(t(B_estimNew),t(B))$rv
  
  return(list(B_estimNew=B_estimNew,moy_estim=moy_estim,diff=diff,cor=cor))
}
  
