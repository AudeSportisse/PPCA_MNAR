#Libraries:
library("FactoMineR")
library("softImpute")
library("MASS")
library("gtools")


source("Function_main_algebraic.R")
source("Function_main_graphical.R")
source("Functions_supp_algebraic.R")
source("Functions_supp_graphical.R")
source("OtherFunctions.R")


######
## Name: ComparMethods_PPCA_iteration_MNARgeneral
## Description: #it returns a list containing the estimations for the mean, the variance and the covariances associated to each MNAR missing variable of the data matrix generated under the PPCA model, for different methods. 
#The missing values are introduced in the variables according to a general MNAR mechanism using a logistic regression. 
## Arguments: 
# seed.num: to fix the random number generator. 
# n: number of observations. 
# p: number of variables. 
# r: number of latent variables. (If r is not well specified, please provide information to r_theo.)
# B: (theoretical) loading matrix of size r_theo*p.
# mean_theo: vector of size p for the theoretical mean of the data matrix. 
# sigma: noise level (standart deviation). 
# indMissVar: indexes of the missing variables. 
# indRegVar: indexes of the pivot variables. If NULL, they are chosen as the complementary of the indexes of the missing variables. 
# r_theo: rank (number of latent variables) with which B is generated. If NULL, r_theo=r. 
###### 

ComparMethods_PPCA_iteration_MNARgeneral <-
  function(seed_num,
           n,
           p,
           r,
           B,
           sigma,
           mean_theo,
           indMissVar,
           indRegVar=NULL,
           r_theo=NULL
  ){
    
    ##Model: Probabilistic PCA
    
    if(is.null(r_theo)){
      r_theo=r
    }
    
    set.seed(seed_num)
    W <- matrix(rnorm(n * r_theo), nrow = n, ncol = r_theo)
    Noise <- matrix(rnorm(p * n, sd = sigma), nrow = n, ncol = p)
    Y <- rep(mean_theo, n) + W %*% B  + Noise
    
    #Covariance matrix (theoretical)
    CovTheo <- t(B) %*% B + sigma ^ 2 * diag(1, ncol = p, nrow = p)
    
    ##Introduction of missing values
    YNA <- Y
    for (j in indMissVar) {
      ## Example: MNAR general
      a <- 3
      b <- 0
      c <- 2
      select_prob <-
        function(x) {
          #probability of selecting coordinate Xij
          #res = pnorm(x)
          x1 = x[1]
          x2 = x[2]
          x3 = x[3]
          res = 1 / (1 + exp(-(a * x1) - (b* x2) - (c* x3)))
          return(res)
        }
      Two_othervariables <- sample(setdiff(indMissVar,j),2,replace=FALSE)
      prob <- apply(Y[, c(j,Two_othervariables)],1, select_prob)
      
      missing = c()
      for (k in 1:n) {
        u <- runif(1)
        if (prob[k] > u) {
          missing = c(missing, k)
        }
      }
      YNA[missing, j] <- NA
    }
    M = 1 - is.na(YNA)
    
    indRegVarAll <- setdiff(1:p,indMissVar)
    
    if (is.null(indRegVar)){
      indRegVar=setdiff(1:p,indMissVar) 
    }
    
    ## MNAR algebraic
    
    res_estim_MNAR <- Mean_covariances_estimations_algebraic(YNA,indMissVar,indRegVar,r,opt="agg")
    MeanMNAR <- res_estim_MNAR$mean
    CovMNAR <- res_estim_MNAR$cov
    
    res_PPCA_MNAR <- Results_PPCA_imputation(CovMNAR,MeanMNAR,YNA,Y,B,indMissVar,M,r,sigma)
    CorrelationMNAR <- res_PPCA_MNAR$Corr
    MSEMNAR <- res_PPCA_MNAR$MSE
    
    ## CC case
    
    Meancc <- apply(YNA,2,mean,na.rm = TRUE)
    Covcc <- matrix(1000,ncol=p,nrow=p)
    for (l in 1:p){
      Covcc[l,] <- apply(YNA,2,cov,YNA[,l],use="complete.obs")
    }
    
    ## MAR case
    
    res_estim_MAR <- Mean_covariances_estimations_MAR(YNA,indMissVar,indRegVarAll,r)
    MeanMAR <- res_estim_MAR$mean
    CovMAR <- res_estim_MAR$cov
    
    res_PPCA_MAR <- Results_PPCA_imputation(CovMAR,MeanMAR,YNA,Y,B,indMissVar,M,r,sigma)
    CorrelationMAR <- res_PPCA_MAR$Corr
    MSEMAR <- res_PPCA_MAR$MSE
    
    ## SoftImpute
    
    res_estim_soft <- Mean_covariances_estimations_imputation_softImpute(YNA,Y,M,indMissVar)
    MeanSoft <- res_estim_soft$mean
    CovSoft <- res_estim_soft$cov
    YSoft <- res_estim_soft$Yimp
    
    res_PPCA_soft <- Results_PPCA(CovSoft,YSoft,Y,B,M,r,sigma)
    CorrelationSoft <- res_PPCA_soft$Corr
    MSESoft <- res_PPCA_soft$MSE
    
    ## Mean Imputation
    
    YMean <- ImputeMean0(YNA) 
    MeanMean <- apply(YMean,2,mean)
    CovMean <- var(YMean)
    res_PPCA_Mean <- Results_PPCA(CovMean,YMean,Y,B,M,r,sigma)
    CorrelationMean <- res_PPCA_Mean$Corr
    MSEMean <- res_PPCA_Mean$MSE
    
    ## MAR EM
    
    res_MAREM <- Estimations_PPCA_MAREM(YNA,Y,M,B,indMissVar,r,sigma)
    MeanMAREM <- res_MAREM$mean
    CovMAREM <- res_MAREM$cov
    CorrelationMAREM <- res_MAREM$Corr
    MSEMAREM <- res_MAREM$MSE
    
    
    ##General results

    Mean <- list(MNAR=MeanMNAR,CC=Meancc,MAR=MeanMAR,Soft=MeanSoft,Mean=MeanMean,MAREM=MeanMAREM)
    Cov <- list(MNAR=CovMNAR,CC=Covcc,MAR=CovMAR,Soft=CovSoft,Mean=CovMean,MAREM=CovMAREM)
    Correlation <- list(MNAR=CorrelationMNAR,MAR=CorrelationMAR,Soft=CorrelationSoft,Mean=CorrelationMean,MAREM=CorrelationMAREM)
    MSEres <-list(MNAR=MSEMNAR,MAR=MSEMAR,Soft=MSESoft,Mean=MSEMean,MAREM=MSEMAREM)
    
    
    result = list(Mean=Mean, Cov=Cov, Correlation=Correlation, MSEres=MSEres)
    return(result)
    
  }

