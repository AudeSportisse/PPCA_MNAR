#Libraries:
library("FactoMineR")
library("softImpute")
library("MASS")
library("gtools")


source("Function_main_algebraic.R")
source("Functions_supp_algebraic.R")
source("OtherFunctions.R")
source("OtherFunctions_fixedeffect.R")



######
## Name: ComparMethods_PPCA_iteration_fixedeffect
## Description: #it returns a list containing the estimations for the mean, the variance and the covariances associated to each MNAR missing variable of the data matrix generated under a low rank model, for different methods. 
                #The missing values are introduced in the variables according to a self-masked MNAR mechanism using a logistic regression. 
## Arguments: 
# seed.num: to fix the random number generator. 
# n: number of observations. 
# p: number of variables. 
# r: rank of the low-rank matrix.
# Xlowrank: low-rank matrix of size n*p and rank r. 
# sigma: noise level (standart deviation). 
# indMissVar: indexes of the missing variables. 
# indRegVar: indexes of the pivot variables. If NULL, they are chosen as the complementary of the indexes of the missing variables. 
# param_logistic: parameters of the logistic regression, vector of 2 elements. (By default, it leads to about 35% missing values in total for n=1000, p=10, r=2 and 7 missing variables.)
###### 

ComparMethods_PPCA_iteration_fixedeffect <-
  function(seed_num,
           n,
           p,
           r,
           Xlowrank,
           sigma,
           indMissVar,
           indRegVar=NULL,
           param_logistic=c(3,0)
  ){
    
    ##Model: Probabilistic PCA
    set.seed(seed_num)
    Noise <- matrix(rnorm(p * n, sd = sigma), nrow = n, ncol = p)
    Y <- Xlowrank + Noise
   
    ##Introduction of missing values
    YNA <- Y
    missingreg <- c()
    a <- param_logistic[1]
    b <- param_logistic[2]
    for (j in indMissVar) {
      ##Logistic regression
      select_prob <-
        function(x, modmecha) {
          #probability of selecting coordinate Xij
          res = 1 / (1 + exp(-a * (x - b)))
          return(res)
        }
      prob <- sapply(Y[, j], select_prob, modmecha)
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
    
    res_PPCA_MNAR <- Results_imputation_fixedeffect(CovMNAR,MeanMNAR,YNA,Y,indMissVar,M,r,sigma)
    MSEMNAR <- res_PPCA_MNAR$MSE
    
    
    ## CC case
    Meancc <- apply(YNA,2,mean,na.rm = TRUE)
    Covcc <- matrix(1000,ncol=p,nrow=p)
    for (l in 1:p){
      Covcc[l,] <- apply(YNA,2,cov,YNA[,l],use="complete.obs")
    }
    
    ## SoftImpute
    
    res_estim_soft <- Mean_covariances_estimations_imputation_softImpute(YNA,Y,M,indMissVar)
    MeanSoft <- res_estim_soft$mean
    CovSoft <- res_estim_soft$cov
    YSoft <- res_estim_soft$Yimp
    
    res_PPCA_soft <- Results_fixedeffect(CovSoft,YSoft,Y,M,r,sigma)
    MSESoft <- res_PPCA_soft$MSE
    
    ## Mean Imputation
    
    YMean <- ImputeMean0(YNA) 
    MeanMean <- apply(YMean,2,mean)
    CovMean <- var(YMean)
    res_PPCA_Mean <- Results_fixedeffect(CovMean,YMean,Y,M,r,sigma)
    MSEMean <- res_PPCA_Mean$MSE
    
    ##General results

    Mean <- list(MNAR=MeanMNAR,CC=Meancc,Soft=MeanSoft,Mean=MeanMean)
    Cov <- list(MNAR=CovMNAR,CC=Covcc,Soft=CovSoft,Mean=CovMean)
    MSEres <-list(MNAR=MSEMNAR,Soft=MSESoft,Mean=MSEMean)
    
    result = list(Mean=Mean, Cov=Cov, MSEres=MSEres)
    return(result)
    
  }

