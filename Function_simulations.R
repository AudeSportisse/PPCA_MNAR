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
## Name: ComparMethods_PPCA_iteration
## Description: #it returns a list containing the estimations for the mean, the variance and the covariances associated to each MNAR missing variable of the data matrix generated under the PPCA model, for different methods. 
               #The missing values are introduced in the variables according to a self-masked MNAR mechanism using a logistic regression. 
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
# param_logistic: parameters of the logistic regression, vector of 2 elements. (By default, it leads to about 35% missing values in total for n=1000, p=10, r=2 and 7 missing variables.)
###### 

ComparMethods_PPCA_iteration <-
  function(seed_num,
           n,
           p,
           r,
           B,
           sigma,
           mean_theo,
           indMissVar,
           indRegVar=NULL,
           r_theo=NULL,
           param_logistic=c(3,0)
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
    a <- param_logistic[1]
    b <- param_logistic[2]
    for (j in indMissVar) {
      #Logistic regression
      select_prob <-
        function(x, modmecha) {
          #probability of selecting coordinate Xij
          res = 1 / (1 + exp(-a * (x - b)))
          return(res)
        }
      prob <- sapply(Y[, j], select_prob, modmecha)
      
      ## Example: MNAR general
      # a <- 3
      # b <- 0
      # c <- 2
      # select_prob <-
      #   function(x) {
      #     #probability of selecting coordinate Xij
      #     #res = pnorm(x)
      #     x1 = x[1]
      #     x2 = x[2]
      #     x3 = x[3]
      #     res = 1 / (1 + exp(-(a * x1) - (b* x2) - (c* x3)))
      #     return(res)
      #   }
      # Two_othervariables <- sample(setdiff(indMissVar,j),2,replace=FALSE)
      # prob <- apply(Y[, c(j,Two_othervariables)],1, select_prob)
      
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
    
    ## MNAR graphical
    
    #res_estim_MNAR_graph <- Mean_covariances_estimations_graphical(YNA,indMissVar,indRegVarAll,opt="agg")
    #MeanMNAR_graph <- res_estim_MNAR_graph$mean
    #CovMNAR_graph <- res_estim_MNAR_graph$cov
    
    #res_PPCA_MNAR_graph <- Results_PPCA_imputation(CovMNAR_graph,MeanMNAR_graph,YNA,Y,B,indMissVar,M,r,sigma)
    #CorrelationMNAR_graph <- res_PPCA_MNAR_graph$Corr
    #MSEMNAR_graph <- res_PPCA_MNAR_graph$MSE
    
    ## MNAR algebraic
    
    res_estim_MNAR <- Mean_covariances_estimations_algebraic(YNA,indMissVar,indRegVar,r,opt="agg")
    MeanMNAR <- res_estim_MNAR$mean
    CovMNAR <- res_estim_MNAR$cov
    
    res_PPCA_MNAR <- Results_PPCA_imputation(CovMNAR,MeanMNAR,YNA,Y,B,indMissVar,M,r,sigma)
    CorrelationMNAR <- res_PPCA_MNAR$Corr
    MSEMNAR <- res_PPCA_MNAR$MSE
    
    ## MNAR graphical no aggregation
    
    #res_estim_MNAR_graph_noagg <- Mean_covariances_estimations_graphical(YNA,indMissVar,indRegVar,opt="noagg")
    #MeanMNAR_graph_noagg <- res_estim_MNAR_graph_noagg$mean
    #CovMNAR_graph_noagg <- res_estim_MNAR_graph_noagg$cov
    
    #res_PPCA_MNAR_graph_noagg <- Results_PPCA_imputation(CovMNAR_graph_noagg,MeanMNAR_graph_noagg,YNA,Y,B,indMissVar,M,r,sigma)
    #CorrelationMNAR_graph_noagg <- res_PPCA_MNAR_graph_noagg$Corr
    #MSEMNAR_graph_noagg <- res_PPCA_MNAR_graph_noagg$MSE
    
    ## MNAR algebraic no aggregation
    
    #res_estim_MNAR_noagg <- Mean_covariances_estimations_algebraic(YNA,indMissVar,indRegVar,r,sigma,opt="noagg")
    #MeanMNAR_noagg <- res_estim_MNAR_noagg$mean
    #CovMNAR_noagg <- res_estim_MNAR_noagg$cov
    
    #res_PPCA_MNAR_noagg <- Results_PPCA_imputation(CovMNAR_noagg,MeanMNAR_noagg,YNA,Y,B,indMissVar,M,r,sigma)
    #CorrelationMNAR_noagg <- res_PPCA_MNAR_noagg$Corr
    #MSEMNAR_noagg <- res_PPCA_MNAR_noagg$MSE
    
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
    
    #Mean <- list(MNARgraph=MeanMNAR_graph,MNAR=MeanMNAR,MNARgraph_noagg=MeanMNAR_graph_noagg,MNARnoagg=MeanMNAR_noagg,CC=Meancc,MAR=MeanMAR,Soft=MeanSoft,Mean=MeanMean,MAREM=MeanMAREM)
    #Cov <- list(MNARgraph=CovMNAR_graph,MNAR=CovMNAR,MNARgraph_noagg=CovMNAR_graph_noagg,MNARnoagg=CovMNAR_noagg,CC=Covcc,MAR=CovMAR,Soft=CovSoft,Mean=CovMean,MAREM=CovMAREM)
    #Correlation <- list(MNARgraph=CorrelationMNAR_graph,MNAR=CorrelationMNAR,MNARgraph_noagg=CorrelationMNAR_graph_noagg,MNARnoagg=CorrelationMNAR_noagg,MAR=CorrelationMAR,Soft=CorrelationSoft,Mean=CorrelationMean,MAREM=CorrelationMAREM)
    #MSEres <-list(MNARgraph=MSEMNAR_graph,MNAR=MSEMNAR,MNARgraph_noagg=MSEMNAR_graph_noagg,MNARnoagg=MSEMNAR_noagg,MAR=MSEMAR,Soft=MSESoft,Mean=MSEMean,MAREM=MSEMAREM)
    
    Mean <- list(MNAR=MeanMNAR,CC=Meancc,MAR=MeanMAR,Soft=MeanSoft,Mean=MeanMean,MAREM=MeanMAREM)
    Cov <- list(MNAR=CovMNAR,CC=Covcc,MAR=CovMAR,Soft=CovSoft,Mean=CovMean,MAREM=CovMAREM)
    Correlation <- list(MNAR=CorrelationMNAR,MAR=CorrelationMAR,Soft=CorrelationSoft,Mean=CorrelationMean,MAREM=CorrelationMAREM)
    MSEres <-list(MNAR=MSEMNAR,MAR=MSEMAR,Soft=MSESoft,Mean=MSEMean,MAREM=MSEMAREM)
    
    
    result = list(Mean=Mean, Cov=Cov, Correlation=Correlation, MSEres=MSEres)
    return(result)
    
  }

