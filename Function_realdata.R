#Libraries:
library("FactoMineR")
library("softImpute")
library("MASS")
library("gtools")


source("Function_main_algebraic.R")
source("Functions_supp_algebraic.R")
source("OtherFunctions_realdata.R")


######
## Name: ComparMethods_PPCA_realdata
## Description: #it returns a list containing the prediction error relative to the mean imputation introducing additional MNAR values in one missing variables and considering the other variables as MCAR for different method.
                #The missing values are introduced in the variables according to a self-masked MNAR mechanism using a logistic regression. 
## Arguments: 
# seed.num: to fix the random number generator. 
# YNA: incomplete data matrix of size n*p
# YNAimp: complete data matrix of size n*p imputed in advance which serves to estimate the mean and covariances using the empirical quantities (in our method). 
# r: number of latent variables. 
# sigma: noise level (standart deviation). 
# indMissVar: indexes of the missing variables. 
# indRegVar: indexes of the pivot variables. If NULL, they are chosen as the complementary of the indexes of the missing variables. 
# param_logistic: parameters of the logistic regression, vector of 2 elements. (By default, it leads to about 35% missing values in total for n=1000, p=10, r=2 and 7 missing variables.)
# scale: if scale=TRUE, the variables are scaled to give the same weight to each variable before applying the algorithm. By default, scale=FALSE.
###### 


ComparMethods_PPCA_realdata <- function(seed_num,YNA,YNAimp,r,sigma,indMissVar,indRegVar,param_logistic,scale=FALSE){
  
  set.seed(seed_num)
  
  if(scale==TRUE){
    meanX <- apply(YNA, 2,mean,na.rm=TRUE)
    YNA <- t(t(YNA) - meanX)
    etX <- apply(YNA, 2,sd,na.rm=TRUE)
    YNA <- t(t(YNA)/etX)
  }
  
  ## Introduction of additionnal MNAR values in YNA_addmiss. 
  YNA_addmiss=YNA
  a <- param_logistic[1]
  b <- param_logistic[2]
  M2=matrix(1,nrow=nrow(YNA),ncol=ncol(YNA))
  for (j in indMissVar){
    select_prob <- function(x){
      res=(1/(1+exp(-a*(x-b))))
      return(res)
    }
    
    prob <- sapply(YNA[,j],select_prob)
    missing=c()
    for (i in 1:nrow(YNA)){
      if(!is.na(YNA[i,j])){
        u<-runif(1)
        if(prob[i]>u){missing=c(missing,i)}
      }
    }
    YNA_addmiss[missing,j]=NA
    M2[missing,j]=0
  }
  
  res_estim_MNAR <- Mean_covariances_estimations_algebraic(YNA_addmiss,indMissVar,indRegVar,r,YNAimp=YNAimp,opt_data="MCAR")
  Mean_estimated <- res_estim_MNAR$mean
  Covmatrix_estimated <- res_estim_MNAR$cov
  
  YNA_MNAR <- Results_PPCA_imputation_realdata(Covmatrix_estimated,Mean_estimated,YNA_addmiss,indMissVar,r,sigma)
  
  ## MAR EM

  res_MAREM <- Estimations_PPCA_MAREM_realdata(YNA_addmiss,indMissVar,r,sigma)
  Y_MAREM <- res_MAREM$Y_estim
  
  if(scale==TRUE){
    YNA_MNAR <- t(t(YNA_MNAR) * etX )
    YNA_MNAR <- t(t(YNA_MNAR) + meanX )
    YNA_addmiss <- t(t(YNA_addmiss) * etX )
    YNA_addmiss <- t(t(YNA_addmiss) + meanX )
    YNA <- t(t(YNA) * etX )
    YNA <- t(t(YNA) + meanX )
    Y_MAREM <- t(t(Y_MAREM) * etX )
    Y_MAREM <- t(t(Y_MAREM) + meanX )
  }
  
  YNA_withoutNA=YNA
  for (i in 1:nrow(YNA)){
    for (j in 1:ncol(YNA)){
      if (is.na(YNA[i,j]))
        YNA_withoutNA[i,j]=1
    }
  }
  
  ## SoftMAR
  
  res_estim_soft <- Mean_covariances_estimations_imputation_softImpute_realdata(YNA_addmiss,M,indMissVar)
  Y_Soft <- res_estim_soft$Yimp
  
  
  ## Results
  Y_Mean <- ImputeMean0(YNA_addmiss) 
  MSEvsmean <- function(Ychap, Y, Y_Mean){ return(sum((Y-Ychap)^2)/sum((Y-Y_Mean[,indMissVar]*(1-M2[,indMissVar]))^2))}
  MSE <- list(MNAR=MSEvsmean(YNA_MNAR[,indMissVar]*(1-M2[,indMissVar]),YNA_withoutNA[,indMissVar]*(1-M2[,indMissVar]),Y_Mean),MAREM=MSEvsmean(Y_MAREM[,indMissVar]*(1-M2[,indMissVar]),YNA_withoutNA[,indMissVar]*(1-M2[,indMissVar]),Y_Mean),Soft=MSEvsmean(Y_Soft[,indMissVar]*(1-M2[,indMissVar]),YNA_withoutNA[,indMissVar]*(1-M2[,indMissVar]),Y_Mean))
  return(list(MSE=MSE))
}