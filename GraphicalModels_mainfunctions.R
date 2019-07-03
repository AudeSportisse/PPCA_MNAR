#Libraries:
library("FactoMineR")
library("softImpute")
library("MASS")
library("gtools")

######
# Name: ComparMethods_PPCA_iteration
# Date: 25/06/2019
# Description: it gives the estimations for the mean, the variance and the covariances associated to each MNAR missing variable of the data matrix generated under the PPCA model using the coefficient matrix and the noise level. 
# Arguments: 
#seed.num: to fix the random number generator. 
#n: number of observations. 
#p: number of variables. 
#r: number of latent variables.
#B: coefficient matrix. 
#mean_theo: vector for the theoretical mean of the data matrix. 
#sigma: noise level (standart deviation). 
#indcolNA: indexes of the missing variables. 
#nbNA: number of missing values (if modmecha=="Censor")
#modmecha: model for the MNAR mechanism: "Logistic" or "Censor". 
#best: 1 to return the results for the graphical and algebraic MNAR methods by choosing the best result of the combinations of observable variables for the regressions, 0 (default) otherwise.
#agg: 1 (default) to return the results for the graphical and algebraic MNAR methods by choosing the aggregation of the combinations of observable variables for the regressions, 0 otherwise.
#noagg: 1 to return the results for the graphical and algebraic MNAR methods by randomly choosing a combination of observable variables for the regressions, 0 (default) otherwise.
#shrink: 1 to return the results for the algebraic MNAR method by using a shrinkage for the matrix inversions, 0 (default) otherwise.
#simplify: FALSE if all the the combinations of observable variables shoud be considered for the graphical and algebraic MNAR methods, TRUE (default) if only a part of them is considered. 
#####


ComparMethods_PPCA_iteration <-
  function(seed_num,
           n,
           p,
           r,
           B,
           mean_theo,
           sigma,
           indcolNA,
           nbNA,
           modmecha,
           best = 0,
           agg = 1,
           noagg = 0,
           shrink = 0,
           simplify = TRUE
           ){
    
    ##Model: Probabilistic PCA
    set.seed(seed_num)
    W <- matrix(rnorm(n * r), nrow = n, ncol = r)
    Noise <- matrix(rnorm(p * n, sd = sigma), nrow = n, ncol = p)
    Y <- rep(mean_theo, n) + W %*% B + Noise
    
    #Covariance matrix (theoretical)
    CovTheo <- t(B) %*% B + sigma ^ 2 * diag(1, ncol = p, nrow = p)
    
    
    ##Introduction of missing values
    YNA <- Y
    missingreg <- c()
    for (j in indcolNA) {
      if (modmecha == "Logistic") {
        ##Logistic regression
        a <- 3
        b <- 0
        select_prob <-
          function(x, modmecha) {
            #probability of selecting coordinate Xij
            res = 1 / (1 + exp(-a * (x - b)))
            return(res)
          }
        prob <- sapply(Y[, j], select_prob, modmecha)
        compt = 0
        missing = c()
        for (k in 1:n) {
          u <- runif(1)
          compt = compt + (prob[k] > u)
          if (prob[k] > u) {
            missing = c(missing, k)
          }
        }
        YNA[missing, j] <- NA
      } else{
        ##Censorship
        YNA[YNA[, j] > sort(YNA[, j], decreasing = TRUE)[nbNA + 1], j] <-
          NA
      }
    }
    M = 1 - is.na(YNA)
    
    #Initialisation
    ##MNAR
    MeanMNAR <- numeric(length(indcolNA))
    MeanMNAR_agg <- numeric(length(indcolNA))
    MeanMNAR_best <- numeric(length(indcolNA))
    
    #Variance graphical
    VarMNARg <- numeric(length(indcolNA))
    VarMNARg_agg <- numeric(length(indcolNA))
    VarMNARg_best <- numeric(length(indcolNA))
    
    #Variance algebraic
    VarMNAR <- numeric(length(indcolNA))
    VarMNAR_agg <- numeric(length(indcolNA))
    VarMNAR_best <- numeric(length(indcolNA))
    VarMNAR_shrink <- numeric(length(indcolNA))
    VarMNAR_agg_shrink <- numeric(length(indcolNA))
    VarMNAR_best_shrink <- numeric(length(indcolNA))
    
    ##MAR
    MeanMAR <- numeric(length(indcolNA))
    VarMAR <- numeric(length(indcolNA))
    
    ##CC
    Meancc <- numeric(length(indcolNA))
    Varcc <- numeric(length(indcolNA))
    
    ##Soft
    MeanSoft <- numeric(length(indcolNA))
    VarSoft <- numeric(length(indcolNA))
    
    ##Mean
    MeanMean <- numeric(length(indcolNA))
    VarMean <- numeric(length(indcolNA))
    
    ##MAR EM
    MeanMAREM <- numeric(length(indcolNA))
    VarMAREM <- numeric(length(indcolNA))
    
    ##Estimation MNAR

    #Computation of the covariances and variances related to the fully observed variables (with empirical quantities).
    CovMNARg <- var(YNA)
    CovMNARg_agg <- var(YNA)
    CovMNARg_best <- var(YNA)
    CovMNAR <- var(YNA)
    CovMNAR_shrink <- var(YNA)
    CovMNAR_agg <- var(YNA)
    CovMNAR_agg_shrink <- var(YNA)
    CovMNAR_best <- var(YNA)
    CovMNAR_best_shrink <- var(YNA)
    
    #Observed variables which are used for the regressions. 
    indcolObs <- 1:p
    for (jbis in indcolNA){
      indcolObs <- indcolObs[indcolObs!=jbis]
    }
    if (simplify == TRUE){
      indcolReg <- sample(indcolObs,r+2,replace=FALSE)
    }else{
      indcolReg <- indcolObs
    }
  
    #MNAR methods
    
    for (j in 1:length(indcolNA)){
      indvalNA <- which(is.na(YNA[,indcolNA[j]])) #indexes for missing values in the considered column. 
      
      YNAcc <- YNA[-indvalNA,] #the matrix containing rows where the considered column is fully observed
      
      ##Graphical method
      
      #Mean and variance estimation
      
      resMean <- Mean_graphical_algebraic(indcolNA,indcolReg,j,YNA,YNAcc,r)
      resVar <- Variance_graphical(indcolNA,indcolReg,j,YNA,YNAcc,r)
      
      if (best==1){
        MeanMNAR_best[j] <- resMean$Mean[which.max(abs(resMean$coeff_mean))]
        VarMNARg_best[j] <- resVar$Var[which.max(resVar$coeff_var)]
      }
      
      if(agg==1){
        MeanMNAR_agg[j] <- median(resMean$Mean)
        VarMNARg_agg[j] <- median(resVar$Var)
      }
      
      if(noagg==1){
        MeanMNAR[j] <- sample(resMean$Mean,1)
        VarMNARg[j] <- sample(resVar$Var,1)
      }
      
      #Covariance estimations between the missing variable and the observable variables
      
      for (jother in indcolObs){
        
        resCov <- Covariance_obsmiss_graphical(indcolNA,indcolReg,j,jother,YNA,YNAcc,r)
        
        if(best==1){
          CovMNARg_best[indcolNA[j],indcolNA[j]]=VarMNARg_best[j]
          CovMNARg_best[jother,indcolNA[j]] <- resCov$cov[which.max(abs(resCov$coeff_cov))]
          CovMNARg_best[indcolNA[j],jother] <- CovMNARg_best[jother,indcolNA[j]]
        }
        if(agg==1){
          CovMNARg_agg[indcolNA[j],indcolNA[j]]=VarMNARg_agg[j]
          CovMNARg_agg[jother,indcolNA[j]] <- median(resCov$cov)
          CovMNARg_agg[indcolNA[j],jother] <- CovMNARg_agg[jother,indcolNA[j]]
        }
        if(noagg==1){
          CovMNARg[indcolNA[j],indcolNA[j]]=VarMNARg[j]
          CovMNARg[jother,indcolNA[j]] <- sample(resCov$cov,1)
          CovMNARg[indcolNA[j],jother] <- CovMNARg[jother,indcolNA[j]]
        }
      }
      
      #Covariance estimations between the missing variable and other missing variables
      
      if(best==1){
        for (jother in 1:(length(indcolNA)-1)){
          if(jother < j){
            CovMNARg_best[indcolNA[j],indcolNA[jother]] <- Covariance_2missing_graphical(VarMNARg_best,CovMNARg_best,"best",indcolNA,indcolReg,j,jother,YNA,YNAcc,r)
            CovMNARg_best[indcolNA[jother],indcolNA[j]] <- CovMNARg_best[indcolNA[j],indcolNA[jother]]
          }
        }
      }
      
      if(agg==1){
        for (jother in 1:(length(indcolNA)-1)){
          if(jother < j){
            CovMNARg_agg[indcolNA[j],indcolNA[jother]] <- Covariance_2missing_graphical(VarMNARg_agg,CovMNARg_agg,"agg",indcolNA,indcolReg,j,jother,YNA,YNAcc,r)
            CovMNARg_agg[indcolNA[jother],indcolNA[j]] <- CovMNARg_agg[indcolNA[j],indcolNA[jother]]
          }
        }
      }
      
      if(noagg==1){
        for (jother in 1:(length(indcolNA)-1)){
          if(jother < j){
            CovMNARg[indcolNA[j],indcolNA[jother]] <- Covariance_2missing_graphical(VarMNARg,CovMNARg,"noagg",indcolNA,indcolReg,j,jother,YNA,YNAcc,r)
            CovMNARg[indcolNA[jother],indcolNA[j]] <- CovMNARg[indcolNA[j],indcolNA[jother]]
          }
        }
      }
      
      ##Algebraic method
      
      #Variance estimation
      
      if (best==1){
        resVarAlg <- Variance_algebraic(MeanMNAR_best,"best",shrink,indcolNA,indcolReg,j,YNA,YNAcc,CovTheo,r)
        if(shrink==1){
          VarMNAR_best_shrink[j] <- resVarAlg[2]
          CovMNAR_best_shrink[indcolNA[j],indcolNA[j]] <- VarMNAR_best_shrink[j]
        }
        VarMNAR_best[j] <- resVarAlg[1]
        CovMNAR_best[indcolNA[j],indcolNA[j]] <-  VarMNAR_best[j]
      }
      
      if (agg==1){
        resVarAlg <- Variance_algebraic(MeanMNAR_agg,"agg",shrink,indcolNA,indcolReg,j,YNA,YNAcc,CovTheo,r)
        if(shrink==1){
          VarMNAR_agg_shrink[j] <- resVarAlg[2]
          CovMNAR_agg_shrink[indcolNA[j],indcolNA[j]] <- VarMNAR_agg_shrink[j]
        }
        VarMNAR_agg[j] <- resVarAlg[1]
        CovMNAR_agg[indcolNA[j],indcolNA[j]] <-  VarMNAR_agg[j]
      }
      
      if (noagg==1){
        resVarAlg <- Variance_algebraic(MeanMNAR,"best",shrink,indcolNA,indcolReg,j,YNA,YNAcc,CovTheo,r)
        if(shrink==1){
          VarMNAR_shrink[j] <- resVarAlg[2]
          CovMNAR_shrink[indcolNA[j],indcolNA[j]] <- VarMNAR_shrink[j]
        }
        VarMNAR[j] <- resVarAlg[1]
        CovMNAR[indcolNA[j],indcolNA[j]] <-  VarMNAR[j]
      }
      
      #Covariance estimations between the missing variable and the observable variables
      
      indcol_cov <- indcolObs
      while(length(indcol_cov)>0){
        
        if (length(indcol_cov)>1){indselec <- sample(indcol_cov,1)}else{indselec <- indcol_cov}
        indcol_cov <- indcol_cov[indcol_cov!=indselec]
        
        if (best==1){
          resCovAlg <- Covariance_obsmiss_algebraic(MeanMNAR_best,"best",shrink,indselec,indcolNA,indcolReg,j,YNA,YNAcc,CovTheo,r)
          if(shrink==1){
            CovMNAR_best_shrink[indselec,indcolNA[j]] <- resCovAlg[2]
            CovMNAR_best_shrink[indcolNA[j],indselec] <- CovMNAR_best_shrink[indselec,indcolNA[j]] 
          }
          CovMNAR_best[indselec,indcolNA[j]] <- resCovAlg[1]
          CovMNAR_best[indcolNA[j],indselec] <-  CovMNAR_best[indselec,indcolNA[j]]
        }
        
        if (agg==1){
          resCovAlg <- Covariance_obsmiss_algebraic(MeanMNAR_agg,"agg",shrink,indselec,indcolNA,indcolReg,j,YNA,YNAcc,CovTheo,r)
          if(shrink==1){
            CovMNAR_agg_shrink[indselec,indcolNA[j]] <- resCovAlg[2]
            CovMNAR_agg_shrink[indcolNA[j],indselec] <- CovMNAR_agg_shrink[indselec,indcolNA[j]] 
          }
          CovMNAR_agg[indselec,indcolNA[j]] <- resCovAlg[1]
          CovMNAR_agg[indcolNA[j],indselec] <-  CovMNAR_agg[indselec,indcolNA[j]]
        }
        
        if (noagg==1){
          resCovAlg <- Covariance_obsmiss_algebraic(MeanMNAR,"noagg",shrink,indselec,indcolNA,indcolReg,j,YNA,YNAcc,CovTheo,r)
          if(shrink==1){
            CovMNAR_shrink[indselec,indcolNA[j]] <- resCovAlg[2]
            CovMNAR_shrink[indcolNA[j],indselec] <- CovMNAR_shrink[indselec,indcolNA[j]] 
          }
          CovMNAR[indselec,indcolNA[j]] <- resCovAlg[1]
          CovMNAR[indcolNA[j],indselec] <-  CovMNAR[indselec,indcolNA[j]]
        }
        
      }
          
      #Covariance estimations between the missing variable and other missing variables
      
      
      if(best==1){
        for (jother in 1:(length(indcolNA)-1)){
          if(jother < j){
            resCovMissAlg <- Covariance_2missing_algebraic(VarMNAR_best,CovMNAR_best,"best",shrink,indcolNA,indcolReg,j,jother,YNA,YNAcc,CovTheo,r)
            if(shrink==1){
              CovMNAR_best_shrink[indselec,indcolNA[j]] <- resCovMissAlg[2]
              CovMNAR_best_shrink[indcolNA[j],indselec] <- CovMNAR_shrink[indselec,indcolNA[j]] 
            }
            CovMNAR_best[indcolNA[j],indcolNA[jother]] <- resCovMissAlg[1]
            CovMNAR_best[indcolNA[jother],indcolNA[j]] <-  CovMNAR_best[indcolNA[j],indcolNA[jother]]
          }
        }
      }
      
      if(agg==1){
        for (jother in 1:(length(indcolNA)-1)){
          if(jother < j){
            resCovMissAlg <- Covariance_2missing_algebraic(VarMNAR_agg,CovMNAR_agg,"agg",shrink,indcolNA,indcolReg,j,jother,YNA,YNAcc,CovTheo,r)
            if(shrink==1){
              CovMNAR_agg_shrink[indselec,indcolNA[j]] <- resCovMissAlg[2]
              CovMNAR_agg_shrink[indcolNA[j],indselec] <- CovMNAR_shrink[indselec,indcolNA[j]] 
            }
            CovMNAR_agg[indcolNA[j],indcolNA[jother]] <- resCovMissAlg[1]
            CovMNAR_agg[indcolNA[jother],indcolNA[j]] <- CovMNAR_agg[indcolNA[j],indcolNA[jother]] 
          }
        }
      }
      
      
      if(noagg==1){
        for (jother in 1:(length(indcolNA)-1)){
          if(jother < j){
            resCovMissAlg <- Covariance_2missing_algebraic(VarMNAR,CovMNAR,"noagg",shrink,indcolNA,indcolReg,j,jother,YNA,YNAcc,CovTheo,r)
            if(shrink==1){
              CovMNAR_shrink[indcolNA[j],indcolNA[j]] <- resCovMissAlg[2]
              CovMNAR_shrink[indcolNA[j],indcolNA[j]] <- CovMNAR_shrink[indcolNA[j],indcolNA[j]] 
            }
            CovMNAR[indcolNA[j],indcolNA[jother]] <- resCovMissAlg[1]
            CovMNAR[indcolNA[jother],indcolNA[j]] <- CovMNAR[indcolNA[j],indcolNA[jother]] 
          }
        }
      }
    
    }
    
    
    #Results
    
    if(best==1){
      res_best <- Results_graph(CovMNARg_best,MeanMNAR_best,YNA,Y,B,indcolNA,M,r,p)
      res2_best <- Results_graph(CovMNAR_best,MeanMNAR_best,YNA,Y,B,indcolNA,M,r,p)
      CorrelationMNARg_best <- res_best[[1]]
      CorrelationMNAR_best <- res2_best[[1]]
      MSEMNARg_best <- res_best[[2]]
      MSEMNAR_best <- res2_best[[2]]
      if(shrink==1){
        res2_best_shrink <- Results_graph(CovMNAR_best_shrink,MeanMNAR_best,YNA,Y,B,indcolNA,M,r,p)
        CorrelationMNAR_best_shrink <- res2_best_shrink[[1]]
        MSEMNAR_best_shrink <- res2_best_shrink[[2]]
      }
    }
    
    if(agg==1){
      res_agg <- Results_graph(CovMNARg_agg,MeanMNAR_agg,YNA,Y,B,indcolNA,M,r,p)
      res2_agg <- Results_graph(CovMNAR_agg,MeanMNAR_agg,YNA,Y,B,indcolNA,M,r,p)
      CorrelationMNARg_agg <- res_agg[[1]]
      CorrelationMNAR_agg <- res2_agg[[1]]
      MSEMNARg_agg <- res_agg[[2]]
      MSEMNAR_agg <- res2_agg[[2]]
      if(shrink==1){
        res2_agg_shrink <- Results_graph(CovMNAR_agg_shrink,MeanMNAR_agg,YNA,Y,B,indcolNA,M,r,p)
        CorrelationMNAR_agg_shrink <- res2_agg_shrink[[1]]
        MSEMNAR_agg_shrink <- res2_agg_shrink[[2]] 
      }
    }
    
    if(noagg==1){
      res_noagg <- Results_graph(CovMNARg,MeanMNAR,YNA,Y,B,indcolNA,M,r,p)
      res2_noagg <- Results_graph(CovMNAR,MeanMNAR,YNA,Y,B,indcolNA,M,r,p)
      CorrelationMNARg <- res_noagg[[1]]
      CorrelationMNAR <- res2_noagg[[1]]
      MSEMNARg <- res_noagg[[2]]
      MSEMNAR <- res2_noagg[[2]]
      if(shrink==1){
        res2_noagg_shrink <- Results_graph(CovMNAR_shrink,MeanMNAR,YNA,Y,B,indcolNA,M,r,p)
        CorrelationMNAR_shrink <- res2_noagg_shrink[[1]]
        MSEMNAR_shrink <- res2_noagg_shrink[[2]] 
      }
    }
    
    #CC case
    
    Covcc <- var(YNA,na.rm = TRUE)
    for (j in 1:length(indcolNA)){
      indvalNA <- which(is.na(YNA[,indcolNA[j]]))
      YNAcctot <- YNA[-indvalNA,]
      Meancc[j] <- mean(YNAcctot[,indcolNA[j]])
      Varcc[j] <- var(YNAcctot[,indcolNA[j]])
    }
    
    ## MAR case
    CovMAR <- var(YNA)
    
    for (j in 1:length(indcolNA)){
      
      indvalNA <- which(is.na(YNA[,indcolNA[j]])) #indexes for missing values in the considered column. 
      
      YNAcc <- YNA[-indvalNA,] #the matrix containing rows where the considered column is fully observed
    
      indcolMAR <- sample(indcolReg,r,replace=FALSE)
      indselec <- sample(indcolMAR,1)
      indcolMAR <- indcolMAR[indcolMAR!=indselec]
      indcolMAR2 <- c(indcolMAR,indselec)
      regcc <- lm(YNAcc[, indcolNA[j]] ~  YNAcc[,indcolMAR2])
      sumMean <- c()
      for (ind in 1:length(indcolMAR2)){ sumMean <- c(sumMean,coefficients(regcc)[ind+1] * mean(YNA[,indcolMAR2[ind]]))}
      MeanMAR[j] <- coefficients(regcc)[1] + sum(sumMean)
      
      #Variance estimation 
      if(length(indcolMAR)==1){indselec2=indcolMAR}else{indselec2 <- sample(indcolMAR,1)}
      regccMAR <- lm(YNAcc[, indcolNA[j]] ~  YNAcc[, c(indselec,indselec2)])
      M21 <- c(cov(YNAcc[, indselec], YNAcc[, indcolNA[j]]), cov(YNAcc[, indselec2], YNAcc[, indcolNA[j]])) 
      MMAR <-solve(var(YNAcc[, c(indselec, indselec2)])) 
      add1 <- coefficients(regccMAR)[2] ^ 2 * var(YNA[, indselec]) + coefficients(regccMAR)[3] ^ 2 * var(YNA[, indselec2])
      add2 <- 2 * coefficients(regccMAR)[2] * coefficients(regccMAR)[3] * cov(YNA[, indselec], YNA[, indselec2])
      VarMAR[j] <- var(YNAcc[, indcolNA[j]]) - t(M21) %*% MMAR %*% M21 + add1 + add2  
      
      
      CovMAR[indcolNA[j],indcolNA[j]]=VarMAR[j]
      for (jother in indcolObs){
        indcolRegMAR <- indcolReg[indcolReg!=jother]
        if (length(indcolRegMAR)==1){indselec=indcolRegMAR}else{indselec <- sample(indcolRegMAR,1,replace=FALSE)}
        regccMAR <- lm(YNAcc[, indcolNA[j]] ~  YNAcc[, c(jother,indselec)])
        Covind <- coefficients(regccMAR)[2]*(var(YNA[,jother])+mean(YNA[,jother])^2)+coefficients(regccMAR)[1]*mean(YNA[,jother])+coefficients(regccMAR)[3]*(cov(YNA[,jother],YNA[,indselec])+mean(YNA[,jother])*mean(YNA[,indselec]))-mean(YNA[,jother])*MeanMAR[j] 
        CovMAR[jother,indcolNA[j]] <- Covind
        CovMAR[indcolNA[j],jother] <- CovMAR[jother,indcolNA[j]]
      } 
      for (jother in 1:(length(indcolNA)-1)){
        if(jother < j){
          if (length(indcolReg)==1){indselec=indcolReg}else{indselec <- sample(indcolReg,1,replace=FALSE)}
          regccMAR <- lm(YNAcc[, indcolNA[j]] ~  YNAcc[, c(indcolNA[jother],indselec)])
          Covind <- coefficients(regccMAR)[2]*(VarMAR[jother]+MeanMAR[jother]^2)+coefficients(regccMAR)[1]*MeanMAR[jother]+coefficients(regccMAR)[3]*(CovMAR[indselec,indcolNA[jother]]+MeanMAR[jother]*mean(YNA[,indselec]))-MeanMAR[jother]*MeanMAR[j]
          CovMAR[indcolNA[j],indcolNA[jother]] <- Covind
          CovMAR[indcolNA[jother],indcolNA[j]] <- CovMAR[indcolNA[j],indcolNA[jother]]
        } 
      } 
      
    }
    
    resMAR <- Results_graph(CovMAR,MeanMAR,YNA,Y,B,indcolNA,M,r,p)
    CorrelationMAR <- resMAR[[1]]
    MSEMAR <- resMAR[[2]]
    
    
    ## SoftImpute
    
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
    for (j in 1:length(indcolNA)){
      missing <- which(is.na(YNA[,indcolNA[j]]))
      YSoft <- X1.soft
      MeanSoft[j] <- mean(YSoft[, indcolNA[j]])
      VarSoft[j] <- var(YSoft[,indcolNA[j]])
    }
    CovSoft <- var(YSoft)
    ressvd <- svd(CovSoft - sigma ^ 2 * diag(1, ncol = p, nrow = p))
    BSoft <- sqrt(diag(ressvd$d)[1:r, 1:p]) %*% t(ressvd$u)
    CorrelationSoft <- coeffRV(t(BSoft),t(B))$rv
    MSESoft <- MSE(YSoft*(1-M),Y*(1-M))
    
    
    ## Mean Imputation
    
    YMean <- ImputeMean0(YNA)
    for (j in 1:length(indcolNA)){
      MeanMean[j] <- mean(YMean[, indcolNA[j]])
      VarMean[j] <- var(YMean[,indcolNA[j]])
    }
    CovMean <- var(YMean)
    ressvd <- svd(var(YMean) - sigma ^ 2 * diag(1, ncol = p, nrow = p))
    BMean <- sqrt(diag(ressvd$d)[1:r, 1:p]) %*% t(ressvd$u)
    CorrelationMean <- coeffRV(t(BMean),t(B))$rv
    MSEMean <- MSE(YMean*(1-M),Y*(1-M))
    
    ## PPCA EM
    
    B_estimNew <- B + matrix(rnorm(p*r,sd=1), nrow = r, ncol = p)
    moy_estim <- rnorm(p,sd=1)
    seuilEM <- 10^-4
    diff <- 100
    
    # Algo EM
    while (diff > seuilEM){
      resPPCA_MAREM_it <- PPCA_MAREM_it(B_estimNew,B,moy_estim,YNA,sigma)
      B_estimNew <- resPPCA_MAREM_it$B_estimNew
      moy_estim <- resPPCA_MAREM_it$moy_estim
      print(resPPCA_MAREM_it$cor)
      diff <- resPPCA_MAREM_it$diff
      #print(diff)
    }
    
    # Imputation 
    Y_MAREM <- YNA
    Cov <- t(B_estimNew)%*%B_estimNew + sigma ^ 2 * diag(1, ncol = p, nrow = p)
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
      remp3 <- t(YNA[, enscol] - remp3bis)
      remp <- moy_estim[j] + remp1 %*% remp2 %*% remp3
      missing <- which(is.na(YNA[,indcolNA[j]]))
      Y_MAREM[missing, indcolNA[j]] <- remp[missing]
    }
    MeanMAREM <- moy_estim[indcolNA]
    for (j in 1:length(indcolNA)){
      VarMAREM[j] <- var(Y_MAREM[,indcolNA[j]])
    }
    CovMAREM <- var(Y_MAREM)
    CorrelationMAREM <- coeffRV(t(B_estimNew),t(B))$rv
    MSEMAREM <- MSE(Y_MAREM*(1-M),Y*(1-M))
    
    
    ##General results
    
    if (best==1 & agg==1 & noagg==1){
      if (shrink==1){
        Mean <- c(MeanMNAR_best,MeanMNAR_agg,MeanMNAR,Meancc,MeanMAR,MeanSoft,MeanMean,MeanMAREM)
        Var <- c(VarMNARg_best,VarMNAR_best,VarMNAR_best_shrink,VarMNARg_agg,VarMNAR_agg,VarMNAR_agg_shrink,VarMNARg,VarMNAR,VarMNAR_shrink,Varcc,VarMAR,VarSoft,VarMean,VarMAREM)
        Cov <- list(CovMNARg_best,CovMNAR_best,CovMNAR_best_shrink,CovMNARg_agg,CovMNAR_agg,CovMNAR_agg_shrink,CovMNARg,CovMNAR,CovMNAR_shrink,Covcc,CovMAR,CovSoft,CovMean,CovMAREM)
        Correlation <- c(CorrelationMNARg_best,CorrelationMNAR_best,CorrelationMNAR_best_shrink,CorrelationMNARg_agg,CorrelationMNAR_agg,CorrelationMNAR_agg_shrink,CorrelationMNARg,VarMNAR,CorrelationMNAR_shrink,CorrelationMAR,CorrelationSoft,CorrelationMean,CorrelationMAREM)
        MSEres <-c(MSEMNARg_best,MSEMNAR_best,MSEMNAR_best_shrink,MSEMNARg_agg,MSEMNAR_agg,MSEMNAR_agg_shrink,MSEMNARg,MSEMNAR,MSEMNAR_shrink,MSEMAR,MSESoft,MSEMean,MSEMAREM)
      }else{
        Mean <- c(MeanMNAR_best,MeanMNAR_agg,MeanMNAR,Meancc,MeanMAR,MeanSoft,MeanMean,MeanMAREM)
        Var <- c(VarMNARg_best,VarMNAR_best,VarMNARg_agg,VarMNAR_agg,VarMNARg,VarMNAR,Varcc,VarMAR,VarSoft,VarMean,VarMAREM)
        Cov <- list(CovMNARg_best,CovMNAR_best,CovMNARg_agg,CovMNAR_agg,CovMNARg,CovMNAR,Covcc,CovMAR,CovSoft,CovMean,CovMAREM)
        Correlation <- c(CorrelationMNARg_best,CorrelationMNAR_best,CorrelationMNARg_agg,CorrelationMNAR_agg,CorrelationMNARg,VarMNAR,CorrelationMAR,CorrelationSoft,CorrelationMean,CorrelationMAREM)
        MSEres <-c(MSEMNARg_best,MSEMNAR_best,MSEMNARg_agg,MSEMNAR_agg,MSEMNARg,MSEMNAR,MSEMAR,MSESoft,MSEMean,MSEMAREM)
      }
    }
    
    if (best==1 & noagg==1 & agg==0){
      if (shrink==1){
        Mean <- c(MeanMNAR_best,MeanMNAR,Meancc,MeanMAR,MeanSoft,MeanMean,MeanMAREM)
        Var <- c(VarMNARg_best,VarMNAR_best,VarMNAR_best_shrink,VarMNARg,VarMNAR,VarMNAR_shrink,Varcc,VarMAR,VarSoft,VarMean,VarMAREM)
        Cov <- list(CovMNARg_best,CovMNAR_best,CovMNAR_best_shrink,CovMNARg,CovMNAR,CovMNAR_shrink,Covcc,CovMAR,CovSoft,CovMean,CovMAREM)
        Correlation <- c(CorrelationMNARg_best,CorrelationMNAR_best,CorrelationMNAR_best_shrink,CorrelationMNARg,CorrelationMNAR,CorrelationMNAR_shrink,CorrelationMAR,CorrelationSoft,CorrelationMean,CorrelationMAREM)
        MSEres <-c(MSEMNARg_best,MSEMNAR_best,MSEMNAR_best_shrink,MSEMNARg,MSEMNAR,MSEMNAR_shrink,MSEMAR,MSESoft,MSEMean,MSEMAREM)
      }else{
        Mean <- c(MeanMNAR_best,MeanMNAR,Meancc,MeanMAR,MeanSoft,MeanMean,MeanMAREM)
        Var <- c(VarMNARg_best,VarMNAR_best,VarMNARg,VarMNAR,Varcc,VarMAR,VarSoft,VarMean,VarMAREM)
        Cov <- list(CovMNARg_best,CovMNAR_best,CovMNARg,CovMNAR,Covcc,CovMAR,CovSoft,CovMean,CovMAREM)
        Correlation <- c(CorrelationMNARg_best,CorrelationMNAR_best,CorrelationMNARg,CorrelationMNAR,CorrelationMAR,CorrelationSoft,CorrelationMean,CorrelationMAREM)
        MSEres <-c(MSEMNARg_best,MSEMNAR_best,MSEMNARg,MSEMNAR,MSEMAR,MSESoft,MSEMean,MSEMAREM)
      }
    }
    
    if (best==1 & agg==1 & noagg==0){
      if (shrink==1){
        Mean <- c(MeanMNAR_best,MeanMNAR_agg,Meancc,MeanMAR,MeanSoft,MeanMean,MeanMAREM)
        Var <- c(VarMNARg_best,VarMNAR_best,VarMNAR_best_shrink,VarMNARg_agg,VarMNAR_agg,VarMNAR_agg_shrink,Varcc,VarMAR,VarSoft,VarMean,VarMAREM)
        Cov <- list(CovMNARg_best,CovMNAR_best,CovMNAR_best_shrink,CovMNARg_agg,CovMNAR_agg,CovMNAR_agg_shrink,Covcc,CovMAR,CovSoft,CovMean,CovMAREM)
        Correlation <- c(CorrelationMNARg_best,CorrelationMNAR_best,CorrelationMNAR_best_shrink,CorrelationMNARg_agg,CorrelationMNAR_agg,CorrelationMNAR_agg_shrink,CorrelationMAR,CorrelationSoft,CorrelationMean,CorrelationMAREM)
        MSEres <-c(MSEMNARg_best,MSEMNAR_best,MSEMNAR_best_shrink,MSEMNARg_agg,MSEMNAR_agg,MSEMNAR_agg_shrink,MSEMAR,MSESoft,MSEMean,MSEMAREM)
      }else{
        Mean <- c(MeanMNAR_best,MeanMNAR_agg,Meancc,MeanMAR,MeanSoft,MeanMean,MeanMAREM)
        Var <- c(VarMNARg_best,VarMNAR_best,VarMNARg_agg,VarMNAR_agg,Varcc,VarMAR,VarSoft,VarMean,VarMAREM)
        Cov <- list(CovMNARg_best,CovMNAR_best,CovMNARg_agg,CovMNAR_agg,Covcc,CovMAR,CovSoft,CovMean,CovMAREM)
        Correlation <- c(CorrelationMNARg_best,CorrelationMNAR_best,CorrelationMNARg_agg,CorrelationMNAR_agg,CorrelationMAR,CorrelationSoft,CorrelationMean,CorrelationMAREM)
        MSEres <-c(MSEMNARg_best,MSEMNAR_best,MSEMNARg_agg,MSEMNAR_agg,MSEMAR,MSESoft,MSEMean,MSEMAREM)
      }
    }
    
    if (noagg==1 & agg==1 & best==0){
      if (shrink==1){
        Mean <- c(MeanMNAR_agg,MeanMNAR,Meancc,MeanMAR,MeanSoft,MeanMean,MeanMAREM)
        Var <- c(VarMNARg_agg,VarMNAR_agg,VarMNAR_agg_shrink,VarMNARg,VarMNAR,VarMNAR_shrink,Varcc,VarMAR,VarSoft,VarMean,VarMAREM)
        Cov <- list(CovMNARg_agg,CovMNAR_agg,CovMNAR_agg_shrink,CovMNARg,CovMNAR,CovMNAR_shrink,Covcc,CovMAR,CovSoft,CovMean,CovMAREM)
        Correlation <- c(CorrelationMNARg_agg,CorrelationMNAR_agg,CorrelationMNAR_agg_shrink,CorrelationMNARg,CorrelationMNAR,CorrelationMNAR_shrink,CorrelationMAR,CorrelationSoft,CorrelationMean,CorrelationMAREM)
        MSEres <-c(MSEMNARg_agg,MSEMNAR_agg,MSEMNAR_agg_shrink,MSEMNARg,MSEMNAR,MSEMNAR_shrink,MSEMAR,MSESoft,MSEMean,MSEMAREM)
      }else{
        Mean <- c(MeanMNAR_agg,MeanMNAR,Meancc,MeanMAR,MeanSoft,MeanMean,MeanMAREM)
        Var <- c(VarMNARg_agg,VarMNAR_agg,VarMNARg,VarMNAR,Varcc,VarMAR,VarSoft,VarMean,VarMAREM)
        Cov <- list(CovMNARg_agg,CovMNAR_agg,CovMNARg,CovMNAR,Covcc,CovMAR,CovSoft,CovMean,CovMAREM)
        Correlation <- c(CorrelationMNARg_agg,CorrelationMNAR_agg,CorrelationMNARg,CorrelationMNAR,CorrelationMAR,CorrelationSoft,CorrelationMean,CorrelationMAREM)
        MSEres <-c(MSEMNARg_agg,MSEMNAR_agg,MSEMNARg,MSEMNAR,MSEMAR,MSESoft,MSEMean,MSEMAREM)
      }
    }
    
    if (noagg==1 & agg==0 & best==0){
      if (shrink==1){
        Mean <- c(MeanMNAR,Meancc,MeanMAR,MeanSoft,MeanMean,MeanMAREM)
        Var <- c(VarMNARg,VarMNAR,VarMNAR_shrink,Varcc,VarMAR,VarSoft,VarMean,VarMAREM)
        Cov <- list(CovMNARg,CovMNAR,CovMNAR_shrink,Covcc,CovMAR,CovSoft,CovMean,CovMAREM)
        Correlation <- c(CorrelationMNARg,CorrelationMNAR,CorrelationMNAR_shrink,CorrelationMAR,CorrelationSoft,CorrelationMean,CorrelationMAREM)
        MSEres <-c(MSEMNARg,MSEMNAR,MSEMNAR_shrink,MSEMAR,MSESoft,MSEMean,MSEMAREM)
      }else{
        Mean <- c(MeanMNAR,Meancc,MeanMAR,MeanSoft,MeanMean,MeanMAREM)
        Var <- c(VarMNARg,VarMNAR,Varcc,VarMAR,VarSoft,VarMean,VarMAREM)
        Cov <- list(CovMNARg,CovMNAR,Covcc,CovMAR,CovSoft,CovMean,CovMAREM)
        Correlation <- c(CorrelationMNARg,CorrelationMNAR,CorrelationMAR,CorrelationSoft,CorrelationMean,CorrelationMAREM)
        MSEres <-c(MSEMNARg,MSEMNAR,MSEMAR,MSESoft,MSEMean,MSEMAREM)
      }
    }
    
    if (agg==1 & noagg==0 & best==0){
      if (shrink==1){
        Mean <- c(MeanMNAR_agg,Meancc,MeanMAR,MeanSoft,MeanMean,MeanMAREM)
        Var <- c(VarMNARg_agg,VarMNAR_agg,VarMNAR_agg_shrink,Varcc,VarMAR,VarSoft,VarMean,VarMAREM)
        Cov <- list(CovMNARg_agg,CovMNAR_agg,CovMNAR_agg_shrink,Covcc,CovMAR,CovSoft,CovMean,CovMAREM)
        Correlation <- c(CorrelationMNARg_agg,CorrelationMNAR_agg,CorrelationMNAR_agg_shrink,CorrelationMAR,CorrelationSoft,CorrelationMean,CorrelationMAREM)
        MSEres <-c(MSEMNARg_agg,MSEMNAR_agg,MSEMNAR_agg_shrink,MSEMAR,MSESoft,MSEMean,MSEMAREM)
      }else{
        Mean <- c(MeanMNAR_agg,Meancc,MeanMAR,MeanSoft,MeanMean,MeanMAREM)
        Var <- c(VarMNARg_agg,VarMNAR_agg,Varcc,VarMAR,VarSoft,VarMean,VarMAREM)
        Cov <- list(CovMNARg_agg,CovMNAR_agg,Covcc,CovMAR,CovSoft,CovMean,CovMAREM)
        Correlation <- c(CorrelationMNARg_agg,CorrelationMNAR_agg,CorrelationMAR,CorrelationSoft,CorrelationMean,CorrelationMAREM)
        MSEres <-c(MSEMNARg_agg,MSEMNAR_agg,MSEMAR,MSESoft,MSEMean,MSEMAREM)
      }
    }
    
    if (best==1 & noagg==0 & agg==0){
      if (shrink==1){
        Mean <- c(MeanMNAR_best,Meancc,MeanMAR,MeanSoft,MeanMean,MeanMAREM)
        Var <- c(VarMNARg_best,VarMNAR_best,VarMNAR_best_shrink,Varcc,VarMAR,VarSoft,VarMean,VarMAREM)
        Cov <- list(CovMNARg_best,CovMNAR_best,CovMNAR_best_shrink,Covcc,CovMAR,CovSoft,CovMean,CovMAREM)
        Correlation <- c(CorrelationMNARg_best,CorrelationMNAR_best,CorrelationMNAR_best_shrink,CorrelationMAR,CorrelationSoft,CorrelationMean,CorrelationMAREM)
        MSEres <-c(MSEMNARg_best,MSEMNAR_best,MSEMNAR_best_shrink,MSEMAR,MSESoft,MSEMean,MSEMAREM)
      }else{
        Mean <- c(MeanMNAR_best,Meancc,MeanMAR,MeanSoft,MeanMean,MeanMAREM)
        Var <- c(VarMNARg_best,VarMNAR_best,Varcc,VarMAR,VarSoft,VarMean,VarMAREM)
        Cov <- list(CovMNARg_best,CovMNAR_best,Covcc,CovMAR,CovSoft,CovMean,CovMAREM)
        Correlation <- c(CorrelationMNARg_best,CorrelationMNAR_best,CorrelationMAR,CorrelationSoft,CorrelationMean,CorrelationMAREM)
        MSEres <-c(MSEMNARg_best,MSEMNAR_best,MSEMAR,MSESoft,MSEMean,MSEMAREM)
      }
    }
    
    result = list(Mean,Var, Cov, Correlation,MSEres)
    return(result)
     
  }

