library(softImpute)
library(randomForest)
library(missForest)
source("FISTA_tools.R", local = TRUE)
source("General_tools.R", local = TRUE)


#######
#Initialization of the EM algorithm
#######


Initialize_theta <- function(X, r) {
  res = svd(X)
  result = (res$u[, 1:r] * res$d[r]) %*% (t(res$v[, 1:r]))
  return(result)
}


#######
#Univariate and multivariate case 
#######

#The mechanism parameters are the same ones for all the variables.
#The missing variables are the first nbcol variables.


IterEM <- function(Xtrue,
                   X,
                   XNA,
                   Thet,
                   a,
                   b,
                   M,
                   Ns,
                   bruit,
                   algo,
                   lam,
                   nbcol) {
  a_initOld = a
  b_initOld = b
  ThetaOld = Thet
  
  SIR <- function(Theta, XNA, a, b, sigma, i, jind) {
    Nn = Ns * 1000
    sim = rnorm(Nn, mean = Theta[i, jind], sd = sigma)
    
    RapportRej <- function(x) {
      g <- function(x) {
        res = ((1 / (1 + exp(-a * (
          x - b
        )))) ^ (1 - M[i, jind])) * ((1 - (1 / (
          1 + exp(-a * (x - b))
        ))) ^ M[i, jind])
        return(res)
      }
      
      num = g(x)
      return(num)
    }
    poids = RapportRej(sim)
    res <- sample(sim, Ns, prob = poids, replace = TRUE)
    return(res)
  }
  
  sigma = sqrt(bruit)
  XNA_remp = XNA
  data_remp <- list()
  nbNA = 0
  for (i in 1:nrow(XNA)) {
    nbNAcol = 0
    for (j in 1:ncol(XNA)) {
      if (is.na(XNA[i, j])) {
        ech = SIR(ThetaOld, XNA, a_initOld, b_initOld, sigma, i, j)
        XNA_remp[i, j] = mean(ech)
        if (nbNAcol == 0) {
          data_remp[[i]] = ech
        } else{
          data_remp[[i]] = append(data_remp[[i]], ech)
        }
        nbNAcol = nbNAcol + 1
      }
    }
    nbNA = nbNA + nbNAcol
  }
  
  if (algo == "soft") {
    RES = NULL
    gridlambda1 = seq(0, lambda0(X) * 1.1, length = 300)
    for (i in 1:length(gridlambda1)) {
      fit1 = softImpute(
        as.matrix(XNA_remp),
        rank = min(dim(XNA_remp)) - 1,
        lambda = gridlambda1[i],
        maxit = 10000,
        type = "svd"
      )
      if (fit1$d[1] == 0) {
        ThetaNew = as.matrix(ImputeMean(XNA_remp))
      } else if (length(fit1$d) == 1) {
        ThetaNew = (fit1$u * fit1$d) %*% t(fit1$v)
      } else{
        ThetaNew = (fit1$u %*% diag(fit1$d)) %*% t(fit1$v)
      }
      if (lam == "Tot") {
        RES[i] = MSE(ThetaNew, Xtrue)
      } else{
        RES[i] = MSE(ThetaNew * (1 - M), X * (1 - M))
      }
    }
    fit1 = softImpute(
      as.matrix(XNA_remp),
      rank = min(dim(XNA_remp)) - 1,
      lambda = gridlambda1[which.min(RES)],
      maxit = 10000,
      type = "svd"
    )
    if (fit1$d[1] == 0) {
      ThetaNew = as.matrix(ImputeMean(XNA_remp))
    } else if (length(fit1$d) == 1) {
      ThetaNew = (fit1$u * fit1$d) %*% t(fit1$v)
    } else{
      ThetaNew = (fit1$u %*% diag(fit1$d)) %*% t(fit1$v)
    }
  } else{
    RES = NULL
    gridParam = seq(0, lambda0(X) * 1.1, length = 100)
    for (i in 1:length(gridParam)) {
      X.FISTA2 = FISTA0(XNA_remp, gridParam[i], bruit)
      if (lam == "Tot") {
        RES[i] = MSE(X.FISTA2, Xtrue)
      } else{
        RES[i] = MSE(X.FISTA2 * (1 - M), X * (1 - M))
      }
    }
    ThetaNew = FISTA0(XNA_remp, gridParam[which.min(RES)], bruit)
  }
  
  M_concat = rep(as.vector(M[,1]), Ns )
  X_concat = matrix(0, nrow = Ns * n , ncol = 1)
  for (i in 1:nrow(XNA)) {
    for (k in (1:Ns)) {
      if (is.na(XNA[i, 1]) == TRUE) {
        X_concat[n * (k - 1) + i, 1] = data_remp[[i]][k]
      } else{
        X_concat[n * (k - 1) + i, 1] = XNA[i, 1]
      }
    }
  }
  
  GLM = glm(M_concat ~ X_concat, family = binomial)
  a_initNew = -GLM$coeff[2]
  b_initNew = -GLM$coeff[1] / GLM$coeff[2]
  conv = GLM$converged
  
  diff = norm(ThetaNew - ThetaOld, type = "F") / (norm(ThetaOld, type =
                                                         "F") + 10 ^ -3)
  
  return(
    list(
      diff = diff,
      ThetaNew = ThetaNew,
      a_initNew = a_initNew,
      b_initNew = b_initNew,
      conv = conv
    )
  )
  
}

