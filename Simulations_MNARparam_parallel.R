#en parallèle
library(randomForest)
library(missForest)
library(doParallel)
library(doSNOW)
cl <- makeCluster(2, outfile="")
registerDoSNOW(cl)


#Parameters
Nbit=10
p=10
r=2
sigma=0.1
indMissVar=1:7
moy=0
set.seed(2002)
B <- matrix(rnorm(p * r), nrow = r, ncol = p)
modmecha<-"Logistic"
MatParam <- B
nbNA <- 350

mod="PPCA"

results.list = foreach (i = 1:Nbit, .combine = "rbind")  %dopar% {
    source("EM_tools.R",local=TRUE)
    source("General_tools.R",local=TRUE)

  set.seed(i)
  print("iteration")
  print(i)
  
  #Model
  ##Probabilistic PCA
  Noise <- matrix(rnorm(p * n, sd = sigma), nrow = n, ncol = p)
  if(mod=="PPCA"){
    W <- matrix(rnorm(n * r), nrow = n, ncol = r)
    Y <- rep(moy, n) + W %*% MatParam + Noise
  }else{
    Y=MatParam+Noise
  }
  
  ##Introduction of missing values
  YNA <- Y
  missingreg <- c()
  for (j in indMissVar) {
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
  
  #Paramétrique
  a=3 
  b=0
  Ns=500
  ThetaNew=Initialize_theta(ImputeMean0(YNA),r)
  aNew=a-1
  bNew=b-1
  M=1-is.na(YNA)
  
  diff=100
  ccompt<-0
  while(ccompt<80){
    Noise <- matrix(rnorm(p * n, sd = sigma), nrow = n, ncol = p)
    if(mod=="ProbaPCA"){
      ParamNew <- IterEM(W %*% t(MatParam),Y,YNA,ThetaNew,aNew,bNew,M,Ns,sigma^2,algo="soft",lam="Pred",nbcol=length(indMissVar))
    }else{
      ParamNew <- IterEM(MatParam,Y,YNA,ThetaNew,aNew,bNew,M,Ns,sigma^2,algo="soft",lam="Pred",nbcol=length(indMissVar))
    }
    diff=ParamNew$diff
    print(MSE(ThetaNew*(1-M),Y*(1-M)))
    ThetaNew=ParamNew$ThetaNew
    aNew=ParamNew$a_initNew
    bNew=ParamNew$b_initNew
    conv1=ParamNew$conv
    ccompt=ccompt+1
  }
  
  YParam <- YNA
  MeanParam <- apply(YNA,2,mean)
  for (j in 1:length(indMissVar)){
    missing <- which(is.na(YNA[,indMissVar[j]]))
    YParam[missing,indMissVar[j]] <- ThetaNew[missing,indMissVar[j]]
    MeanParam[indMissVar[j]] <- mean(YParam[,indMissVar[j]])
  }
  CovParam <- var(YParam)
  MSEParam <- MSE(YParam*(1-M),Y*(1-M)) 
  
  MeanParam <- list(MeanParam)
  CovParam <- list(CovParam)
  MSEParam <- list(MSEParam)
  
  cbind(MeanParam,CovParam,MSEParam)
    
}


