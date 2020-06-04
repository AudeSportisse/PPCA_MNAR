library(corpcor)
source("General_tools.R",local=TRUE)


prox_l1 <- function(x,alg,tau,ncp){ #computes the projection on the set M.*X = X_na
  if (alg=="hard"){res=(abs(x)>=sort(abs(x),decreasing=TRUE)[ncp])*x}else{res=max(0,1-tau/max(1e-15,abs(x)))*x}
  return(res)
}


soft_thresh_svd <- function(X,alg,gamma,ncp){ #performs soft thresholding in the value decomposition
  res=fast.svd(X)
  U=res$u
  V=res$v
  S=diag(res$d)
  if(alg=="hard"){
    s=prox_l1(diag(S),alg,gamma,ncp) 
  }else{
    s=sapply(diag(S),prox_l1,alg,gamma,ncp) 
  }
  if(alg=="hard"){
    if(ncp==0){
      X_res=as.matrix(ImputeMean(X))
    }else{
      S[1:length(s),1:length(s)] = diag(s)
      X_res=U%*%S%*%t(V) 
    }
  }else{
    S[1:length(s),1:length(s)] = diag(s)
    X_res=U%*%S%*%t(V) 
  }
  return(X_res)
}

FISTA0 <- function(X,param,bruit){ #algorithm FISTA without missing values
  X=as.matrix(X) 
  niter=200
  lambdaNew=0.1
  xNew=yNew=matrix(0,ncol=ncol(X),nrow=nrow(X))
  diff=1
  while (diff>10^(-6)){
    xOld=xNew
    yOld=yNew
    lambdaOld=lambdaNew
    xNew=soft_thresh_svd(as.matrix(yOld-(yOld-X)),alg="soft",param,ncp=NULL)
    lambdaNew=(1+sqrt(1+4*lambdaOld**2))/2
    yNew=xNew+((lambdaOld-1)/lambdaNew)*(xNew-xOld)
    diff=norm(xNew-xOld,type="F")
  }
  return(xNew)
}


FISTANA <- function(X,M,param,bruit,alg,ncp=NULL){ #algorithm FISTA with missing values
  
  #Initialization:
  X=as.matrix(X) 
  missing=which(is.na(X))
  X=ImputeMean0(X) #imputes the missing values by the mean
  
  lambdaNew=0.1
  xNew=yNew=matrix(0,ncol=ncol(X),nrow=nrow(X))
  diff=1
  while(diff>10^(-6)){
    xOld=xNew
    yOld=yNew
    lambdaOld=lambdaNew
    xNew=soft_thresh_svd(as.matrix(yOld-(M*(yOld-X))),alg,param,ncp)
    if(alg=="thresh"){
      thresh=max(xNew[missing]) 
      for( ind in missing){
        xNew[ind]=max(thresh,xNew[ind])
      }
    }
    lambdaNew=(1+sqrt(1+4*lambdaOld**2))/2
    yNew=xNew+((lambdaOld-1)/lambdaNew)*(xNew-xOld)
    diff=norm(xNew-xOld,type="F")
  }
  
  return(xNew)
}

FISTA <- function(X,M,param,a,b,bruit){ #Algorithms FISTA by modelling the mechanism with the intuitive idea
  
  #Initialization:
  X=as.matrix(X)
  X=ImputeMean0(X) #imputes the missing values by the mean
  
  lambdaNew=0.1
  xNew=yNew=matrix(0,ncol=ncol(X),nrow=nrow(X))
  diff=1
  while(diff>10^(-6)){
    xOld=xNew
    yOld=yNew
    lambdaOld=lambdaNew
    grad=matrix(0,nrow=nrow(X),ncol=ncol(X))
    for(i in 1:nrow(X)){
      for(j in 1:ncol(X)){
        grad[i,j]=(a*exp(a*yOld[i,j])/(exp(a*b)+exp(a*yOld[i,j])))
      }
    }
    xNew=soft_thresh_svd(as.matrix(yOld-(1/((1/bruit)+a^2/4))*((M*(1/bruit)*(yOld-X)+M*grad))),param)
    lambdaNew=(1+sqrt(1+4*lambdaOld**2))/2
    yNew=xNew+((lambdaOld-1)/lambdaNew)*(xNew-xOld)
    diff=norm(xNew-xOld,type="F")
  }
  return(xNew)
}

