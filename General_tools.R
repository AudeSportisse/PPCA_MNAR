library(FactoMineR)

#MSE
MSE2 <- function(X1, X2){ return(sum((X1 - X2)^2))}
MSE <- function(X1, X2){ return(norm(as.matrix(X1 - X2),type="F")/norm(as.matrix(X2),type="F"))}

#SNR:
SNR <- function(X1, X2){ return(-20*log10(norm(as.matrix(X1 - X2), type = "F")/norm(as.matrix(X2), type = "F")))}

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

#Data matrix and parameter matrix simulations
simu <- function(n=100,p=10,S=4,sig=0.25){
  Xtilde= matrix(rnorm(n*p), nrow=n, ncol=p)
  decomp = svd(Xtilde)
  d = diag(decomp$d[1:S])
  Xw = decomp$u[,1:S]%*%d%*%t(decomp$v[,1:S])
  Xs=Xw
  don = Xs  + matrix(rnorm(n*p, 0, sig), n, p)
  return(list(don=don,mu=Xs))
}


