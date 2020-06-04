source("Function_simulations_fixedeffect.R")
library(ggplot2)
library(gridExtra)


simu_lowrank <- function(n,p,S,sig){
  Xtilde= matrix(rnorm(n*p), nrow=n, ncol=p)
  decomp = svd(Xtilde)
  d = diag(decomp$d[1:S])
  Xw = decomp$u[,1:S]%*%d%*%t(decomp$v[,1:S])
  Xs=Xw
  don = Xs  + matrix(rnorm(n*p, 0, sig), n, p)
  return(list(don=don,mu=Xs))
}


set.seed(2002)
Nbit <- 40
n <- 200
p <- 10
r <- 2
sigma <- 0.1
indMissVar <- c(1,2,3,4,5,9,10)
Xlowrank <- simu_lowrank(n,p,r,sigma)


result_fixedeffect <- lapply(1:Nbit,ComparMethods_PPCA_iteration_fixedeffect,
                                     n,
                                     p,
                                     r,
                                     Xlowrank,
                                     sigma,
                                     indMissVar)

result <- result_fixedeffect

##Plot 

#Mean, Variance
j <- 1 #the missing variable for which the mean, variance and covariances is ploted.
Meanglob <- c()
Varglob <- c()
for (k in 1:Nbit){
  Meanglob <- c(Meanglob,result[[k]]$Mean$MNAR[j],result[[k]]$Mean$CC[j],result[[k]]$Mean$Soft[j],result[[k]]$Mean$Mean[j])
  Varglob <- c(Varglob,result[[k]]$Cov$MNAR[j,j],result[[k]]$Cov$CC[j,j],result[[k]]$Cov$Soft[j,j],result[[k]]$Cov$Mean[j,j])
}
data_plot_Mean <- data.frame(result=Meanglob,meth=rep(c("2.MNAR","8.CC","6.Soft","7.Mean"),Nbit))
data_plot_Var <- data.frame(result=Varglob,meth=rep(c("2.MNAR","8.CC","6.Soft","7.Mean"),Nbit))

#use ylim(c(...,...)) if outliers. 
plot1 <- ggplot(data=data_plot_Mean,aes(x=meth,y=result))+geom_boxplot()+ylab("")+xlab("")+scale_x_discrete(labels=c("MNAR","SoftMAR","Mean","Del")) +theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x=element_text(size=16,angle=60),axis.text.y=element_text(size=16),title=element_text(hjust = 0.5,size=16)) + theme(legend.position='none')
plot2 <- ggplot(data=data_plot_Var,aes(x=meth,y=result))+geom_boxplot()+ylab("")+xlab("")+scale_x_discrete(labels=c("MNAR","SoftMAR","Mean","Del")) +theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x=element_text(size=16,angle=60),axis.text.y=element_text(size=16),title=element_text(hjust = 0.5,size=16)) + theme(legend.position='none')

grid.arrange(plot1,plot2,ncol=2) 


#MSE
MSEglob <- c()
for (k in 1:Nbit){
  MSEglob <- c(MSEglob,result[[k]]$MSE$MNAR,result[[k]]$MSE$Soft,result[[k]]$MSE$Mean)
}
data_plot_MSE <- data.frame(result=MSEglob,meth=rep(c("2.MNAR","6.Soft","7.Mean"),Nbit))
#use ylim(c(...,...)) if outliers. 
plot3 <- ggplot(data=data_plot_MSE,aes(x=meth,y=result))+geom_boxplot()+ylab("")+xlab("")+scale_x_discrete(labels=c("MNAR","SoftMAR","Mean"))+theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x=element_text(size=16,angle=60),axis.text.y=element_text(size=16),title=element_text(size=16)) + theme(legend.position='none')
grid.arrange(plot1,plot2,plot3,ncol=3)

