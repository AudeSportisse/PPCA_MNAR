source("Function_simulations.R")
library(ggplot2)
library(gridExtra)

## Toy example 

set.seed(2002)
Nbit <- 20
p <- 10
r <- 2
sigma <- 0.1
indMissVar <- 1:7
mean_theo <- 0
B <- matrix(rnorm(p * r), nrow = r, ncol = p)
n <- 1000


result_01 <- lapply(1:Nbit,ComparMethods_PPCA_iteration,
                      n=n,
                      p=p,
                      r=r,
                      B=B,
                      mean_theo=mean_theo,
                      sigma=sigma,
                      indMissVar=indMissVar
)


CovTheo <- t(B) %*% B + sigma ^ 2 * diag(1, ncol = p, nrow = p)

##Plot 

#Mean, Variance, Covariances
j <- 1 #the missing variable for which the mean, variance and covariances is ploted.
jmiss <- 2 #the missing variable for which the covariance is ploted in plot3.
jobs <- 8 #the observed variable for which the covariance is ploted in plot4.
Meanglob <- c()
Varglob <- c()
Covglob1 <- c()
Covglob2 <- c()
for (k in 1:Nbit){
  Meanglob <- c(Meanglob,result[[k]]$Mean$MNAR[j],result[[k]]$Mean$CC[j],result[[k]]$Mean$MAR[j],result[[k]]$Mean$Soft[j],result[[k]]$Mean$Mean[j],result[[k]]$Mean$MAREM[j])
  Varglob <- c(Varglob,result[[k]]$Cov$MNAR[j,j],result[[k]]$Cov$CC[j,j],result[[k]]$Cov$MAR[j,j],result[[k]]$Cov$Soft[j,j],result[[k]]$Cov$Mean[j,j],result[[k]]$Cov$MAREM[j,j])
  Covglob1 <- c(Covglob1,result[[k]]$Cov$MNAR[j,jmiss],result[[k]]$Cov$CC[j,jmiss],result[[k]]$Cov$MAR[j,jmiss],result[[k]]$Cov$Soft[j,jmiss],result[[k]]$Cov$Mean[j,jmiss],result[[k]]$Cov$MAREM[j,jmiss])
  Covglob2 <- c(Covglob2,result[[k]]$Cov$MNAR[j,jobs],result[[k]]$Cov$CC[j,jobs],result[[k]]$Cov$MAR[j,jobs],result[[k]]$Cov$Soft[j,jobs],result[[k]]$Cov$Mean[j,jobs],result[[k]]$Cov$MAREM[j,jobs])
}
data_plot_Mean <- data.frame(result=Meanglob,meth=rep(c("2.MNAR","8.CC","4.MAR","6.Soft","7.Mean","5.MAREM"),Nbit))
data_plot_Var <- data.frame(result=Varglob,meth=rep(c("2.MNAR","8.CC","4.MAR","6.Soft","7.Mean","5.MAREM"),Nbit))
data_plot_Cov1 <- data.frame(result=Covglob1,meth=rep(c("2.MNAR","8.CC","4.MAR","6.Soft","7.Mean","5.MAREM"),Nbit))
data_plot_Cov2 <- data.frame(result=Covglob2,meth=rep(c("2.MNAR","8.CC","4.MAR","6.Soft","7.Mean","5.MAREM"),Nbit))

#use ylim(c(...,...)) if outliers. 
plot1 <- ggplot(data=data_plot_Mean,aes(x=meth,y=result))+geom_boxplot()+geom_hline(yintercept = mean_theo,color = "red")+ylab("")+xlab("")+scale_x_discrete(labels=c("MNAR","MAR","EMMAR","SoftMAR","Mean","Del")) +theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x=element_text(size=16,angle=60),axis.text.y=element_text(size=16),title=element_text(hjust = 0.5,size=16)) + theme(legend.position='none')
plot2 <- ggplot(data=data_plot_Var,aes(x=meth,y=result))+geom_boxplot()+geom_hline(yintercept = CovTheo[j,j],color = "red")+ylab("")+xlab("")+scale_x_discrete(labels=c("MNAR","MAR","EMMAR","SoftMAR","Mean","Del")) +theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x=element_text(size=16,angle=60),axis.text.y=element_text(size=16),title=element_text(hjust = 0.5,size=16)) + theme(legend.position='none')
plot3 <- ggplot(data=data_plot_Cov1,aes(x=meth,y=result))+geom_boxplot()+geom_hline(yintercept = CovTheo[j,jmiss],color = "red")+ylab("")+xlab("")+scale_x_discrete(labels=c("MNAR","MAR","EMMAR","SoftMAR","Mean","Del")) +theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x=element_text(size=16,angle=60),axis.text.y=element_text(size=16),title=element_text(hjust = 0.5,size=16)) + theme(legend.position='none')
plot4 <- ggplot(data=data_plot_Cov2,aes(x=meth,y=result))+geom_boxplot()+geom_hline(yintercept = CovTheo[j,jobs],color = "red")+ylab("")+xlab("")+scale_x_discrete(labels=c("MNAR","MAR","EMMAR","SoftMAR","Mean","Del")) +theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x=element_text(size=16,angle=60),axis.text.y=element_text(size=16),title=element_text(hjust = 0.5,size=16)) + theme(legend.position='none')

grid.arrange(plot1,plot2,ncol=2) 
grid.arrange(plot3,plot4,ncol=2) 



#Correlation and MSE
Correlationglob <- c()
MSEglob <- c()
for (k in 1:Nbit){
  Correlationglob <- c(Correlationglob,result[[k]]$Correlation$MNAR,result[[k]]$Correlation$MAR,result[[k]]$Correlation$Soft,result[[k]]$Correlation$Mean,result[[k]]$Correlation$MAREM)
  MSEglob <- c(MSEglob,result[[k]]$MSE$MNAR,result[[k]]$MSE$MAR,result[[k]]$MSE$Soft,result[[k]]$MSE$Mean,result[[k]]$MSE$MAREM)
}
data_plot_Correlation <- data.frame(result=Correlationglob,meth=rep(c("2.MNAR","4.MAR","6.Soft","7.Mean","5.MAREM"),Nbit))
data_plot_MSE <- data.frame(result=MSEglob,meth=rep(c("2.MNAR","4.MAR","6.Soft","7.Mean","5.MAREM"),Nbit))

#use ylim(c(...,...)) if outliers. 
plot5 <- ggplot(data=data_plot_Correlation,aes(x=meth,y=result))+geom_boxplot()+ylab("")+xlab("")+scale_x_discrete(labels=c("MNAR","MAR","EMMAR","SoftMAR","Mean"))+theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x=element_text(size=16,angle=60),axis.text.y=element_text(size=16),title=element_text(size=16)) + theme(legend.position='none') 
plot6 <- ggplot(data=data_plot_MSE,aes(x=meth,y=result))+geom_boxplot()+ylab("")+xlab("")+scale_x_discrete(labels=c("MNAR","MAR","EMMAR","SoftMAR","Mean"))+theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x=element_text(size=16,angle=60),axis.text.y=element_text(size=16),title=element_text(size=16)) + theme(legend.position='none')
grid.arrange(plot5,plot6,ncol=2)

