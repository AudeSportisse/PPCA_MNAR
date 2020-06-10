source("Function_simulations_MNARgeneral.R")
library(ggplot2)
library(gridExtra)

set.seed(2002)
Nbit <- 20
p <- 20
r <- 3
sigma <- 0.8
indMissVar <- 1:10
mean_theo <- 0
B <- matrix(rnorm(p * r), nrow = r, ncol = p)
n <- 1000

#Covariance matrix (theoretical)
CovTheo <- t(B) %*% B + sigma ^ 2 * diag(1, ncol = p, nrow = p)


result_MNARgeneral <- lapply(1:Nbit,ComparMethods_PPCA_iteration_MNARgeneral,
                             n=n,
                             p=p,
                             r=r,
                             B=B,
                             mean_theo=mean_theo,
                             sigma=sigma,
                             indMissVar=indMissVar
)

result <- result_MNARgeneral

##Plot 

#Mean, Variance
j <- 3 #the missing variable for which the mean, variance and covariances is ploted.
Meanglob <- c()
Varglob <- c()
for (k in 1:Nbit){
  Meanglob <- c(Meanglob,result[[k]]$Mean$MNAR[j],result[[k]]$Mean$CC[j],result[[k]]$Mean$MAR[j],result[[k]]$Mean$Soft[j],result[[k]]$Mean$Mean[j],result[[k]]$Mean$MAREM[j])
  Varglob <- c(Varglob,result[[k]]$Cov$MNAR[j,j],result[[k]]$Cov$CC[j,j],result[[k]]$Cov$MAR[j,j],result[[k]]$Cov$Soft[j,j],result[[k]]$Cov$Mean[j,j],result[[k]]$Cov$MAREM[j,j])
}
data_plot_Mean <- data.frame(result=Meanglob,meth=rep(c("2.MNAR","8.CC","4.MAR","6.Soft","7.Mean","5.MAREM"),Nbit))
data_plot_Var <- data.frame(result=Varglob,meth=rep(c("2.MNAR","8.CC","4.MAR","6.Soft","7.Mean","5.MAREM"),Nbit))

#use ylim(c(...,...)) if outliers. 
plot1 <- ggplot(data=data_plot_Mean,aes(x=meth,y=result))+geom_boxplot()+geom_hline(yintercept = mean_theo,color = "red")+ylab("")+xlab("")+scale_x_discrete(labels=c("MNAR","MAR","EMMAR","SoftMAR","Mean","Del")) +theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x=element_text(size=16,angle=60),axis.text.y=element_text(size=16),title=element_text(hjust = 0.5,size=16)) + theme(legend.position='none')
plot2 <- ggplot(data=data_plot_Var,aes(x=meth,y=result))+geom_boxplot()+geom_hline(yintercept = CovTheo[j,j],color = "red")+ylab("")+xlab("")+scale_x_discrete(labels=c("MNAR","MAR","EMMAR","SoftMAR","Mean","Del")) +theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x=element_text(size=16,angle=60),axis.text.y=element_text(size=16),title=element_text(hjust = 0.5,size=16)) + theme(legend.position='none')



#MSE
MSEglob <- c()
for (k in 1:Nbit){
  MSEglob <- c(MSEglob,result[[k]]$MSE$MNAR,result[[k]]$MSE$MAR,result[[k]]$MSE$Soft,result[[k]]$MSE$Mean,result[[k]]$MSE$MAREM)
}
data_plot_MSE <- data.frame(result=MSEglob,meth=rep(c("2.MNAR","4.MAR","6.Soft","7.Mean","5.MAREM"),Nbit))
#use ylim(c(...,...)) if outliers. 
plot3 <- ggplot(data=data_plot_MSE,aes(x=meth,y=result))+geom_boxplot()+ylab("")+xlab("")+scale_x_discrete(labels=c("MNAR","MAR","EMMAR","SoftMAR","Mean"))+theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x=element_text(size=16,angle=60),axis.text.y=element_text(size=16),title=element_text(size=16)) + theme(legend.position='none')
grid.arrange(plot1,plot2,plot3,ncol=3)

