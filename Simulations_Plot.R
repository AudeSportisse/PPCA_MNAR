#Libraries:
library("FactoMineR")
library("softImpute")
library("MASS")
library("gtools")

#Parameters
n <- 1000
p <- 50 #50
r <- 10 ##No rank 1. #10
sigma <- 10^-4
mean_theo <- 0
indcolNA <- 30:50 #c(1,  2 , 3,  4, 5, 9,10)

set.seed(23)
B <- matrix(rnorm(p * r,sd=3), nrow = r, ncol = p)

#Mechanism choice
nbNA <- 500 #0.2*n
mecha <- "MNAR"
modmecha <- "Logistic" #or "Censor"

#Number of iterations
Nbit <- 10

#Aggregation
best <- 0 #choice of the best coefficient for the expectation
agg <- 1 #aggregation with median
noagg <- 0 #at random
shrink <- 0 #regularization for the algebraic method

#simplify <- FALSE #FALSE if all combinations of observed variables are used for the regressions. 


result_r2=lapply(1:Nbit,ComparMethods_PPCA_iteration,n=n,
                  p=p,
                  r=r,
                  B=B,
                  mean_theo=mean_theo,
                  sigma=sigma,
                  indcolNA=indcolNA,
                  nbNA=nbNA,
                  modmecha=modmecha,
                  best = 0,
                  agg = 1,
                  noagg = 0,
                  shrink = 0,
                  simplify = TRUE)


#Plot
library(ggplot2)
library(gridExtra)

#Covariance matrix (theoretical)
CovTheo <- t(B) %*% B + sigma ^ 2 * diag(1, ncol = p, nrow = p)

#Plot mean and variance estimations
for (j in 1:length(indcolNA)){
  Meanglob <- c()
  Varglob <- c()
  for (k in 1:Nbit){
    Meanglob <- c(Meanglob,result_r2[[k]][[1]][seq(j,length(result_r2[[k]][[1]]),by=length(indcolNA))])
    Varglob <- c(Varglob,result_r2[[k]][[2]][seq(j,length(result_r2[[k]][[2]]),by=length(indcolNA))])
  }
  data_plot_Mean <- data.frame(result=Meanglob,meth=rep(c("1.MNAR algebraic","6.CC","2.MAR","4.Soft","5.Mean","3.MAREM"),Nbit))
  data_plot_Var <- data.frame(result=Varglob,meth=rep(c("1.MNAR graphical","2.MNAR algebraic","7.CC","3.MAR","5.Soft","6.Mean","4.MAREM"),Nbit))
  plot1=ggplot(data=data_plot_Mean,aes(x=data_plot_Mean$meth,y=data_plot_Mean$result))+geom_boxplot()+geom_hline(yintercept = mean_theo,color = "red")+ylab(paste("Mean",indcolNA[j]))+xlab("")+scale_x_discrete(labels=c("MNAR","MAR","EMMAR","SoftMAR","Mean","Del")) +theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=22,angle=45),axis.text.y=element_text(size=22),title=element_text(size=22)) + theme(legend.position='none')
  plot2=ggplot(data=data_plot_Var,aes(x=data_plot_Var$meth,y=data_plot_Var$result))+geom_boxplot()+geom_hline(yintercept = CovTheo[indcolNA[j],indcolNA[j]],color = "red")+ylab(paste("Variance",indcolNA[j]))+xlab("")+scale_x_discrete(labels=c("MNAR graph","MNAR alg","MAR","EMMAR","SoftMAR","Mean","Del")) +theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=22,angle=45),axis.text.y=element_text(size=22),title=element_text(size=22)) + theme(legend.position='none')
  grid.arrange(plot1,plot2,ncol=2)
}

#Plot covariances
j=32
data_plot_Cov=NULL
for (l in 29:34){
  matrixres <- matrix(0,nrow=7,ncol=Nbit)
  k=1
  for (ibis in 1:Nbit){
    matrixres[1,k]<-result_r2[[k]][[3]][[1]][j,l]
    matrixres[2,k]<-result_r2[[k]][[3]][[2]][j,l]
    matrixres[3,k]<-result_r2[[k]][[3]][[3]][j,l]
    matrixres[4,k]<-result_r2[[k]][[3]][[4]][j,l]
    matrixres[5,k]<-result_r2[[k]][[3]][[5]][j,l]
    matrixres[6,k]<-result_r2[[k]][[3]][[6]][j,l]
    matrixres[7,k]<-result_r2[[k]][[3]][[7]][j,l]
    k <- k+1
  }
  data_plot_Cov[[l]]=data.frame(result33=c(matrixres[1,],matrixres[2,],matrixres[3,],matrixres[4,],matrixres[5,],matrixres[6,],matrixres[7,]),meth=c(rep("1.MNAR graphical",Nbit),rep("2.MNAR algebraic",Nbit),rep("7.CC",Nbit),rep("3.MAR",Nbit),rep("5.Soft",Nbit),rep("6.Mean",Nbit),rep("4.MAREM",Nbit)))
}
plot2=ggplot(data=data_plot_Cov[[29]],aes(x=data_plot_Cov[[29]]$meth,y=data_plot_Cov[[29]]$result33)) + geom_boxplot()  + ylab(paste("Covariance",j,29))+xlab("")+scale_x_discrete(labels=c("MNAR graph","MNAR alg","MAR","EMMAR","SoftMAR","Mean","Del")) + geom_hline(yintercept = CovTheo[j,29],color = "red") +theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=22,angle=45),axis.text.y=element_text(size=22),title=element_text(size=22)) + theme(legend.position='none')#+scale_y_continuous(limits=c(lim1,lim2)) 
#plot3=ggplot(data=data_plot_Cov[[34]],aes(x=data_plot_Cov[[34]]$meth,y=data_plot_Cov[[34]]$result33)) + geom_boxplot()  + ylab(paste("Covariance",j,30))+xlab("")+scale_x_discrete(labels=c("MNAR graph","MNAR alg","MAR","EMMAR","SoftMAR","Mean","Del")) + geom_hline(yintercept = CovTheo[j,30],color = "red") +theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=22,angle=45),axis.text.y=element_text(size=22),title=element_text(size=22)) + theme(legend.position='none')#+scale_y_continuous(limits=c(lim1,lim2)) 
plot4=ggplot(data=data_plot_Cov[[31]],aes(x=data_plot_Cov[[31]]$meth,y=data_plot_Cov[[31]]$result33)) + geom_boxplot()  + ylab(paste("Covariance",j,31))+xlab("")+scale_x_discrete(labels=c("MNAR graph","MNAR alg","MAR","EMMAR","SoftMAR","Mean","Del")) + geom_hline(yintercept = CovTheo[j,31],color = "red") +theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=22,angle=45),axis.text.y=element_text(size=22),title=element_text(size=22)) + theme(legend.position='none')#+scale_y_continuous(limits=c(lim1,lim2)) 
grid.arrange(plot2,plot4,ncol=2)




#Correlation and MSE
Correlationglob <- c()
MSEglob <- c()
for (k in 1:Nbit){
  Correlationglob <- c(Correlationglob,result_r2[[k]][[4]])
  MSEglob <- c(MSEglob,result_r2[[k]][[5]])
}
data_plot_Correlation <- data.frame(result3=Correlationglob,meth=rep(c("1.MNAR graphical","2.MNAR algebraic","3.MAR","5.Soft","6.Mean","4.MAREM"),Nbit))
data_plot_MSE <- data.frame(result3=MSEglob,meth=rep(c("1.MNAR graphical","2.MNAR algebraic","3.MAR","5.Soft","6.Mean","4.MAREM"),Nbit))
plot1=ggplot(data=data_plot_Correlation,aes(x=data_plot_Correlation$meth,y=data_plot_Correlation$result3))+geom_boxplot()+ylab("Correlation")+xlab("")+scale_x_discrete(labels=c("MNAR graph","MNAR alg","MAR","EMMAR","SoftMAR","Mean"))+theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=22,angle=45),axis.text.y=element_text(size=22),title=element_text(size=22)) + theme(legend.position='none')
plot2=ggplot(data=data_plot_MSE,aes(x=data_plot_MSE$meth,y=data_plot_MSE$result3))+geom_boxplot()+ylab("Prediction error")+xlab("")+scale_x_discrete(labels=c("MNAR graph","MNAR alg","MAR","EMMAR","SoftMAR","Mean"))+scale_y_continuous(limits=c(0,2)) +theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=22,angle=45),axis.text.y=element_text(size=22),title=element_text(size=22)) + theme(legend.position='none')
grid.arrange(plot1,plot2,ncol=2)

    
  


