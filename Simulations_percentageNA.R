source("Function_simulations.R")
library(ggplot2)
library(gridExtra)

get_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

## Varying the percentage of missing values

set.seed(2002)
Nbit <- 2
p <- 10
r <- 2
sigma <- 0.1
indMissVar <- 1:7
mean_theo <- 0
B <- matrix(rnorm(p * r), nrow = r, ncol = p)
n <- 1000

result_10percNA <- lapply(1:Nbit,ComparMethods_PPCA_iteration,
                       n=n,
                       p=p,
                       r=r,
                       B=B,
                       mean_theo=mean_theo,
                       sigma=sigma,
                       indMissVar=indMissVar,
                       param_logistic=c(-1,-2.5)
)

result_35percNA <- lapply(1:Nbit,ComparMethods_PPCA_iteration,
                 n=n,
                 p=p,
                 r=r,
                 B=B,
                 mean_theo=mean_theo,
                 sigma=sigma,
                 indMissVar=indMissVar,
                 param_logistic=c(3,0)
)

result_50percNA <- lapply(1:Nbit,ComparMethods_PPCA_iteration,
                 n=n,
                 p=p,
                 r=r,
                 B=B,
                 mean_theo=mean_theo,
                 sigma=sigma,
                 indMissVar=indMissVar,
                 param_logistic=c(1,-1.3)
)


##Plot 

#Mean, Variance, Covariances

j <- 3 #the missing variable for which the mean, variance and covariances is ploted.
jmiss <- 4 #the missing variable for which the covariance is ploted in plot3.
jobs <- 9 #the observed variable for which the covariance is ploted in plot4.

Meanglob <- c()
Varglob <- c()
Covglob1 <- c()
Covglob2 <- c()
Correlationglob <- c()
MSEglob <- c()

result <- result_10percNA
for (k in 1:Nbit){
  Meanglob <- c(Meanglob,result[[k]]$Mean$MNAR[j],result[[k]]$Mean$Soft[j],result[[k]]$Mean$Mean[j],result[[k]]$Mean$MAREM[j])
  Varglob <- c(Varglob,result[[k]]$Cov$MNAR[j,j],result[[k]]$Cov$Soft[j,j],result[[k]]$Cov$Mean[j,j],result[[k]]$Cov$MAREM[j,j])
  Covglob1 <- c(Covglob1,result[[k]]$Cov$MNAR[j,jmiss],result[[k]]$Cov$Soft[j,jmiss],result[[k]]$Cov$Mean[j,jmiss],result[[k]]$Cov$MAREM[j,jmiss])
  Covglob2 <- c(Covglob2,result[[k]]$Cov$MNAR[j,jobs],result[[k]]$Cov$Soft[j,jobs],result[[k]]$Cov$Mean[j,jobs],result[[k]]$Cov$MAREM[j,jobs])
  Correlationglob <- c(Correlationglob,result[[k]]$Correlation$MNAR,result[[k]]$Correlation$Soft,result[[k]]$Correlation$Mean,result[[k]]$Correlation$MAREM)
  MSEglob <- c(MSEglob,result[[k]]$MSE$MNAR,result[[k]]$MSE$Soft,result[[k]]$MSE$Mean,result[[k]]$MSE$MAREM)
}
result <- result_35percNA
for (k in 1:Nbit){
  Meanglob <- c(Meanglob,result[[k]]$Mean$MNAR[j],result[[k]]$Mean$Soft[j],result[[k]]$Mean$Mean[j],result[[k]]$Mean$MAREM[j])
  Varglob <- c(Varglob,result[[k]]$Cov$MNAR[j,j],result[[k]]$Cov$Soft[j,j],result[[k]]$Cov$Mean[j,j],result[[k]]$Cov$MAREM[j,j])
  Covglob1 <- c(Covglob1,result[[k]]$Cov$MNAR[j,jmiss],result[[k]]$Cov$Soft[j,jmiss],result[[k]]$Cov$Mean[j,jmiss],result[[k]]$Cov$MAREM[j,jmiss])
  Covglob2 <- c(Covglob2,result[[k]]$Cov$MNAR[j,jobs],result[[k]]$Cov$Soft[j,jobs],result[[k]]$Cov$Mean[j,jobs],result[[k]]$Cov$MAREM[j,jobs])
  Correlationglob <- c(Correlationglob,result[[k]]$Correlation$MNAR,result[[k]]$Correlation$Soft,result[[k]]$Correlation$Mean,result[[k]]$Correlation$MAREM)
  MSEglob <- c(MSEglob,result[[k]]$MSE$MNAR,result[[k]]$MSE$Soft,result[[k]]$MSE$Mean,result[[k]]$MSE$MAREM)
}
result <- result_50percNA
for (k in 1:Nbit){
  Meanglob <- c(Meanglob,result[[k]]$Mean$MNAR[j],result[[k]]$Mean$Soft[j],result[[k]]$Mean$Mean[j],result[[k]]$Mean$MAREM[j])
  Varglob <- c(Varglob,result[[k]]$Cov$MNAR[j,j],result[[k]]$Cov$Soft[j,j],result[[k]]$Cov$Mean[j,j],result[[k]]$Cov$MAREM[j,j])
  Covglob1 <- c(Covglob1,result[[k]]$Cov$MNAR[j,jmiss],result[[k]]$Cov$Soft[j,jmiss],result[[k]]$Cov$Mean[j,jmiss],result[[k]]$Cov$MAREM[j,jmiss])
  Covglob2 <- c(Covglob2,result[[k]]$Cov$MNAR[j,jobs],result[[k]]$Cov$Soft[j,jobs],result[[k]]$Cov$Mean[j,jobs],result[[k]]$Cov$MAREM[j,jobs])
  Correlationglob <- c(Correlationglob,result[[k]]$Correlation$MNAR,result[[k]]$Correlation$Soft,result[[k]]$Correlation$Mean,result[[k]]$Correlation$MAREM)
  MSEglob <- c(MSEglob,result[[k]]$MSE$MNAR,result[[k]]$MSE$Soft,result[[k]]$MSE$Mean,result[[k]]$MSE$MAREM)
}


meth <- rep(c("1.MNAR","3.Soft","4.Mean",'2.MAREM'),Nbit*3)
noise <- noise=c(rep("0.1",Nbit*4),rep("0.35",Nbit*4),rep("0.5",Nbit*4))

Mean <- data.frame(result=Meanglob,meth=meth,noise=noise)
#use ylim(c(...,...)) if outliers. 
plotMean <- ggplot(Mean, aes(x=noise, y=result,color=meth)) +geom_boxplot(size=1.3) + xlab("") + ylab("")+theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),title=element_text(size=24),legend.text=element_text(size=20))+ scale_color_manual("",labels=c('MNAR','EMMAR','SoftMAR','Mean'),values=c("#669900","#FF9933","#993300","black"))+geom_hline(yintercept=0,color="red")
l <- get_legend(plotMean)
plotMean <- plotMean + theme(legend.position = "none")

theo <- c()
sigChoice <- c(0.1,0.25,0.5,0.75,1)
for (i in 1:5){
  theo <- c(theo,(t(B) %*% B + sigChoice[i] ^ 2 * diag(1, ncol = p, nrow = p))[j,j])
}
Var <- data.frame(result=Varglob,meth=meth,noise=noise)
#use ylim(c(...,...)) if outliers. 
plotVar <- ggplot(Var, aes(x=noise, y=result,color=meth)) +geom_boxplot(size=1.3) + xlab("") + ylab("")+theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),title=element_text(size=24),legend.text=element_text(size=26))+ scale_color_manual("",labels=c('MNAR','EMMAR','SoftMAR','Mean'),values=c("#669900","#FF9933","#993300","black"))+ theme(legend.position = "none")+geom_hline(yintercept = CovTheo[j,j],color = "red")

Mse <- data.frame(result=MSEglob,meth=meth,noise=noise)
#use ylim(c(...,...)) if outliers. 
plotMSE <- ggplot(Mse, aes(x=noise, y=result,color=meth)) +geom_boxplot(size=1.3) + xlab("") + ylab("")+theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),title=element_text(size=24),legend.text=element_text(size=26))+ scale_color_manual("",labels=c('MNAR','EMMAR','SoftMAR','Mean'),values=c("#669900","#FF9933","#993300","black"))+ theme(legend.position = "none")


grid.arrange(plotMean,plotVar,plotMSE,l,ncol=4) 
