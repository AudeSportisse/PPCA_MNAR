source("Function_simulations.R")
library(ggplot2)
library(gridExtra)

get_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

## Varying the noise level

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

sigma <- 0.3
result_03 <- lapply(1:Nbit,ComparMethods_PPCA_iteration,
                 n=n,
                 p=p,
                 r=r,
                 B=B,
                 mean_theo=mean_theo,
                 sigma=sigma,
                 indMissVar=indMissVar
)

sigma <- 0.5
result_05 <- lapply(1:Nbit,ComparMethods_PPCA_iteration,
                  n=n,
                  p=p,
                  r=r,
                  B=B,
                  mean_theo=mean_theo,
                  sigma=sigma,
                  indMissVar=indMissVar
)

sigma <- 0.7
result_07 <- lapply(1:Nbit,ComparMethods_PPCA_iteration,
                 n=n,
                 p=p,
                 r=r,
                 B=B,
                 mean_theo=mean_theo,
                 sigma=sigma,
                 indMissVar=indMissVar
)

sigma <- 1
result_1 <- lapply(1:Nbit,ComparMethods_PPCA_iteration,
                  n=n,
                  p=p,
                  r=r,
                  B=B,
                  mean_theo=mean_theo,
                  sigma=sigma,
                  indMissVar=indMissVar
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

result <- result_01
for (k in 1:Nbit){
  Meanglob <- c(Meanglob,result[[k]]$Mean$MNAR[j],result[[k]]$Mean$Soft[j],result[[k]]$Mean$Mean[j],result[[k]]$Mean$MAREM[j])
  Varglob <- c(Varglob,result[[k]]$Cov$MNAR[j,j],result[[k]]$Cov$Soft[j,j],result[[k]]$Cov$Mean[j,j],result[[k]]$Cov$MAREM[j,j])
  Covglob1 <- c(Covglob1,result[[k]]$Cov$MNAR[j,jmiss],result[[k]]$Cov$Soft[j,jmiss],result[[k]]$Cov$Mean[j,jmiss],result[[k]]$Cov$MAREM[j,jmiss])
  Covglob2 <- c(Covglob2,result[[k]]$Cov$MNAR[j,jobs],result[[k]]$Cov$Soft[j,jobs],result[[k]]$Cov$Mean[j,jobs],result[[k]]$Cov$MAREM[j,jobs])
  Correlationglob <- c(Correlationglob,result[[k]]$Correlation$MNAR,result[[k]]$Correlation$Soft,result[[k]]$Correlation$Mean,result[[k]]$Correlation$MAREM)
  MSEglob <- c(MSEglob,result[[k]]$MSE$MNAR,result[[k]]$MSE$Soft,result[[k]]$MSE$Mean,result[[k]]$MSE$MAREM)
}
result <- result_03
for (k in 1:Nbit){
  Meanglob <- c(Meanglob,result[[k]]$Mean$MNAR[j],result[[k]]$Mean$Soft[j],result[[k]]$Mean$Mean[j],result[[k]]$Mean$MAREM[j])
  Varglob <- c(Varglob,result[[k]]$Cov$MNAR[j,j],result[[k]]$Cov$Soft[j,j],result[[k]]$Cov$Mean[j,j],result[[k]]$Cov$MAREM[j,j])
  Covglob1 <- c(Covglob1,result[[k]]$Cov$MNAR[j,jmiss],result[[k]]$Cov$Soft[j,jmiss],result[[k]]$Cov$Mean[j,jmiss],result[[k]]$Cov$MAREM[j,jmiss])
  Covglob2 <- c(Covglob2,result[[k]]$Cov$MNAR[j,jobs],result[[k]]$Cov$Soft[j,jobs],result[[k]]$Cov$Mean[j,jobs],result[[k]]$Cov$MAREM[j,jobs])
  Correlationglob <- c(Correlationglob,result[[k]]$Correlation$MNAR,result[[k]]$Correlation$Soft,result[[k]]$Correlation$Mean,result[[k]]$Correlation$MAREM)
  MSEglob <- c(MSEglob,result[[k]]$MSE$MNAR,result[[k]]$MSE$Soft,result[[k]]$MSE$Mean,result[[k]]$MSE$MAREM)
}
result <- result_05
for (k in 1:Nbit){
  Meanglob <- c(Meanglob,result[[k]]$Mean$MNAR[j],result[[k]]$Mean$Soft[j],result[[k]]$Mean$Mean[j],result[[k]]$Mean$MAREM[j])
  Varglob <- c(Varglob,result[[k]]$Cov$MNAR[j,j],result[[k]]$Cov$Soft[j,j],result[[k]]$Cov$Mean[j,j],result[[k]]$Cov$MAREM[j,j])
  Covglob1 <- c(Covglob1,result[[k]]$Cov$MNAR[j,jmiss],result[[k]]$Cov$Soft[j,jmiss],result[[k]]$Cov$Mean[j,jmiss],result[[k]]$Cov$MAREM[j,jmiss])
  Covglob2 <- c(Covglob2,result[[k]]$Cov$MNAR[j,jobs],result[[k]]$Cov$Soft[j,jobs],result[[k]]$Cov$Mean[j,jobs],result[[k]]$Cov$MAREM[j,jobs])
  Correlationglob <- c(Correlationglob,result[[k]]$Correlation$MNAR,result[[k]]$Correlation$Soft,result[[k]]$Correlation$Mean,result[[k]]$Correlation$MAREM)
  MSEglob <- c(MSEglob,result[[k]]$MSE$MNAR,result[[k]]$MSE$Soft,result[[k]]$MSE$Mean,result[[k]]$MSE$MAREM)
}
result <- result_07
for (k in 1:Nbit){
  Meanglob <- c(Meanglob,result[[k]]$Mean$MNAR[j],result[[k]]$Mean$Soft[j],result[[k]]$Mean$Mean[j],result[[k]]$Mean$MAREM[j])
  Varglob <- c(Varglob,result[[k]]$Cov$MNAR[j,j],result[[k]]$Cov$Soft[j,j],result[[k]]$Cov$Mean[j,j],result[[k]]$Cov$MAREM[j,j])
  Covglob1 <- c(Covglob1,result[[k]]$Cov$MNAR[j,jmiss],result[[k]]$Cov$Soft[j,jmiss],result[[k]]$Cov$Mean[j,jmiss],result[[k]]$Cov$MAREM[j,jmiss])
  Covglob2 <- c(Covglob2,result[[k]]$Cov$MNAR[j,jobs],result[[k]]$Cov$Soft[j,jobs],result[[k]]$Cov$Mean[j,jobs],result[[k]]$Cov$MAREM[j,jobs])
  Correlationglob <- c(Correlationglob,result[[k]]$Correlation$MNAR,result[[k]]$Correlation$Soft,result[[k]]$Correlation$Mean,result[[k]]$Correlation$MAREM)
  MSEglob <- c(MSEglob,result[[k]]$MSE$MNAR,result[[k]]$MSE$Soft,result[[k]]$MSE$Mean,result[[k]]$MSE$MAREM)
}
result <- result_1
for (k in 1:Nbit){
  Meanglob <- c(Meanglob,result[[k]]$Mean$MNAR[j],result[[k]]$Mean$Soft[j],result[[k]]$Mean$Mean[j],result[[k]]$Mean$MAREM[j])
  Varglob <- c(Varglob,result[[k]]$Cov$MNAR[j,j],result[[k]]$Cov$Soft[j,j],result[[k]]$Cov$Mean[j,j],result[[k]]$Cov$MAREM[j,j])
  Covglob1 <- c(Covglob1,result[[k]]$Cov$MNAR[j,jmiss],result[[k]]$Cov$Soft[j,jmiss],result[[k]]$Cov$Mean[j,jmiss],result[[k]]$Cov$MAREM[j,jmiss])
  Covglob2 <- c(Covglob2,result[[k]]$Cov$MNAR[j,jobs],result[[k]]$Cov$Soft[j,jobs],result[[k]]$Cov$Mean[j,jobs],result[[k]]$Cov$MAREM[j,jobs])
  Correlationglob <- c(Correlationglob,result[[k]]$Correlation$MNAR,result[[k]]$Correlation$Soft,result[[k]]$Correlation$Mean,result[[k]]$Correlation$MAREM)
  MSEglob <- c(MSEglob,result[[k]]$MSE$MNAR,result[[k]]$MSE$Soft,result[[k]]$MSE$Mean,result[[k]]$MSE$MAREM)
}

meth <- rep(c("1.MNAR","3.Soft","4.Mean",'2.MAREM'),Nbit*5)
noise <- c(rep("0.1",Nbit*4),rep("0.25",Nbit*4),rep("0.5",Nbit*4),rep("0.75",Nbit*4),rep("1",Nbit*4))

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
Var=data.frame(result=Varglob,meth=meth,noise=noise)
plotVar <- ggplot(Var, aes(x=noise, y=result,color=meth)) +geom_boxplot(size=1.3) + xlab("") + ylab("")+theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),title=element_text(size=24),legend.text=element_text(size=26))+ scale_color_manual("",labels=c('MNAR','EMMAR','SoftMAR','Mean'),values=c("#669900","#FF9933","#993300","black"))
plotVar <- plotVar + geom_segment(aes(x=0.6,y=theo[1],yend=theo[1],xend=1.4),color="red")
plotVar <- plotVar + geom_segment(aes(x=1.6,y=theo[2],yend=theo[2],xend=2.4),color="red")
plotVar <- plotVar + geom_segment(aes(x=2.6,y=theo[3],yend=theo[3],xend=3.4),color="red")
plotVar <- plotVar + geom_segment(aes(x=3.6,y=theo[4],yend=theo[4],xend=4.4),color="red")
#use ylim(c(...,...)) if outliers. 
plotVar <- plotVar + geom_segment(aes(x=4.6,y=theo[5],yend=theo[5],xend=5.4),color="red") + theme(legend.position = "none")

grid.arrange(plotMean,plotVar,l,ncol=3) 

theo <- c()
sigChoice <- c(0.1,0.25,0.5,0.75,1)
for (i in 1:5){
  theo <- c(theo,(t(B) %*% B + sigChoice[i] ^ 2 * diag(1, ncol = p, nrow = p))[j,jmiss])
}

CovMiss <- data.frame(result=Covglob1,meth=meth,noise=noise)
plotMiss <- ggplot(CovMiss, aes(x=noise, y=result,color=meth)) +geom_boxplot(size=1.3) + xlab("") + ylab("")+theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),title=element_text(size=24),legend.text=element_text(size=26))+ scale_color_manual("",labels=c('MNAR','EMMAR','SoftMAR','Mean'),values=c("#669900","#FF9933","#993300","black"))
plotMiss <- plotMiss + geom_segment(aes(x=0.6,y=theo[1],yend=theo[1],xend=1.4),color="red")
plotMiss <- plotMiss + geom_segment(aes(x=1.6,y=theo[2],yend=theo[2],xend=2.4),color="red")
plotMiss <- plotMiss + geom_segment(aes(x=2.6,y=theo[3],yend=theo[3],xend=3.4),color="red")
plotMiss <- plotMiss + geom_segment(aes(x=3.6,y=theo[4],yend=theo[4],xend=4.4),color="red")
#use ylim(c(...,...)) if outliers. 
plotMiss <- plotMiss + geom_segment(aes(x=4.6,y=theo[5],yend=theo[5],xend=5.4),color="red") + theme(legend.position = "none")

theo2 <- c()
sigChoice <- c(0.1,0.25,0.5,0.75,1)
for (i in 1:5){
  theo2 <- c(theo2,(t(B) %*% B + sigChoice[i] ^ 2 * diag(1, ncol = p, nrow = p))[j,9])
}

CovObs <- data.frame(result=Covglob2,meth=meth,noise=noise)
plotObs <- ggplot(CovObs, aes(x=noise, y=result,color=meth)) +geom_boxplot(size=1.3) + xlab("") + ylab("")+theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),title=element_text(size=24),legend.text=element_text(size=26))+ scale_color_manual("",labels=c('MNAR','EMMAR','SoftMAR','Mean'),values=c("#669900","#FF9933","#993300","black"))
plotObs <- plotObs + geom_segment(aes(x=0.6,y=theo2[1],yend=theo2[1],xend=1.4),color="red")
plotObs <- plotObs + geom_segment(aes(x=1.6,y=theo2[2],yend=theo2[2],xend=2.4),color="red")
plotObs <- plotObs + geom_segment(aes(x=2.6,y=theo2[3],yend=theo2[3],xend=3.4),color="red")
plotObs <- plotObs + geom_segment(aes(x=3.6,y=theo2[4],yend=theo2[4],xend=4.4),color="red")
#use ylim(c(...,...)) if outliers. 
plotObs <- plotObs + geom_segment(aes(x=4.6,y=theo2[5],yend=theo2[5],xend=5.4),color="red") + theme(legend.position = "none")

grid.arrange(plotMiss,plotObs,l,ncol=3) 


#Correlation and MSE

Mse <- data.frame(result=MSEglob,meth=meth,noise=noise)
#use ylim(c(...,...)) if outliers. 
plotMSE <- ggplot(Mse, aes(x=noise, y=result,color=meth)) +geom_boxplot(size=1.3) + xlab("") + ylab("")+theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),title=element_text(size=24),legend.text=element_text(size=26))+ scale_color_manual("",labels=c('MNAR','EMMAR','SoftMAR','Mean'),values=c("#669900","#FF9933","#993300","black")) + theme(legend.position = "none")
Corr <- data.frame(result=Correlationglob,meth=meth,noise=noise)
#use ylim(c(...,...)) if outliers. 
plotCorr <- ggplot(Corr, aes(x=noise, y=result,color=meth)) +geom_boxplot(size=1.3) + xlab("") + ylab("")+theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),title=element_text(size=24),legend.text=element_text(size=26))+ scale_color_manual("",labels=c('MNAR','EMMAR','SoftMAR','Mean'),values=c("#669900","#FF9933","#993300","black")) + theme(legend.position = "none")

grid.arrange(plotCorr,plotMSE,l,ncol=3) 

