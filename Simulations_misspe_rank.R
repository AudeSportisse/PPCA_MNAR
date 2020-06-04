source("Function_simulations.R")
library(ggplot2)
library(gridExtra)

set.seed(2002)
Nbit <- 20
p <- 20
r_theo <- 3
sigma <- 0.8
indMissVar <- 1:10
mean_theo <- 0
B <- matrix(rnorm(p * r_theo), nrow = r_theo, ncol = p)
n <- 1000

result_restim2 <- lapply(1:Nbit,ComparMethods_PPCA_iteration,
                    n=n,
                    p=p,
                    r=2,
                    B=B,
                    mean_theo=mean_theo,
                    sigma=sigma,
                    indMissVar=indMissVar,
                    r_theo=r_theo
)

result_restim3 <- lapply(1:Nbit,ComparMethods_PPCA_iteration,
                         n=n,
                         p=p,
                         r=3,
                         B=B,
                         mean_theo=mean_theo,
                         sigma=sigma,
                         indMissVar=indMissVar,
                         r_theo=r_theo
)

result_restim4 <- lapply(1:Nbit,ComparMethods_PPCA_iteration,
                         n=n,
                         p=p,
                         r=4,
                         B=B,
                         mean_theo=mean_theo,
                         sigma=sigma,
                         indMissVar=indMissVar,
                         r_theo=r_theo
)

Correlationglob <- c()
MSEglob <- c()
for (k in 1:Nbit){
  Correlationglob <- c(Correlationglob,result_restim2[[k]]$Correlation$MNAR,result_restim3[[k]]$Correlation$MNAR,result_restim4[[k]]$Correlation$MNAR)
  MSEglob <- c(MSEglob,result_restim2[[k]]$MSE$MNAR,result_restim3[[k]]$MSE$MNAR,result_restim4[[k]]$MSE$MNAR)
}
data_plot_Correlation <- data.frame(result=Correlationglob,meth=rep(c("r=2","r=3","r=4"),Nbit))
data_plot_MSE <- data.frame(result=MSEglob,meth=c(rep(c("r=2","r=3","r=4"),Nbit)))
plotCorr=ggplot(data=data_plot_Correlation,aes(x=meth,y=result))+geom_boxplot()+ylab("")+xlab("")+scale_x_discrete(labels=c("r=2","r=3","r=4"))+theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x=element_text(vjust=0.98,hjust=0.8,size=20,angle=60),axis.text.y=element_text(size=16),title=element_text(size=16)) + theme(legend.position='none') 
plotMSE=ggplot(data=data_plot_MSE,aes(x=meth,y=result))+geom_boxplot()+ylab("")+xlab("")+scale_x_discrete(labels=c("r=2","r=3","r=4"))+theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),axis.text.x=element_text(vjust=0.98,hjust=0.8,size=20,angle=60),axis.text.y=element_text(size=16),title=element_text(size=16)) + theme(legend.position='none')
grid.arrange(plotCorr,plotMSE,ncol=2)