# Find the cases where pairwise comparison fails but simulatanous confidence interval approach succeed.
rm(list=ls())

library(ggplot2)

setwd("//fda.gov/wodc/CDER/Users05/Donglei.Yin/Result")

filename='n10_99_101_100_equal_var_sdR1_1.33_14_rp1000'
#filename='n10_90_110_100_equal_var_sdR1_2_14_rp1000_0706'
#filename='n10_95_105_100_equal_var_sdR1_2_14_rp1000_0709'

sd_margin=1.33

arcd_result<-read.csv(paste0("./Arcd_",filename, ".csv"),header=T)[,-1]

dim(arcd_result)

SD=c(1.33,2,4,6,8,10,12,14)

mu.R1 <- 99
mu.R2 <- 101
mu.T <- 100
n=10
p.level=0.1

out=c()


k=4
  
  
  true12 <- abs(mu.R1-mu.R2)/SD[k]
  true1T <- abs(mu.R1-mu.T)/SD[k]
  true2T <- abs(mu.T-mu.R2)/SD[k]
  
  arcd.data<-arcd_result[,1:400]
  arcd.mean<-arcd_result[,401:424]
  arcd.sd<-arcd_result[,425:464]
  arcd.sig<- arcd_result[,465:488]
  arcd.simuCI1 <- arcd_result[,489:560]
  
  rcd.data<-arcd.data[,(n*5*(k-1)+1):(n*5*k)]
  rcd.mean<-arcd.mean[,(3*(k-1)+1):(3*k)]
  rcd.std<-arcd.sd[,((3+2)*(k-1)+1):((3+2)*k)]
  rcd.sig<-arcd.sig[,(3*(k-1)+1):(3*k)]
  rcd.simuCI1<-arcd.simuCI1[,(9*(k-1)+1):(9*k)]
  
  rcd.sig.pairwise<-cbind(rcd.sig,apply(rcd.sig[2:3],1,prod),apply(rcd.sig,1,prod))
  
  summ.simuCI1 <- cbind(rcd.simuCI1[,1],
                        (rcd.simuCI1[,1]>=p.level),
                        (rcd.simuCI1[,2]/rcd.std[,4]),
                        (rcd.simuCI1[,2]/rcd.std[,4]>=true1T),
                        (rcd.simuCI1[,3]/rcd.std[,4]),
                        (rcd.simuCI1[,3]/rcd.std[,4]>=max(true1T,true2T,true12)),
                        (rcd.simuCI1[,4]),
                        (rcd.simuCI1[,4]>=p.level),
                        (rcd.simuCI1[,5]*1.5),
                        (rcd.simuCI1[,5]*1.5>=true1T),
                        (rcd.simuCI1[,6]*1.5),
                        (rcd.simuCI1[,6]*1.5>=max(true1T,true2T,true12)),
                        (rcd.simuCI1[,7]),
                        (rcd.simuCI1[,7]>=p.level),
                        (rcd.simuCI1[,8]/rcd.std[,4]),
                        (rcd.simuCI1[,8]/rcd.std[,4]>=true1T),
                        (rcd.simuCI1[,9]/rcd.std[,4]),
                        (rcd.simuCI1[,9]/rcd.std[,4]>=max(true1T,true2T,true12)))
  
  
  temp<-cbind(rcd.sig.pairwise,summ.simuCI1)
  
  colnames(temp)<-c("Pair_R1_R2","Pair_R1_T","Pair_R2_T","Pair_23","Pair_123",
                    "Orig_FP","Orig_Power","Orig_CI1_D","Orig_CI1_CR","Orig_CI2_D","Orig_CI2_CR",
                    "Inte_FP","Inte_Power","Inte_CI1_D","Inte_CI1_CR","Inte_CI2_D","Inte_CI2_CR",
                    "LF_FP","LF_Power","LF_CI1_D","LF_CI1_CR","LF_CI2_D","LF_CI2_CR")
  result=as.data.frame(temp[,c("Pair_123","Orig_Power","Inte_Power","LF_Power")])
  
  result$diff_FN<-ifelse((result$Pair_123==0) & (result$Orig_Power==1 | (result$Inte_Power==1) | (result$LF_Power==1)),TRUE,FALSE)
  result$diff_FP<-ifelse((result$Pair_123==1) & (result$Orig_Power==0 | (result$Inte_Power==0) | (result$LF_Power==0)),TRUE,FALSE)
  result$diff_FN_all<-ifelse((result$Pair_123==0) & (result$Orig_Power==1 & (result$Inte_Power==1) & (result$LF_Power==1)),TRUE,FALSE)
  
  index.1=which(complete.cases(temp), arr.ind=TRUE)
  
  index.2=which(result$diff_FN_all==TRUE, arr.ind=TRUE)
  
  index=intersect(index.1,index.2)
  
  
  example.summ<-temp[index[1],]
  example.summ.CI<-rcd.simuCI1[index[1],]
  
  example.data<-matrix(rcd.data[index[1],],nrow=5,byrow = T)
  example.mean<-rcd.mean[index[1],]
  example.sd<-rcd.std[index[1],]
  
  write.csv(example.data,file = "./example_FN_data.csv")
  
  # scatterplot of the samples
  
  value<-unlist(c(example.data[1,],example.data[2,],example.data[3,]))
  data<-data.frame(Group=c(rep("US",10),rep("EU",10),rep("T",10)),value=value,Lot=rep(seq(1:10),3))
  data$Group<-factor(data$Group,levels = c("US","EU","T"))
  
  plot <- ggplot(data=data, aes(x=Lot, y=value, group = Group, colour=Group,shape=Group))+
    geom_point(size=3)+
    theme(legend.position="right")+
    labs(x="Lot",y="Value")+
    scale_x_continuous(breaks=data$Lot)+
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  
  ggsave(plot,filename = "./Scatterplot_FN.png", height=5, width=7, units='in', dpi=600)
  

  
  pairwise.CI <- function(n1,mu1,n2,mu2,refstd,alpha){
    mar <- 1.5*refstd
    Za <- qnorm(1-alpha)
    sig.low <- ( (mu2-mu1)+Za*refstd*sqrt(1/n1+1/n2) < mar )
    sig.up <- ( (mu2-mu1)-Za*refstd*sqrt(1/n1+1/n2) > (-mar) )
    U=(mu2-mu1)+Za*refstd*sqrt(1/n1+1/n2)
    L=(mu2-mu1)-Za*refstd*sqrt(1/n1+1/n2) 
    sig <- sig.low*sig.up
    result<-c(L,U,sig,mar)
    return(result)
  }
  
  pairwise.12<-pairwise.CI(n1=10,mu1=example.mean[1],n2=10,mu2=example.mean[2],refstd=example.sd[4],alpha=0.05)
  pairwise.13<-pairwise.CI(n1=10,mu1=example.mean[1],n2=10,mu2=example.mean[3],refstd=example.sd[4],alpha=0.05)
  pairwise.23<-pairwise.CI(n1=10,mu1=example.mean[2],n2=10,mu2=example.mean[3],refstd=example.sd[5],alpha=0.05)

  
  Pairwise.result<-cbind(c(example.mean[1]-example.mean[2],example.mean[1]-example.mean[3],example.mean[2]-example.mean[3]),
                         (matrix(c(pairwise.12,pairwise.13,pairwise.23),byrow = T,nrow=3)))
  
  Simu.CI<-matrix(c(example.summ$Orig_CI1_D*example.sd[4]*1.5,example.summ$Orig_CI2_D*example.sd[4]*1.5,example.summ$Orig_FP,
             example.summ$Inte_CI1_D*example.sd[4]*1.5,example.summ$Inte_CI2_D*example.sd[4]*1.5,example.summ$Inte_FP,
             example.summ$LF_CI1_D*example.sd[4]*1.5,example.summ$LF_CI2_D*example.sd[4]*1.5,example.summ$LF_FP),byrow = T,nrow = 3)
  row.names(Simu.CI)<-c("Original Version","Integrated Version","Least Favaroble Version")
  colnames(Simu.CI)<-c("Type 1 Confidence Interval","Type 2 Confidence Interval","Fiducial probability")
  
  write.csv(Pairwise.result,file = "./example_FN_pair.csv")
  write.csv(Simu.CI,file = "./example_FN_simu.csv")
  