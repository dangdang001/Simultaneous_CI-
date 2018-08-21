# Purpose: Find the cases where pairwise comparison fails but simulatanous confidence interval approach succeed.
# Name: Donglei Yin
# Date: 2018.07

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
alpha <- 0.05
p.level <- 1 - 2*alpha

out=c()


for (k in 1:length(SD)){
  

true12 <- abs(mu.R1-mu.R2)/SD[k]
true1T <- abs(mu.R1-mu.T)/SD[k]
true2T <- abs(mu.T-mu.R2)/SD[k]
true2T.2 <- abs(mu.T-mu.R2)/SD[k]

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

rcd.sig.pairwise<-apply(rcd.sig,1,prod)

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

colnames(temp)<-c("Pair_123",
                         "Orig_FP","Orig_Power","Orig_CI1_D","Orig_CI1_CR","Orig_CI2_D","Orig_CI2_CR",
                         "Inte_FP","Inte_Power","Inte_CI1_D","Inte_CI1_CR","Inte_CI2_D","Inte_CI2_CR",
                         "LF_FP","LF_Power","LF_CI1_D","LF_CI1_CR","LF_CI2_D","LF_CI2_CR")
result=as.data.frame(temp[,c("Pair_123","Orig_Power","Inte_Power","LF_Power")])

result$diff_FN<-ifelse((result$Pair_123==0) & (result$Orig_Power==1 | (result$Inte_Power==1) | (result$LF_Power==1)),TRUE,FALSE)
result$diff_FP<-ifelse((result$Pair_123==1) & (result$Orig_Power==0 | (result$Inte_Power==0) | (result$LF_Power==0)),TRUE,FALSE)
result$diff_FN_all<-ifelse((result$Pair_123==0) & (result$Orig_Power==1 & (result$Inte_Power==1) & (result$LF_Power==1)),TRUE,FALSE)

out=c(out,sum(result$diff_FN),sum(result$diff_FP),sum(result$diff_FN_all))

}

matrix(out,ncol=3,byrow = TRUE)


