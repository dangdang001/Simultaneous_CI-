rm(list=ls())

library(ggplot2)
library(gridExtra)

setwd("//fda.gov/wodc/CDER/Users05/Donglei.Yin/Result")
#setwd("C:/Users/Donglei/Documents/2018 summer intern at FDA/Simu_CI")

filename='n10_0_2_1_equal_var__sdR1_1.33_12_rp1000_0627_v2'

simu_result<-read.csv(paste0("./summ_",filename, ".csv"),header=T)[,-1]

dim(simu_result)

colnames(simu_result)<-c("Pair_R1_R2","Pair_R1_T","Pair_R2_T","Pair_23","Pair_123",
                         "Orig_FP","Orig_Power","Orig_CI1_D","Orig_CI1_CR","Orig_CI2_D","Orig_CI2_CR",
                         "Inte_FP","Inte_Power","Inte_CI1_D","Inte_CI1_CR","Inte_CI2_D","Inte_CI2_CR",
                         "LF_FP","LF_Power","LF_CI1_D","LF_CI1_CR","LF_CI2_D","LF_CI2_CR")

simu_result$SD=c(1.33,2,4,6,8,10,12)

Power=simu_result[,c("SD","Pair_123","Orig_Power","Inte_Power","LF_Power")]
colnames(Power)=c("SD","Pairwise comparison","Original version","Integrated version","Least favorable version")

CR1=simu_result[,c("SD","Orig_CI1_CR","Inte_CI1_CR","LF_CI1_CR")]
colnames(CR1)=c("SD","Original version","Integrated version","Least favorable version")

CR2=simu_result[,c("SD","Orig_CI2_CR","Inte_CI2_CR","LF_CI2_CR")]
colnames(CR2)=c("SD","Original version","Integrated version","Least favorable version")

## Wide to long

Power.t <- reshape(Power, 
                   varying = c("Pairwise comparison","Original version","Integrated version","Least favorable version"), 
                   v.names = "Value",
                   timevar = "Method", 
                   times = c("Pairwise comparison","Original version","Integrated version","Least favorable version"), 
                   new.row.names = 1:1000,
                   direction = "long")

CR1.t <- reshape(CR1, 
                 varying = c("Original version","Integrated version","Least favorable version"), 
                 v.names = "Value",
                 timevar = "Method", 
                 times = c("Original version","Integrated version","Least favorable version"), 
                 new.row.names = 1:1000,
                 direction = "long")

CR2.t <- reshape(CR2, 
                 varying = c("Original version","Integrated version","Least favorable version"), 
                 v.names = "Value",
                 timevar = "Method", 
                 times = c("Original version","Integrated version","Least favorable version"), 
                 new.row.names = 1:1000,
                 direction = "long")

Power.t$Group<-"Power"
CR1.t$Group<-"Coverage Rate 1"
CR2.t$Group<-"Coverage Rate 2"


out<-rbind(Power.t,CR1.t,CR2.t)

out$Group_f = factor(out$Group, levels=c('Power','Coverage Rate 1','Coverage Rate 2'))
out$Method = factor(out$Method, levels=c('Original version','Integrated version','Least favorable version','Pairwise comparison'))

plot <- ggplot(data=out, aes(x=SD, y=Value, group = Method, colour=Method,shape=Method))+
  geom_point()+
  geom_line()+
  facet_grid(.~Group_f,scales="free_y")+
  theme(legend.position="bottom")+
  labs(x="SD",y="Value")


ggsave(paste0("./Plot_",filename,".png"), height=6, width=12, units='in', dpi=2000)

