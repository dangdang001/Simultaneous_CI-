# Purpose: Power and CR plots comparing pairwise comparison with simulatanous confidence interval approach.
# Name: Donglei Yin
# Date: 2018.07

# Note:
# change the filename for different mu settings, here mu_R1=99, mu_R2=101, mu_T=100
# remember to change the sd_margin and simu_result$SD also

rm(list=ls())

library(ggplot2)

setwd("//fda.gov/wodc/CDER/Users05/Donglei.Yin/Result")


filename='n10_99_101_100_equal_var_sdR1_1.33_14_rp1000'

sd_margin=1.33

simu_result<-read.csv(paste0("./summ_",filename, ".csv"),header=T)[,-1]

dim(simu_result)

colnames(simu_result)<-c("Pair_R1_R2","Pair_R1_T","Pair_R2_T","Pair_23","Pair_123",
                         "Orig_FP","Orig_Power","Orig_CI1_D","Orig_CI1_CR","Orig_CI2_D","Orig_CI2_CR",
                         "Inte_FP","Inte_Power","Inte_CI1_D","Inte_CI1_CR","Inte_CI2_D","Inte_CI2_CR",
                         "LF_FP","LF_Power","LF_CI1_D","LF_CI1_CR","LF_CI2_D","LF_CI2_CR")

simu_result$SD=c(1.33,2,4,6,8,10,12,14)

Power=simu_result[,c("SD","Pair_123","Orig_Power","Inte_Power","LF_Power")]
colnames(Power)=c("SD","Pairwise comparison","Original version","Integrated version","Least favorable version")

CR1=simu_result[,c("SD","Orig_CI1_CR","Inte_CI1_CR","LF_CI1_CR")]
colnames(CR1)=c("SD","Original version","Integrated version","Least favorable version")

CR2=simu_result[,c("SD","Orig_CI2_CR","Inte_CI2_CR","LF_CI2_CR")]
colnames(CR2)=c("SD","Original version","Integrated version","Least favorable version")

## Wide to long

Power.t <- reshape(Power, 
                   varying = c("Pairwise comparison","Original version","Integrated version","Least favorable version"), 
                   v.names = "Power",
                   timevar = "Method", 
                   times = c("Pairwise comparison","Original version","Integrated version","Least favorable version"), 
                   new.row.names = 1:1000,
                   direction = "long")

CR1.t <- reshape(CR1, 
                 varying = c("Original version","Integrated version","Least favorable version"), 
                 v.names = "Coverage_Rate_1",
                 timevar = "Method", 
                 times = c("Original version","Integrated version","Least favorable version"), 
                 new.row.names = 1:1000,
                 direction = "long")

CR2.t <- reshape(CR2, 
                 varying = c("Original version","Integrated version","Least favorable version"), 
                 v.names = "Coverage_Rate_2",
                 timevar = "Method", 
                 times = c("Original version","Integrated version","Least favorable version"), 
                 new.row.names = 1:1000,
                 direction = "long")


Power.t$Method = factor(Power.t$Method, levels=c('Original version','Integrated version','Least favorable version','Pairwise comparison'))
CR1.t$Method = factor(CR1.t$Method, levels=c('Original version','Integrated version','Least favorable version','Pairwise comparison'))
CR2.t$Method = factor(CR2.t$Method, levels=c('Original version','Integrated version','Least favorable version','Pairwise comparison'))


power.plot <- ggplot(data=Power.t, aes(x=SD, y=Power, group = Method, colour=Method,linetype=Method))+
  geom_point(aes(shape=Method),size=1.5)+
  geom_line(aes(linetype=Method))+
  theme(legend.position="bottom")+
  ylim(0,1)+
  labs(x="SD",y="Power")+
  # labs(x="SD",y=paste("Power/ Type )", as.roman(1), " Error"))+
  scale_x_continuous(breaks=simu_result$SD, labels=simu_result$SD)+
  geom_vline(xintercept = sd_margin, linetype="dotted", 
               color = "red")

ggsave(paste0("./Power_",filename,".png"), height=5, width=7, units='in', dpi=600)


CR1.plot <- ggplot(data=CR1.t, aes(x=SD, y=Coverage_Rate_1, group = Method, colour=Method))+
  geom_point()+
  geom_line()+
  theme(legend.position="bottom")+
  ylim(0,1)+
  labs(x="SD",y="Coverage Rate (CI1)")+
  scale_x_continuous(breaks=simu_result$SD, labels=simu_result$SD)+
  geom_vline(xintercept = sd_margin, linetype="dotted", 
             color = "red")

ggsave(paste0("./CR1_",filename,".png"), height=5, width=7, units='in', dpi=600)

CR2.plot <- ggplot(data=CR2.t, aes(x=SD, y=Coverage_Rate_2, group = Method, colour=Method))+
  geom_point()+
  geom_line()+
  theme(legend.position="bottom")+
  ylim(0,1)+
  labs(x="SD",y="Coverage Rate (CI2)")+
  scale_x_continuous(breaks=simu_result$SD, labels=simu_result$SD)+
  geom_vline(xintercept = sd_margin, linetype="dotted", 
             color = "red")

ggsave(paste0("./CR2_",filename,".png"), height=5, width=7,units='in', dpi=600)

