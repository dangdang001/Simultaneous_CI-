rm(list=ls())

library(ggplot2)
library(gridExtra)

setwd("//fda.gov/wodc/CDER/Users05/Donglei.Yin/Result")

filename='n10_99_101_100_equal_var_sdR1_1.33_12_rp200_0628_v2'

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


out<-Reduce(function(x, y) merge(x, y, by=c("SD","Method","id")), list(Power.t,CR1.t,CR2.t))

power.plot <- ggplot(data=Power.t, aes(x=SD, y=Power, group = Method, colour=Method))+
  geom_point()+
  geom_line()+
  theme(legend.position="bottom")+
  labs(x="SD",y="Power")


CR1.plot <- ggplot(data=CR1.t, aes(x=SD, y=Coverage_Rate_1, group = Method, colour=Method))+
  geom_point()+
  geom_line()+
  theme(legend.position="bottom")+
  labs(x="SD",y="Coverage Rate (CI1)")

CR2.plot <- ggplot(data=CR2.t, aes(x=SD, y=Coverage_Rate_2, group = Method, colour=Method))+
  geom_point()+
  geom_line()+
  theme(legend.position="bottom")+
  labs(x="SD",y="Coverage Rate (CI2)")



g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(power.plot)


p <- grid.arrange(arrangeGrob(power.plot + theme(legend.position="none"),
                               CR1.plot + theme(legend.position="none"),
                               CR2.plot + theme(legend.position="none"),
                               nrow=1),
                   mylegend, nrow=2,heights=c(10, 1))

ggsave(plot=p,paste0("./Plot_",filename,".png"), height=6, width=12, units='in', dpi=2000)
