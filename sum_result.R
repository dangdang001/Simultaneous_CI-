rm(list=ls())

library(ggplot2)
library(gridExtra)

#setwd("//fda.gov/wodc/CDER/Users05/Donglei.Yin/Result")
setwd("C:/Users/Donglei/Documents/2018 summer intern at FDA/Simultaneous_CI--master")

#filename='n10_99_100_101_equal_var_sdR1_1.33_14_rp1000'
filename='n10_99_101_100_equal_var_sdR1_1.33_14_rp1000'



simu_result<-read.csv(paste0("./summ_",filename, ".csv"),header=T)[,-1]

dim(simu_result)

colnames(simu_result)<-c("Pair_R1_R2","Pair_R1_T","Pair_R2_T","Pair_23","Pair_123",
                         "Orig_FP","Orig_Power","Orig_CI1_D","Orig_CI1_CR","Orig_CI2_D","Orig_CI2_CR",
                         "Inte_FP","Inte_Power","Inte_CI1_D","Inte_CI1_CR","Inte_CI2_D","Inte_CI2_CR",
                         "LF_FP","LF_Power","LF_CI1_D","LF_CI1_CR","LF_CI2_D","LF_CI2_CR")

simu_result$SD=c(1.33,2,4,6,8,10,12,14)

simu_out=simu_result[,c("Pair_R1_R2","Pair_R1_T","Pair_R2_T","Pair_123",
                        "Orig_Power","Orig_CI1_CR","Orig_CI2_CR",
                        "Inte_Power","Inte_CI1_CR","Inte_CI2_CR",
                        "LF_Power","LF_CI1_CR","LF_CI2_CR")]


write.csv(simu_out,file = "./table2_99_101_100.csv")
