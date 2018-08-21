rm(list=ls())

library(ggplot2)
library(gridExtra)

setwd("//fda.gov/wodc/CDER/Users05/Donglei.Yin/Result")
#setwd("C:/Users/Donglei/Documents/2018 summer intern at FDA/Simu_CI")

filename='n10_0_2_1_equal_var__sdR1_1.33_12_rp1000_0627_v2'

simu_result<-read.csv(paste0("./Arcd_",filename, ".csv"),header=T)[,-1]

dim(simu_result)

temp<-simu_result[,482:490]

sum(is.na(temp[,2]))
sum(is.na(temp[,3]))

sum(is.na(temp[,5]))
sum(is.na(temp[,6]))

sum(is.na(temp[,8]))
sum(is.na(temp[,9]))

