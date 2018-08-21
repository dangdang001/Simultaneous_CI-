################ Figure 5_2
rm(list=ls())
setwd("C:/Sophie/project/Rcode")
library(parallel)
#require(utils)
require(ggplot2)
library(dplyr)
library(parallel)
library(RColorBrewer) 
source("./all3_rej_pow.R")
source("multiplot.R")
load("data_5_1.RData")

eq_power=function(p,n,data){##### p=c(pt,pr), n=c(nt,nr)
  return(sum(dbinom((data$pt*n[1]),n[1],p[1])*dbinom((data$pr*n[2]),n[2],p[2])))
}


ns=c(100,200,300)
rr_margins=c("2.0","2.5","3.0","One Step")
pts=0:80/100
prs=0:60/100
ptrs=expand.grid(pts,prs)
nnode=3
cl = makeCluster(nnode)
power=rep(0,5)
for(i in 1:length(rr_margins)){
  for(j in 1:length(ns)){
    temp=subset(result,each_n==ns[j]&rr_margin==rr_margins[i]&TR==1,select=c("pt","pr"))
    p_tr=parApply(cl,ptrs,1,eq_power,n=c(ns[j],ns[j]),data=temp)
    p_temp=data.frame(pt=ptrs[,1],pr=ptrs[,2],rr_margin=rr_margins[i],each_n=ns[j],power=p_tr)
    power=rbind(power,p_temp)
  }
}
power=power[-1,]
stopCluster(cl)


brks <- cut(power$power,breaks=c(-0.01,0.20,0.40,0.60,0.80,0.90,1.01))
brks <- gsub(","," - ",brks,fixed=TRUE)
brks=gsub("-0.01","0.0",brks,fixed=TRUE)
brks=gsub("1.01","1.0",brks,fixed=TRUE)
power$brks <- gsub("\\(|\\]","",brks)  # reformat guide labels
ggplot(power,aes(pr,pt)) + 
  geom_tile(aes(fill=brks))+
  scale_fill_manual("Pass Rate",values=brewer.pal(6,"YlOrRd"))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+facet_grid(each_n ~ rr_margin)+
  xlab(expression(p[A]))+ylab(expression(p[B]))

save(power,file="data_5_2.RData")

