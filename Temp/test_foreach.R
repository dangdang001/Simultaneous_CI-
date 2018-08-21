# install.packages('doParallel')
library(doParallel)

##################################################################
##################################################################
# install.packages("R2Cuba")


rm(list = ls())
setwd("M:/")
library(R2Cuba)
library(parallel)

time.start<-proc.time()

# Calculate the number of cores
no_cores <- detectCores() - 1

mu.R1 <- 0
mu.R2 <- 2
mu.T <- 1
rsd21 <- 1.25  #single 0.5  2   1.25 0.8 # 0.8  0.8 1.25 1.25 # 0.8  0.8 1.25 1.25
rsdT1 <- 2     #single 1.25 0.8 0.5  2   # 1.25 0.8 0.8 1.25  # 0.5  2   0.5  2
# sd.R1 <- 4
# sd.R2 <- sd.R1*rsd21
# sd.T <- sd.R1*rsdT1
n <- 10
alpha <- 0.05
p.level <- 1 - 2*alpha

NDIM <- 3
NCOMP <- 1
p.int <- 0.999
p.int1 <- 0.005

sample.1product <- function(n,mu,std){
  return(rnorm(n,mean=mu,sd=std))
}

biosim.test <- function(n1,mu1,n2,mu2,refstd,alpha){
  mar <- 1.5*refstd
  Za <- qnorm(1-alpha)
  sig.low <- ( (mu2-mu1)+Za*refstd*sqrt(1/n1+1/n2) < mar )
  sig.up <- ( (mu2-mu1)-Za*refstd*sqrt(1/n1+1/n2) > (-mar) )
  sig <- sig.low*sig.up
  return(sig)
}

sd.all <- c((max(mu.R1,mu.R2,mu.T)-min(mu.R1,mu.R2,mu.T))/1.5,seq(2,12,2))
# sd.all <- c(seq(8,10,2))
# sd.all <- c((max(mu.R1,mu.R2,mu.T)-min(mu.R1,mu.R2,mu.T))/1.5,2)
# sd.all <- 8
ksce <- length(sd.all)
repet <- 1000

Arcd.data <- matrix(NA,nrow=repet,ncol=n*5*ksce)
Arcd.mean <- matrix(NA,nrow=repet,ncol=3*ksce)
Arcd.std <- matrix(NA,nrow=repet,ncol=(3+2)*ksce)
Arcd.sig <- matrix(NA,nrow=repet,ncol=3*ksce)
Arcd.simuCI1 <- matrix(NA,nrow=repet,ncol=9*ksce)


summ.sig <- matrix(NA,nrow=ksce,ncol=5)
summ.simuCI1 <- matrix(NA,nrow=ksce,ncol=18)


for(k in 1:ksce){
  sd.R1 <- sd.all[k]
  sd.R2 <- sd.R1*rsd21 # ?
  sd.T <- sd.R1*rsdT1  # ?
  
  true12 <- abs(mu.R1-mu.R2)/sd.R1
  true1T <- abs(mu.R1-mu.T)/sd.R1
  true2T <- abs(mu.T-mu.R2)/sd.R1
  true2T.2 <- abs(mu.T-mu.R2)/sd.R2
  
  rcd.data <- matrix(NA,nrow=repet,ncol=n*5)
  rcd.mean <- matrix(NA,nrow=repet,ncol=3)
  rcd.std <- matrix(NA,nrow=repet,ncol=3+2)
  rcd.sig <- matrix(NA,nrow=repet,ncol=3)
  rcd.simuCI1 <- matrix(NA,nrow=repet,ncol=9)
  rcd.simuCI2 <- matrix(NA,nrow=repet,ncol=9)
  rcd.simuCI3 <- matrix(NA,nrow=repet,ncol=9)
  
  simu.one.PAR.mpi <- function(i=1,k,rsd21,rsdT1,mu.R1,mu.R2,mu.T,sd.all,n){
    
    library(R2Cuba)
    
    ksce <- length(sd.all)
    
    sd.R1 <- sd.all[k]
    sd.R2 <- sd.R1*rsd21
    sd.T <- sd.R1*rsdT1
    
    alpha <- 0.05
    p.level <- 1 - 2*alpha
    NDIM <- 3
    NCOMP <- 1
    p.int <- 0.999
    p.int1 <- 0.005
    
    sample.1product <- function(n,mu,std){
      return(rnorm(n,mean=mu,sd=std))
    }
    
    biosim.test <- function(n1,mu1,n2,mu2,refstd,alpha){
      mar <- 1.5*refstd
      Za <- qnorm(1-alpha)
      sig.low <- ( (mu2-mu1)+Za*refstd*sqrt(1/n1+1/n2) < mar )
      sig.up <- ( (mu2-mu1)-Za*refstd*sqrt(1/n1+1/n2) > (-mar) )
      sig <- sig.low*sig.up
      return(sig)
    }
    
    data.R1refv <- sample.1product(n,mu.R1,sd.R1) # data.R1refv and data.R1?
    data.R2refv <- sample.1product(n,mu.R2,sd.R2)
    data.R1 <- sample.1product(n,mu.R1,sd.R1)
    data.R2 <- sample.1product(n,mu.R2,sd.R2)
    data.T <- sample.1product(n,mu.T,sd.T)
    
    rcd.data.tmp <- c(data.R1, data.R2, data.T, data.R1refv, data.R2refv)
    
    sd.ref <- sd(data.R1refv)
    sd.ref2 <- sd(data.R2refv)
    sd.ref.lf <- sqrt((n-1)*sd.ref^2/qchisq(1 - alpha,df=n-1))
    sd.ref.lf1 <- sqrt((n-1)*sd.ref^2/qchisq(1-(1-sqrt(1-2*alpha))/2,df=n-1)) # 2*alpha?
    sd.ref.lf2 <- sqrt((n-1)*sd.ref2^2/qchisq(1-(1-sqrt(1-2*alpha))/2,df=n-1))
    mean.data <- c(mean(data.R1),mean(data.R2),mean(data.T))
    sd.data <- c(sd(data.R1),sd(data.R2),sd(data.T))
    
    intmar.l <- mean.data - qnorm(p.int)*sd.ref/sqrt(n)
    intmar.u <- mean.data + qnorm(p.int)*sd.ref/sqrt(n)
    intmarv <- sqrt((n-1)*sd.ref^2/qchisq(1-p.int,df=n-1))
    tmp1 <- sqrt((n-1)*sd.ref^2/qchisq(p.int1,df=n-1))
    intmar.l.int  <- mean.data - qnorm(p.int)*tmp1/sqrt(n) # tmp1 seems not correct?
    intmar.u.int  <- mean.data + qnorm(p.int)*tmp1/sqrt(n)
    intmar.l.lf <- mean.data - qnorm(p.int)*sd.ref.lf/sqrt(n)
    intmar.u.lf <- mean.data + qnorm(p.int)*sd.ref.lf/sqrt(n)
    
    intmar.l.ie <- mean.data - qt(p.int,df=n-1)*sd.data/sqrt(n)
    intmar.u.ie <- mean.data + qt(p.int,df=n-1)*sd.data/sqrt(n)
    
    intmar.l2 <- mean.data - qnorm(p.int)*c(sd.ref,sd.ref2,sd.ref)/sqrt(n)
    intmar.u2 <- mean.data + qnorm(p.int)*c(sd.ref,sd.ref2,sd.ref)/sqrt(n)
    intmar.l2.lf <- mean.data - qnorm(p.int)*c(sd.ref.lf1,sd.ref.lf2,sd.ref.lf)/sqrt(n)
    intmar.u2.lf <- mean.data + qnorm(p.int)*c(sd.ref.lf1,sd.ref.lf2,sd.ref.lf)/sqrt(n)
    
    intmarv1t <- sqrt((n-1)*sd.ref^2/qchisq(1-p.int,df=n-1))
    intmarv2t <- sqrt((n-1)*sd.ref2^2/qchisq(1-p.int,df=n-1))
    tmp2 <- sqrt((n-1)*sd.ref2^2/qchisq(p.int1,df=n-1))
    intmar.l.int2  <- mean.data - qnorm(p.int)*c(tmp1,tmp2,tmp1)/sqrt(n)
    intmar.u.int2  <- mean.data + qnorm(p.int)*c(tmp1,tmp2,tmp1)/sqrt(n)
    
    rcd.mean.tmp <- mean.data
    rcd.std.tmp <- c(sd.data, sd.ref, sd.ref2)
    rcd.sig.tmp <- c(biosim.test(n,mean.data[1],n,mean.data[2],sd.ref,alpha),
                     biosim.test(n,mean.data[1],n,mean.data[3],sd.ref,alpha),
                     biosim.test(n,mean.data[2],n,mean.data[3],sd.ref2,alpha))
    
    ### fiducial probability, simultanuous CI
    #
    ## with the assumption of equal variances--------------------------------------- 
    #
    fun.eFP1 <- function(x){
      refstd <- 1.5*sd.ref
      indt <- ( abs(x[1]-x[2])<=refstd )*( abs(x[1]-x[3])<=refstd )*( abs(x[2]-x[3])<=refstd )
      # mean and sd correct?
      fxyz <- dnorm(x[1], mean = mean.data[1], sd = sd.ref/sqrt(n))*dnorm(x[2], mean = mean.data[2], sd = sd.ref/sqrt(n))*dnorm(x[3], mean = mean.data[3], sd = sd.ref/sqrt(n))*indt
      return(fxyz)
    }
    eFP1.res <- cuhre(NDIM, NCOMP, fun.eFP1, lower=intmar.l, upper=intmar.u, flags= list(verbose=0, final=0))
    rcd.simuCI1.tmp1 <- eFP1.res$value
    
    rcd.simuCI1.tmp2 <- NA
    # Why the if condition?
    if(eFP1.res$value>p.level){
      eFP1.int1 <- function(delta){
        fun.eFP1t1 <- function(x){
          refstd <- 1.5*sd.ref
          indt <- ( abs(x[1]-x[2])<=refstd )*( abs(x[1]-x[3])<=delta )*( abs(x[2]-x[3])<=refstd )
          fxyz <- dnorm(x[1], mean = mean.data[1], sd = sd.ref/sqrt(n))*dnorm(x[2], mean = mean.data[2], sd = sd.ref/sqrt(n))*dnorm(x[3], mean = mean.data[3], sd = sd.ref/sqrt(n))*indt
          return(fxyz)
        }
        return(cuhre(NDIM, NCOMP, fun.eFP1t1, lower=intmar.l, upper=intmar.u, flags= list(verbose=0, final=0))$value)
      }
      # uniroot function
      rcd.simuCI1.tmp2 <- uniroot(function(delta) eFP1.int1(delta)-p.level, lower = 0, upper = 1.5*sd.ref)$root
    }
    eFP1.int2 <- function(delta){
      fun.eFP1t2 <- function(x){
        indt <- ( abs(x[1]-x[2])<=delta )*( abs(x[1]-x[3])<=delta )*( abs(x[2]-x[3])<=delta )
        fxyz <- dnorm(x[1], mean = mean.data[1], sd = sd.ref/sqrt(n))*dnorm(x[2], mean = mean.data[2], sd = sd.ref/sqrt(n))*dnorm(x[3], mean = mean.data[3], sd = sd.ref/sqrt(n))*indt
        return(fxyz)
      }
      return(cuhre(NDIM, NCOMP, fun.eFP1t2, lower=intmar.l, upper=intmar.u, flags= list(verbose=0, final=0))$value)
    }
    rcd.simuCI1.tmp3 <- uniroot(function(delta) eFP1.int2(delta)-p.level, lower = 0, upper = 20*sd.ref)$root
    
    # the integrated version
    fun.eFP1.int <- function(x){
      refstd <- 1.5*x[4]
      indt <- ( abs(x[1]-x[2])<=refstd )*( abs(x[1]-x[3])<=refstd )*( abs(x[2]-x[3])<=refstd )
      # derive the f.d of sd(r1)?
      fxyzu <- dnorm(x[1], mean = mean.data[1], sd = x[4]/sqrt(n))*dnorm(x[2], mean = mean.data[2], sd = x[4]/sqrt(n))*dnorm(x[3], mean = mean.data[3], sd = x[4]/sqrt(n))*indt* (2*x[4]^(-3)*(n-1)*sd.ref^2)*dchisq((n-1)*sd.ref^2/(x[4]^2), df=n-1)
      return(fxyzu)
    }
    eFP1.res.mod <- cuhre(NDIM+1, NCOMP, fun.eFP1.int, lower=c(intmar.l.int,0), upper=c(intmar.u.int,intmarv), flags= list(verbose=0, final=0))
    rcd.simuCI1.tmp4 <- eFP1.res.mod$value
    
    rcd.simuCI1.tmp5 <- NA
    if(eFP1.res.mod$value>p.level){
      eFP1.int1.mod <- function(delta){
        fun.eFP1t1.mod <- function(x){
          refstd <- 1.5*x[4]
          indt <- ( abs(x[1]-x[2])<=refstd )*( abs(x[1]-x[3])<=delta*refstd )*( abs(x[2]-x[3])<=refstd )
          fxyzu <- dnorm(x[1], mean = mean.data[1], sd = x[4]/sqrt(n))*dnorm(x[2], mean = mean.data[2], sd = x[4]/sqrt(n))*dnorm(x[3], mean = mean.data[3], sd = x[4]/sqrt(n))*indt* (2*x[4]^(-3)*(n-1)*sd.ref^2)*dchisq((n-1)*sd.ref^2/(x[4]^2), df=n-1)
          return(fxyzu)
        }
        return(cuhre(NDIM+1, NCOMP, fun.eFP1t1.mod, lower=c(intmar.l.int,0), upper=c(intmar.u.int,intmarv), flags= list(verbose=0, final=0))$value)
      }
      rcd.simuCI1.tmp5 <- uniroot(function(delta) eFP1.int1.mod(delta)-p.level, lower = 0, upper = 1)$root
    }
    eFP1.int2.mod <- function(delta){
      fun.eFP1t2.mod <- function(x){
        refstd <- 1.5*x[4]
        indt <- ( abs(x[1]-x[2])<=delta*refstd )*( abs(x[1]-x[3])<=delta*refstd )*( abs(x[2]-x[3])<=delta*refstd )
        fxyzu <- dnorm(x[1], mean = mean.data[1], sd = x[4]/sqrt(n))*dnorm(x[2], mean = mean.data[2], sd = x[4]/sqrt(n))*dnorm(x[3], mean = mean.data[3], sd = x[4]/sqrt(n))*indt* (2*x[4]^(-3)*(n-1)*sd.ref^2)*dchisq((n-1)*sd.ref^2/(x[4]^2), df=n-1)
        return(fxyzu)
      }
      return(cuhre(NDIM+1, NCOMP, fun.eFP1t2.mod, lower=c(intmar.l.int,0), upper=c(intmar.u.int,intmarv), flags= list(verbose=0, final=0))$value)
    }
    rcd.simuCI1.tmp6 <- uniroot(function(delta) eFP1.int2.mod(delta)-p.level, lower = 0, upper = 20)$root
    
    # the least favorable version
    fun.eFP1.lf <- function(x){
      refstd <- 1.5*sd.ref.lf
      indt <- ( abs(x[1]-x[2])<=refstd )*( abs(x[1]-x[3])<=refstd )*( abs(x[2]-x[3])<=refstd )
      fxyz <- dnorm(x[1], mean = mean.data[1], sd = sd.ref.lf/sqrt(n))*dnorm(x[2], mean = mean.data[2], sd = sd.ref.lf/sqrt(n))*dnorm(x[3], mean = mean.data[3], sd = sd.ref.lf/sqrt(n))*indt
      return(fxyz)
    }
    eFP1.res.lf <- cuhre(NDIM, NCOMP, fun.eFP1.lf, lower=intmar.l.lf, upper=intmar.u.lf, flags= list(verbose=0, final=0))
    rcd.simuCI1.tmp7 <- eFP1.res.lf$value
    
    rcd.simuCI1.tmp8 <- NA
    if(eFP1.res.lf$value>p.level){
      eFP1.int1.lf <- function(delta){
        fun.eFP1t1.lf <- function(x){
          refstd <- 1.5*sd.ref.lf
          indt <- ( abs(x[1]-x[2])<=refstd )*( abs(x[1]-x[3])<=delta )*( abs(x[2]-x[3])<=refstd )
          fxyz <- dnorm(x[1], mean = mean.data[1], sd = sd.ref.lf/sqrt(n))*dnorm(x[2], mean = mean.data[2], sd = sd.ref.lf/sqrt(n))*dnorm(x[3], mean = mean.data[3], sd = sd.ref.lf/sqrt(n))*indt
          return(fxyz)
        }
        return(cuhre(NDIM, NCOMP, fun.eFP1t1.lf, lower=intmar.l.lf, upper=intmar.u.lf, flags= list(verbose=0, final=0))$value)
      }
      rcd.simuCI1.tmp8 <- uniroot(function(delta) eFP1.int1.lf(delta)-p.level, lower = 0, upper = 1.5*sd.ref.lf)$root
    }
    eFP1.int2.lf <- function(delta){
      fun.eFP1t2.lf <- function(x){
        indt <- ( abs(x[1]-x[2])<=delta )*( abs(x[1]-x[3])<=delta )*( abs(x[2]-x[3])<=delta )
        fxyz <- dnorm(x[1], mean = mean.data[1], sd = sd.ref.lf/sqrt(n))*dnorm(x[2], mean = mean.data[2], sd = sd.ref.lf/sqrt(n))*dnorm(x[3], mean = mean.data[3], sd = sd.ref.lf/sqrt(n))*indt
        return(fxyz)
      }
      return(cuhre(NDIM, NCOMP, fun.eFP1t2.lf, lower=intmar.l.lf, upper=intmar.u.lf, flags= list(verbose=0, final=0))$value)
    }
    rcd.simuCI1.tmp9 <- uniroot(function(delta) eFP1.int2.lf(delta)-p.level, lower = 0, upper = 20*sd.ref.lf)$root
    
    
    # ---------------- combine data and results: n*5, 3+5+3, 9 ---------------
    results <- c(rcd.data.tmp,
                 rcd.mean.tmp,
                 rcd.std.tmp,
                 rcd.sig.tmp,
                 rcd.simuCI1.tmp1,
                 rcd.simuCI1.tmp2,
                 rcd.simuCI1.tmp3,
                 rcd.simuCI1.tmp4,
                 rcd.simuCI1.tmp5,
                 rcd.simuCI1.tmp6,
                 rcd.simuCI1.tmp7,
                 rcd.simuCI1.tmp8,
                 rcd.simuCI1.tmp9
                 
    )
    
    return(results)
  }
  
  # system.time(out.tmp.cont <- sapply(1:repet, simu.one.PAR.mpi, k=k,mu.R1=mu.R1,mu.R2=mu.R2,mu.T=mu.T,rsd21=rsd21,rsdT1=rsdT1,sd.all=sd.all,n=n))
  
  cl <- makeCluster(no_cores)
  
  clusterExport(cl, c("k","repet","rsd21","rsdT1","mu.R1","mu.R2","mu.T","sd.all","n"), envir=environment())
  
  clusterEvalQ(cl, {library(R2Cuba)})
  
  system.time(out.tmp.cont <- parSapply(cl, 1:repet, simu.one.PAR.mpi, k=k,rsd21=rsd21,rsdT1=rsdT1,mu.R1=mu.R1,mu.R2=mu.R2,mu.T=mu.T,sd.all=sd.all,n=n))
  
  stopCluster(cl)
  
  restmp <- t(out.tmp.cont)
  
  rcd.data <- restmp[,1:(n*5)]
  rcd.mean <- restmp[,(n*5+1):(n*5+3)]
  rcd.std <- restmp[,(n*5+3+1):(n*5+3+5)]
  rcd.sig <- restmp[,(n*5+3+5+1):(n*5+3+5+3)]
  rcd.simuCI1 <- restmp[,(n*5+3+5+3+1):(n*5+3+5+3+9)]
  
  
  Arcd.data[,(n*5*(k-1)+1):(n*5*k)] <- rcd.data
  Arcd.mean[,(3*(k-1)+1):(3*k)] <- rcd.mean
  Arcd.std[,((3+2)*(k-1)+1):((3+2)*k)] <- rcd.std
  Arcd.sig[,(3*(k-1)+1):(3*k)] <- rcd.sig
  Arcd.simuCI1[,(9*(k-1)+1):(9*k)] <- rcd.simuCI1
  
  
  summ.sig[k,] <- c(apply(rcd.sig,2,mean), sum(apply(rcd.sig[,2:3],1,prod))/repet, sum(apply(rcd.sig,1,prod))/repet)
  
  summ.simuCI1[k,] <- c(mean(rcd.simuCI1[,1]),
                        mean(rcd.simuCI1[,1]>=p.level),
                        mean(rcd.simuCI1[,2]/rcd.std[,4], na.rm=T),
                        mean(rcd.simuCI1[,2]/rcd.std[,4]>=true1T, na.rm=T),
                        mean(rcd.simuCI1[,3]/rcd.std[,4]),
                        mean(rcd.simuCI1[,3]/rcd.std[,4]>=max(true1T,true2T,true12)),
                        mean(rcd.simuCI1[,4]),
                        mean(rcd.simuCI1[,4]>=p.level),
                        mean(rcd.simuCI1[,5]*1.5, na.rm=T),
                        mean(rcd.simuCI1[,5]*1.5>=true1T, na.rm=T),
                        mean(rcd.simuCI1[,6]*1.5),
                        mean(rcd.simuCI1[,6]*1.5>=max(true1T,true2T,true12)),
                        mean(rcd.simuCI1[,7]),
                        mean(rcd.simuCI1[,7]>=p.level),
                        mean(rcd.simuCI1[,8]/rcd.std[,4], na.rm=T),
                        mean(rcd.simuCI1[,8]/rcd.std[,4]>=true1T, na.rm=T),
                        mean(rcd.simuCI1[,9]/rcd.std[,4]),
                        mean(rcd.simuCI1[,9]/rcd.std[,4]>=max(true1T,true2T,true12))
  )
  
  
}

Arcd <- cbind(Arcd.data, Arcd.mean, Arcd.std, Arcd.sig, Arcd.simuCI1)
summ <- round(cbind(summ.sig, summ.simuCI1),3)


rundate <- c("0627")


write.csv(Arcd, paste("Arcd_n",n,"_",mu.R1,"_",mu.R2,"_",mu.T,"_equal_var_","_sdR1_",round(sd.all[1],2),"_",round(sd.all[ksce],2),"_rp",repet,"_",rundate,"_v2.csv",sep=""))
write.csv(summ, paste("summ_n",n,"_",mu.R1,"_",mu.R2,"_",mu.T,"_equal_var_","_sdR1_",round(sd.all[1],2),"_",round(sd.all[ksce],2),"_rp",repet,"_",rundate,"_v2.csv",sep=""))


time.end<-proc.time()

time.end-time.start