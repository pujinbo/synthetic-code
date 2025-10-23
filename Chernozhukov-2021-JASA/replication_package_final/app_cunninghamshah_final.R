###########################################################################################################
# Application: Cunningham and Shah (2018, Review of Economic Studies)
# Paper: An exact and robust conformal inference approach for counterfactual and synthetic controls
# Authors: V. Chernozhukov, K. Wuthrich, Y. Zhu
# DISCLAIMER: This software is provided "as is" without warranty of any kind, expressed or implied. 
# Questions/error reports: kwuthrich@ucsd.edu
###########################################################################################################

###########################################################################################################
# Preliminaries
###########################################################################################################

rm(list = ls())

setwd("/Users/kasparwuthrich/Dropbox/research/SC/SC with Victor and Yinchu/Conformal Paper/Code/replication_package_final")

set.seed(12345)

### Packages
library(xtable)
library(plotrix)

### Functions
source("functions_conformal_final.R")

###########################################################################################################
# Analysis
###########################################################################################################

### Specify norm for test statistic

q_norm <- 1 # L_1 norm

### Data (obtained from Scott Cunningham and Manisha Shah)

logfemrate <- read.delim("logfemrate.txt", header=T)

# First column: RI; columns 2-51: controls

Y1go <- as.matrix(logfemrate[,1])
Y0go <- as.matrix(logfemrate[,2:ncol(logfemrate)])

# Check I: comparison to Table 1 in Cunningham&Shah (2018)
# data19992009 <- as.matrix(logfemrate[15:25,])
# c(mean(data19992009),sd(data19992009))
# dim(data19992009)
# Check II: comparison to Table 9 in Cunningham&Shah (2018)
# cbind(1985:2009,Y1go)

# Time periods

T0go  <- 19
T1go  <- 6
T01go <- T0go+T1go

# Raw data plots

time <- c(seq(1985,2009,1))

pdf("graphics/gonorrhoea_data_raw.pdf",pointsize=14,width=8.0,height=6.0)
plot(range(time),c(0,10), ylab="Log Female Gonorrhea per 100,000", xlab="Time", main="", type="n")
for (j in 1:dim(Y0go)[2]) lines(time,Y0go[,j],col="darkgrey",lwd=0.5,lty=1)
lines(time,Y1go,col="black",lwd=5)
abline(v=time[T0go]+0.5,col="darkgrey",lty=2,lwd=1.5)
legend("topleft",legend=c("Other U.S. States","Rhode Island"), seg.len=2, col=c("darkgrey","black"),fill=NA,border=NA, lty=c(1,1),lwd=c(0.5,5), merge=T,bty="n")
graphics.off()

### Placebo specification tests

sens.classo.mb <- sens.sc.mb <- sens.did.mb <- sens.classo.all <- sens.sc.all <- sens.did.all <- matrix(NA,3,1)

Y1pre <- Y1go[1:T0go]
Y0pre <- Y0go[1:T0go,]

for (t in 1:3){
  sens.classo.mb[t,1]   <- moving.block.q(Y1pre,Y0pre,(T0go-t),t,"classo",q_norm)
  sens.classo.all[t,1]  <- all.q(Y1pre,Y0pre,(T0go-t),t,"classo",10000,q_norm)
  sens.sc.mb[t,1]       <- moving.block.q(Y1pre,Y0pre,(T0go-t),t,"sc",q_norm)
  sens.sc.all[t,1]      <- all.q(Y1pre,Y0pre,(T0go-t),t,"sc",10000,q_norm)
  sens.did.mb[t,1]      <- moving.block.q(Y1pre,Y0pre,(T0go-t),t,"did",q_norm)
  sens.did.all[t,1]     <- all.q(Y1pre,Y0pre,(T0go-t),t,"did",10000,q_norm)
}

xtable(cbind(sens.did.mb[,1],sens.sc.mb[,1],sens.classo.mb[,1],sens.did.all[,1],sens.sc.all[,1],sens.classo.all[,1]))

### Residual plots pre-treatment period

r.pre.classo      <- classo(Y1pre,Y0pre,1)
r.pre.sc          <- sc(Y1pre,Y0pre)

u.hat.go.pre.did    <- did(Y1pre,Y0pre)
u.hat.go.pre.sc     <- r.pre.sc$u.hat
u.hat.go.pre.classo <- r.pre.classo$u.hat

pdf("graphics/gonorrhoea_resid_pre_did.pdf",pointsize=16,width=6.0,height=6.0)
plot(range(seq(1985,2003,1)),c(-1,1), ylab="Residuals Pre-Treatment Period", xlab="Time", main="Difference-in-Differences",type="n")
points(seq(1985,2003,1),u.hat.go.pre.did,col="black",pch=20, lwd=2)
abline(h=0,col="grey",lty=2,lwd=1)
graphics.off()

pdf("graphics/gonorrhoea_resid_pre_sc.pdf",pointsize=16,width=6.0,height=6.0)
plot(range(seq(1985,2003,1)),c(-1,1), ylab="Residuals Pre-Treatment Period", xlab="Time", main="Synthetic Control",type="n")
points(seq(1985,2003,1),u.hat.go.pre.sc,col="black",pch=20, lwd=2)
abline(h=0,col="grey",lty=2,lwd=1)
graphics.off()

pdf("graphics/gonorrhoea_resid_pre_classo.pdf",pointsize=16,width=6.0,height=6.0)
plot(range(seq(1985,2003,1)),c(-1,1), ylab="Residuals Pre-Treatment Period", xlab="Time", main="Constrained Lasso",type="n")
points(seq(1985,2003,1),u.hat.go.pre.classo,col="black",pch=20, lwd=2)
abline(h=0,col="grey",lty=2,lwd=1)
graphics.off()

### No effect

p.noeff.did.mb      <- moving.block.q(Y1go,Y0go,T0go,T1go,"did",q_norm)
p.noeff.sc.mb       <- moving.block.q(Y1go,Y0go,T0go,T1go,"sc",q_norm)
p.noeff.classo.mb   <- moving.block.q(Y1go,Y0go,T0go,T1go,"classo",q_norm)

p.noeff.did.all     <- all.q(Y1go,Y0go,T0go,T1go,"did",10000,q_norm)
p.noeff.sc.all      <- all.q(Y1go,Y0go,T0go,T1go,"sc",10000,q_norm)
p.noeff.classo.all  <- all.q(Y1go,Y0go,T0go,T1go,"classo",10000,q_norm)

xtable(cbind(p.noeff.did.mb,p.noeff.sc.mb,p.noeff.classo.mb,p.noeff.did.all,p.noeff.sc.all,p.noeff.classo.all))

### Pointwise CI

alpha <- 0.1
grid <- c(seq(-5,2,0.01))

vec.ci.classo <- vec.ci.sc <- vec.ci.did <- m.ci.classo <- m.ci.sc <- m.ci.did <- NULL

for (t in 1:T1go){
  indices       <- c(1:T0go,T0go+t)
  Y1ci          <- Y1go[indices]
  Y0ci          <- Y0go[indices,]
  ci.sc         <- ci(Y1ci,Y0ci,"sc",alpha,grid)
  vec.ci.sc     <- rbind(vec.ci.sc,cbind((t+2003),ci.sc))
  m.ci.sc       <- cbind(m.ci.sc,mean(ci.sc))
  ci.classo     <- ci(Y1ci,Y0ci,"classo",alpha,grid)
  vec.ci.classo <- rbind(vec.ci.classo,cbind((t+2003),ci.classo))
  m.ci.classo   <- cbind(m.ci.classo,mean(ci.classo))
  ci.did        <- ci(Y1ci,Y0ci,"did",alpha,grid)
  vec.ci.did    <- rbind(vec.ci.did,cbind((t+2003),ci.did))
  m.ci.did      <- cbind(m.ci.did,mean(ci.did))
}

time <- c(seq(2004,2009,1))

pdf("graphics/ci_gonorrhoea_classo.pdf", pointsize=16,width=6.0,height=6.0)
plot(vec.ci.classo[,1],vec.ci.classo[,2], ylab="Gap in Log Female Gonorrhea per 100,000", xlab="Years", main="", col="black", pch=".", xlim = c(2003,2010), ylim=c(-2,2))
title("Constrained Lasso")
abline(h=0,col="grey",lty=2,lwd=2)
points(time, m.ci.classo,pch=16,cex=0.8,col="black")
legend(2003,2,legend=c("90% Confidence Interval"), col="black", cex=1,lty=c(1), lwd=c(2.25), bty=("n"))
legend(2003+0.3,2,legend=c(""), cex=1, col=c("black"), pch=16, bty=("n"))
graphics.off()

pdf("graphics/ci_gonorrhoea_sc.pdf", pointsize=16,width=6.0,height=6.0)
plot(vec.ci.sc[,1],vec.ci.sc[,2], ylab="Gap in Log Female Gonorrhea per 100,000", xlab="Years", main="", col="black", pch=".", xlim = c(2003,2010), ylim=c(-2,2))
title("Synthetic Control")
abline(h=0,col="grey",lty=2,lwd=2)
points(time, m.ci.sc,pch=16,cex=0.8,col="black")
legend(2003,2,legend=c("90% Confidence Interval"), col="black", cex=1,lty=c(1), lwd=c(2.25), bty=("n"))
legend(2003+0.3,2,legend=c(""), cex=1, col=c("black"), pch=16, bty=("n"))
graphics.off()

pdf("graphics/ci_gonorrhoea_did.pdf", pointsize=16,width=6.0,height=6.0)
plot(vec.ci.did[,1],vec.ci.did[,2], ylab="Gap in Log Female Gonorrhea per 100,000", xlab="Years", main="", col="black", pch=".", xlim = c(2003,2010), ylim=c(-2,2))
title("Difference-in-Differences")
abline(h=0,col="grey",lty=2,lwd=2)
points(time, m.ci.did,pch=16,cex=0.8,col="black")
legend(2003,2,legend=c("90% Confidence Interval"), col="black", cex=1,lty=c(1), lwd=c(2.25), bty=("n"))
legend(2003+0.3,2,legend=c(""), cex=1, col=c("black"), pch=16, bty=("n"))
graphics.off()

### Robustness checks

# Generate indicators for units which get non-zero weight in pre-treatment period

ind.sc      <- abs(r.pre.sc$w.hat)>1e-5
ind.classo  <- abs(r.pre.classo$w.hat[-1])>1e-5
ind.joint   <- ind.sc + ind.classo > 0
id.nonzero  <- which(ind.joint==1)

# Obtain p-values excluding one of the important control units at the time

p.vec.mb <- p.vec.all <- matrix(NA,length(id.nonzero),3)
for (i in 1:length(id.nonzero)) {
  ind             <- id.nonzero[i]
  p.vec.mb[i,1]   <- moving.block.q(Y1go,Y0go[,-ind],T0go,T1go,"did",q_norm)
  p.vec.mb[i,2]   <- moving.block.q(Y1go,Y0go[,-ind],T0go,T1go,"sc",q_norm)
  p.vec.mb[i,3]   <- moving.block.q(Y1go,Y0go[,-ind],T0go,T1go,"classo",q_norm)
  p.vec.all[i,1]  <- all.q(Y1go,Y0go[,-ind],T0go,T1go,"did",10000,q_norm)
  p.vec.all[i,2]  <- all.q(Y1go,Y0go[,-ind],T0go,T1go,"sc",10000,q_norm)
  p.vec.all[i,3]  <- all.q(Y1go,Y0go[,-ind],T0go,T1go,"classo",10000,q_norm)
}

# Plots

x <- c(rep(1,length(id.nonzero)),rep(2,length(id.nonzero)),rep(3,length(id.nonzero)))

pdf("graphics/robustness_mb.pdf", pointsize=14,width=6.0,height=6.0)
sizeplot(x,  c(p.vec.mb) , main="P-values, Moving Block Permutations", xlab="", ylab="p-value",xaxt='n',ylim=c(0,0.2))
axis(1, at=1:3, labels=c("DID","SC","CL"))
graphics.off()

pdf("graphics/robustness_iid.pdf", pointsize=14,width=6.0,height=6.0)
sizeplot(x,  c(p.vec.all) , main="P-values, iid Permutations", xlab="", ylab="p-value",xaxt='n',ylim=c(0,0.2))
axis(1, at=1:3, labels=c("DID","SC","CL"))
graphics.off()


