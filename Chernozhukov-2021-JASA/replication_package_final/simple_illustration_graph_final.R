###########################################################################################################
# Stylized simulation to illustrate the importance of imposing the null
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

### Functions
source("functions_conformal_final.R")

###########################################################################################################
# Additional functions
###########################################################################################################

generate.AR.series <- function(T01,rho){
  u           <- rep(NA,T01)
  epsl        <- rnorm(T01)*sqrt(1-rho^2) 
  startvalue  <- rnorm(1)
  u[1]        <- rho*startvalue + epsl[1]
  for (t in 2:T01){
    u[t] <- rho*u[(t-1)]+epsl[t]
  }
  return(u)
}

sc.nonull <- function(Y1,Y0){
  T0    <- length(Y1)-1
  Y1est <- Y1[1:T0]
  Y0est <- Y0[1:T0,]
  w.hat <- sc(Y1est,Y0est)$w.hat
  u.hat <- Y1 - Y0 %*% w.hat
  return(u.hat)
}

sc.null <- function(Y1,Y0){
  w.hat <- sc(Y1,Y0)$w.hat
  u.hat <- Y1 - Y0 %*% w.hat
  return(u.hat)
}

sim <- function(T0,J,alpha.sig,rho){
  T01     <- T0 + 1
  Y0      <- matrix(rnorm((J)*T01),T01,J)
  w       <- matrix(0,J,1)
  w[1:3]  <- 1/3

  Y1 <- Y0%*%w + generate.AR.series(T01,rho)
  
  u.nonull  <- sc.nonull(Y1,Y0)
  p.nonull  <- mean(abs(u.nonull)>=abs(u.nonull[T0+1]))
  u.null    <- sc.null(Y1,Y0)
  p.null    <- mean(abs(u.null)>=abs(u.null[T0+1]))
  
  rej.null    <- (p.null<=alpha.sig)
  rej.nonull  <- (p.nonull<=alpha.sig)

  return(cbind(rej.null,rej.nonull))
}

###########################################################################################################
# Simulations
###########################################################################################################

nreps     <- 100000
alpha.sig <- 0.1
T0        <- 19
J         <- 50
rho.vec   <- c(seq(0.0,0.9,0.1),0.95)

results.overall <- matrix(NA,length(rho.vec),2)
for (r in 1:length(rho.vec)){
  rho     <- rho.vec[r]
  results <- matrix(NA,nreps,2)
  for (rep in 1:nreps){
    results[rep,] <- sim(T0,J,alpha.sig,rho)
  }
  results.overall[r,] <- colMeans(results)
}

T0 <- 99
results.overall.ls <- matrix(NA,length(rho.vec),2)
for (r in 1:length(rho.vec)){
  rho         <- rho.vec[r]
  results.ls  <- matrix(NA,nreps,2)
  for (rep in 1:nreps){
    results.ls[rep,] <- sim(T0,J,alpha.sig,rho)
  }
  results.overall.ls[r,] <- colMeans(results.ls)
}


###########################################################################################################
# Figures
###########################################################################################################

graphics.off()
pdf("graphics/imposing_null_intro.pdf",pointsize=14,width=8.0,height=6.0)
plot(1, ylab="Empirical Rejection Probability", xlab="AR(1) Coefficient", main="", col="black", pch=".", xlim = range(0,1), ylim=c(0,0.5))
lines(rho.vec,results.overall[,1],col="black",lwd=3,lty=1,type="b",pch=1)
lines(rho.vec,results.overall[,2],col="black",lwd=3,lty=1,type="b",pch=2)
abline(h=0.1,col="black",lty=2,lwd=1)
legend("topleft",legend=c("Imposing the null hypothesis","Not imposing the null hypothesis"), seg.len=2, col=c("black","black"),fill=NA,border=NA, lty=c(1,1),lwd=c(3,3),pch=c(1,2), merge=T,bty="n")
graphics.off()

pdf("graphics/imposing_null_comparison.pdf",pointsize=14,width=8.0,height=6.0)
plot(1, ylab="Empirical Rejection Probability", xlab="AR(1) Coefficient", main="", col="black", pch=".", xlim = range(0,1), ylim=c(0,0.5))
lines(rho.vec,results.overall[,1],col="black",lwd=3,lty=1,type="b",pch=1)
lines(rho.vec,results.overall[,2],col="black",lwd=3,lty=1,type="b",pch=2)
lines(rho.vec,results.overall.ls[,1],col="black",lwd=3,lty=1,type="b",pch=3)
lines(rho.vec,results.overall.ls[,2],col="black",lwd=3,lty=1,type="b",pch=4)
abline(h=0.1,col="black",lty=2,lwd=1)
legend("topleft",legend=c("Imposing the null hypothesis (T=20)","Not imposing the null hypothesis (T=20)","Imposing the null hypothesis (T=100)","Not imposing the null hypothesis (T=100)"), seg.len=2, col=c("black","black","black","black"),fill=NA,border=NA, lty=c(1,1,1,1),lwd=c(3,3,3,3),pch=c(1,2,3,4), merge=T,bty="n")
graphics.off()

