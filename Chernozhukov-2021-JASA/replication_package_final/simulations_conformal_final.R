########################################################################################################
# Simulation study
# Paper: An exact and robust conformal inference approach for counterfactual and synthetic controls
# Authors: V. Chernozhukov, K. Wuthrich, Y. Zhu
# DISCLAIMER: This software is provided "as is" without warranty of any kind, expressed or implied. 
# Questions/error reports: kwuthrich@ucsd.edu
########################################################################################################

########################################################################################################
# Preliminaries
########################################################################################################

rm(list = ls())

setwd("/Users/kasparwuthrich/Dropbox/research/SC/SC with Victor and Yinchu/Conformal Paper/Code/replication_package_final")

set.seed(12345)

### Packages

library(xtable)
library(VGAM)

source("functions_conformal_final.R")

########################################################################################################
# Additional functions
########################################################################################################

# Generate AR series
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

# Simulation function
sim <- function(stationary,T0,T1,J,rho1,rho2,alpha.sig,DGP,alternative){
  T01 <- T0 + T1
  
  F1 <- rnorm(T01)
  if (stationary==1){
    F2 <- rnorm(T01)
  }
  if (stationary==0){
    F2 <- rnorm(T01)+c(1:T01)
  }
  
  lambda1 <- c(seq(1,J,1)/J)
  lambda2 <- c(seq(1,J,1)/J)
  
  w1  <- cbind(rep(1/J,J))
  w2  <- cbind(c(1/3,1/3,1/3,rep(0,J-3)))
  w3  <- (-1)*cbind(rep(1/J,J))
  w4  <- cbind(c(1,-1,rep(0,J-2)))
  w   <- w1*(DGP==1)+w2*(DGP==2)+w3*(DGP==3)+w4*(DGP==4)
  
  y <- kronecker(lambda1,matrix(1,T01,1)) + kronecker(matrix(1,J,1),F1) + kronecker(cbind(lambda2),cbind(F2)) 
  
  eps <- matrix(NA,T01,J)
  for (j in 1:J){
    eps[,j] <- generate.AR.series(T01,rho1)
  }
  
  Y0  <- matrix(y,T01,J)+eps
  u   <- generate.AR.series(T01,rho2)
  Y1  <- Y0 %*% w + u
  Y1[(T0+1):T01]  <- Y1[(T0+1):T01] + alternative

  p.did     <- moving.block.q(Y1,Y0,T0,T1,"did",1)
  p.sc      <- moving.block.q(Y1,Y0,T0,T1,"sc",1)
  p.classo  <- moving.block.q(Y1,Y0,T0,T1,"classo",1)
  
  rej.did     <- (p.did<=alpha.sig)
  rej.sc      <- (p.sc<=alpha.sig)
  rej.classo  <- (p.classo<=alpha.sig)

  return(cbind(rej.did,rej.sc,rej.classo))
}

# Wrapper
sim.wrapper <- function(stationary,T0vec,Jvec,rho1,rho2,alpha.sig,DGP,alternative,nreps){
  
  results.sc      <- matrix(NA,length(T0vec),length(Jvec))
  results.did     <- matrix(NA,length(T0vec),length(Jvec))
  results.classo  <- matrix(NA,length(T0vec),length(Jvec))
  
  for(j in 1:length(Jvec)){
    J <- Jvec[j]
    for (t in 1:length(T0vec)){
      
      T0  <- T0vec[t]
      T1  <- 1
      T01 <- T0+T1
      
      #print(c(T0,J))
      simulation.result <- matrix(NA,nreps,3)
      for (r in 1:nreps){
        simulation.result[r,] <- sim(stationary,T0,T1,J,rho1,rho2,alpha.sig,DGP,alternative)
      }
      
      results <- colMeans(simulation.result)
      
      results.did[t,j]    <- results[1]
      results.sc[t,j]     <- results[2]
      results.classo[t,j] <- results[3]
      
    }
  }
  return(cbind(results.did,results.sc,results.classo))
}

# Power curve
power.curve <- function(T0,J,rho1,rho2,alpha.sig,DGP,alternative.vec,nreps){
  stationary <- 1
  power.curve.results <- matrix(NA,length(alternative.vec),3)
  for (i in 1:length(alternative.vec)){
    alternative <- alternative.vec[i]
    power.curve.results[i,] <- sim.wrapper(stationary,T0,J,rho1,rho2,alpha.sig,DGP,alternative,nreps) 
  }
  return(power.curve.results)
}

# Oracle power
oracle.power <- function(alt){
  power <- pfoldnorm(qfoldnorm(1-alpha.sig), mean = alt, lower.tail = FALSE)
  return(power)
}


########################################################################################################
# Simulations
########################################################################################################

### Setting

set.seed(12345)

nreps     <- 5000
alpha.sig <- 0.1

### Tables

T0vec <- c(20,50,100)
Jvec  <- c(20,50,100)

# Size stationary design

size.iid <- NULL 
for (dgp in 1:4){
  res       <- sim.wrapper(1,T0vec,Jvec,0,0,alpha.sig,dgp,0,nreps)
  size.iid  <- rbind(size.iid,res)
  print(dgp)
  print(xtable(cbind(T0vec,res),digits=c(0,0,rep(2,3*length(Jvec)))),include.rownames=F)
}

size.wd <- NULL 
for (dgp in 1:4){
  res       <- sim.wrapper(1,T0vec,Jvec,0.6,0.6,alpha.sig,dgp,0,nreps)
  size.wd   <- rbind(size.wd,res)
  print(dgp)
  print(xtable(cbind(T0vec,res),digits=c(0,0,rep(2,3*length(Jvec)))),include.rownames=F)
}

# Size nonstationary design

size.iid.nonstat <- NULL 
for (dgp in 1:4){
  res               <- sim.wrapper(0,T0vec,Jvec,0,0,alpha.sig,dgp,0,nreps)
  size.iid.nonstat  <- rbind(size.iid.nonstat,res)
  print(dgp)
  print(xtable(cbind(T0vec,res),digits=c(0,0,rep(2,3*length(Jvec)))),include.rownames=F)
}

size.wd.nonstat <- NULL 
for (dgp in 1:4){
  res             <- sim.wrapper(0,T0vec,Jvec,0.6,0.6,alpha.sig,dgp,0,nreps)
  size.wd.nonstat <- rbind(size.wd.nonstat,res)
  print(dgp)
  print(xtable(cbind(T0vec,res),digits=c(0,0,rep(2,3*length(Jvec)))),include.rownames=F)
}

### Power curve

T0              <- 19
J               <- 50
alternative.vec <- seq(0,6,1)

oracle.power.vec <- unlist(lapply(alternative.vec,oracle.power))

results.pc.wd <- NULL
for (dgp in 1:4){
  pc.wd         <- power.curve(T0,J,0.6,0.6,alpha.sig,dgp,alternative.vec,nreps)
  results.pc.wd <- cbind(results.pc.wd,pc.wd)
  
  pdf(paste("graphics/pc",dgp,".pdf",sep=""), pointsize=13,width=6.0,height=6.0)
  plot(range(alternative.vec),c(0,1), type="n", ylab="Empirical Rejection Rates", xlab="Alternative", main="")
  title(main= paste("DGP",dgp,sep=""))
  par(lty = 1, col = 1, lwd = 3)
  lines(alternative.vec, pc.wd[,1],type="b",pch=1)
  par(lty = 1, col = 1, lwd = 3)
  lines(alternative.vec, pc.wd[,2],type="b",pch=2)
  par(lty = 1, col = 1, lwd = 3)
  lines(alternative.vec, pc.wd[,3],type="b",pch=4)
  par(lty = 1, col = 1, lwd = 3)
  lines(alternative.vec, oracle.power.vec,type="b",pch=8)
  abline(h=0.1,col="darkgrey",lty=3,lwd=0.5)
  par(lty = 1, lwd = 1.5, col = 1)
  legend("topleft", c("Difference-in-Differences","Synthetic Control", "Constrained Lasso","Oracle Power"), lty = c(1,1,1,1), col=c(1,1,1,1), pch=c(1,2,4,8), lwd=c(3,3,3,3), bty = "n");
  dev.off()
}

filename <- paste("RData/simulations",format(Sys.time(),"%Y-%m-%d_%T"),".RData",sep="")
save.image(file=filename)

