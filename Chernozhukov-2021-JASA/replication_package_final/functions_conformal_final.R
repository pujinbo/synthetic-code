########################################################################################################
# Generic functions for conformal inference
# Paper: An exact and robust conformal inference approach for counterfactual and synthetic controls
# Authors: V. Chernozhukov, K. Wuthrich, Y. Zhu
# DISCLAIMER: This software is provided "as is" without warranty of any kind, expressed or implied. 
# Questions/error reports: kwuthrich@ucsd.edu
########################################################################################################

### Required packages

library(limSolve)
library(CVXR)

### Functions

# Difference-in-differences
did <- function(Y1,Y0){
  u.hat <- Y1 - mean(Y1-rowMeans(Y0)) - rowMeans(Y0)
  return(u.hat)
}

# Synthetic control (requires R-package limSolve)
sc <- function(Y1,Y0){
  J <- dim(Y0)[2]
  e <- matrix(1,1,J)
  f <- 1
  g <- diag(x=1,J,J)
  h <- matrix(0,J,1)  
  w.hat <- lsei(A=Y0,B=Y1,E=e,F=f,G=g,H=h,type=2)$X
  u.hat <- Y1-Y0%*%w.hat
  return(list(u.hat=u.hat,w.hat=w.hat))
}

# Constrained Lasso (requires R-package CVXR)
classo <- function(Y1,Y0,K){
  J       <- dim(Y0)[2]
  w       <- Variable((J+1))
  loss    <- mean((Y1-cbind(1,Y0)%*%w)^2)
  constr  <- list(sum(abs(w[2:(J+1)])) <= K)
  prob    <- Problem(Minimize(loss),constr)
  result  <- solve(prob)
  w.hat   <- result$getValue(w)
  u.hat   <- Y1 - cbind(1,Y0)%*%w.hat
  return(list(u.hat=u.hat,w.hat=w.hat))
}

# Moving block permutations
moving.block.q <- function(Y1,Y0,T0,T1,M,q){
  T01 <- T0+T1
  if (M=="classo"){
    u.hat <- classo(Y1,Y0,1)$u.hat
  } 
  if (M=="sc"){
    u.hat <- sc(Y1,Y0)$u.hat
  }
  if (M=="did"){
    u.hat <- did(Y1,Y0)
  }
  sub.size  <- T1
  u.hat.c   <- c(u.hat,u.hat)
  S.vec     <- matrix(NA,T01,1)
  if (q==1){
    for (s in 1:(T01)){
      S.vec[s,1]  <- sum(abs(u.hat.c[s:(s+sub.size-1)])) 
    }
  }
  if (q==2){
    for (s in 1:(T01)){
      S.vec[s,1]  <- sqrt(sum((u.hat.c[s:(s+sub.size-1)])^2))
    }
  }
  p <- mean(S.vec>=S.vec[T0+1])
  return(p)
}

# All/iid permutations (use random subsample of all permutations)
all.q <- function(Y1,Y0,T0,T1,M,nperm,q){
  T01 <- T0+T1
  if (M=="classo"){
    u.hat <- classo(Y1,Y0,1)$u.hat
  } 
  if (M=="sc"){
    u.hat <- sc(Y1,Y0)$u.hat
  }
  if (M=="did"){
    u.hat <- did(Y1,Y0)
  }
  post.ind  <- ((T0+1):T01)
  pre.ind   <- (1:T0)
  S.vec     <- matrix(NA,nperm,1)
  if (q==1){
    Sq <- sum(abs(u.hat[post.ind]))
    for (r in 1:nperm){
      u.hat.p     <- u.hat[sample(1:T01,replace=F)]
      S.vec[r,1]  <- sum(abs(u.hat.p[post.ind])) 
    }
  }
  if (q==2){
    Sq <- sqrt(sum((u.hat[post.ind])^2))
    for (r in 1:nperm){
      u.hat.p     <- u.hat[sample(1:T01,replace=F)]
      S.vec[r,1]  <- sqrt(sum((u.hat.p[post.ind])^2)) 
    }
  }
  p <- 1/(nperm+1)*(1+sum(S.vec>=Sq))
  return(p)
}

# Confidence interval via test inversion 
ci <- function(Y1,Y0,M,alpha,thetagrid){
  T01 <- dim(Y0)[1]
  T0 <- T01-1
  p.vec <- rep(NA,length(thetagrid))
  for (ind in 1:length(thetagrid)){
    Y1_0 <- Y1
    Y1_0[T01] <- Y1[T01] - thetagrid[ind]
    if (M=="classo"){
      u.hat <- classo(Y1_0,Y0,1)$u.hat
    } 
    if (M=="sc"){
      u.hat <- sc(Y1_0,Y0)$u.hat
    }
    if (M=="did"){
      u.hat <- did(Y1_0,Y0)
    }
    p.vec[ind] <- mean(abs(u.hat)>=abs(u.hat[T01]))
  }
  ci <- thetagrid[p.vec>alpha]
  return(ci)
}