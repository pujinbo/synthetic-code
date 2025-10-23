# Empirical application for "Low-Rank Approximations  of Nonseparable Panel Models," by Ivan Fernandez-Val, Hugo Freeman and Martin Weidner

# Data source: Xu, Y. (2017) "Generalized synthetic control method: Causal inference with interactive fixed effects model," Political Analysis 25(1), 57--76.
# url: https://yiqingxu.org/papers/english/2016_Xu_gsynth/Xu_PA_replication.zip

# Main variables:
#     1 - abb: state id
#     2 - year
#     3 - turnout: election turnout
#     4 - policy_edr: indicator for election day registration policy
#     5 - policy_mail_in: indicator for universal mail in registration (not used in the analysis)
#     6 - policy_motor: indicator for motor voter registration (not used in the analysis)

rm(list=ls(all=TRUE)) ## eliminating everything in  memory

# Set working directory to replication folder
setwd("/Replication Package")

# General two-way matching
gtw<-function(Y,X,lambda = 0, R = 0){
  I = dim(Y)[1]; J = dim(Y)[2]
  Y.null = array(NA, c(I,J));
  joint.dist<-vector(mode = "list", length = 1);joint.feasible<-joint.dist; joint.set = joint.dist
  mu = array(NA, dim = c(I,J)); joint.nn <- array(NA, dim = c(I,J,2))
  Y.null = Y; Y.null[X!=0] = NA
  if(R == 0){Y.null. <- softImpute(Y.null, lambda = lambda, rank.max = (min(I,J)-1))
  }else{Y.null. = softImpute(Y.null, lambda = 0, rank.max = R)}
  Y.null = Y.null.$u%*%diag(Y.null.$d, nrow = length(Y.null.$d), ncol = length(Y.null.$d))%*%t(Y.null.$v)
  phi.alpha = Y.null.$u%*%sqrt(diag(Y.null.$d, nrow = length(Y.null.$d), ncol = length(Y.null.$d))); psi.gamma = Y.null.$v%*%sqrt(diag(Y.null.$d, nrow = length(Y.null.$d), ncol = length(Y.null.$d)))
  joint.set = which(X==0, arr.ind = TRUE)
  for (t in 1:J){
    for (i in 1:I){
      if (X[i,t] == 0){mu[i,t]<-Y[i,t]
      }else{
        i.nn.set = setdiff(which(X[,t]==0),i); t.nn.set = setdiff(which(X[i,]==0),t); it.nn.set = as.matrix(expand.grid(i.nn.set, t.nn.set))
        feasible.set = as.matrix(joint.set[which(apply(it.nn.set, MARGIN = 1, function(x){t(apply(1*(x == t(joint.set)), MARGIN = 2,FUN = sum))}) == 2,arr.ind = TRUE)[,1],])
        if(dim(feasible.set)[1]==1){# deal with vector valued feasible set
          euclid.i = sum((phi.alpha[i,] - phi.alpha[feasible.set[1,1],])^2)
          euclid.t = sum((psi.gamma[t,] - psi.gamma[feasible.set[1,2],])^2)
          euclid = list(euclid.i + euclid.t)
          joint.feasible[[i+(t-1)*I]] <- t(feasible.set)
        }else{
          dist.i = data.frame((phi.alpha[feasible.set[,1],] - t(array(phi.alpha[i,], dim = c(length(phi.alpha[i,]),length(feasible.set[,1])))))^2); 
          euclid.i = apply(dist.i, MARGIN = 1, sum)
          dist.t = data.frame((psi.gamma[feasible.set[,2],] - t(array(psi.gamma[t,], dim = c(length(psi.gamma[t,]),length(feasible.set[,1])))))^2); 
          euclid.t = apply(dist.t, MARGIN = 1, sum)
          # dist.i = t(sapply(phi.alpha[i,], function(x){x - phi.alpha[feasible.set[,1]]}));euclid.i = apply(dist.i,MARGIN = 2, function(x)(sum(x)^2))
          # dist.t = t(sapply(psi.gamma[t,], function(x){x - psi.gamma[feasible.set[,2]]}));euclid.t = apply(dist.t,MARGIN = 2, function(x)(sum(x)^2))
          euclid = as.vector(euclid.i + euclid.t)
          joint.feasible[[i+(t-1)*I]] <- as.matrix(feasible.set)
        }
        joint.dist[[i+(t-1)*I]] <- euclid
        joint.nn[i,t,] <- joint.feasible[[i+(t-1)*I]][which.min(joint.dist[[i+(t-1)*I]]),]
        mu[i,t] <- Y[i,joint.nn[i,t,2]] + Y[joint.nn[i,t,1],t] - Y[joint.nn[i,t,1],joint.nn[i,t,2]]
      }}
  }
  return(list(mu = mu, Y.null = Y.null,svd = Y.null., joint.nn = joint.nn, joint.dist = joint.dist, joint.feasible = joint.feasible))
}

nn <-function(Y, R = 3){
  s <- svd(Y); N <- dim(Y)[1]; T <- dim(Y)[2]
  alpha <- s$u[,1:R]%*%diag(s$d, nrow = R,ncol = R); gamma <- s$v[,1:R]
  alpha.nn <- array(NA, dim = N); gamma.nn <- array(NA, dim = T)
  for (i in 1:N){
    dist = t(sapply(alpha[i], function(x){x - alpha}))
    euclid = sapply(dist, function(x)(sum(x)^2))
    alpha.nn[i] = which.min(rank(euclid[-i])) + 1*(which.min(rank(euclid[-i]))+1>i)
  }
  for (i in 1:N){
    dist = t(sapply(gamma[i], function(x){x - gamma}))
    euclid = sapply(dist, function(x)(sum(x)^2))
    gamma.nn[i] = which.min(rank(euclid[-i])) + 1*(which.min(rank(euclid[-i]))+1>i)
  }
  joint.nn <- array(NA, dim = c(N,T,2))
  joint.nn[,,1] <- replicate(T,alpha.nn);joint.nn[,,2] <- t(replicate(N,gamma.nn))
  return(list(joint.nn = joint.nn))
}

# Matrix completion using soft imputation
mcsi<-function(Y,X, lambda = 0, tolerance = 1e-5, max.iterations = 5000){
  PYobs <- Y*(1-X); metric = 10e15; counter <- 1;check.improve = TRUE
  Z.old <- Z.new<-Y*(1-X); 
  while (metric > tolerance & counter<max.iterations & check.improve){
    PYnobs <- Z.old*X
    R <- max(which(svd(PYobs + PYnobs)$d > lambda))
    svd1 <- svd(PYobs + PYnobs, nu = R, nv = R); svd1$d<-sapply(svd1$d[1:R] - lambda, function(x){max(0, x)})
    Z.new <- svd1$u%*%diag(svd1$d, nrow = length(svd1$d), ncol = length(svd1$d))%*%t(svd1$v)
    check.improve<-(norm(Z.new - Z.old, type = "F") < metric); metric <- norm(Z.new - Z.old, type = "F"); 
    Z.old <- Z.new; counter <-counter+1 ; 
  }
  return(svd.complete = svd1)
}


# Fixed-Effects with matrix completion iteration scheme 
femc<-function(Y,X, lambda = 0, tolerance = 1e-5, max.iterations = 5000){
  PYobs <- Y*(1-X); metric = 10e15; counter <- 1;check.improve = TRUE
  Z.old <- Z.new<-Y*(1-X);  CY0 <- 0
  i<-row(Y); t = col(Y); i.flat <- c(i); t.flat <-c(t); X.flat <- c(X);
  while (metric > tolerance & counter<max.iterations){
    # PYnobs <- (Z.old + CY0)*X lambda <- lambda*sum(1-X)/2
    Y.tilde <- Y - CY0; 
    PYobs <- (Y.tilde)*(1-X); PYnobs <- (Z.old)*X
    R <- min(max(which(svd(PYobs + PYnobs)$d > lambda)), min(dim(Y))-1)
    svd1 <- svd(PYobs + PYnobs, nu = R, nv = R); svd1$d<-sapply(svd1$d[1:R] - lambda, function(x){max(0, x)})
    Z.new <- svd1$u%*%diag(svd1$d, nrow = length(svd1$d), ncol = length(svd1$d))%*%t(svd1$v)
    
    fit1 <- lm(c((Y -  Z.new)) ~  factor(i.flat) + factor(t.flat), subset = (X.flat == 0));
    CY0  <- matrix(predict(fit1, data.frame(i.flat,t.flat)), ncol = ncol(Y), nrow = nrow(Y), byrow = TRUE);
    
    
    check.improve<-(norm(Z.new - Z.old, type = "F") < metric); metric <- norm(Z.new - Z.old, type = "F"); 
    Z.old <- Z.new; counter <-counter+1 ; 
  }
  return(list(svd.FEout = svd1, CY = CY0, iterations = counter))
}


# Main program;

load("Data/election.RData") # read election data
state = unique(turnout$abb)
year = unique(turnout$year)
K = 3

# Pretrends analysis;

ever.treated <- function(data)
{
  return(max(data$policy_edr));
}
et <- by(turnout, turnout$abb, ever.treated);
et.states <- names(et[et == 1]);
mean.et <- by(turnout$turnout[turnout$abb %in% et.states], turnout$year[turnout$abb %in% et.states], mean);
mean.nt <- by(turnout$turnout[!(turnout$abb %in% et.states)], turnout$year[!(turnout$abb %in% et.states)], mean);
year    <- names(mean.et);

pdf("Results/ptrends.pdf", height=6, width=6);
plot(year[year<1976], mean.et[year<1976], ylim = c(44,76), type = "l", xlab = "year", ylab = "turnout (%)", lty = 1);
points(year[year<1976], mean.et[year<1976], ylim = c(45,75), pch = 16);
lines(year[year<1976], mean.nt[year<1976], col =2, lty = 1);
points(year[year<1976], mean.nt[year<1976], ylim = c(45,75), pch = 15, col =2);
legend("topleft", c("Ever Treated","Never Treated"), lty = c(1,1), lwd = c(1,1), col = c(1,2), pch = c(16,15), bty = "n");
dev.off();



# Create N*T arrays with the outcome (turnout) and treatment (EDR)
Y = array(NA, dim = c(length(state), length(year)))
for (s in 1:length(state)){
  Y[s,] = turnout$turnout[turnout$abb == state[s]]
}
X = array(NA, dim = c(length(state), length(year),K))
for (k in 1:K){
  for (s in 1:length(state)){
    X[s,,k] = turnout[turnout$abb == state[s], k +3]
  } 
}
X <- X[,,1]

I <- length(state); J <- length(year)

# Model with additive individual and time effects

# Matric completion method
# Setting lambda for soft-impute method so roughly cut-off past two factors. This is a hyperparameter we need to choose
mcY0.0 <- femc(Y,X,lambda = 0, max.iterations = 250) # This procedure takes a long time if initialised with lambda = 0
plot(log(mcY0.0$svd.FEout$d)) # looks like rank 6 or 8 is optimal from "elbow" point
lambda <- mcY0.0$svd.FEout$d[6]; 

# Estimate matrix completion counter-factual for treated 
mcY0 <- femc(Y,X,lambda = lambda); mcY <- mcY0$svd.FEout$u%*%diag(mcY0$svd.FEout$d, nrow = length(mcY0$svd.FEout$d),ncol = length(mcY0$svd.FEout$d))%*%t(mcY0$svd.FEout$v)
CY0 <- mcY0$CY; mcY0 <- mcY0$CY + mcY;
YC <- Y - CY0; 

# Crude search for optimal rank cut-off point for matching 
plot(log(svd(YC)$d)) # looks like rank 4/6 is optimal from "elbow" point

# Two-way matching estimator
# Estimate general two-way counter-factual for treated. The "est" step takes a long time
# install.packages("softImpute")
library(softImpute)
est  <- gtw(YC,X, R=0, lambda = lambda); 
twY0 <- YC; 
twY0[X==0]<-NA; one.way<- twY0
for (i in which(X==1)){
  index<-est$joint.nn[i - I*(ceiling(i/I)-1),ceiling(i/I),];
  twY0[i]<-YC[i - I*(ceiling(i/I)-1),index[2]] + YC[index[1],ceiling(i/I)] - YC[index[1],index[2]];
  one.way[i]<-YC[index[1],index[2]];
}
twY0 <- CY0 + twY0;
one.way <- CY0 + one.way

# Two-way matching estimator with multiple matches
# Estimate weighted many match counter-factual for treated
# Matching only on residuals
mwY0 <- YC; 
mwY0[X==0]<-NA; 
matches <- 10;
for (i in which(X==1)){
  index <-cbind(est$joint.feasible[[i]],est$joint.dist[[i]]);
  index <- (index[order(index[,3],decreasing = FALSE),])[1:min(matches,dim(index)[1]),];
  mwY0[i]<-diag(YC[i - I*(ceiling(i/I)-1),index[,2]] + YC[index[,1],ceiling(i/I)] - YC[index[,1],index[,2]])%*%as.matrix(1/(index[,3]*sum(1/index[,3])));
}
mwY0 <- CY0 + mwY0;

# One-way matching estimator with multiple matches
# Estimate weighted many match counter-factual for treated
# Matching only on residuals
one.way.many <- mwY0
matches2 <- 5;
for (i in which(X==1)){
  index <-cbind(est$joint.feasible[[i]],est$joint.dist[[i]]);
  index <- (index[order(index[,3],decreasing = FALSE),])[1:min(matches2,dim(index)[1]),];
  one.way.many[i]<-diag(YC[index[,1],index[,2]])%*%as.matrix(1/(index[,3]*sum(1/index[,3])));
}
one.way.many <- CY0 + one.way.many

# Estimate simple averaging counter-factual for treated
saY1 <- Y; saY1[X==0] <- NA

# Differences in means
Dmeans <- apply(X*Y,2,sum) / apply(X,2,sum) - apply((1-X)*Y,2,sum) / apply(1-X,2,sum); # Comparison of means for treated and nontreated by year
Dmeans   <- Dmeans[!is.na(Dmeans)];

# Difference-in-difference estimators;
fit_did1 <- lm(turnout ~ policy_edr + factor(year) + factor(abb), data = turnout); # Diff-in-diff estimates with time varying ATT
did1      <- array(coef(fit_did1)[2], dim = length(Dmeans));
fit_did2 <- lm(turnout ~ policy_edr:factor(year) + factor(abb), data = turnout); # Diff-in-diff estimates with time varying ATT
did2      <- coef(fit_did2)[62:71];


# Take time specific averages, ignoring NA's, which are the non-treated
gtwm  <- apply(twY0,MARGIN = 2,function(x){if(length(na.omit(x))>0){sum(na.omit(x))/length(na.omit(x))}else{NA}});gtwm<-gtwm[!is.na(gtwm)]
gmwm  <- apply(mwY0,MARGIN = 2,function(x){if(length(na.omit(x))>0){sum(na.omit(x))/length(na.omit(x))}else{NA}});gmwm<-gmwm[!is.na(gmwm)]
g1w.single  <- apply(one.way,MARGIN = 2,function(x){if(length(na.omit(x))>0){sum(na.omit(x))/length(na.omit(x))}else{NA}});g1w.single<-g1w.single[!is.na(g1w.single)]
g1w.many  <- apply(one.way.many,MARGIN = 2,function(x){if(length(na.omit(x))>0){sum(na.omit(x))/length(na.omit(x))}else{NA}});g1w.many<-g1w.many[!is.na(g1w.many)]
mcY0[X==0]<-NA
mcse  <- apply(mcY0,MARGIN = 2,function(x){if(length(na.omit(x))>0){sum(na.omit(x))/length(na.omit(x))}else{NA}});mcse<-mcse[!is.na(mcse)]
y1bar <- apply(saY1,MARGIN = 2,function(x){if(length(na.omit(x))>0){sum(na.omit(x))/length(na.omit(x))}else{NA}});y1bar<-y1bar[!is.na(y1bar)]

# Compute the ATT
ATT     <- cbind(Dmeans, did1,did2, array(y1bar, dim = c(length(y1bar),5)) - cbind(gtwm,gmwm,g1w.single,g1w.many,mcse)); row.names(ATT)<-year[15:24]

# Compute time-averaged ATT
treated.sum <- apply(X,2,sum);
treated.sum <- treated.sum[treated.sum > 0];
ATT.avg  <- apply(ATT, 2, weighted.mean, treated.sum);

# Figure of results
ATT.plot <- cbind(Dmeans, did2, array(y1bar, dim = c(length(y1bar),3)) - cbind(gmwm,g1w.many,mcse));
pdf("Results/ATT_with_additive_effects.pdf", height=6, width=6);
matplot(ATT.plot, type = 'l', lty = c(1,1,1,1,1,1,1),xlab='year', xaxt = 'n', ylab="turnout (%)"); axis(1, at=1:10, labels=rownames(ATT),las=0)
matpoints(ATT.plot, type = 'p', pch = c(16,15,18,17,19,20)); 
abline(h=0, col = "grey", lty=2)
legend("bottomleft", c("Dmeans","DiD",paste("TWM-",matches, sep = ''), paste("SM-",matches2, sep = ''),"MC"), pch = c(16,15,18,17,19), lty = c(1,1,1,1,1), lwd = c(1,1,1,1,1), col = c(1,2,3,4,5), bty = "n");
dev.off()


# Distribution effects;

ys <- sort(turnout$turnout[turnout$policy_edr==1]);
ys <- quantile(turnout$turnout, c(10:98)/100);
DTT0.mc.fe   <- rep(0,length(ys));
DTT0.twm.fe  <- rep(0,length(ys));
DTT0.sm.fe   <- rep(0,length(ys));
newdata   <- turnout;
newdata$policy_edr <- 0;
approx.factors.fe  <- 3;
matches <- 10;
matches2 <- 5;


Iy = array(NA, dim = c(length(state), length(year)));

trim01 <- function(x) { ifelse(x < 0, 0, ifelse(x > 1, 1, x))};


for (y in 1:length(ys)){
  for (s in 1:length(state)){
    Iy[s,] = 1*(turnout$turnout[turnout$abb == state[s]] <= ys[y]);
  }

  # Matrix completion with additive effects
  mcY0.0 <- femc(Iy,X,lambda = 0, max.iterations = 250) # This procedure takes a long time if initialised with lambda = 0
  lambda <- mcY0.0$svd.FEout$d[approx.factors.fe];
  mcY0 <- femc(Iy,X,lambda = lambda); mcY <- mcY0$svd.FEout$u%*%diag(mcY0$svd.FEout$d, nrow = length(mcY0$svd.FEout$d),ncol = length(mcY0$svd.FEout$d))%*%t(mcY0$svd.FEout$v)
  CY0 <- mcY0$CY; mcY0 <- mcY0$CY + mcY; mcY0[X==0]<-NA
  DTT0.mc.fe[y] <- mean(mcY0, na.rm = TRUE);
  DTT0.mc.fe[y] <- ifelse(DTT0.mc.fe[y]<0,0,ifelse(DTT0.mc.fe[y]>1,1,DTT0.mc.fe[y]));
  IyC <- Iy - CY0;

  # Two-way matching estimator with multiple matches
  est  <- gtw(IyC,X, R=0, lambda = lambda); 
  mwY0 <- IyC; 
  mwY0[X==0]<-NA; 
  for (i in which(X==1)){
    index <-cbind(est$joint.feasible[[i]],est$joint.dist[[i]]);
    index <- (index[order(index[,3],decreasing = FALSE),])[1:min(matches,dim(index)[1]),];
    mwY0[i]<-diag(IyC[i - I*(ceiling(i/I)-1),index[,2]] + IyC[index[,1],ceiling(i/I)] - IyC[index[,1],index[,2]])%*%as.matrix(1/(index[,3]*sum(1/index[,3])));
  }
  DTT0.twm.fe[y] <- mean(apply(CY0 + mwY0, c(1,2), trim01), na.rm = TRUE);

  # One-way matching estimator with multiple matches
  one.way.many <- mwY0;
  one.way.many[X==0]<-NA; 
  for (i in which(X==1)){
    index <-cbind(est$joint.feasible[[i]],est$joint.dist[[i]]);
    index <- (index[order(index[,3],decreasing = FALSE),])[1:min(matches2,dim(index)[1]),];
    one.way.many[i]<-diag(IyC[index[,1],index[,2]])%*%as.matrix(1/(index[,3]*sum(1/index[,3])));
  }
  DTT0.sm.fe[y] <- mean(apply(CY0 + one.way.many, c(1,2), trim01), na.rm = TRUE);
  
}

taus <- c(10:90)/100;
us   <- c(0:100)/100;

QY1    <- approxfun(us, quantile(turnout$turnout[turnout$policy_edr==1],us));
QY0.dm <- approxfun(us, quantile(turnout$turnout[turnout$policy_edr==0],us));
QY0.mc.fe <- approxfun(sort(DTT0.mc.fe), ys);
QY0.twm.fe <- approxfun(sort(DTT0.twm.fe), ys);
QY0.sm.fe <- approxfun(sort(DTT0.sm.fe), ys);


pdf("Results/QTT_with_additive_effects.pdf", height=6, width=6);
plot(taus, QY1(taus) - QY0.dm(taus), ylim = c(-0.5,17.5), type = "l", xlab = "quantile index", ylab = "turnout (%)", lty = 3);
lines(taus, QY1(taus) - QY0.mc.fe(taus), col = 3, lty =1);
lines(taus, QY1(taus) - QY0.twm.fe(taus), col = 4, lty =2);
lines(taus, QY1(taus) - QY0.sm.fe(taus), col = 2, lty = 4);
abline(h=0, col = "grey", lty =2);
legend("topright", c("Dquantiles",paste("TWM-",matches, sep = ''), paste("SM-",matches2, sep = ''),"MC"), lty = c(3,1,2,4), lwd = c(1,1,1,1), col = c(1,3,4,2), bty = "n");
dev.off();


save.image(file = "Results/election-results.RData");

