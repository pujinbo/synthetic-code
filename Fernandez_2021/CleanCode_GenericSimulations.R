# Simulations for "Low-Rank Approximations  of Nonseparable Panel Models," by Ivan Fernandez-Val, Hugo Freeman and Martin Weidner

## Clear environment
rm(list=ls(all=TRUE)) ## eliminating everything in  memory
## softImpute package installation
# install.packages("softImpute")

# Set working directory to replication folder
setwd("/Replication Package")

# General Two way function will return, for each i and t, the list of feasible indices to match with along with associated distance metric
gtw<-function(Y,X,lambda = 0, R = 0){
  #Define dimensions 
  I = dim(Y)[1]; J = dim(Y)[2] 
  #Create empty matrix completion counterfactual array; NA for softImpute function 
  Y.null = Y; Y.null[X!=0] = NA
  #Create empty singular vector distance, feasible set and joint feasible set arrays
  joint.dist<-vector(mode = "list", length = 1)
  joint.feasible<-joint.dist
  joint.set = joint.dist 
  # Create empty completed matrix. This is not used when I call the functions in the simulations below but included for completeness
  mu = array(NA, dim = c(I,J)); 
  # Create empty nearest neighbour array for i and t, hence third dimension is "2"
  joint.nn <- array(NA, dim = c(I,J,2))
  # Call softImpute function depending on input for R
  # R = 0 uses soft imputation for the chosen lambda
  if(R == 0){Y.null. <- softImpute(Y.null, lambda = lambda, rank.max = (min(I,J)-1))
  # R > 0 uses hard impute with output rank = R 
  }else{Y.null. = softImpute(Y.null, lambda = 0, rank.max = R)}
  # Use softImpute output to construct matrix completion solution
  Y.null = Y.null.$u%*%diag(Y.null.$d, nrow = length(Y.null.$d), ncol = length(Y.null.$d))%*%t(Y.null.$v)
  # Create proxy phi and psi arrays for nearest neightbour calculation
  phi.alpha = Y.null.$u%*%sqrt(diag(Y.null.$d, nrow = length(Y.null.$d), ncol = length(Y.null.$d))); 
  psi.gamma = Y.null.$v%*%sqrt(diag(Y.null.$d, nrow = length(Y.null.$d), ncol = length(Y.null.$d)))
  # Create index list where X = 0, i.e. indices to use for matching (this gets more complicated below)
  joint.set = which(X==0, arr.ind = TRUE)
  # A for loop must be used to perform our nearest neightbour matching
  for (t in 1:J){
    for (i in 1:I){
      # When X = 0 we don't apply any matching, just use actual values
      if (X[i,t] == 0){mu[i,t]<-Y[i,t]
      }else{
        # Finding feasible matching set: 
        # For i and t start with indices not equal i or t but with X = 0
        i.nn.set = setdiff(which(X[,t]==0),i); 
        t.nn.set = setdiff(which(X[i,]==0),t); 
        # Join these sets for all combinations 
        it.nn.set = as.matrix(expand.grid(i.nn.set, t.nn.set))
        # Include only where X_is = X_jt = X_js = 0
        feasible.set = as.matrix(joint.set[which(apply(it.nn.set, MARGIN = 1, function(x){t(apply(1*(x == t(joint.set)), MARGIN = 2,FUN = sum))}) == 2,arr.ind = TRUE)[,1],])
        # If feasible set is just one then need vector calculation (this could probably be folded into the dim != 1 part)
        if(dim(feasible.set)[1]==1){# deal with vector valued feasible set
          # Calculate distance
          euclid.i = sum((phi.alpha[i,] - phi.alpha[feasible.set[1,1],])^2)
          euclid.t = sum((psi.gamma[t,] - psi.gamma[feasible.set[1,2],])^2)
          euclid = list(euclid.i + euclid.t)
          # since only one feasible observation this becomes the joint feasible set 
          joint.feasible[[i+(t-1)*I]] <- t(feasible.set)
        }else{ # matrix valued feasible set
          # calculate distance
          dist.i = data.frame((phi.alpha[feasible.set[,1],] - t(array(phi.alpha[i,], dim = c(length(phi.alpha[i,]),length(feasible.set[,1])))))^2); 
          euclid.i = apply(dist.i, MARGIN = 1, sum)
          dist.t = data.frame((psi.gamma[feasible.set[,2],] - t(array(psi.gamma[t,], dim = c(length(psi.gamma[t,]),length(feasible.set[,1])))))^2); 
          euclid.t = apply(dist.t, MARGIN = 1, sum)
          euclid = as.vector(euclid.i + euclid.t)
          # create array of which indices in feasible set these distances are associated with
          joint.feasible[[i+(t-1)*I]] <- as.matrix(feasible.set)
        }
        # assign distance to list of feasible indices
        joint.dist[[i+(t-1)*I]] <- euclid
        # indice with smallest join distance. Agian, I do not use this in the simulations but include this as an automatic single match nearest neighbour estimator
        joint.nn[i,t,] <- joint.feasible[[i+(t-1)*I]][which.min(joint.dist[[i+(t-1)*I]]),]
        mu[i,t] <- Y[i,joint.nn[i,t,2]] + Y[joint.nn[i,t,1],t] - Y[joint.nn[i,t,1],joint.nn[i,t,2]]
      }}
  }
  return(list(mu = mu, Y.null = Y.null,svd = Y.null., joint.nn = joint.nn, joint.dist = joint.dist, joint.feasible = joint.feasible))
}
# Matrix completion method based on softImpute R function. This returns very close output to softImpute and is written as stepping stone to femc below. This can be skipped if you wish
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
# Fixed-Effects with matrix completion iteration scheme. Returns matrix completion with additive FE
femc<-function(Y,X, lambda = 0, tolerance = 1e-5, max.iterations = 5000){
  # Initialise observed matrix 
  PYobs <- Y*(1-X); 
  # Initialise convergence metric
  metric = 10e15; 
  # Initialise loop counter
  counter <- 1;
  # Flag to check if iteration has improved convergence metric. Not used here
  check.improve = TRUE
  # Initialise FE matrix and unobserved matrix entries
  Z.old <- Z.new<-Y*(1-X);  CY0 <- 0
  # Vectorise indices and X
  i<-row(Y); t = col(Y); i.flat <- c(i); t.flat <-c(t); X.flat <- c(X);
  while (metric > tolerance & counter<max.iterations){
    # Ignore: PYnobs <- (Z.old + CY0)*X lambda <- lambda*sum(1-X)/2
    
    # Project out FE
    Y.tilde <- Y - CY0; 
    # Define observed and unobserved outcome matrix
    PYobs <- (Y.tilde)*(1-X); PYnobs <- (Z.old)*X
    # Set rank according to value for lambda
    R <- min(max(which(svd(PYobs + PYnobs)$d > lambda)), min(dim(Y))-1)
    # Singular value decomposition after shrinking singular values by lambda
    svd1 <- svd(PYobs + PYnobs, nu = R, nv = R); svd1$d<-sapply(svd1$d[1:R] - lambda, function(x){max(0, x)})
    # Update unobserved entries
    Z.new <- svd1$u%*%diag(svd1$d, nrow = length(svd1$d), ncol = length(svd1$d))%*%t(svd1$v)
    # Calculate matrix of FE (note here we don't actually get each FE, just the matrix of FE_i + FE_t)
    fit1 <- lm(c((Y -  Z.new)) ~  factor(i.flat) + factor(t.flat), subset = (X.flat == 0));
    CY0  <- matrix(predict(fit1, data.frame(i.flat,t.flat)), ncol = ncol(Y), nrow = nrow(Y), byrow = TRUE);
    
    # Check if improved
    check.improve<-(norm(Z.new - Z.old, type = "F") < metric); 
    # Update convergence metric
    metric <- norm(Z.new - Z.old, type = "F"); 
    # Update unoberved entries 
    Z.old <- Z.new; 
    # Update counter
    counter <-counter+1 ; 
  }
  return(list(svd.FEout = svd1, CY = CY0))
}
# Function to calculate the Gaussian kernel function 
g.func<- function(a,b,theta = 1/2){(1/(theta*sqrt(2*pi)))*exp(-(a-b)^2/theta^2)}

library(softImpute)
library(foreach)
library(doParallel)
library(parallel)

# Linear generic setting
## Uncomment to run simulation##
# Set seed and initialise parallel clusters
set.seed(1);
# Set number of cores to use 
cl <- makeCluster(7);
registerDoParallel(cl)
# Set simulation count and matrix dimensions. These simulations are slow, N = T = 15 to 20 is reasonably fast, but gets exponentially slow larger than that
MC = 1000; 
N = 30; 
T = N
# Set values for alpha and gamma
alpha <- runif(N); 
alpha<- (ecdf(alpha)(alpha)-1/N)/(1-1/N); 
gamma <- runif(T);
gamma<- (ecdf(gamma)(gamma)-1/T)/(1-1/T)
# create matrix for g(alpha, gamma)
g<- g.func(array(alpha, dim = c(N,T)), t(array(gamma, dim = c(T,N))), theta = 1/4);
# take values of alpha and gamma to define values for X
X<-array(0, dim = c(N,T))
## Missing at random not used in paper
# xu <- array(runif(N*T), dim = c(N,T))
# X[1:round(N/2),round(T/2):T] <- 1*(array(ecdf(xu)(xu[1:round(N/2),round(T/2):T]), dim = c(round(N/2),length(round(T/2):T)))>1/2);
# Missingness used in paper: 
X[1:round(N/2),round(T/2):T] <- 1*(array(ecdf(g)(g[1:round(N/2),round(T/2):T]), dim = c(round(N/2),length(round(T/2):T)))>1/2);
# g matrix for g(X = 1)
g1<-g;
g1[X==0]<-NA;
# set rank for matrix completion and nearest neighbour calculations
approx.factors <- 5
# set of matches to consider for nearest neighbour matching
match.grid <- c(3,5,10,30,50)
# run parallel simulations
rr <- foreach(m = 1:MC) %dopar% {
  library(softImpute)
  U <- array(rnorm(N*T), dim = c(N,T))/2
  Y <- g + X + U ; YNA <- Y; YNA[X==1]<-NA
  mc <- mcsi(Y,X, lambda = 0); lambda = mc$d[approx.factors]
  fe <- femc(Y,X, lambda = 0, max.iterations = 250); lambda1 <- fe$svd.FEout$d[approx.factors]
  est1<-gtw(Y,X, R = approx.factors);
  mc <- mcsi(Y,X, lambda = lambda)
  si <- softImpute(YNA, rank.max = min(dim(YNA) - 1), lambda = lambda);
  fe <- femc(Y,X, lambda = lambda1)
  # est<-gtw(Y,M, lambda = lambda, allx = 0)
  Y0.mc <- mc$u%*%diag(mc$d, nrow = length(mc$d), ncol =length(mc$d))%*%t(mc$v); Y0.mc[X == 0] <- NA
  Y0.si <- si$u%*%diag(si$d, nrow = length(si$d), ncol =length(si$d))%*%t(si$v); Y0.si[X == 0] <- NA
  Y0.fe <- fe$svd.FEout$u%*%diag(fe$svd.FEout$d, nrow = length(fe$svd.FEout$d), ncol =length(fe$svd.FEout$d))%*%t(fe$svd.FEout$v); Y0.fe[X == 0] <- NA
  Y0.gt1 <- Y; Y0.gt1[X==0]<-NA; Y0.g11 <- Y0.gt1;
  
  Y0.gt3<-Y0.g13<-array(Y0.gt1, dim = c(dim(Y0.gt1)[1],dim(Y0.gt1)[1],length(match.grid)))
  for (m in 1:length(match.grid)){
    matches <-match.grid[m]
    gt3.placeholder<-g13.placeholder<-Y0.gt1
    for (i in which(X==1)){index<-cbind(est1$joint.feasible[[i]],est1$joint.dist[[i]]);index <- (index[order(index[,3],decreasing = FALSE),])[1:min(matches,dim(index)[1]),];gt3.placeholder[i]<-diag(Y[i - N*(ceiling(i/N)-1),index[,2]] + Y[index[,1],ceiling(i/N)] - Y[index[,1],index[,2]])%*%as.matrix(1/(index[,3]*sum(1/index[,3])))}
    for (i in which(X==1)){index<-cbind(est1$joint.feasible[[i]],est1$joint.dist[[i]]);index <- (index[order(index[,3],decreasing = FALSE),])[1:min(matches,dim(index)[1]),];g13.placeholder[i]<-diag(Y[index[,1],index[,2]])%*%as.matrix(1/(index[,3]*sum(1/index[,3])))}
    Y0.gt3[,,m]<-gt3.placeholder
    Y0.g13[,,m]<-g13.placeholder
  }
  for (i in which(X==1)){index<-est1$joint.nn[i - N*(ceiling(i/N)-1),ceiling(i/N),];Y0.gt1[i]<-Y[i - N*(ceiling(i/N)-1),index[2]] + Y[index[1],ceiling(i/N)] - Y[index[1],index[2]]}
  for (i in which(X==1)){index<-est1$joint.nn[i - N*(ceiling(i/N)-1),ceiling(i/N),];Y0.g11[i]<-Y[index[1],index[2]]}
  year<-c(col(Y)); state<-c(row(Y)); Y.sim.flat<-c(Y); X.flat<-c(X)
  fit_did1 <- lm(Y.sim.flat ~ X.flat + factor(year) + factor(state)); # Diff-in-diff estimates with time varying ATT
  Y.did1 <- Y; Y.did1[X==0] <- NA
  for (i in which(X==1)){Y.did1[i]<-Y[i] - fit_did1$coefficients["X.flat"]};
  fit_did2 <- lm(Y.sim.flat ~ X.flat:factor(year) + factor(state)); # Diff-in-diff estimates with time varying ATT
  fit_did2$coefficients[2:N]<-fit_did2$coefficients[2:N] + fit_did2$coefficients[1]
  Y.did2 <- Y; Y.did2[X==0] <- NA
  for (i in which(X==1)){Y.did2[i]<-fit_did2$coefficients[i - N*(ceiling(i/N)-1)]};
  Y0 <- Y; Y0[X == 1] <- NA;Y1 <- Y; Y1[X == 0] <- NA
  r <- array(c(Y0,U, Y1, Y0.fe, Y0.si, Y0.gt1, Y0.g11, Y.did1, Y.did2, Y0.g13,Y0.gt3) , dim = c(N,T,19))
  dimnames(r)[[3]]<-c("Y0","U","Y1","Matrix Completion FE","MC(Soft--Impute R Package)","Two-Way", "One-Way",  'DiD',  'Interactive DiD',sapply(match.grid, function(x){paste('One-Way',x, sep = '')}),sapply(match.grid, function(x){paste('Two-Way',x, sep = '')}))
  r
}
stopCluster(cl);
closeAllConnections()
r <- list(array(NA, dim = c(MC,N,T,19)), g1, g, match.grid)
for (i in 1:length(rr)){r[[1]][i,,,]<-rr[[i]]}
save(r , file = "Simulation Results/sim.Rdata")
load("Simulation Results/sim.Rdata"); 
match.grid <- c(3,5,10,30,50)
r$Y0<-r[[1]][,,,1]; 
r[[1]]<-r[[1]][,,,-1]
r$U <-r[[1]][,,,1]; 
r[[1]]<-r[[1]][,,,-1]
N <- dim(r[[2]])[1]; 
T <- dim(r[[2]])[2]; 
L <- dim(r[[1]])[4]
# Names of estimators 
estimators<-c("Mean--Differencing","MCFE","MC","TWM", "SM",  'Simple DiD',  'DiD',sapply(match.grid, function(x){paste('SM-',x, sep = '')}),sapply(match.grid, function(x){paste('TWM-',x, sep = '')}))

stdev.it.estimator<-array(NA, dim = c(N*T,L)); 
for (e in 1:dim(r[[1]])[4]){stdev.it.estimator[,e]<-c(apply(r[[1]][,,,e], MARGIN = c(2,3), sd)) };
stdev.it.estimator<-na.omit(stdev.it.estimator)
stdev.t.estimator<-array(NA, dim = c(dim(r[[1]])[1],T,L)); 
for (e in 1:dim(r[[1]])[4]){stdev.t.estimator[,,e]<-apply(r[[1]][,,,e], MARGIN = c(1,3), mean,na.rm = TRUE) };
stdev.t.estimator[,!is.na(stdev.t.estimator[1,,2]),1]<- apply(r$Y0[,,!is.na(stdev.t.estimator[1,,2])], MARGIN = c(1,3), mean,na.rm = TRUE)
stdev.t.estimator<-apply(stdev.t.estimator,MARGIN = c(2,3), sd)
stdev.t.estimator<-na.omit(stdev.t.estimator)
bias.t<- replicate(L , apply(r[[2]][,round(T/2):T], MARGIN = 2, mean, na.rm = TRUE)) - na.omit(apply(r[[1]], MARGIN = c(3,4), mean, na.rm = TRUE))
bias.t[,1]<- na.omit(apply(r[[2]],MARGIN = 2, mean,na.rm = TRUE) - apply(r$Y0,MARGIN = 2, mean,na.rm = TRUE))
bias.it<- array(apply(apply(r[[1]], MARGIN = c(2,3,4), mean), MARGIN = 3, function(x){(r[[2]] - x)}),dim = c(N,T,L))
stdev.estimator <- array(apply(apply(r[[1]], MARGIN = c(1,4), mean, na.rm = TRUE),MARGIN = 2, sd))
stdev.estimator[1] <- sd(apply(r$Y0[,,], MARGIN = c(1), mean, na.rm = TRUE))
bias<-array(apply(bias.it, MARGIN = 3, mean, na.rm = TRUE))
bias[1]<- mean(r[[2]],na.rm = TRUE) - mean(r$Y0, na.rm = TRUE) 
rmse.t.estimator<- sqrt(bias.t^2 + stdev.t.estimator^2)
rmse<-sqrt(bias^2 + stdev.estimator^2)

dimnames(rmse)[[1]]<-dimnames(stdev.estimator)[[1]]<-dimnames(bias)[[1]]<-estimators
print(xtable::xtable(as.table(cbind(bias,stdev.estimator,rmse)), type = "latex"), file = paste("Simulation Results/outputSynthetic.tex", sep =''))

## Plot bias, stdev and rmse
jpeg(file=paste("Simulation Results/generic_biasOt.jpeg", sep = ''))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
est.plot<- c(3,7,4,15,5,9)
plotbias.t <- bias.t[,est.plot]
matplot(c(round(T/2):T),plotbias.t , type = 'l', xlab = "Year",ylab = '', lwd = 2)
title("Bias")
legend(x = (T + 1.75), y =  max(plotbias.t), legend = estimators[est.plot], cex = 0.8, col = 1:length(est.plot), lty = 1:length(est.plot))
dev.off()
jpeg(file=paste("Simulation Results/generic_stdevOt.jpeg", sep = ''))
par(mar=c(5.1, 4.1, 4.1, 5.1), xpd=TRUE)
est.plot<- c(3,7,4,15,5,9)
plotsd.t <- stdev.t.estimator[,est.plot] 
matplot(c(round(T/2):T),plotsd.t, type = 'l', xlab = "Year",ylab = '', lwd = 2)
title("Standard Deviation")
y.legend <- max(plotsd.t)
dev.off()
jpeg(file=paste("Simulation Results/generic_RMSEOt.jpeg", sep = ''))
par(mar=c(5.1, 4.1, 4.1, 5.1), xpd=TRUE)
est.plot<- c(3,7,4,15,5,9)
plotrmse.t <- rmse.t.estimator[,est.plot] 
matplot(c(round(T/2):T),plotrmse.t, type = 'l', xlab = "Year",ylab = '', lwd = 2
        ,ylim = c(min(rmse.t.estimator[,c(3,4,5,16,11,6)]), max(rmse.t.estimator[,c(3,4,5,16,11,6)])))
title("RMSE")
dev.off()








