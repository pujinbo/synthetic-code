## Replication file: Liu, Wang, Xu (2022)
# Compare IFEct and MC

sink("log/log_sim_ife_mc.txt") # start log
begin.time<-Sys.time()

library(fect)
library(ggplot2)
library(grid)
library(gridExtra)
library(fastplm)
library(future)
library(doParallel)
library(abind)
cores <- detectCores()
cl2 <- makeClusterPSOCK(cores,verbose = FALSE)
registerDoParallel(cl2)

# load DGP
source("code/simulateDID.R")

## IFE vs. MC
num_factors <- c(1:9)
times <- 20
K <- 25
nboots <- times * K

## note: av for att.avg; t for att_t; else for individual effects
ifeo_mse <-ife_mse <- mc_mse <- matrix(NA, nboots, 9)
ifeo_mse_t <-ife_mse_t <- mc_mse_t <- matrix(NA, nboots, 9)

## other parameters
ifeo_iter <- ife_iter <- mc_iter <- ife_r <- matrix(NA, nboots, 9)

  
oneboot <- function(){

  ifeo_iter <- ifeo_mse <- ifeo_mse_t <- rep(NA,9)
  ife_r <- ife_iter <- ife_mse <- ife_mse_t <- rep(NA,9) 
  mc_iter <- mc_mse <- mc_mse_t <- rep(NA,9) 
  
  for (i in 1:9) { 
    
    num_f <- num_factors[i]
    f_size <- sqrt(3/num_f)

    ## same timing
    data <- simulateData(N = 200, TT = 30, r = num_f, 
     tr.threshold = 0.5, tr.start = 21,
     p = 0, force = 3, mu = 5, 
     Rtype = "n", Ftype = "ar1", 
     Fsize = f_size, eff.size = 0, eff.noise = 0)

    # library(panelView)
    # panelview(Y ~ D, index = c("id","time"), data = data, by.timing = TRUE)        
    
    ## oracle ife
    out_ifeo <- fect(Y~D, data = data, index = c("id", "time"),
      method = "ife", force = 3, r = num_f, CV = 0, parallel = 0,
      tol = 1e-4, se = 0)

    ifeo_eff <- as.vector(out_ifeo$eff[out_ifeo$D.dat == 1])
    ifeo_mse[i] <- mean((ifeo_eff)^2)    
    ifeo_mse_t[i] <- out_ifeo$att[30]^2
    ifeo_iter[i] <- out_ifeo$niter
    
    ## ife cv
    out_ife <- fect(Y~D, data = data, index = c("id", "time"),
                    method = "ife",
                    force = 3, r = c(0, 10), CV = 1, parallel = 0,
                    tol = 1e-4, se = 0, 
                    cv.prop = 0.2, cv.nobs = 1, cv.treat = 0)

    ife_eff <- as.vector(out_ife$eff[out_ife$D.dat == 1])
    ife_mse[i] <- mean((ife_eff)^2)    
    ife_mse_t[i] <- out_ife$att[30]^2
    ife_iter[i] <- out_ife$niter
    ife_r[i] <- out_ife$r.cv    
    
    ## mc
    out_mc <- fect(Y~D, data = data, index = c("id", "time"),
                   method = "mc",
                   force = 3, CV = 1, parallel = 0,
                   tol = 1e-4, se = 0, nboots = 200, 
                   cv.prop = 0.2, cv.nobs = 1, cv.treat = 0)

    mc_eff <- as.vector(out_mc$eff[out_mc$D.dat == 1])
    mc_mse[i] <- mean((mc_eff)^2)    
    mc_mse_t[i] <- out_mc$att[30]^2
    mc_iter[i] <- out_mc$niter

  }

  out <- cbind(
    ifeo_mse, ife_mse, mc_mse, 
    ifeo_mse_t, ife_mse_t, mc_mse_t,
    ifeo_iter, ife_iter, mc_iter, ife_r) 

}

## start ----------------

f <- function(){
  function(...) abind(...,along=3)  
}

cores <- detectCores()
para.clusters <- makeCluster(cores)
registerDoParallel(para.clusters)

## parallel in batches
for (k in 1:K) {

	boot.out <- foreach(j = 1:times, .combine=f(), .inorder = FALSE,
		.export = c("fect"), .packages = c("fect")
		) %dopar% {
		return(oneboot())
	}

	for (j in 1:times) { 
		i <- (k - 1)*times + j
		ifeo_mse[i,] <- boot.out[,1,j]
    ife_mse[i,] <- boot.out[,2,j]
    mc_mse[i,] <- boot.out[,3,j]
    ifeo_mse_t[i,] <- boot.out[,4,j]
		ife_mse_t[i,] <- boot.out[,5,j]
    mc_mse_t[i,] <- boot.out[,6,j]
    ifeo_iter[i,] <- boot.out[,7,j]
		ife_iter[i,] <- boot.out[,8,j]
    mc_iter[i,] <- boot.out[,9,j]
		ife_r[i,] <- boot.out[,10,j]
	}

  ## summary results
  mean_ife_mse_t <- apply(ife_mse_t, 2, mean, na.rm = TRUE)
  mean_ifeo_mse_t <- apply(ifeo_mse_t, 2, mean, na.rm = TRUE)
  mean_mc_mse_t <- apply(mc_mse_t, 2, mean, na.rm = TRUE)

  mean_ife_mse <- apply(ife_mse, 2, mean, na.rm = TRUE)
  mean_ifeo_mse <- apply(ifeo_mse, 2, mean, na.rm = TRUE)
  mean_mc_mse <- apply(mc_mse, 2, mean, na.rm = TRUE)


  save(ifeo_mse, ife_mse, mc_mse,
    ifeo_mse_t, ife_mse_t, mc_mse_t, 
    ifeo_iter, ife_iter, mc_iter, ife_r,
    file = "results/sim_ife_mc.RData")

  # check results
  pdf("results/est_sim_compare_check.pdf")
  plot(1:9, mean_ife_mse, type = "l", main = paste("k =",k), ylim = c(0,8))
  lines(1:9, mean_ifeo_mse, col = 2)
  lines(1:9, mean_mc_mse, col = 4)
  graphics.off()

  cat(k,"\n")
  print(Sys.time()-begin.time) 

}


## end -------------------------------------

stopCluster(cl2)
sink()
