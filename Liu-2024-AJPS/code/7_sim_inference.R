## Replication file: Liu, Wang, Xu (2022)
# Inferential methods

sink("log/log_inference.txt") # start log
begin.time <- Sys.time()

library(fect)
library(ggplot2)
library(grid)
library(gridExtra)
library(fastplm)
library(extraDistr)
library(doParallel)
library(doRNG)
library(future)


# load DGP: staggered adoption 
source("code/simulateDID.R")
seed <- 12345
cores <- detectCores()
para.clusters <- future::makeClusterPSOCK(cores)
registerDoParallel(para.clusters)
if (is.null(seed) == FALSE) {
    registerDoRNG(seed)
}

## ---------------- N = 50 case ---------------- 
cat("Parallel computing ...\n")

nboots <- 1000
TT <- 20 
N <- 50

## save results 
pdata_lmboot <- pdata_boot <- pdata_jack <- rep(NA, nboots)


oneboot <- function() {

	data <- simulateData(N = N, TT = TT, r = 0, 
		tr.threshold = seq(0.5,0.9,0.1), tr.start = c(19, 17, 15, 13, 11), 
		p = 0, force = 3, mu = 5, Rtype = "u", Ftype = "arma", 
		Fsize = 0, eff.size = 2, eff.noise = 1)

	# library(panelView)
	# panelview(Y ~ D, data = data, index = c("id","time"), by.timing = TRUE)
		
    ## true att 
    att <- mean(data$eff[which(data$D == 1)])

    ## two-way fe model 
    out_fe <- fastplm(Y~D, data = data, index = c("id", "time"),
			                se = 1, vce = "bootstrap", cluster = "id",
                      wild = FALSE, nboots = 200, parallel = 0)

    sub.pdata_lmboot <- (out_fe$est.coefficients[1] - att) / out_fe$est.coefficients[2]
    
	## block bootstrap
	out_fect1 <- fect(Y~D, data = data, index = c("id", "time"),
		method = "fe", vartype = "bootstrap",
		force = 3, r = 0, CV = 0, parallel = 0,
		tol = 1e-4, se = 1, nboots = 200)

	sub.pdata_boot <- (out_fect1$est.avg[1] - att) / out_fect1$est.avg[2]
    
	## jackknife
	out_fect2 <- fect(Y~D, data = data, index = c("id", "time"),
		method = "fe", vartype = "jackknife",
		force = 3, r = 0, CV = 0, parallel = 0,
		tol = 1e-4, se = 1, nboots = 200)

	sub.pdata_jack <- (out_fect2$est.avg[1] - att) / out_fect2$est.avg[2]
    

    output <- list(sub.pdata_lmboot = sub.pdata_lmboot,
    	             sub.pdata_boot = sub.pdata_boot, 
    	             sub.pdata_jack = sub.pdata_jack)

    return(output)
}

boot.out <- foreach(j=1:nboots, 
                    .inorder = FALSE,
                    .export = c("fastplm", "fect", "rbern"),
                    .packages = c("fastplm", "fect", "base")
                    ) %dopar% {
                        return(oneboot())
                    }


for (j in 1:nboots) { 
    pdata_lmboot[j] <- boot.out[[j]]$sub.pdata_lmboot
    pdata_boot[j] <- boot.out[[j]]$sub.pdata_boot
    pdata_jack[j] <- boot.out[[j]]$sub.pdata_jack
}
	
	
save(pdata_lmboot, pdata_boot, pdata_jack, file = "results/qqplots_N50.RData")
print(Sys.time()-begin.time) 

## ---------------- N = 100 case ---------------- 
cat("Parallel computing ...\n")

nboots <- 1000
TT <- 20 
N <- 100

## save results 
pdata_lmboot <- pdata_boot <- pdata_jack <- rep(NA, nboots)


oneboot <- function() {

	data <- simulateData(N = N, TT = TT, r = 0, 
		tr.threshold = seq(0.5,0.9,0.1), tr.start = c(19, 17, 15, 13, 11), 
		p = 0, force = 3, mu = 5, Rtype = "u", Ftype = "arma", 
		Fsize = 0, eff.size = 2, eff.noise = 1)

	# library(panelView)
	# panelview(Y ~ D, data = data, index = c("id","time"), by.timing = TRUE)
		
    ## true att 
    att <- mean(data$eff[which(data$D == 1)])

    ## two-way fe model 
    out_fe <- fastplm(Y~D, data = data, index = c("id", "time"),
			                se = 1, vce = "bootstrap", cluster = "id",
                      wild = FALSE, nboots = 200, parallel = 0)

    sub.pdata_lmboot <- (out_fe$est.coefficients[1] - att) / out_fe$est.coefficients[2]
    
	## block bootstrap
	out_fect1 <- fect(Y~D, data = data, index = c("id", "time"),
		method = "fe", vartype = "bootstrap",
		force = 3, r = 0, CV = 0, parallel = 0,
		tol = 1e-4, se = 1, nboots = 200)

	sub.pdata_boot <- (out_fect1$est.avg[1] - att) / out_fect1$est.avg[2]
    
	## jackknife
	out_fect2 <- fect(Y~D, data = data, index = c("id", "time"),
		method = "fe", vartype = "jackknife",
		force = 3, r = 0, CV = 0, parallel = 0,
		tol = 1e-4, se = 1, nboots = 200)

	sub.pdata_jack <- (out_fect2$est.avg[1] - att) / out_fect2$est.avg[2]
    

    output <- list(sub.pdata_lmboot = sub.pdata_lmboot,
    	             sub.pdata_boot = sub.pdata_boot, 
    	             sub.pdata_jack = sub.pdata_jack)

    return(output)
}

boot.out <- foreach(j=1:nboots, 
                    .inorder = FALSE,
                    .export = c("fastplm", "fect", "rbern"),
                    .packages = c("fastplm", "fect", "base")
                    ) %dopar% {
                        return(oneboot())
                    }


for (j in 1:nboots) { 
    pdata_lmboot[j] <- boot.out[[j]]$sub.pdata_lmboot
    pdata_boot[j] <- boot.out[[j]]$sub.pdata_boot
    pdata_jack[j] <- boot.out[[j]]$sub.pdata_jack
}
	
	
save(pdata_lmboot, pdata_boot, pdata_jack, file = "results/qqplots_N100.RData")
print(Sys.time()-begin.time) 

stopCluster(para.clusters)
sink() # close log

