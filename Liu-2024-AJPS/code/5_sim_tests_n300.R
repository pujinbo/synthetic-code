## Replication file: Liu, Wang, Xu (2022)
# Compare performances between F test and the equivalence test (N = 300)

sink("log/log_sim_tests_n300.txt") # start log


begin.time<-Sys.time()

library(fect)
library(ggplot2)
library(grid)
library(gridExtra)
library(fastplm)
library(future)
library(doParallel)
cores <- detectCores()
cl2 <- makeClusterPSOCK(cores,verbose = FALSE)
registerDoParallel(cl2)


# load DGP
source("code/simulateDID.R")


## ---------------------------- ##
## fixed units, vary bias.      ##
## ---------------------------- ##

set.seed(1234)
TT <- 31
ts <- arima.sim(list(order = c(1,0,0), ar = 0.5), n = (TT+20))[-c(1:20)] + 0.3 * c(1:TT)
Factor <- matrix(scale(ts), TT, 1)
#plot(Factor, type = "l")

## Example

# f_size <- 0.5
# data <- simulateData(N = n, TT = TT, r = 1, 
# 	tr.threshold = 0.5, tr.start = 21,
# 	p = 0, force = 3, mu = 5, Rtype = "u", factor = Factor,
# 	Fsize = f_size, eff.size = 0, eff.noise = 1)			

# out_fect <- fect(Y~D, data = data, index = c("id", "time"),
# 	method = "fe", force = 3, r = 0, CV = 0, parallel = 0,
# 	tol = 1e-4, se = 1, nboots = 200)
#plot(out_fect)


library(abind)
f <- function(){
  function(...) abind(...,along=3)  
}


## repeated samples
K <- 20
nsims <- 30
ntot <- K * nsims

n <- 300
f_size_seq <- c(seq(0, 0.3, by = 0.05), seq(0.4, 0.8, by = 0.1))
num_f <- length(f_size_seq)


sim.att.avg <- matrix(NA, ntot, num_f)
sim.F <- matrix(NA, ntot, num_f)
sim.F.p <- matrix(NA, ntot, num_f)
sim.equiv.p <- matrix(NA, ntot, num_f)
sim.y.sd <- matrix(NA, ntot, num_f)
sim.res.sd <- matrix(NA, ntot, num_f)


for (k in 1:K) {

	result <- foreach(i = 1:nsims,.combine=f(),.inorder=FALSE,
	    .export = c("fect"), .packages = c("fect")) %dopar% {

		out <- matrix(NA, 6, num_f) # att.avg, F stat, F p-value, y.sd, res.sd

		for (j in 1:num_f){

			f_size <- f_size_seq[j]	
			data <- simulateData(N = n, TT = TT, r = 1, 
				tr.threshold = 0.5, tr.start = 21,
				p = 0, force = 3, mu = 5, Rtype = "u", factor = Factor,
				Fsize = f_size, eff.size = 0, eff.noise = 1)			
			
			out_fect <- fect(Y~D, data = data, index = c("id", "time"),
				method = "fe", force = 3, r = 0, CV = 0, parallel = 0,
				tol = 1e-4, se = 1, nboots = 200)

			# att
			out[1,j] <- out_fect$att.avg

			## F
			out[2,j] <- out_fect$test.out$f.stat

			## F p-values 
			out[3,j] <- out_fect$test.out$f.p

			## equivalence p-value
			out[4,j] <- out_fect$test.out$tost.equiv.p

	    ## y sd /res sd
	    out[5,j] <- sd(data$Y[which(data$D == 0)])
	    out[6,j] <- sqrt(out_fect$sigma2)	
			
		}
		return(out)		
	}

	# save
	part <- 1:nsims + (k-1)*nsims
	sim.att.avg[part,] <- t(result[1,,])
	sim.F[part,] <- t(result[2,,])
	sim.F.p[part,] <- t(result[3,,])
	sim.equiv.p[part,] <- t(result[4,,])
	sim.y.sd[part,] <- t(result[5,,])
	sim.res.sd[part,] <- t(result[6,,])

	save(f_size_seq, sim.att.avg, sim.F, sim.F.p, 
		sim.equiv.p, sim.y.sd, sim.res.sd,
		file = paste0("results/sim_tests_n",n,".RData"))
	print(Sys.time()-begin.time) 

}

stopCluster(cl2)

sink()





