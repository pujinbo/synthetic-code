## Replication file: Liu, Wang, Xu (2022)
# Simulated Sample

sink("log/log_simdata.txt") # start log
begin.time<-Sys.time()

library(fect)
library(foreach)
library(ggplot2)
library(grid)
library(gridExtra)
library(fastplm)
library(panelView)

nboots <- 1000
set.seed(12345)
source("code/simulateData.R")
simdata <- simulateData(N = 200, TT = 35, r = 2, 
  p = 2, beta = c(1, 3), force = 3, mu = 5, Rtype = "n", Ftype = c("trend","white"), 
  Fsize = c(1.5,1), eff.size = 0.4, eff.noise = 0.2, seed = 12345)

## Factors
plot(simdata$F1[1:30], type = "l", ylim = c(-40,40))
lines(simdata$F2[1:30],  col = 2)

## Treated probability
panelview(tr_prob ~ D, data = simdata, index = c("id","time"), type = "outcome",
          axis.lab = "time", xlab = "Time", ylab = "Unit", main = "Treated Probability")

panelview(Y ~ D, data = simdata, index = c("id","time"), by.timing = TRUE,
          axis.lab = "time", xlab = "Time", ylab = "Unit")

panelview(Y ~ D, data = simdata, index = c("id","time"), 
          axis.lab = "time", xlab = "Time", ylab = "Unit")

panelview(eff ~ D, data = simdata, index = c("id","time"), type = "outcome",
  axis.lab = "time", xlab = "Time", ylab = "Unit", main = "Treated Probability")


##################
## look at data
##################


## treatment
panelview(Y ~ D, data = simdata, index = c("id","time"), by.timing = TRUE,
  axis.lab = "time", xlab = "Time", ylab = "Unit", show.id = c(1:40), 
  color = c("gray80","gray50"), background= "white", theme.bw = TRUE)

## outcome
panelview(Y ~ D, data = simdata, index = c("id","time"), 
  axis.lab = "time", xlab = "Time", ylab = "Outcome", type = "raw", by.group = 0)



# fect
out.fe <- fect(Y~D+X1+X2, data = simdata, index = c("id", "time"),
  method = "fe", force = 3, parallel = 1, loo = 1,
  tol = 1e-4, se = 1, nboots = nboots)

plot(out.fe, ylim = c(-5,5), xlim = c(-20, 6), loo = 1)

out.fe.p <- fect(Y~D+X1+X2, data = simdata, index = c("id", "time"),
  method = "fe", force = 3, parallel = 1,
  tol = 1e-4, se = 1, nboots = nboots, 
  placeboTest = 1, placebo.period = c(-2, 0))

plot(out.fe.p, ylim = c(-5,5), xlim = c(-20, 6))

out.fe.co <- fect(Y~D+X1+X2, data = simdata, index = c("id", "time"),
  method = "fe", force = 3, parallel = 1,
  tol = 1e-4, se = 1, nboots = nboots,
  carryoverTest = TRUE, carryover.period = c(1,3))

plot(out.fe.co, ylim = c(-5,5), xlim = c(-20, 6))

# ife
out.ife <- fect(Y~D+X1+X2, data = simdata, index = c("id", "time"),
  method = "ife", force = 3, r = 2, CV = 0, parallel = 1, loo = 1,
  tol = 1e-4, se = 1, nboots = nboots)

plot(out.ife, ylim = c(-5,5), xlim = c(-20, 6), loo = 1)


out.ife.p <- fect(Y~D+X1+X2, data = simdata, index = c("id", "time"),
  method = "ife",force = 3, r = 2, CV = 0, parallel = 1,
  tol = 1e-4, se = 1, nboots = nboots, 
  placeboTest = 1, placebo.period = c(-2, 0))

plot(out.ife.p, ylim = c(-10,10), xlim = c(-20, 5))

out.ife.co <- fect(Y~D+X1+X2, data = simdata, index = c("id", "time"),
  method = "ife", force = 3, parallel = 1,
  tol = 1e-4, se = 1, nboots = nboots, r = 2, CV = FALSE,
  carryoverTest = TRUE, carryover.period = c(1,3))

plot(out.ife.co, ylim = c(-5,5), xlim = c(-20, 6))

# mc
out.mc <- fect(Y~D+X1+X2, data = simdata, index = c("id", "time"),
  method = "mc", force = 3, CV = 1, parallel = 1, loo = 1,
  tol = 1e-4, se = 1, nboots = nboots)

plot(out.mc, ylim = c(-10,10), xlim = c(-20, 5), loo = 1)

out.mc.p <- fect(Y~D+X1+X2, data = simdata, index = c("id", "time"),
  method = "mc", force = 3, r = 1, CV = 0, parallel = 1, 
  lambda = out.mc$lambda.cv,
  tol = 1e-4, se = 1, nboots = nboots, 
  placeboTest = 1, placebo.period = c(-2, 0))

plot(out.mc.p, ylim = c(-10,10), xlim = c(-20, 5))

out.mc.co <- fect(Y~D+X1+X2, data = simdata, index = c("id", "time"),
  method = "mc", force = 3, parallel = 1,
  tol = 1e-4, se = 1, nboots = nboots, lambda = out.mc$lambda.cv, 
  CV = FALSE, carryoverTest = TRUE, carryover.period = c(1,3))

plot(out.mc.co, ylim = c(-5,5), xlim = c(-20, 6))

save(simdata, out.fe, out.fe.p, out.fe.co, 
  out.ife, out.ife.p, out.ife.co, 
  out.mc, out.mc.p, out.mc.co,
  file = "results/simexample.RData")

##############################
## Save Figures
##############################

library(fect)
library(ggplot2)
library(grid)
library(gridExtra)
library(panelView)
load("results/simexample.RData")
d_ylim <- c(-4, 6)
e_ylim <- c(-4, 4)
xlim <- c(-15,5)
text.pos <- c(-15,5)
bound <- 0.36 * sqrt(out.fe$sigma2)
eff <- aggregate(simdata$eff, list(simdata$tr_cum), FUN=mean)[-1,2]
x1 <- seq(xlim[1], xlim[2], by = 1)
y1 <- c(rep(0, (1-xlim[1])),eff[1:xlim[2]])
att <- cbind.data.frame(x1, y1)

## gap plot
xlab <- "Time Since the Treatment Started"
ylab <- "Effect of D on Y"
p.gap1 <- plot(out.fe, type = "gap", ylab = ylab, xlab = xlab, ylim = d_ylim, xlim = xlim,
  main = "FEct", theme.bw = TRUE, p.value = 1,  cex.text = 0.8)  + 
	geom_line(aes(x1, y1), colour = "red", att, linetype = "dashed")
p.gap2 <- plot(out.ife, type = "gap", ylab = ylab, xlab = xlab, ylim = d_ylim, xlim = xlim,
  main = "IFEct", theme.bw = TRUE, p.value = 1,  cex.text = 0.8) +
	geom_line(aes(x1, y1), colour = "red", att, linetype = "dashed")
p.gap3 <- plot(out.mc, type = "gap", ylab = ylab, xlab = xlab, ylim = d_ylim, xlim = xlim,
  main = "MC", theme.bw = TRUE,  cex.text = 0.8) + 
	geom_line(aes(x1, y1), colour = "red", att, linetype = "dashed")
margin = theme(plot.margin = unit(c(1,1,1,1), "line"))
p.gap <- grid.arrange(grobs = lapply(list(p.gap1, p.gap2, p.gap3), "+", margin),
  ncol = 3, widths = c (1, 1, 1))
ggsave('graph/sim0_gap.pdf', p.gap, width = 17, height = 6)

## placebo test
ylab <- "Effect of D on Y"
stats.labs <- c("t test p-value","TOST p-value")
p.placebo1 <- plot(out.fe.p, type = "gap", ylab = ylab, xlab = xlab, placeboTest = 1, 
  ylim = d_ylim, xlim = xlim, main = "FEct", text.pos = NULL, 
  theme.bw = TRUE, stats = c("placebo.p","equiv.p"), stats.labs = stats.labs)
p.placebo1 <- p.placebo1 +  geom_line(aes(x1, y1), colour = "red", att, linetype = "dashed")
# IFE
p.placebo2 <- plot(out.ife.p, type = "gap", ylab = ylab, xlab = xlab, placeboTest = 1, 
  ylim = d_ylim, xlim = xlim, main = "IFEct", theme.bw = TRUE, 
  stats = c("placebo.p","equiv.p"), stats.labs = stats.labs)
p.placebo2 <- p.placebo2 + geom_line(aes(x1, y1), colour = "red", att, linetype = "dashed")
# MC
p.placebo3 <- plot(out.mc.p, type = "gap", ylab = ylab, xlab = xlab, placeboTest = 1, 
  ylim = d_ylim, xlim = xlim, main = "MC", theme.bw = TRUE, 
  stats = c("placebo.p","equiv.p"), stats.labs = stats.labs)
p.placebo3 <- p.placebo3 + geom_line(aes(x1, y1), colour = "red", att, linetype = "dashed")
# merge
margin = theme(plot.margin = unit(c(1,1,1,1), "line"))
p.placebo <- grid.arrange(grobs = lapply(list(p.placebo1, p.placebo2, p.placebo3), "+", margin),
  ncol = 3, widths = c (1, 1, 1))
ggsave('graph/sim0_placebo.pdf', p.placebo, width = 17, height = 6)

## carryover test
xlab <- "Time Since the Treatment Ended"
ylab <- "Effect of D on Y"
stats.labs <- c("t test p-value","TOST p-value")
stats.pos <- c(0.5,6)
p1.carryover <- plot(out.fe.co, type = "exit", ylab = ylab, xlab = xlab,  
  ylim = d_ylim, xlim = xlim, main = "FEct", theme.bw = TRUE, stats.pos = stats.pos,
  stats = c("carryover.p","equiv.p"), stats.labs = stats.labs)
# IFE
p2.carryover <- plot(out.ife.co, type = "exit", ylab = ylab, xlab = xlab, 
  ylim = d_ylim, xlim = xlim, main = "IFEct", theme.bw = TRUE, stats.pos = stats.pos, 
  stats = c("carryover.p","equiv.p"), stats.labs = stats.labs)
# MC
p3.carryover <- plot(out.mc.co, type = "exit", ylab = ylab, xlab = xlab, 
  ylim = d_ylim, xlim = xlim, main = "MC", theme.bw = TRUE, stats.pos = stats.pos, 
  stats = c("carryover.p","equiv.p"), stats.labs = stats.labs)
# merge
margin = theme(plot.margin = unit(c(1,1,1,1), "line"))
p.carryover <- grid.arrange(grobs = lapply(list(p1.carryover, p2.carryover, p3.carryover), "+", margin),
  ncol = 3, widths = c (1, 1, 1))
ggsave('graph/sim0_carryover.pdf', p.carryover, width = 17, height = 6)


# equivalence test
xlab <- "Time Since the Treatment Started"
ylab <- "Average Prediction Error"
stats.labs <- c("F test p-value","TOST max p-value")
p.equiv1 <- plot(out.fe, type = "equiv", ylab = ylab, xlab = xlab, 
  ylim = e_ylim, xlim = c(-20, 0), main = "FEct", theme.bw = TRUE, 
  legendOff = TRUE, stats.labs = stats.labs, loo = 1)
p.equiv2 <- plot(out.ife, type = "equiv", ylab = ylab, xlab = xlab, 
  ylim = e_ylim, xlim = c(-20, 0), main = "IFEct", theme.bw = TRUE, 
  legendOff = TRUE, stats.labs = stats.labs, loo = 1)
p.equiv3 <- plot(out.mc, type = "equiv", ylab = ylab, xlab = xlab, 
  ylim = e_ylim, xlim = c(-20, 0), main = "MC", theme.bw = TRUE, 
  legendOff = TRUE, stats.labs = stats.labs, loo = 1)
margin = theme(plot.margin = unit(c(1,1,1,1), "line"))
p.equiv <- grid.arrange(grobs = lapply(list(p.equiv1, p.equiv2, p.equiv3), "+", margin),
  ncol = 3, widths = c (1, 1, 1))
ggsave('graph/sim0_equiv.pdf', p.equiv, width = 17, height = 6)


## Outcome and treatment
p.outcome <- panelview(Y ~ D, data = simdata, index = c("id","time"), theme.bw = TRUE,
  axis.lab = "time", xlab = "Time", ylab = "Outcome", type = "raw", by.group = 0)
ggsave('graph/sim0_outcome.pdf', p.outcome, width = 8, height = 6)


p.treat <- panelview(Y ~ D, data = simdata, index = c("id","time"), 
  axis.lab = "time", xlab = "Time", ylab = "Unit", show.id = c(1:100), 
  background= "white", theme.bw = TRUE)
ggsave('graph/sim0_treat.pdf', p.treat, width = 8, height = 5)

print(Sys.time()-begin.time) 
sink() # close log


