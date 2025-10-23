## Replication file: Liu, Wang, Xu (2022)
# Empirical example: Hainmueller and Hangartner (2015)

sink("log/log_hh2015.txt") # start log
begin.time<-Sys.time()

require(readstata13)
require(fect)
library(grid)
library(gridExtra)
library(ggplot2)
library(fastplm)
library(panelView)

o.data <- read.dta13("data/hh2015.dta")
dim <- nrow(o.data)
dim(o.data)
Y <- "nat_rate_ord"
D <- "indirect"
unit <- "bfs"
period <- "year"
FE <- c(unit, period)
cl <- unit
seed <- 1453
d <- o.data[complete.cases(o.data[,c(Y,D,FE)]), c(Y,D,FE)]

# statistics for the table
FEct.formula <- as.formula(paste0(Y, "~", D))
ylimit <- c(-3, 6)
elimit <- c(-0.4, 0.4)
nboots <- 1000

## Fixed effects
out.fe <- fastplm(FEct.formula, data = d, index=FE, se = 1, nboots = nboots)
summary(out.fe)


## FEct/IFEct
out.fect <- fect(formula = as.formula(FEct.formula), data = d, method = "ife",
  index=FE, tol = 1e-5, se = 1, nboots = nboots, loo = 1,
  r = 0, CV = FALSE, force = "two-way", parallel = 1)

out.fect.p <- fect(formula = as.formula(FEct.formula), data = d, method = "fe",
 index=FE, tol = 1e-5, se = 1, nboots = nboots, r = 0,
 CV = FALSE, force = "two-way", parallel = 1,
 placeboTest = 1, placebo.period = c(-2, 0))


## MC give essentially the same result as FEct (lambda is bigger than the first eigenvalue)
out.mc <- fect(as.formula(FEct.formula), data = d, method = "mc",
 index=FE, tol = 1e-3, se = TRUE, nboots = nboots,
 cv.prop = 0.1, cv.nobs = 3, cv.treat = 1, loo = 0,
 CV = TRUE, force = "two-way", parallel = 1)

save(d, FEct.formula, out.fe, out.fect, out.fect.p, out.mc, 
  file = "results/ex_HH2015.RData")


## Plotting
library(panelview)
library(ggplot2)
require(fect)
load("results/ex_HH2015.RData")
d_ylim <- c(-4, 4)


## Outcome
p.outcome <- panelview(FEct.formula, data = d, index = c("bfs","year"), type = "raw", 
	axis.lab = "time", xlab = "Year", ylab = "Naturalization Rate (%)", main = "", by.group = TRUE,
  legend.labs = c("Direct Democracy", "Indirect Democracy"))
ggsave('graph/ex_HH2015_outcome.pdf', p.outcome, width = 8, height = 8)

## treatment
p.treat <- panelview(FEct.formula, data = d, 
  index = c("bfs","year"), by.timing = TRUE, 
  axis.lab = "time", xlab = "Year", ylab = "Municipality", main = "", gridOff = TRUE,
  background= "white", 
  legend.labs = c("Direct Democracy", "Indirect Democracy"))
ggsave('graph/ex_HH2015_treat.pdf', p.treat, width = 8, height = 6)

## FEct 
xlab <- "Time Since the Treatment Started"
p.gap <- plot(out.fect, type = "gap", ylab = "Effect of Indirect Democracy on Naturalization Rate (%)", ylim = d_ylim, 
  xlim = c(-14,5), main = "", theme.bw = TRUE, xlab = xlab,
  p.value = TRUE, text.pos = c(-14.5,3)) 
ggsave('graph/ex_HH2015_gap.pdf', p.gap, width = 6, height = 6)


p.equiv <- plot(out.fect, type = "equiv", ylab = "Average Prediction Error", ylim = d_ylim, 
  xlim = c(-14,0), main = "", theme.bw = TRUE, xlab = xlab,
  stats.labs = c("F test p-value","TOST max p-value"),
  legendOff = TRUE, loo = 1) 
ggsave('graph/ex_HH2015_equiv.pdf', p.equiv, width = 6, height = 6)
 

p.placebo <- plot(out.fect.p, type = "gap", ylab = "Effect of Indirect Democracy on Naturalization Rate (%)", 
  placeboTest = 1, placebo.period = c(-1, 0), xlab = xlab,
  stats.labs = c("t test p-value","TOST p-value"),
  xlim = c(-14,5), ylim = d_ylim,  label.pos = c(-12.5, 3.8), 
  main = "", theme.bw = TRUE, show.stats = TRUE, stats = c("placebo.p","equiv.p"))
ggsave('graph/ex_HH2015_placebo.pdf', p.placebo, width = 6, height = 6)


## Coefficients
summary(out.fe)
out.fect$est.avg
#out.ife$est.avg
out.mc$est.avg
out.mc$lambda.cv
out.mc$CV.out.mc

print(Sys.time()-begin.time) 
sink()






