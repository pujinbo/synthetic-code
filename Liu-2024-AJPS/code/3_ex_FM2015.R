## Replication file: Liu, Wang, Xu (2022)
# Empirical example: Fouirnaies and Mutlu-Eren (2015)

sink("log/log_fm2015.txt") # start log
begin.time<-Sys.time()

require(readstata13)
require(fect)
library(grid)
library(gridExtra)
library(ggplot2)
library(fastplm)
library(panelView)

o.data <- read.dta13("data/FM2015.dta")
Y <- "logSgwaPercap"
D <- "treat"
X <- NULL
unit <- "councilnumber"
period <- "year"
FE <- c(unit, period)

d <- o.data[complete.cases(o.data[,c(Y,D,X,FE)]),]
d <- d[order(d[,unit], d[,period]),]
fect.formula <- as.formula(paste0(Y,"~",D))
dim(d)

panelview(fect.formula, data = d, index = FE,
  axis.lab = "time", by.timing = 1, show.id = 1:100)

panelview(fect.formula, data = d, index = FE, type = "outcome")

nboots <- 1000

## fixed-effects
set.seed(1234)
out.fe <- fastplm(formula = as.formula(fect.formula), data = d,
                  index=FE, se = 1, vce = "boot", nboots = nboots, cluster = unit)
summary(out.fe)

######################################
# Do not consider carry over effects
######################################

## Main Estimation
out.fect <- fect(formula = as.formula(fect.formula), data = d, method = "fe",
 index=FE, tol = 1e-3, se = 1, nboots = nboots, r = 0, min.T0 = 4,
 loo = 1, pre.period = c(-6,0), group = 'cohort', 
 CV = FALSE, force = "two-way", parallel = 1, cores = 4)

out.ife <- fect(formula = as.formula(fect.formula), data = d, method = "ife",
  index=FE, tol = 1e-3, se = 1, nboots = nboots, min.T0 = 4,
  loo = 1, pre.period = c(-6,0), group = 'cohort', 
  CV = FALSE, r =2, force = "two-way", parallel = 1, cores = 4)

out.mc <- fect(as.formula(fect.formula), data = d, method = "mc",
 index=FE, tol = 1e-4, se = 1, nboots = nboots, min.T0 = 4, 
 loo = 1, pre.period = c(-6,0), group = 'cohort', 
 CV = 1, force = "two-way", parallel = 1, cores = 4)

out.mc$lambda.cv/out.mc$eigen.all[1] # 0.0749

## Placebo tests

out.fect.p <- fect(formula = as.formula(fect.formula), data = d, method = "fe",
 index=FE, tol = 1e-3, se = 1, nboots = nboots, r = 0, min.T0 = 1,
 CV = FALSE, force = "two-way", parallel = 1, cores = 4, placeboTest = 1, 
 placebo.period = c(-2, 0))

out.ife.p <- fect(formula = as.formula(fect.formula), data = d, method = "ife",
  index=FE, tol = 1e-3, se = 1, nboots = nboots, r = 2, min.T0 = 1, 
  CV = 0, force = "two-way", parallel = 1, cores = 4, 
  placeboTest = 1, placebo.period = c(-2, 0))

out.mc.p <- fect(formula = as.formula(fect.formula), data = d, method = "mc",
 index=FE, tol = 1e-4, se = 1, nboots = nboots, lambda = out.mc$lambda.cv, 
 CV = 0, force = "two-way", parallel = 1, cores = 4, min.T0 = 1, 
 placeboTest = 1, placebo.period = c(-2, 0))


## Tests for carry-over effects

out.fect.c <- fect(formula = as.formula(fect.formula), data = d, method = "fe",
 index=FE, tol = 1e-3, se = 1, nboots = nboots, force = "two-way", min.T0 = 1, 
 parallel = 1, cores = 4, carryoverTest = 1, carryover.period = c(1, 5))

out.ife.c <- fect(formula = as.formula(fect.formula), data = d, method = "ife",
  index=FE, tol = 1e-3, se = 1, nboots = nboots, r = 2, min.T0 = 1,
  CV = 0, force = "two-way", parallel = 1, cores = 4, 
  carryoverTest = TRUE, carryover.period = c(1,5))

out.mc.c <- fect(formula = as.formula(fect.formula), data = d, method = "mc",
 index=FE, tol = 1e-4, se = 1, nboots = nboots, lambda = out.mc$lambda.cv, 
 CV = 0, force = "two-way", parallel = 1, cores = 4, min.T0 = 1,
 carryoverTest = 1, carryover.period = c(1, 5))


##########################################################
# Analysis after removing three post-treatment periods
##########################################################


## Main Estimation
out2.fect <- fect(formula = as.formula(fect.formula), data = d, method = "fe",
 index=FE, tol = 1e-3, se = 1, nboots = nboots, r = 0, min.T0 = 4,
 loo = 0, pre.period = c(-6,0), group = 'cohort', carryover.rm = 3,
 CV = FALSE, force = "two-way", parallel = 1, cores = 4)

out2.ife <- fect(formula = as.formula(fect.formula), data = d, method = "ife",
  index=FE, tol = 1e-3, se = 1, nboots = nboots, min.T0 = 4,
  loo = 0, pre.period = c(-6,0), group = 'cohort', carryover.rm = 3,
  CV = FALSE, r =2, force = "two-way", parallel = 1, cores = 4)

out2.mc <- fect(as.formula(fect.formula), data = d, method = "mc",
 index=FE, tol = 1e-4, se = 1, nboots = nboots, min.T0 = 4, 
 loo = 0, pre.period = c(-6,0), group = 'cohort', carryover.rm = 3,
 CV = 1, force = "two-way", parallel = 1, cores = 4)

out2.mc$lambda.cv/out2.mc$eigen.all[1] # 0.0749

## Placebo tests

out2.fect.p <- fect(formula = as.formula(fect.formula), data = d, method = "fe",
 index=FE, tol = 1e-3, se = 1, nboots = nboots, r = 0, carryover.rm = 3, min.T0 = 1,
 CV = FALSE, force = "two-way", parallel = 1, cores = 4, placeboTest = 1, 
 placebo.period = c(-2, 0))

out2.ife.p <- fect(formula = as.formula(fect.formula), data = d, method = "ife",
  index=FE, tol = 1e-3, se = 1, nboots = nboots, r = 2, carryover.rm = 3, min.T0 = 1,
  CV = 0, force = "two-way", parallel = 1, cores = 4, 
  placeboTest = 1, placebo.period = c(-2, 0))

out2.mc.p <- fect(formula = as.formula(fect.formula), data = d, method = "mc",
 index=FE, tol = 1e-4, se = 1, nboots = nboots, lambda = out2.mc$lambda.cv, 
 CV = 0, force = "two-way", parallel = 1, cores = 4, carryover.rm = 3, min.T0 = 1,
 placeboTest = 1, placebo.period = c(-2, 0))


## Tests for carry-over effects

out2.fect.c <- fect(formula = as.formula(fect.formula), data = d, method = "fe",
 index=FE, tol = 1e-3, se = 1, nboots = nboots, force = "two-way",  carryover.rm = 3, min.T0 = 1,
 parallel = 1, cores = 4, carryoverTest = 1, carryover.period = c(1, 3))

out2.ife.c <- fect(formula = as.formula(fect.formula), data = d, method = "ife",
  index=FE, tol = 1e-3, se = 1, nboots = nboots, r = 2,  carryover.rm = 3, min.T0 = 1,
  CV = 0, force = "two-way", parallel = 1, cores = 4, 
  carryoverTest = TRUE, carryover.period = c(1,3))

out2.mc.c <- fect(formula = as.formula(fect.formula), data = d, method = "mc",
 index=FE, tol = 1e-4, se = 1, nboots = nboots, lambda = out2.mc$lambda.cv, 
 CV = 0, force = "two-way", parallel = 1, cores = 4,  carryover.rm = 3, min.T0 = 1,
 carryoverTest = 1, carryover.period = c(1, 3))


# Save
save(d, fect.formula, out.fe, 
  out.fect, out.fect.p, out.fect.c,
  out.ife, out.ife.p, out.ife.c,
  out.mc, out.mc.p, out.mc.c,
  out2.fect, out2.fect.p, out2.fect.c,
  out2.ife, out2.ife.p, out2.ife.c,
  out2.mc, out2.mc.p, out2.mc.c,
  file = "results/ex_FM2015.RData")

###############################
## Plotting
###############################

load("results/ex_FM2015.RData")

library(panelview)
library(ggplot2)
require(fect)
library(grid)
library(gridExtra)


## Outcome
unit <- "councilnumber"
period <- "year"
FE <- c(unit, period)
p.outcome <- panelview(as.formula(fect.formula), data = d, index = FE, theme.bw = TRUE,
          type = "outcome", main = "", xlab = "Year", ylab = "Specific Grants (log)", 
          legend.labs = c("w/o Partisan Alignment", "w/ Partisan Alignment"), by.group = FALSE)
ggsave('graph/ex_FM2015_outcome.pdf', p.outcome, width = 8, height = 5)


## Treatment
p.treat <- panelview(as.formula(fect.formula), data = d, index = FE, 
          pre.post = FALSE, by.timing = TRUE, xlab = "Year", ylab = "Local Council", 
          legend.labs = c("w/o Partisan Alignment", "w/ Partisan Alignment"),
          background= "white", axis.lab = "time",
          gridOff = TRUE) 
p.treat <- p.treat + theme(axis.text.x = element_text(angle = 90))
ggsave('graph/ex_FM2015_treat.pdf', p.treat, width = 8, height = 8)

######################################
# Do not consider carry over effects
######################################


## Gap plots
d_ylim <- c(-0.4, 0.4)
xlab <- "Year(s) Since Partisan Alignment"
p1.gap <- plot(out.fect, type = "gap", ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  ylim = d_ylim, proportion = 0.1, xlab = xlab,
  xlim = c(-6,4), main = "FEct", theme.bw = TRUE, cex.text = 0.8)

p2.gap <- plot(out.ife, type = "gap", ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  ylim = d_ylim, proportion = 0.1, xlab = xlab,
  xlim = c(-6,4), main = "IFEct", theme.bw = TRUE, cex.text = 0.8)

p3.gap <- plot(out.mc, type = "gap", ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  ylim = d_ylim, proportion = 0.1, xlab = xlab,
  xlim = c(-6,4), main = "MC", theme.bw = TRUE, cex.text = 0.8)
margin = theme(plot.margin = unit(c(1,1,1,1), "line"))
p.gap <- grid.arrange(grobs = lapply(list(p1.gap, p2.gap, p3.gap), "+", margin),
  ncol = 3, widths = c (1, 1, 1))
ggsave('graph/ex_FM2015_gap.pdf', p.gap, width = 17, height = 6)

## Placebo plots
stats.labs <- c("t test p-value","TOST p-value")
xlab <- "Year(s) Since Partisan Alignment"
p1.placebo <- plot(out.fect.p, ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  placeboTest = 1, proportion = 0.1, xlab = xlab,
  xlim = c(-6,4), ylim = d_ylim, main = "FEct", theme.bw = TRUE, 
  stats = c("placebo.p","equiv.p"), stats.labs = stats.labs) 

p2.placebo <- plot(out.ife.p, ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  placeboTest = 1, proportion = 0.1, xlab = xlab,
  xlim = c(-6,4), ylim = d_ylim, main = "IFEct", theme.bw = TRUE, 
  stats = c("placebo.p","equiv.p"), stats.labs = stats.labs) 

p3.placebo <- plot(out.mc.p, ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  placeboTest = 1, proportion = 0.1, xlab = xlab,
  xlim = c(-6,4), ylim = d_ylim, main = "MC", theme.bw = TRUE, 
  stats = c("placebo.p","equiv.p"), stats.labs = stats.labs) 

margin = theme(plot.margin = unit(c(1,1,1,1), "line"))
p.placebo <- grid.arrange(grobs = lapply(list(p1.placebo, p2.placebo, p3.placebo), "+", margin),
  ncol = 3, widths = c (1, 1, 1))
ggsave('graph/ex_FM2015_placebo.pdf', p.placebo, width = 17, height = 6)


## Equivalence plots
stats.labs <- c("F test p-value","TOST max p-value")
xlab <- "Year(s) Since Partisan Alignment"
p1.equiv <- plot(out.fect, type = "equiv", ylab = "Average Prediction Error", ylim = d_ylim, 
  xlim = c(-6,0), main = "FEct", theme.bw = TRUE, legendOff = TRUE, loo = 1, xlab = xlab,
  stats.labs = stats.labs) 

p2.equiv <- plot(out.ife, type = "equiv", ylab = "Average Prediction Error", ylim = d_ylim, 
  xlim = c(-6,0), main = "IFEct", theme.bw = TRUE, legendOff = TRUE, loo = 1, xlab = xlab,
  stats.labs = stats.labs) 

p3.equiv <- plot(out.mc, type = "equiv", ylab = "Average Prediction Error", ylim = d_ylim, 
  xlim = c(-6,0), main = "MC", theme.bw = TRUE, legendOff = TRUE, loo = 1, xlab = xlab,
  stats.labs = stats.labs) 

margin = theme(plot.margin = unit(c(1,1,1,1), "line"))
p.equiv <- grid.arrange(grobs = lapply(list(p1.equiv, p2.equiv, p3.equiv), "+", margin),
  ncol = 3, widths = c (1, 1, 1))
ggsave('graph/ex_FM2015_equiv.pdf', p.equiv, width = 17, height = 6)


## Carryover effects
c_ylim <- c(-0.4,0.4)
xlab <- "Year(s) Since Partisan Disalignment"
stats.labs <- c("  t test p-value","TOST p-value")
stats.pos <- c(2.8, 0.42)
p1 <- plot(out.fect.c, type = "exit", 
  ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  carryoverTest = 1, proportion = 0.1, xlab = xlab, stats.pos = stats.pos,
  xlim = c(-3,9), ylim = c_ylim, main = "FEct", theme.bw = TRUE, 
  stats = c("carryover.p","equiv.p"), stats.labs = stats.labs) 

p2 <- plot(out.ife.c, type = "exit",
  ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  carryoverTest = 1, proportion = 0.1, xlab = xlab, stats.pos = stats.pos,
  xlim = c(-3,9), ylim = c_ylim, main = "IFEct", theme.bw = TRUE, 
  stats = c("carryover.p","equiv.p"), stats.labs = stats.labs) 

p3 <- plot(out.mc.c, type = "exit",
  ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  carryoverTest = 1, proportion = 0.1, xlab = xlab, stats.pos = stats.pos,
  xlim = c(-3,9), ylim = c_ylim, main = "MC", theme.bw = TRUE, 
  stats = c("carryover.p","equiv.p"), stats.labs = stats.labs) 

margin = theme(plot.margin = unit(c(1,1,1,1), "line"))
p.carryover <- grid.arrange(grobs = lapply(list(p1, p2, p3), "+", margin),
  ncol = 3, widths = c (1, 1, 1))
ggsave('graph/ex_FM2015_carryover.pdf', p.carryover, width = 17, height = 6)


d_ylim <- c(-0.4, 0.4)
xlab <- "Year(s) Since Partisan Alignment"
ylab <- "The Effect of Partisan Alignment on Specific Grants (log)"


p0 <- plot(out.ife,main = "All Cohorts",ylim=c(-0.4,0.4),
  xlab = xlab, ylab = ylab)

p1 <- plot(out.ife,show.group = '1',main = "Cohort: First Treated in [1992,1996]",ylim=c(-0.4,0.4),
  xlab = xlab, ylab = ylab)

p2 <- plot(out.ife,show.group = '2',main = "Cohort: First Treated in [1997,2009]",ylim=c(-0.4,0.4),
  xlab = xlab, ylab = ylab)

p3 <- plot(out.ife,show.group = '3',main = "Cohort: First Treated after 2010",ylim=c(-0.4,0.4),
  xlab = xlab, ylab = ylab)

margin = theme(plot.margin = unit(c(1,1,1,1), "line"))
p <- grid.arrange(grobs = lapply(list(p1, p2, p3), "+", margin),
  ncol = 3, widths = c (1, 1, 1))
ggsave('graph/ex_FM2015_cohorts.pdf', p, width = 17, height = 6)

d_ylim <- c(-0.4, 0.4)
xlab <- "Year(s) Since Partisan Alignment"
ylab <- "The Effect of Partisan Alignment on Specific Grants (log)"
p0 <- plot(out.fect, main = "All Cohorts",ylim=c(-0.4,0.4),
  xlab = xlab, ylab = ylab)

p1 <- plot(out.fect,show.group = '1',main = "Cohort: First Treated in [1992,1996]",ylim=c(-0.4,0.4),
  xlab = xlab, ylab = ylab)

p2 <- plot(out.fect,show.group = '2',main = "Cohort: First Treated in [1997,2009]",ylim=c(-0.4,0.4),
  xlab = xlab, ylab = ylab)

p3 <- plot(out.fect,show.group = '3',main = "Cohort: First Treated after 2010",ylim=c(-0.4,0.4),
  xlab = xlab, ylab = ylab)

margin = theme(plot.margin = unit(c(1,1,1,1), "line"))
p <- grid.arrange(grobs = lapply(list(p1, p2, p3), "+", margin),
  ncol = 3, widths = c (1, 1, 1))
ggsave('graph/ex_FM2015_cohorts0.pdf', p, width = 17, height = 6)

##########################################################
# Analysis after removing three post-treatment periods
##########################################################

## Gap plots
d_ylim <- c(-0.4, 0.4)
xlab <- "Year(s) Since Partisan Alignment"
p1.gap <- plot(out2.fect, type = "gap", ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  ylim = d_ylim, proportion = 0.1, xlab = xlab,
  xlim = c(-6,4), main = "FEct", theme.bw = TRUE, cex.text = 0.8)

p2.gap <- plot(out2.ife, type = "gap", ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  ylim = d_ylim, proportion = 0.1, xlab = xlab,
  xlim = c(-6,4), main = "IFEct", theme.bw = TRUE, cex.text = 0.8)

p3.gap <- plot(out2.mc, type = "gap", ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  ylim = d_ylim, proportion = 0.1, xlab = xlab,
  xlim = c(-6,4), main = "MC", theme.bw = TRUE, cex.text = 0.8)
margin = theme(plot.margin = unit(c(1,1,1,1), "line"))
p.gap <- grid.arrange(grobs = lapply(list(p1.gap, p2.gap, p3.gap), "+", margin),
  ncol = 3, widths = c (1, 1, 1))
ggsave('graph/ex_FM2015b_gap.pdf', p.gap, width = 17, height = 6)

## Placebo plots
stats.labs <- c("t test p-value","TOST p-value")
xlab <- "Year(s) Since Partisan Alignment"
p1.placebo <- plot(out2.fect.p, ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  placeboTest = 1, proportion = 0.1, xlab = xlab,
  xlim = c(-6,4), ylim = d_ylim, main = "FEct", theme.bw = TRUE, 
  stats = c("placebo.p","equiv.p"), stats.labs = stats.labs) 

p2.placebo <- plot(out2.ife.p, ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  placeboTest = 1, proportion = 0.1, xlab = xlab,
  xlim = c(-6,4), ylim = d_ylim, main = "IFEct", theme.bw = TRUE, 
  stats = c("placebo.p","equiv.p"), stats.labs = stats.labs) 

p3.placebo <- plot(out2.mc.p, ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  placeboTest = 1, proportion = 0.1, xlab = xlab,
  xlim = c(-6,4), ylim = d_ylim, main = "MC", theme.bw = TRUE, 
  stats = c("placebo.p","equiv.p"), stats.labs = stats.labs) 

margin = theme(plot.margin = unit(c(1,1,1,1), "line"))
p.placebo <- grid.arrange(grobs = lapply(list(p1.placebo, p2.placebo, p3.placebo), "+", margin),
  ncol = 3, widths = c (1, 1, 1))
ggsave('graph/ex_FM2015b_placebo.pdf', p.placebo, width = 17, height = 6)


## Carryover effects
c_ylim <- c(-0.4,0.4)
xlab <- "Year(s) Since Partisan Disalignment"
stats.labs <- c("  t test p-value","TOST p-value")
stats.pos <- c(2.8, 0.42)
p1 <- plot(out2.fect.c, type = "exit", 
  ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  carryoverTest = 1, proportion = 0.1, xlab = xlab, stats.pos = stats.pos,
  xlim = c(-3,9), ylim = c_ylim, main = "FEct", theme.bw = TRUE, 
  stats = c("carryover.p","equiv.p"), stats.labs = stats.labs) 

p2 <- plot(out2.ife.c, type = "exit",
  ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  carryoverTest = 1, proportion = 0.1, xlab = xlab, stats.pos = stats.pos,
  xlim = c(-3,9), ylim = c_ylim, main = "IFEct", theme.bw = TRUE, 
  stats = c("carryover.p","equiv.p"), stats.labs = stats.labs) 

p3 <- plot(out2.mc.c, type = "exit",
  ylab = "The Effect of Partisan Alignment on Specific Grants (log)", 
  carryoverTest = 1, proportion = 0.1, xlab = xlab, stats.pos = stats.pos,
  xlim = c(-3,9), ylim = c_ylim, main = "MC", theme.bw = TRUE, 
  stats = c("carryover.p","equiv.p"), stats.labs = stats.labs) 

margin = theme(plot.margin = unit(c(1,1,1,1), "line"))
p.carryover <- grid.arrange(grobs = lapply(list(p1, p2, p3), "+", margin),
  ncol = 3, widths = c (1, 1, 1))
ggsave('graph/ex_FM2015b_carryover.pdf', p.carryover, width = 17, height = 6)

d_ylim <- c(-0.4, 0.4)
xlab <- "Year(s) Since Partisan Alignment"
ylab <- "The Effect of Partisan Alignment on Specific Grants (log)"


p0 <- plot(out2.ife,main = "All Cohorts",ylim=c(-0.4,0.4),
  xlab = xlab, ylab = ylab)

p1 <- plot(out2.ife,show.group = '1',main = "Cohort: First Treated in [1992,1996]",ylim=c(-0.4,0.4),
  xlab = xlab, ylab = ylab)

p2 <- plot(out2.ife,show.group = '2',main = "Cohort: First Treated in [1997,2009]",ylim=c(-0.4,0.4),
  xlab = xlab, ylab = ylab)

p3 <- plot(out2.ife,show.group = '3',main = "Cohort: First Treated after 2010",ylim=c(-0.4,0.4),
  xlab = xlab, ylab = ylab)

margin = theme(plot.margin = unit(c(1,1,1,1), "line"))
p <- grid.arrange(grobs = lapply(list(p1, p2, p3), "+", margin),
  ncol = 3, widths = c (1, 1, 1))
ggsave('graph/ex_FM2015b_cohorts.pdf', p, width = 17, height = 6)

## FEct (cohort effects)
d_ylim <- c(-0.4, 0.4)
xlab <- "Year(s) Since Partisan Alignment"
ylab <- "The Effect of Partisan Alignment on Specific Grants (log)"
p0 <- plot(out2.fect, main = "All Cohorts",ylim=c(-0.4,0.4),
  xlab = xlab, ylab = ylab)

p1 <- plot(out2.fect,show.group = '1',main = "Cohort: First Treated in [1992,1996]",ylim=c(-0.4,0.4),
  xlab = xlab, ylab = ylab)

p2 <- plot(out2.fect,show.group = '2',main = "Cohort: First Treated in [1997,2009]",ylim=c(-0.4,0.4),
  xlab = xlab, ylab = ylab)

p3 <- plot(out2.fect,show.group = '3',main = "Cohort: First Treated after 2010",ylim=c(-0.4,0.4),
  xlab = xlab, ylab = ylab)

margin = theme(plot.margin = unit(c(1,1,1,1), "line"))
p <- grid.arrange(grobs = lapply(list(p1, p2, p3), "+", margin),
  ncol = 3, widths = c (1, 1, 1))
ggsave('graph/ex_FM2015b_cohorts0.pdf', p, width = 17, height = 6)



print(Sys.time()-begin.time) 