#-------------------------------------------------#

# Can synthetic controls improve causal inference in 
# interrupted time series evaluations of public health interventions? 

# First author & code developer: Michelle Degli Esposti
# Co-authors: Thees Sprecklesen, Antonio Gasparrini, Carl Bonander, Douglas J. Wiebe, Alexa R. Yakubovich, David K. Humphreys

# 15 October 2019

#-------------------------------------------------#

# Illustrative example for Figure 3:
    # DATA: Real-world data from WONDER, CDC
    # STEP 1. Replicating the original CITS analyses by Humphreys et al. (2017) using four comparison states (New York, New Jersey, Ohio, Virginia) 
    # STEP 2. Extending the original CITS analyses by Humphreys et al. (2017) using synthetic control methods (15 comparison states in the donor pool) 

#-------------------------------------------------#


#----------------------#
# 1. PREPARE ENVIRONMENT
#----------------------#

# 1.1 Clear environment
rm(list=ls())

# 1.2 Load packages
library(foreign); library(tsModel); library("lmtest"); library("Epi");library(dplyr)
library("splines"); library("vcd"); library(sandwich); library("reshape2"); library(tidyr); library(plyr)
library(Synth)


#----------------------#
# 2. LOAD DATA
#----------------------#

# Wide format (real-world data from WONDER, CDC)
flhom <- read.csv("./syg_dat1.csv") 
# Long format (real-world data from WONDER, CDC)
flcom <- read.csv("./syg_dat2.csv") 


#----------------------#
# 3. PREPARE DATA
#----------------------#

# 3.1 Outlier Sept, 2001 = 2753 homicides due to 9/11; remove from the analyses
flhom[33, 7] = NA
 
# 3.2 Derive homicide rates per 100,000 for both dataframes
flhom <- flhom %>% mutate(fl_rates = fl_hom/fl_stdpop*10^5,
                          c_rates = c_hom/c_stdpop*10^5)
flcom <- flcom %>% mutate(rates = hom/stdpop*10^5)



#----------------------#
# 4. SIMPLE ITS MODELS
#----------------------#
# Modelling: segmented log-linked Gaussian generalized linear model for monthly homicide rates, 1999-2014
# Results for these simple ITS models are not explicitly reported in the paper but are used for plotting Figure 3(a)

# 4.1 Florida (treated group)
flhom.m <- glm(fl_rates ~  Effective + time +
                 harmonic(month,2,12), family=gaussian(link = "log"), flhom)
# Summarise model & beta coefficients
summary(flhom.m)
confint(flhom.m, level = 0.95)
# Exponentiate to obtain RRs and 95% CIs
round(ci.lin(flhom.m,Exp=T),2)

# 4.2 Four comparison states (New York, New Jersey, Ohio, Virginia) 
chom.m <- glm(c_rates ~ Effective + time + 
                harmonic(month,2,12), family=gaussian(link = "log"),flhom)
# Summarise model & beta coefficients
summary(chom.m)
confint(chom.m, level = 0.95)
# Exponentiate to obtain RRs and 95% CIs
round(ci.lin(chom.m,Exp=T),2)



#----------------------#
# 5. CONTROLLED ITS MODELS: Using four comparison states
#----------------------#
# Modelling: segmented log-linked Gaussian generalized linear model for monthly homicide rates, 1999-2014
# Results for this replication of the original CITS are explicitly reported in the paper (see 14pg & Figure 3a in the manuscript)

# 5.1 Florida (treated group) vs Four comparison states (New York, New Jersey, Ohio, Virginia) 
int.m1 <- glm(rates ~ Effective*case + time*case
                      + harmonic(month,2,12),family=gaussian(link = "log"),flcom)
summary(int.m1)
confint(int.m1, level = 0.95)
# Exponentiate to obtain RRs and 95% CIs
round(ci.lin(int.m1,Exp=T),2)
# Interpreting Difference-in-difference estimation
    # The effect in the control states is under 'Effective'
    # The difference-in-difference effect is under 'Effective:case'
    # Therefore, to compute the results for Florida multiply 'Effective' x 'Effective:case'
    # OR relabel Florida=0 and Controls=1 in the 'case' column and re-run model to easily get estimate & corresponding CIs

# 5.2 Relabel Florida=0 and Controls=1 for easy interpretation of estimates
flcom <- flcom %>% mutate(case = -(case - 1))
int.m1.Florida <- glm(rates ~ Effective*case + time*case
              + harmonic(month,2,12), family=gaussian(link = "log"), flcom)
summary(int.m1.Florida)
confint(int.m1.Florida, level = 0.95)
# Exponentiate to obtain RRs and 95% CIs
round(ci.lin(int.m1.Florida,Exp=T),2)
# Interpreting Difference-in-difference estimation
    # The effect in Florida is now under 'Effective'

# 5.3 Sensitivity analysis: AR(1) models with Newey-West Standard Errors
# Get AR(1) covariance matrix using the sandwich package
ar1.vcov = NeweyWest(int.m1.Florida, lag=1)
# Find std.err for parameter on 'Effective'
se.effective = sqrt(ar1.vcov[2,2])
# Lower 95% CI
round(exp(coef(int.m1.Florida)[2]-se.effective*1.96), 2)
# Upper 95% CI
round(exp(coef(int.m1.Florida)[2]+se.effective*1.96), 2)
# = A bit wider CI, but almost identical



#----------------------#
# 6. PLOT FIGURE 3 (a) 
#----------------------#

# 6.1 Create data points
# 6.1.1 Florida
flhom.obs1 <- with(flhom, fl_hom/fl_stdpop*10^5)
flhom.datanew <- data.frame(fl_stdpop=mean(flhom$fl_stdpop),Effective=rep(c(0,1),c(819,1101)),
                            time= 1:1920/10,month=rep(1:120/10,16))
flhom.m.pred11 <- predict(flhom.m,type="response",flhom.datanew)
# 6.1.2 Comparison states 
chom.obs1 <- with(flhom, c_hom/c_stdpop*10^5)
chom.datanew <- data.frame(c_stdpop=mean(flhom$c_stdpop),Effective=rep(c(0,1),c(819,1101)),
                           time= 1:1920/10,month=rep(1:120/10,16))
chom.m.pred11 <- predict(chom.m,type="response",chom.datanew)

# 6.2 Plot & save Figure 3(a) in the manuscript
tiff(file="./Figure3a",width=4000,height=3000,res=600)
par (mar=c(4,4,1.5,1.5))
plot(1:192,flhom.obs1,type="n",ylim=c(0.0,1.0),xlab="Year",
     ylab="Rate per 100,000", main= NULL, cex.lab=.75, cex.axis=.65, frame.plot=F,xaxt="n",las=2) 
rect(82,0.0,192,1.0, col=grey(0.9),border=F) 
points(1:192,flhom.obs1,cex=0.7, pch=16, col="blue3") 
points(1:192,chom.obs1 ,cex=0.7, pch=16, col="darkgoldenrod1")
axis(1,at=0:16*12,labels=F, cex.lab=.75, cex.axis=.75) 
axis(1,at=0:15*12+6,tick=F,labels=1999:2014, cex.lab=.75, cex.axis=.65)
lines(1:1920/10,flhom.m.pred11,col="blue3", lwd=2)
lines(1:1920/10,chom.m.pred11,col="darkgoldenrod1", lwd=2)
legend ("topleft", legend=c("Florida", "Comparison states"), col=c("blue3", "darkgoldenrod1"), lty=1:1, cex=0.6, bty="n", inset=0.015)
dev.off()



# --------------------------------------------------------------------------------------------------
# 7. DERIVING A SYNTHETIC CONTROL (Synth Package)
# --------------------------------------------------------------------------------------------------

# 7.1 Load data for synthetic control analyses (includes 15 comparison state in donor pool) -- real-world data from WONDER, CDC
full_data <- read.csv("./syg_dat3.csv")
full_data <- subset(full_data, select = -c(X)) 


# 7.2.1 Convert variables to numeric
data <- subset(full_data, select = -c(State, Year.month)) # drop non-numeric variables
data <- data.frame(lapply(data, function(x) as.numeric(x))) # make all variables numeric
State <- full_data$State
data <- cbind(data, State) # re-add state names
data$State <- as.character(data$State)

# 7.2.2 Derive homicide rates
data$Population <- round(data$Population, digits = -2) # Round population to the nearest 100 to match Humphreys et al (2017) analysis
data$Homicide_count <- data$homicide_total
data$HomicideRates <- as.numeric((data$Homicide_count/data$Population)*100000)


# 7.2.3 Derive suicide rates as otherwise they lie outside the convex hull
data$Firearm.suicides <- as.numeric((data$Firearm.suicides/data$Population)*100000)
data$Suicides <- as.numeric((data$Suicides/data$Population)*100000)

# 7.2.4 Remove NAs
data <- data %>% filter(!is.na(State.Code))

# 7.3 Derive synthetic controls using: 
    # Outcome: homicide rates per 100,000 population
    # State-level characteristics (predictors & special predictors), see Supplementary Table 1
        # NB population size is not included as incoporated in outcome variable as they are entered as homicide rates per 100,000 population
    # 15 comparison states in the donor pool, see Supplementary Table 2
dataprep.out <- dataprep(foo = data, 
                         predictors = c("Unemployment_adj","Firearm.suicides", "Suicides"),
                         predictors.op = "mean", 
                         time.predictors.prior = c(1:81),
                         special.predictors = list(
                           # Convention, separate: list(VARNAME, PERIOD, AVERAGING METHOD)
                           # Yearly data
                           list("Paid.Hunting.License.Holders", seq(1,81, 12),"mean"),
                           list("Annual.Burglary.Rate", seq(1,81, 12),"mean"),
                           list("Personal.income.per.capita..dollars.", seq(1,81, 12),"mean"),
                           list("Annual.robbery.rate", seq(1,81, 12),"mean"),
                           list("Annual.incarceration.rate", seq(1,81, 12),"mean"),
                           list("Num_pop_over15", seq(1,81, 12),"mean"),
                           list("Proportion.of.population.hispanic", seq(1,81, 12),"mean"),
                           list("Proportion.of.population.AA.or.B", seq(1,81, 12),"mean") ,
                           list("Percentage.4.year.bachelors.or.above..25.years.old.and.over", seq(1,81, 12),"mean"),
                           list("Proportion.of.15.24", seq(1,81, 12),"mean"),
                           list("Gallons.of.ethanol", seq(1,81, 12),"mean"),
                           list("Num_pov", seq(1,81, 12),"mean"),
                           # Less than yearly e.g. every 4 years
                           list("Number.of.sworn.officers.per.1.000.U.S..residents", seq(13,81, 48), "mean"), # every 4 years: 2000, 2004
                           list("Percentage.of.republican.voters.in.presidential.election", seq(13,81, 48), "mean"), # every 4 years: 2000, 2004
                           list("Density", seq(13,81, 120), "mean"), # every decade starting 2000
                           list("MSA", seq(13,81, 120), "mean"), # every decade starting 2000
                           list("prop_urban", seq(13,81, 120), "mean") # every decade starting 2000
                         ),
                         dependent = "HomicideRates", 
                         unit.variable = "State.Code", 
                         unit.names.variable = "State",
                         time.variable = "time", 
                         treatment.identifier = 12, # Identifies Florida in State
                         controls.identifier = c(5,9,10,15,19,23,24,25,31,34,36,38,39,44,56), # Identifies comparison states in the donor pool
                         time.optimize.ssr = c(1:81), 
                         time.plot = c(1:192))
# Summarise
synth.out <- synth(dataprep.out)
# Unit weights
round(synth.out$solution.w,2)
# Predictor weights
synth.out$solution.v
# Save treatment and synthetic control for subsequent data analysis & figures
path.case <- dataprep.out$Y1plot
path.synth <- dataprep.out$Y0plot %*% synth.out$solution.w


# 7.4 Extract synthetic control and create a new dataframe to run the ITS models
# Filter dataframe for easier spreading to wide format
data_sub <- data %>%
  dplyr::select(State.Code, Year, Population, Homicide_count, Month.code, time, Case, State)
# Rearrange column order for easier spreading to wide format
data_sub_hom <- data_sub[c("time", "Year","Month.code","State", "Homicide_count")] # homicides
data_sub_hom_wide <- spread(data_sub_hom, State, Homicide_count) # homicides
data_sub_pop <- data_sub[c("time", "State", "Population")] # population
data_sub_pop_wide <- spread(data_sub_pop, State, Population) # population
# Merge where .x represents homicides and .y represents pop
data_sub_wide <- merge(data_sub_hom_wide, data_sub_pop_wide, by="time") 
# Add in effective to denote the enactment of SYG 0 <= time=81; 1>= time=82
data_sub_wide <- data_sub_wide %>% 
  mutate(Effective = ifelse(time <= 81, 0, 1))
# Add results from the synthetic control analyses
path.case <- as.data.frame(path.case)
names(path.case)[1] <- "path.case"
path.synth <- as.data.frame(path.synth)
names(path.synth)[1] <- "path.synth"
# Merge
data_wide <- cbind(data_sub_wide, path.case, path.synth)
# Prepare new dataframe in the long format (ie Florida & Synthetic control appended)
flcom <- data_wide %>% 
  dplyr::select(time, Year, Month.code, Effective, path.case, path.synth)
# Spilt to gather & linking variable & then merge [NB messy code]
flcom <- gather(flcom, case, Homicides, c(path.case, path.synth), factor_key=TRUE) # reshape data frame
# Tidy dataframe: 'case' variable indicates "1" = Florida & "0" = Comparison states
flcom <- flcom %>% 
  mutate(case = ifelse(case == "path.case", 1, 0))



#----------------------#
# 8. CONTROLLED ITS MODELS: Using a synthetic control
#----------------------#
# Modelling: segmented log-linked Gaussian generalized linear model for monthly homicide rates, 1999-2014
# Results are for the extended CITS models using a synthetic control are explicitly reported in the paper (see 14-15pg & Figure 3b in the manuscript)


# 8.1 Simple ITS models -- used to plot Figure 3b
# 8.1.1 Florida (treated group)
flhom.m <- glm(path.case ~  Effective + time +
                 harmonic(Month.code,2,12), family=gaussian(link = "log"), data_wide)
# 8.1.2 Synthetic control 
chom.m <- glm(path.synth ~ Effective + time + 
                harmonic(Month.code,2,12), family=gaussian(link = "log"),data_wide)

# 8.2 Controlled ITS models
# 8.2.1 Florida (treated group) vs synthetic control
int.m1 <- glm(Homicides ~ Effective*case + time*case
              + harmonic(Month.code,2,12), family=gaussian(link = "log"), flcom)
summary(int.m1)
confint(int.m1, level = 0.95)
# Exponentiate to obtain RRs and 95% CIs
round(ci.lin(int.m1,Exp=T),2)
# Interpreting Difference-in-difference estimation
      # The effect in the synthetic control is under 'Effective'
      # The difference-in-difference effect is under 'Effective:case'
      # Therefore, to compute the results for Florida multiply 'Effective' x 'Effective:case'
      # OR relabel Florida=0 and Syntheic control=1 in the 'case' column and re-run model to easily get estimate & corresponding CIs

# 8.2.2 Relabel Florida=0 and Controls=1 for easy interpretation of estimates
flcom <- flcom %>% mutate(case = -(case - 1))
int.m1.Florida <- glm(Homicides ~ Effective*case + time*case
                      + harmonic(Month.code,2,12), family=gaussian(link = "log"), flcom)
summary(int.m1.Florida)
confint(int.m1.Florida, level = 0.95)
# Exponentiate to obtain RRs and 95% CIs
round(ci.lin(int.m1.Florida,Exp=T),2)


# 8.2.3 Sensitivity analysis: AR(1) models with Newey-West Standard Errors
      # Running two models (one for Florida, one for synthetic) and combining the estimates.
flcom.old <- flcom # Keep old flcom
# 8.2.3.1 Florida (treated group) vs synthetic control
int.m1.intharm <- glm(Homicides ~ Effective*case + time*case
                      + case*harmonic(Month.code,2,12), family=gaussian(link = "log"), flcom.old)
# Exponentiate to obtain RRs and 95% CIs
round(ci.lin(int.m1.intharm,Exp=T),2)
# Get AR(1) covariance matrix for Florida using models from wide form data
ar1.vcov.flor = NeweyWest(flhom.m, lag=1)
# Find std.err for parameter on 'Effective'
se.florida = sqrt(ar1.vcov.flor[2,2])
# Get AR(1) covariance matrix for synthetic control from wide form data
ar1.vcov.scm = NeweyWest(chom.m, lag=1)
# Find std.err for parameter on 'Effective'
se.scm = sqrt(ar1.vcov.scm[2,2])
# Combine estimates (this estimate is identical to estimate from int.m1.intharm,which also allows seasonality to differ between the groups.)
est.comb = coef(flhom.m)[2]-coef(chom.m)[2]
se.comb = sqrt(se.florida^2+se.scm^2) # Assume no covariance
# Lower 95% CI
round(exp(est.comb-se.comb*1.96), 2)
# Upper 95% CI
round(exp(est.comb+se.comb*1.96), 2)
# = A bit wider CI, but again almost identical



#----------------------#
# 9. PLOT FIGURE 3 (b) 
#----------------------#

# 9.1 Create data points
# 9.1.1 Florida
flhom.obs <- with(data_wide, path.case)
flhom.datanew <- data.frame(Effective=rep(c(0,1),c(819,1101)),
                            time= 1:1920/10,Month.code=rep(1:120/10,16))
flhom.m.pred1 <- predict(flhom.m,type="response",flhom.datanew)
# 9.1.2 Synthetic control
chom.obs <- with(data_wide, path.synth)
chom.datanew <- data.frame(Effective=rep(c(0,1),c(819,1101)),
                           time= 1:1920/10,Month.code=rep(1:120/10,16))
chom.m.pred1 <- predict(chom.m,type="response",chom.datanew)

# 9.2 Plot & save Figure 3(b) in the manuscript
tiff(file="./Figure3b",width=4000,height=3000,res=600)
par (mar=c(4,4,1.5,1.5))
plot(1:192,flhom.obs,type="n",ylim=c(0.0,1.0),xlab="Year",
     ylab="Rate per 100,000", main= NULL, cex.lab=.75, cex.axis=.65, frame.plot=F,xaxt="n",las=2) 
rect(82,0.0,192,1.0, col=grey(0.9),border=F) 
points(1:192,flhom.obs,cex=0.7, pch=16, col="blue3") 
points(1:192,chom.obs ,cex=0.7, pch=16, col="darkgoldenrod1")
axis(1,at=0:16*12,labels=F, cex.lab=.75, cex.axis=.75) 
axis(1,at=0:15*12+6,tick=F,labels=1999:2014, cex.lab=.75, cex.axis=.65)
lines(1:1920/10,chom.m.pred1,col="darkgoldenrod1", lwd=2)
lines(1:1920/10,flhom.m.pred1,col="blue3", lwd=2)
legend ("topleft", legend=c("Florida", "Synthetic control"), col=c("blue3", "darkgoldenrod1"), lty=1:1, cex=0.6, bty="n", inset=0.015)
dev.off()


#----------------------#
# 10. PLOT FIGURE 3 by combining panels (a) and (b)
#----------------------#

# Plot & save Figure 3 in the manuscript
tiff(file="./Figure3",width=4000,height=6000,res=600)
par (mfrow=c(2,1))
plot(1:192,flhom.obs1,type="n",ylim=c(0.0,1.0),xlab="Year",
     ylab="Rate per 100,000", main= NULL, cex.lab=.75, cex.axis=.65, frame.plot=F,xaxt="n",las=2) 
title(main= "(a)", adj = 0)
rect(82,0.0,192,1.0, col=grey(0.9),border=F) 
points(1:192,flhom.obs1,cex=0.7, pch=16, col="blue3") 
points(1:192,chom.obs1 ,cex=0.7, pch=16, col="darkgoldenrod1")
axis(1,at=0:16*12,labels=F, cex.lab=.75, cex.axis=.75) 
axis(1,at=0:15*12+6,tick=F,labels=1999:2014, cex.lab=.75, cex.axis=.65)
lines(1:1920/10,flhom.m.pred11,col="blue3", lwd=2)
lines(1:1920/10,chom.m.pred11,col="darkgoldenrod1", lwd=2)
legend ("topleft", legend=c("Florida", "Comparison states"), col=c("blue3", "darkgoldenrod1"), lty=1:1, cex=0.6, bty="n", inset=0.015)
plot(1:192,flhom.obs,type="n",ylim=c(0.0,1.0),xlab="Year",
     ylab="Rate per 100,000",  cex.lab=.75, cex.axis=.65, frame.plot=F,xaxt="n",las=2) 
title(main= "(b)", adj = 0)
rect(82,0.0,192,1.0, col=grey(0.9),border=F) 
points(1:192,flhom.obs,cex=0.7, pch=16, col="blue3") 
points(1:192,chom.obs ,cex=0.7, pch=16, col="darkgoldenrod1")
axis(1,at=0:16*12,labels=F, cex.lab=.75, cex.axis=.75) 
axis(1,at=0:15*12+6,tick=F,labels=1999:2014, cex.lab=.75, cex.axis=.65)
lines(1:1920/10,chom.m.pred1,col="darkgoldenrod1", lwd=2)
lines(1:1920/10,flhom.m.pred1,col="blue3", lwd=2)
legend ("topleft", legend=c("Florida", "Synthetic control"), col=c("blue3", "darkgoldenrod1"), lty=1:1, cex=0.6, bty="n", inset=0.015)
dev.off()
