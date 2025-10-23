#Clears previously loaded objects
rm(list=ls())

#Sets working directory
setwd("...")

#Required packages here
for(pack in c("gsynth","forecast","xtable","CADFtest","car","geometry"))
  if(!require(pack, character.only = T))
  {
    install.packages(pack)
    library(pack, character.only = T)
  }

#Auxiliary
source("../aux.R")

#Sets seed to allow for replication
set.seed(1234)

application = "basque"

outcome = "gdpcap"
label.outcome = "real GDP per capita"

time_field = "year"
label.time_field = "Year"
pre_treat_end = 1969

state_field = "regionname"
label.state_field = "Country"
remove_states = c("Syntetic Basque Country","Spain (Espana)")

base_state = read.csv(paste(application,"/data.csv",sep=""),stringsAsFactors = F)

base_state$outcome = base_state[,outcome]
base_state$time = base_state[,time_field]
base_state$state = base_state[,state_field]

base_state = base_state[order(base_state$state,base_state$time),]

sc_original =  base_state[base_state$state == "Syntetic Basque Country",]

base_state = base_state[!(base_state$state%in%remove_states),]


outcome.mat = matrix(base_state$outcome, nrow = length(unique(base_state$time)))
date.mat = matrix(base_state$time, nrow =length(unique(base_state$time)))

state_labels = unique(base_state$state)

j_treated = which(state_labels == "Basque Country (Pais Vasco)")

outcome.mat = cbind(outcome.mat[,j_treated],outcome.mat[,-j_treated])
pret=outcome.mat[date.mat[,1]<=pre_treat_end,]
postt= outcome.mat[date.mat[,1]>pre_treat_end,]

state_labels = state_labels[state_labels!="Basque Country (Pais Vasco)"]

#Conformal inference
conformal_fp(pret, postt)


#Plots

#PPlot 1.A
pdf("data.pdf", width = 10, height =8)
par(mar=c(10, 4, 4, 4))
plot(date.mat[,1],  outcome.mat[,1], type = "l",
     xlab = "Year", ylab = "Thousands of 1986 USD", ylim = c(min(outcome.mat), max(outcome.mat)))
grid(nx = NA, ny = NULL)
for(jjj in 2:ncol(outcome.mat))
  lines(date.mat[,1], outcome.mat[,jjj], col = "gray", lty = 2)
lines(date.mat[,1], outcome.mat[,1])
abline(v= pre_treat_end)
text(x = pre_treat_end-2, y =  max(outcome.mat)-2, labels = "Pre-treatment", pos = 2)
text(x = pre_treat_end+2, y =  max(outcome.mat)-2, labels = "Post-treatment", pos=4)
par(xpd = TRUE)
legend("bottom",inset = c(-0.3), legend = c("Treated unit","Controls"), lty =  c(1,2), ncol =2)
dev.off()


#Constructing counterfactuals for plot
did_eff = (outcome.mat[,1] - rowMeans(outcome.mat[,-1]))- (mean(pret[,1]) - mean(pret[,-1]))

sc_mode = synth_control_est(pret, outcome.mat)
sc_eff =sc_mode$effects


demean_pre = t(t(pret)-colMeans(pret))
demean_post = t(t(postt) - colMeans(pret))
demeaned_sc_mode = tryCatch({synth_control_est_demean(pret, outcome.mat)},
                            error = function(e){
                              print(e)
                              synth_control_est_allowrankdef(demean_pre, rbind(demean_pre,demean_post), add_up = T, nonnegative = T)
                            }
                )

#demeaned_sc_mode$effects  = rbind(demean_pre, demean_post)%*%c(1,-demeaned_sc_mode$solution.w)
demeaned_sc_eff = demeaned_sc_mode$effects

weight_values =  data.frame("SC (all lags)" = round(sc_mode$w,digits=4), "Demeaned SC (all lags)"= round(demeaned_sc_mode$weights,digits=4), "SC (Abadie et al. 2003)" = 0)
rownames(weight_values) = state_labels

weight_values["Cataluna",3] = 0.8508
weight_values["Madrid (Comunidad De)",3] = 0.1492

print(xtable(weight_values, digits=4, caption = "Estimated weights"), file = "weight_matrix.tex", caption.placement = "top")


did_counter  = outcome.mat[,1] - did_eff
sc_counter = outcome.mat[,1] - sc_eff
sc_demean_counter =   outcome.mat[,1] - demeaned_sc_eff

sc_original_counter = sc_original$gdpcap

#PPlot 1.B
pdf("sc_compare.pdf", width = 10, height =8)
par(mar=c(10, 4, 4, 4))
plot(date.mat[,1],  outcome.mat[,1], type = "l",
     xlab = "Year", ylab = "Thousands of 1986 USD", ylim = c(min(outcome.mat, sc_counter, sc_demean_counter, sc_original_counter),max(outcome.mat, sc_counter, sc_demean_counter, sc_original_counter)))
grid(nx = NA, ny = NULL)
lines(date.mat[,1], sc_counter, lty = 2)
lines(date.mat[,1], sc_demean_counter, lty = 3)
lines(date.mat[,1], sc_original_counter, lty = 4)
abline(v= pre_treat_end)
text(x = pre_treat_end-2, y =  max(outcome.mat)-2, labels = "Pre-treatment", pos = 2)
text(x = pre_treat_end+2, y =  max(outcome.mat)-2, labels = "Post-treatment", pos=4)
par(xpd = TRUE)
legend("bottom",inset = c(-0.3), legend = c("Treated unit","Demeaned SC (all lags)", "SC (all lags)", "SC (Abadie et al. 2003)"), lty =  c(1,3,2,4), ncol =2)
dev.off()


outcome.mat.detrend = outcome.mat - rowMeans(outcome.mat[,-1])
sc_counter.detrend = sc_counter -rowMeans(outcome.mat[,-1])
sc_demean_counter.detrend = sc_demean_counter -rowMeans(outcome.mat[,-1])
sc_original_counter.detrend = sc_original_counter -rowMeans(outcome.mat[,-1])


#PPlot 1.B
pdf("sc_compare.pdf", width = 10, height =8)
par(mar=c(10, 4, 4, 4))
plot(date.mat[,1],  outcome.mat[,1], type = "l",
     xlab = "Year", ylab = "Thousands of 1986 USD", ylim = c(min(outcome.mat[,1], sc_counter, sc_demean_counter, sc_original_counter),max(outcome.mat[,1], sc_counter, sc_demean_counter, sc_original_counter)))
grid(nx = NA, ny = NULL)
lines(date.mat[,1], sc_counter, lty = 2)
lines(date.mat[,1], sc_demean_counter, lty = 3)
lines(date.mat[,1], sc_original_counter, lty = 4)
abline(v= pre_treat_end)
text(x = pre_treat_end-2, y =  max(outcome.mat)-2, labels = "Pre-treatment", pos = 2)
text(x = pre_treat_end+2, y =  max(outcome.mat)-2, labels = "Post-treatment", pos=4)
par(xpd = TRUE)
legend("bottom",inset = c(-0.3), legend = c("Treated unit","Demeaned SC (all lags)", "SC (all lags)", "SC (Abadie et al. 2003)"), lty =  c(1,3,2,4), ncol =2)
dev.off()


#PPlot 1.C
pdf("sc_detrend.pdf", width = 10, height =8)
par(mar=c(10, 4, 4, 4))
plot(date.mat[,1],  outcome.mat.detrend[,1], type = "l",
     xlab = "Year", ylab = "Thousands of 1986 USD", ylim = c(0,4))
grid(nx = NA, ny = NULL)
lines(date.mat[,1], sc_counter.detrend, lty = 2)
lines(date.mat[,1], sc_demean_counter.detrend, lty = 3)
lines(date.mat[,1], sc_original_counter.detrend, lty = 4)
abline(v= pre_treat_end)
text(x = pre_treat_end-2, y =  max(outcome.mat.detrend)-1, labels = "Pre-treatment", pos = 2)
text(x = pre_treat_end+2, y =  max(outcome.mat.detrend)-1, labels = "Post-treatment", pos=4)
par(xpd = TRUE)
legend("bottom",inset = c(-0.3), legend = c("Treated unit","Demeaned SC (all lags)", "SC (all lags)", "SC (Abadie et al. 2003)"), lty =  c(1,3,2,4), ncol =2)
dev.off()




#Figure 1.D
pdf("spec_test.pdf", width = 10, height =8)
par(mar=c(10, 4, 4, 4))
plot(date.mat[,1],  rbind(pret,postt)[,1], type = "l",
     xlab = "Year", ylab = "Thousands of 1986 USD", ylim = c(min(did_counter, sc_demean_counter, rbind(pret,postt)[,1]), max(did_counter, sc_demean_counter, rbind(pret,postt)[,1])))
grid(nx = NA, ny = NULL)
lines(date.mat[,1], sc_demean_counter, lty = 2)
lines(date.mat[,1], did_counter, lty = 4)
abline(v= pre_treat_end)
text(x = pre_treat_end-2, y =  max(outcome.mat)-2, labels = "Pre-treatment", pos = 2)
text(x = pre_treat_end+2, y =  max(outcome.mat)-2, labels = "Post-treatment", pos=4)
par(xpd = TRUE)
legend("bottom",inset = c(-0.3), legend = c("Treated unit","Diff-in-Diff", "Demeaned SC (all lags)"), lty =  c(1,4,2), ncol =2)
dev.off()



