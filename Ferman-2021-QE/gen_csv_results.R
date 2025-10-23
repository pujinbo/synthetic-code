#Clears previously loaded objects
rm(list=ls())

#Sets working directory
setwd("~/Dropbox/CPS_Data")

library("xtable")
library("parallel")

application = "state_gdp"

outcome = "employment_per"
label.outcome = "Employment rate (%)"
add.nonstationary = F

load(paste("fixed_env_",outcome,"_added_nonstationary_",add.nonstationary,".Rdata",sep=""))


pdf(paste("res_plot_",outcome, ".pdf",sep=""),width = 20, height=7)
par(mfrow = c(5,11))
par(mar = c(2,2,2,1))
for(j in 1:ncol(residuals))
  plot(as.Date(date.mat[,1],origin = "1970-01-01"),residuals[,j], type = "l",col = "grey", lty = 1,  lwd = 0.5, ylab = "", xlab = "")
dev.off()

table.uroot = c()
for(rr in 1:ncol(residuals))
{
 
  model.small = CADFtest(residuals[,rr],max.lag.y = floor(12*((Tt)/100)^(1/4)), type = "none", criterion = "MAIC")
  table.uroot = rbind(table.uroot, "p-value" = model.small$p.value)
}
state_labels = unique(cps_state$state)

table.uroot = cbind("State FIPS" =  state_labels, table.uroot)

write.csv(table.uroot, paste(outcome,"added_nonstationary_",add.nonstationary, "_results_nonstationary_test.csv",sep=""))


state =   readRDS(paste("results_per_state/fixed_",outcome,"_state_",1,"_added_nonstationary_",add.nonstationary,".RDS",sep=""))

J.list = unique(sapply(state, function(x){x$J}))
T0.list = unique(sapply(state, function(x){x$T0}))
T1.list = unique(sapply(state, function(x){x$T1}))

mat_resultados = c()
for(j in 1:Nn)
{
  print(j)
  full =  readRDS(paste("results_per_state/fixed_",outcome,"_state_",j,"_added_nonstationary_",add.nonstationary,".RDS",sep=""))
  
  for(J in J.list)
  {
    for(T0 in T0.list)
    {
      for(T1 in T1.list)
      {
          state.keep = sapply(full, function(x){(x$J==J)&(x$T0==T0)&(x$T1==T1)})
          state = full[state.keep]
          
          eff.mat.sc = sapply(state, function(x){
            if(length(x$sc)>0)
              sapply(x$sc, function(y){mean(y$effect)}) else rep(NA, ncol(Factors)+1)
          })
          
          
          eff.mat.did = sapply(state, function(x){
            if(length(x$did)>0)
              sapply(x$did, function(y){mean(y$effect)}) else rep(NA, ncol(Factors)+1)
          })
          
          eff.mat.sc_demean = sapply(state, function(x){
            if(length(x$sc_demean)>0)
              sapply(x$sc_demean, function(y){mean(y$effect)}) else rep(NA, ncol(Factors)+1)
          })
          
          
          rec.mat.sc = rowMeans(sapply(state, function(x){
            if(length(x$sc)>0)
              rowMeans(sapply(x$sc, function(y){c(as.vector(x$FE.case[-1]%*%y$w),as.vector(t(x$Mu.case[-1,])%*%y$w))})) else rep(NA, ncol(Factors)+1)
          }),na.rm = T)
          
          rec.mat.sc_demean = rowMeans(sapply(state, function(x){
            if(length(x$sc_demean)>0)
              rowMeans(sapply(x$sc_demean, function(y){c(as.vector(t(x$Mu.case[-1,])%*%y$w))})) else rep(NA, ncol(Factors))
          }),na.rm=T)
          
          rec.mat.did = rowMeans(sapply(state, function(x){
            if(length(x$did)>0)
              rowMeans(sapply(x$did, function(y){c(as.vector(t(x$Mu.case[-1,])%*%y$weights))})) else rep(NA, ncol(Factors))
          }),na.rm=T)
          
          
          weight.stats_sc = t(rowMeans(sapply(state, function(x){
            if(length(x$sc)>0)
              rowMeans(sapply(x$sc, function(y){c(max(y$w),sum(y$w^2))})) else rep(NA, 2)
          }),na.rm = T))
          
          weight.stats_sc_demean = t(rowMeans(sapply(state, function(x){
            if(length(x$sc_demean)>0)
              rowMeans(sapply(x$sc_demean, function(y){c(max(y$w),sum(y$w^2))})) else rep(NA, 2)
          }),na.rm = T))
          
          colnames(weight.stats_sc) = paste("sc_",c("max_w", "sum_w_sqr"),sep="")    
          colnames(weight.stats_sc_demean) = paste("sc_demean_",c("max_w", "sum_w_sqr"),sep="")          
          
         
          names_effs = paste("eff_", c("no_brk", paste("brk_fac",1:(nrow(eff.mat.sc)-1),sep="")),sep="")
          
          factors = t(c(state[[1]]$FE.case[1], state[[1]]$Mu.case[1,]))
          colnames(factors) = c("fe", paste("fac",1:(length(factors)-1),sep= ""))
          
          factors_sc = t(rec.mat.sc)
          colnames(factors_sc) =  paste("rec_sc_",c("fe", paste("fac",1:(length(factors)-1),sep= "")),sep="")
          
          factors_sc_demean = t(rec.mat.sc_demean)
          colnames(factors_sc_demean) =  paste("rec_sc_demean_",c( paste("fac",1:(length(factors)-1),sep= "")),sep="")
          
          factors_did = t(rec.mat.did)
          colnames(factors_did) =  paste("rec_did_",c( paste("fac",1:(length(factors)-1),sep= "")),sep="")
          
          mean_sc = t(rowMeans(eff.mat.sc, na.rm=T))
          sd_sc = t(apply(eff.mat.sc, 1, sd,na.rm=T))
          colnames(mean_sc) = paste("sc_mean_",names_effs,sep="")
          colnames(sd_sc) = paste("sc_sd_",names_effs,sep="")
          
          mean_did = t(rowMeans(eff.mat.did, na.rm = T))
          sd_did = t(apply(eff.mat.did, 1, sd, na.rm = T))
          colnames(mean_did) = paste("did_mean_",names_effs,sep="")
          colnames(sd_did) = paste("did_sd_",names_effs,sep="")
          
          
          mean_sc_demean = t(rowMeans(eff.mat.sc_demean,na.rm=T))
          sd_sc_demean = t(apply(eff.mat.sc_demean, 1, sd, na.rm=T))
          colnames(mean_sc_demean) = paste("sc_demean_mean_",names_effs,sep="")
          colnames(sd_sc_demean) = paste("sc_demean_sd_",names_effs,sep="")
          
          
          share_fail_sc = rowMeans(is.na(eff.mat.sc))[1]
          share_fail_sc_demean = rowMeans(is.na(eff.mat.sc_demean))[1]
          
          
          #Spec tests
          test.q2 = sapply(state, function(x){
            sapply(x$break_test, function(z){z[1]})
          })
          test.q1 = sapply(state, function(x){
            sapply(x$break_test, function(z){z[2]})
          })
          
          test.q2.05 = t(rowMeans(test.q2<0.05))
          test.q1.05 = t(rowMeans(test.q1<0.05))
          
          test.q2.10 = t(rowMeans(test.q2<0.10))
          test.q1.10 = t(rowMeans(test.q1<0.10))
          
          colnames(test.q2.05) = paste("rej_05_q2_", c("no_brk", paste("brk_fac",1:(nrow(eff.mat.sc)-1),sep="")),sep="")
          colnames(test.q2.10) = paste("rej_10_q2_", c("no_brk", paste("brk_fac",1:(nrow(eff.mat.sc)-1),sep="")),sep="")
          colnames(test.q1.05) = paste("rej_05_q1_", c("no_brk", paste("brk_fac",1:(nrow(eff.mat.sc)-1),sep="")),sep="")
          colnames(test.q1.10) = paste("rej_10_q1_", c("no_brk", paste("brk_fac",1:(nrow(eff.mat.sc)-1),sep="")),sep="")
          
          
          
          line_resultados =  data.frame("state_no" = state[[1]]$j, "state_label" = state[[1]]$state_FIPS,
                                  "T0" = T0, "T1" = T1, "J" = J, 
                                  "in_hull"= 1*state[[1]]$in.hull,
                                  "in_hull_incl_FE" = 1*state[[1]]$in.hull.incl.FE,
                                  "share_fail_sc" = share_fail_sc,
                                  "share_fail_sc_demean" = share_fail_sc_demean, stringsAsFactors = F)
          
          line_resultados = cbind(line_resultados, factors, factors_sc, factors_sc_demean, factors_did, weight.stats_sc,
                                  weight.stats_sc_demean, mean_sc, sd_sc, mean_did, sd_did, mean_sc_demean,
                                  sd_sc_demean,test.q2.05,test.q2.10, test.q1.05,test.q1.10)
          
          mat_resultados = rbind(mat_resultados, line_resultados)
          
        }
        
        
      }
    }
}
    
write.csv(mat_resultados, paste(outcome,"added_nonstationary_",add.nonstationary, "_results.csv",sep=""))
    
    
    
    