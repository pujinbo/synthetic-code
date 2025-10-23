#Clears previously loaded objects
rm(list=ls())

#Sets working directory
setwd("~/Dropbox/CPS_Data")

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

outcome = "employment_per"
label.outcome = "Employment rate (%)"

#outcome = "lrwage"
#label.outcome = "Log real weekly earnings (1982-1984 avg. prices)"

cps_state = read.csv("Data/CPS/base_state_final.csv",stringsAsFactors = F)

cps_state = cps_state[cps_state$year>=1982&cps_state$year<=2019,]

cps_state$outcome = cps_state[,outcome]
cps_state$year_month = paste(cps_state$year,sapply(cps_state$month,function(x){if(x<10)return(paste("0",x,sep="")) else return(x)}),sep="-")

cps_state = cps_state[order(cps_state$state,cps_state$year_month),]

#Part 1. Plot
outcome.mat = matrix(cps_state$outcome, nrow = length(unique(cps_state$year_month)))
date.mat = matrix(as.Date(cps_state$year_month_day), nrow =length(unique(cps_state$year_month)))

pdf(paste("plot_",outcome, ".pdf",sep=""),width = 14, height=7)
plot(as.Date(date.mat[,1],origin = "1970-01-01"),rowMeans(outcome.mat), type = "l", ylab = label.outcome, xlab = "Month/Year", ylim = c((min(outcome.mat)),(max(outcome.mat))))
for(j in 1:ncol(outcome.mat))
  lines(as.Date(date.mat[,1],origin = "1970-01-01"), outcome.mat[,j], col = "lightgrey", lty = 2,  lwd = 0.5)
lines(as.Date(date.mat[,1],origin = "1970-01-01"),rowMeans(outcome.mat))
dev.off()

#Part 2. Estimating factor model
unit.fe = as.matrix(model.matrix(~-1+as.factor(state), data = cps_state))
unit.fe = unit.fe[,-1]
colnames(unit.fe) = paste("zzz",1:ncol(unit.fe),sep="")

cps_state = cbind(cps_state, unit.fe)

formula = paste("outcome ~",paste(colnames(unit.fe),collapse = " + "))


Nn = length(unique(cps_state$state))
Tt = length(unique(paste(cps_state$year,cps_state$month)))

run.fac.model <- function(r){
  factor.model = interFE(as.formula(formula),data = cps_state, r= r, se = F, index = c("state","year_month"), force = "time", normalize = F)
  
  icp1 = log(factor.model$sigma2) + r*((Nn+Tt)/(Nn*Tt))*log((Nn*Tt)/(Nn+Tt))
  
  return(list("icp1"=icp1, "model"= factor.model))
}

models = lapply(0:10, run.fac.model)
ind.min = which.min(sapply(models, function(x){x$icp1}))

factor.model =  models[[ind.min]]$model

state.effects = factor.model$beta 
state.effects = c(0, state.effects)

Mu = factor.model$lambda
Factors = factor.model$factor

time.effects = factor.model$xi

time.effects = time.effects+factor.model$mu
residuals = factor.model$residuals

pdf(paste("factor_plot_",outcome, ".pdf",sep=""),width = 14, height=7)
plot(as.Date(date.mat[,1],origin = "1970-01-01"),Factors[,1], type = "l", ylab = "Factors", xlab = "Month/Year", ylim = c((min(Factors)),(max(Factors))))
for(j in 1:ncol(Factors))
  lines(as.Date(date.mat[,1],origin = "1970-01-01"), Factors[,j],  lty = j)
legend("topright", legend= paste("Factor",1:ncol(Factors)), lty = 1:ncol(Factors), bty = "n")
dev.off()


arma.list = list()
#Fitting best ARMA model
for(j in 1:ncol(residuals))
{
  fit <- auto.arima(residuals[,j], seasonal = F,max.p = floor(12*(Tt/100)^(1/4)), max.q = floor(12*(Tt/100)^(1/4)), max.order = floor(12*(Tt/100)^(1/4)), stationary = T, stepwise = F, ic = "bic")
  print(fit)
  arma.list[[j]] = fit
}

#Max T for simulations /
Lim = 1300 

arma.sim <- function(fit){
  if(length(fit$coef)==0)
    return(rnorm(1000+Lim, sd = sqrt(fit$sigma2))[(1001):(1000+Lim)]) else {
      
      arma.spec = list(order = c(fit$arma[1],0,fit$arma[2]))
      if(fit$arma[1]>0)
        arma.spec$ar = fit$coef[grepl("ar", names(fit$coef))]
      
      if(fit$arma[2]>0)
        arma.spec$ma = fit$coef[grepl("ma", names(fit$coef))]
      
    }
  
  if(is.na(fit$coef["intercept"]))
    mean = 0 else
      if(fit$arma[1]>0)
        mean = (1- sum(fit$coef[grepl("ar", names(fit$coef))]))*fit$coef["intercept"] else
          mean = fit$coef["intercept"]
        
  return(arima.sim(n=1000+Lim, model = arma.spec, rand.gen = function(x){rnorm(x, mean=mean, sd = sqrt(fit$sigma2))})[(1001):(1000+Lim)])
  
}

#Testing for factors now

table.uroot = c()
is.stationary = c()
for(ff in 1:ncol(Factors))
{
  model = CADFtest(Factors[,ff],max.lag.y = floor(12*((Tt)/100)^(1/4)), type = "drift", criterion = "MAIC")
  tstat.drift = summary(model$est.model)$coefficients[1,3]
  
  fstat = linearHypothesis(model$est.model, c("(Intercept) = 0", "L(y, 1) = 0"))$F[2]
  
  model.small = CADFtest(Factors[,ff],max.lag.y = floor(12*((Tt)/100)^(1/4)), type = "none", criterion = "MAIC")
  table.uroot = rbind(table.uroot, cbind("P-value (model w/ drift)"=model$p.value,"P-value (model w/o drift)" = model.small$p.value,   "T-stat drift" = tstat.drift, "F-stat drift" =fstat))

  is.stationary =c(is.stationary, ifelse(abs(tstat.drift)<qnorm(0.95), model.small$p.value < 0.1, model$p.value < 0.1))
}

rownames(table.uroot) = paste("Factor",1:ncol(Factors))

print(xtable(table.uroot,caption = paste(label.outcome, "-- Factor unit root tests"), align = rep("c",5)),
      file = paste("uroot_",outcome, ".tex",sep=""),
      include.rownames = T, caption.placement = "top", table.placement = "H")


arma.coefs = lapply(1:ncol(Factors), function(ff){
  if(is.stationary[ff])
  {
    fit <- auto.arima(Factors[,ff], seasonal = F,max.p = floor(12*(Tt/100)^(1/4)), max.q = floor(12*(Tt/100)^(1/4)), max.order = floor(12*(Tt/100)^(1/4)), stationary = T, stepwise = F, ic = "bic")
    print(fit)
    return(fit) 
  } else {
   return(NULL)
  }
})

simFactors <- function()
{
  sapply(1:ncol(Factors), function(ff){
    if(is.stationary[ff])
      return(arma.sim(arma.coefs[[ff]])) else {
        est = c(Factors[,ff])
        
        SS = ceiling(Lim/Tt)
        
        for(ss in 2:SS)
          est = c(est, Factors[,ff]+est[length(est)])
        
        est = est[1:Lim]
        
        return(est)
      }
  })
}

#Here we have the option to include non-linear heterogenous trends. In the main paper, we did not include them. 
add.nonstationary  = F
if(add.nonstationary) 
{
  Factors = cbind(Factors,1:nrow(Factors)) 
  
  vvv = runif(nrow(Mu), min = -1, max = 1)
  vvv[order(vvv)[1:2]] = mean(vvv[order(vvv)][1:2])
  vvv[order(vvv,decreasing = T)[1:2]] = mean(vvv[order(vvv,decreasing = T)][1:2])
  Mu = cbind(Mu,vvv)
  
  arma.coefs = c(arma.coefs,list(NULL))
  
  is.stationary = c(is.stationary, F)
}

Replications = 5000
simulations = list()
fac.list = list()
for(rr in 1:Replications)
{
  fac.rr = simFactors()
  fac.list[[rr]] = fac.rr
  simulations[[rr]] = fac.rr%*%t(Mu) + rep(1,Lim)%*%t(state.effects) + sapply(arma.list, arma.sim)
 
}

state_codes = unique(cps_state$state)

#Loading distance matrix
dist_mat = read.csv("Data/distance_centroid.csv")
rownames(dist_mat) = dist_mat[,1]

dist_mat = as.matrix(dist_mat[,-1])

dist_mat = dist_mat[rownames(dist_mat)%in%state_codes,rownames(dist_mat)%in%state_codes]

#Saves environment
save.image(paste("fixed_env_",outcome,"_added_nonstationary_",add.nonstationary,".Rdata",sep=""))

lista.resultados =  list()

J.vec = c(Nn-1)
T1.vec = c(12)
T0.vec = c(120,240,480,1200)

fips_labels = unique(cps_state$state)
for(j in 1:Nn)
{
  print(j)
  Mu.case = rbind(Mu[j,], Mu[-j,])
  
  state.fe.case = c(state.effects[j],state.effects[-j])
  
  for(J in J.vec){
    print(J)
  
  controls_fixed = order(dist_mat[j,-j])[1:J]
  hull_fixed = inhulln(convhulln(Mu.case[controls_fixed+1,]),t(Mu.case[1,]))
  hull_fixed_incl_fe = inhulln(convhulln(cbind(Mu.case[controls_fixed+1,],state.fe.case[controls_fixed+1])),t(c(Mu.case[1,],state.fe.case[1])))
  lab_controls_fixed = paste(controls_fixed, collapse = ";")
  
   for(rr in 1:Replications){
     for(Jmethod in c("neighbor"))
     {
       if(Jmethod=="random")
       {
         controls = sample(1:(Nn-1),J, replace = F)
         hull = inhulln(convhulln(Mu.case[controls+1,]),t(Mu.case[1,]))
         hull_incl_fe = inhulln(convhulln(cbind(Mu.case[controls+1,],state.fe.case[controls+1])),t(c(Mu.case[1,],state.fe.case[1])))
         lab_controls = paste(controls, collapse = ";")
       } else {
         controls = controls_fixed
         hull = hull_fixed
         hull_incl_fe = hull_fixed_incl_fe
         lab_controls = lab_controls_fixed
       }
       
        Ybase = simulations[[rr]]
        Y = cbind(Ybase[,j], Ybase[,-j])
        Y = Y[,c(1,1+controls)]
        
        for(T0 in T0.vec){
          for(T1 in T1.vec){
            sc_list = list()
            sc_demean_list = list()
            did_list = list()
            bai_list = list()
            
            break.test.list = list()
            
            
            for(structural.break in 0:ncol(Factors))
            {
             
              Y_before = Y[1:T0,]
              Y_after = Y[(T0+1):(T0+T1),]
              
              if(structural.break > 0)
                Y_after = t(t(Y_after) + 2*Mu.case[c(1,controls+1),structural.break]*sd(Factors[,structural.break]))
              
              if(structural.break==0)
              {
                sc = synth_control_est(Y_before,Y_after)
                sc_demean = synth_control_est_demean(Y_before, Y_after)
                
              }
              
              if(structural.break>0)
              {
                sc$effects =  Y_after%*%c(1,-sc$w)
                sc_demean$effects =  -sc_demean$intercept + Y_after%*%c(1,-sc_demean$w)
              }
              
              did = list("effect" = mean(Y_after[,1]) - mean(Y_before[,1]) - (mean(Y_after[,-1]) - mean(Y_before[,-1]) ),
                         "weights" = rep(1/(J), J))
              
              
              
              
              sc_list[[paste(structural.break)]] = sc
              sc_demean_list[[paste(structural.break)]]  = sc_demean
              did_list[[paste(structural.break)]] = did
              break.test.list[[paste(structural.break)]] =  conformal_fp(Y_before, Y_after)
               
              
              
                           
            }
            
            lista.resultados[[paste(j,rr,Jmethod,J,T0,T1,sep=".")]] = list("j"=j, "state_FIPS" = fips_labels[j], "J"= J, "Jmethod" = Jmethod, "in.hull" = hull, "in.hull.incl.FE" = hull_incl_fe,
                                                                                            "lab.controls" = lab_controls, "Mu.case" = Mu.case[c(1,controls+1),],
                                                                                            "FE.case" = state.fe.case[c(1,controls+1)],
                                                                                            "replication"=rr, "T0"=T0, "T1"=T1,
                                                                                            "sc" = sc_list, "sc_demean" = sc_demean_list, "did" = did_list, "bai" = bai_list,
                                                                                             "break_test" = break.test.list)
            
        }
      }
     }
   }
  }
  
  saveRDS(lista.resultados, paste("results_per_state/fixed_",outcome,"_state_",j,"_added_nonstationary_",add.nonstationary,".RDS",sep=""))
  
  rm(list="lista.resultados")
  
  lista.resultados = list()
  
}
