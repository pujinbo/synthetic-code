#Loads required package, installing if needed
for(pack in c("nnls","quadprog", "parallel"))
  if(!require(pack, character.only = T))
  {
    install.packages(pack)
    require(pack, character.only = T)
  }


#Function to simulate AR(1) for T0 periods. If is.na(start) == T, then initial value
#is drawn from stationary distribution. We discard the initial draw after that.
simulate_ar1 <- function(rho, var_shock, T0, intercept = 0, start = NA)
{
  y = ifelse(is.na(start), rnorm(1,mean=intercept/(1-rho),sd = sqrt(var_shock/(1-rho^2))), start)
  
  for(t in 1:T0)
    y = c(y,intercept + rho*y[t]+rnorm(1,mean=0,sd = sqrt(var_shock)))
  
  return(y[-1])
}

#Synthetic control estimator
synth_control_est <-function(y_before, y_after)
{
  y= y_before[,1]
  X =y_before[,-1]
  
  Dmat = t(X)%*%X
  dvec = t(X)%*%y
  
  Amat = t(rbind(rep(1,ncol(X)),diag(ncol(X))))
  bvec = c(1, rep(0,ncol(X)))
  
  synth_model = solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  
  w = synth_model$solution
  effects = y_after%*%c(1,-w)
  
  return(list("w" = w, "effects" = effects, "value_min" = (2*synth_model$value +y%*%y) ))
  
}

synth_control_est_demean <-function(y_before, y_after)
{
  y= y_before[,1]
  X =y_before[,-1]
  
  X = cbind(1,X)
  
  Dmat = t(X)%*%X
  dvec = t(X)%*%y
  
  Amat = t(rbind(c(0,rep(1,ncol(X)-1)),cbind(0,diag(ncol(X)-1))))
  bvec = c(1, rep(0,ncol(X)-1))
  
  synth_model = solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  
  w = synth_model$solution[-1]
  intercept = synth_model$solution[1]
  effects = -intercept + y_after%*%c(1,-w)
  
  return(list("w" = w, "intercept" = intercept, "effects" = effects, "value_min" = (2*synth_model$value +y%*%y) ))
  
}


conformal_fp <- function(pret, postt, q = 2){
  
  full = rbind(pret, postt)
  
  #estimating under the null
  sc_mod = synth_control_est_demean(full,full)
  
  did_mod = (full[,1] - rowMeans(full[,-1])) - (mean(full[,1]) - mean(full[,-1]))
  
  res = sc_mod$effects - did_mod
  
  T0 = nrow(pret) 
  T1 = nrow(postt)
  Tt = T0+T1
  
  conformal_distr <- lapply(0:(Tt-1), function(j){
    
    index = 1:(T0+T1) + j
    
    index = ifelse(index > Tt, index - Tt, index )
    
    res_ord = res[index]
    
    res_post = res_ord[(T0+1):(T0+T1)]
    
    
    return(c((sqrt(1/T1)*sum(abs(res_post)^q))^(1/q),sqrt(1/T1)*abs(sum(res_post))))
  })

  
  conformal_distr = do.call(rbind, conformal_distr)
  
  return(c(mean(conformal_distr[1,1]<=conformal_distr[,1]),mean(conformal_distr[1,2]<=conformal_distr[,2])))
}


conformal_fp_parallel <- function(pret, postt, q = 2){
  
  full = rbind(pret, postt)
  
  #estimating under the null
  sc_mod = synth_control_est_demean(full,full)
  
  did_mod = (full[,1] - rowMeans(full[,-1])) - (mean(full[,1]) - mean(full[,-1]))
  
  res = sc_mod$effects - did_mod
  
  T0 = nrow(pret) 
  T1 = nrow(postt)
  Tt = T0+T1
  
  conformal_distr <- mclapply(0:(Tt-1), function(j){
    
    index = 1:(T0+T1) + j
    
    index = ifelse(index > Tt, index - Tt, index )
    
    res_ord = res[index]
    
    res_post = res_ord[(T0+1):(T0+T1)]
    
    
    return(c((sqrt(1/T1)*sum(abs(res_post)^q))^(1/q),sqrt(1/T1)*abs(sum(res_post))))
  })
  
  
  conformal_distr = do.call(rbind, conformal_distr)
  
  return(c(mean(conformal_distr[1,1]<=conformal_distr[,1]),mean(conformal_distr[1,2]<=conformal_distr[,2])))
}


#Synthetic control estimation which allows for rank-defficient problem
#y_before = (T0) x (J=1) matrix with outcomes before intervention. First column is treated.
#y_after = (T - T0) x (J+1) matrix with outcomes after intervetion. First column is treated.
#constrained = should weights be constrained so optimisation is convex?
synth_control_est_allowrankdef <- function(y_before, y_after, add_up = F, nonnegative = F, start_value = c())
{
  X = y_before[,-1]
  y = y_before[,1]
  if((!add_up)&(!nonnegative))
  {
    w = as.vector(solve(t(X)%*%X)%*%t(X)%*%y) 
    convergence = T
    value = sum((y - X%*%w)^2)
  } else if(add_up&(!nonnegative)){
    y_transf = y - X[,1]
    X_transf = X[,-1] - X[, 1]
    w = as.vector(solve(t(X_transf)%*%X_transf)%*%t(X_transf)%*%y_transf)
    w = c(1-sum(w), w)
    convergence = T
    value = sum((y - X%*%w)^2)
  }else if((!add_up)&nonnegative){
    model = nnls(X,y)
    convergence = T
    w = model$x
    value = model$deviance
  } else{
    #Tries NNLS first
    y_transf = y - X[,1]
    X_transf = X[,-1] - X[, 1]
    model_try = nnls(X_transf, y_transf)
    
    if(sum(model_try$x)<=1)
    {
      w = c(1-sum(model_try$x), model_try$x)
      convergence = T
      value = model_try$deviance
    } else {
      print("Unconstrained nao deu certo")
      ui = rbind(diag(ncol(X)-1), rep(-1, ncol(X)-1))
      ci = c(rep(0, ncol(X)-1),-1)
      
      #BFGS optimization
      w_bfgs = constrOptim(f = function(w,y,X){
        w = c(1- sum(w),w)
        return(sum((y-X%*%w)^2)) }, grad = function(w,y,X){
          w = c(1- sum(w),w)
          return(colSums(2*as.vector((y-X%*%w))*(X[,1]-X[,-1])))
        },
        theta = (model_try$x +1e-5)/(sum(model_try$x)+(length(model_try$x)+1)*1e-5), ui = ui, ci = ci,
        method = "BFGS", y = y, X = X, control = list("maxit"=10000))
      
      #Nelder-Mead optimization
      w_nelder = constrOptim(f = function(w,y,X){
        w = c(1- sum(w),w)
        return(sum((y-X%*%w)^2)) },
        theta = (model_try$x +1e-5)/(sum(model_try$x)+(length(model_try$x)+1)*1e-5), ui = ui, ci = ci,
        method = "Nelder-Mead", y = y, X = X, control = list("maxit"=10000))
      
      if(w_nelder$value<w_bfgs$value)
        w = w_nelder else
          w = w_bfgs
      
      if(length(start_value)>0)
      {
        
        #BFGS optimization
        w_bfgs_1 = constrOptim(f = function(w,y,X){
          w = c(1- sum(w),w)
          return(sum((y-X%*%w)^2)) }, grad = function(w,y,X){
            w = c(1- sum(w),w)
            return(colSums(2*as.vector((y-X%*%w))*(X[,1]-X[,-1])))
          },
          theta = start_value[-1], ui = ui, ci = ci,
          method = "BFGS", y = y, X = X, control = list("maxit"=10000))
        
        #Nelder-Mead optimization
        w_nelder_1 = constrOptim(f = function(w,y,X){
          w = c(1- sum(w),w)
          return(sum((y-X%*%w)^2)) },
          theta = start_value[-1], ui = ui, ci = ci,
          method = "Nelder-Mead", y = y, X = X, control = list("maxit"=10000))
        
        
        if(w_nelder_1$value<w_bfgs_1$value)
          w_1 = w_nelder_1 else
            w_1 = w_bfgs_1
        
        if(w$value>w_1$value)
        {
          print("Valor na ponta eh melhor")
          w= w_1
        }
        
      }
      
      convergence = (w$convergence == 0)
      value = w$value
      w = c(1-sum(w$par),w$par)
    }
  }
  return(list("weights"= w, "effects" = y_after%*%c(1,-w), "convergence" = convergence, "value" = value))
}
