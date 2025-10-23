density.fn=function(t, lambda, pi ) sum(dexp(t,rate=lambda) * pi)
survival.fn=function(t, lambda, pi ) sum(pexp(t,rate=lambda,lower.tail = F) * pi)  

density.fn.lognorm=function(t, mu,s, pi ) sum(dnorm(t,mean=mu, sd=s) * pi)
survival.fn.lognorm=function(t, mu,s, pi ) sum(pnorm(t,mean=mu, sd=s,lower.tail = F) * pi)  

hazard.rate=function(t, lambda, pi )
  density.fn(t, lambda, pi ) / survival.fn(t, lambda, pi )


model_fitting_qqplot=function(result, st1,nu1, st2,nu2, burnin=1e3){
  # library(plyr);library(reshape2); library(ggplot2)
  # st1=dat.clinical$OS.weeks+inc; nu1=dat.clinical$OS.failure
  # st2=dat.synthetic$OS.weeks; nu2=dat.synthetic$OS.failure
  
  n1=length(st1); n2=length(st2)
  ind1=which(nu1); ind2=which(nu2)
  
  alloc_vars=result$`Allocation variables`+1 ###to change from C++ indexing to R
  cluster1= alply( alloc_vars[-(1:burnin),ind1],1, .fun = identity)
  cluster2= alply( alloc_vars[-(1:burnin),n1+ind2],1, .fun = identity)
  lambda1= alply( result$Avg_response_clinical[-(1:burnin),], 1,.fun=identity)
  lambda2= alply( result$Avg_response_control[-(1:burnin),], 1,.fun=identity)
  pi1=alply( result$pimat1[-(1:burnin),] , 1,.fun=identity) 
  pi2=alply( result$pimat2[-(1:burnin),] , 1,.fun=identity) 
  
  
  unifs1=t(mapply(function(clust.id, theta ,st){
    pexp(st, rate = theta[clust.id] )
  },cluster1, lambda1,  MoreArgs = list(st= st1[ind1] ) ))
  
  unifs2=t(mapply(function(clust.id, theta ,st){
    pexp(st, rate = theta[clust.id] )
  },cluster2, lambda2,  MoreArgs = list( st=st2[ind2] ) ))
  unifs=data.frame(unifs=  cbind(unifs1,unifs2)[round(seq(1, nrow(unifs1), length.out = 100 ),0), ] , MCMC=as.factor(1:100) , colval =as.factor(sample.int(9,size=100,replace = T) ))
  
  qplot.dat=melt(unifs,id.vars= c("MCMC","colval"))
  
  ggplot(qplot.dat,aes(sample=value,group=MCMC,colour=colval)) +stat_qq(distribution = stats::qunif,geom = "line",size=.3) + scale_color_brewer(palette="Blues") +
    theme(legend.position = "none") 
}

hazard_plots=function(result, timepoints, burnin=1e3){
  library(plyr); library(reshape2); library(ggplot2)
  # alloc_vars=result$`Allocation variables`+1 ###to change from C++ indexing to R
  # cluster1= alply( alloc_vars[-(1:burnin),1:n1],1, .fun = identity)
  # cluster2= alply( alloc_vars[-(1:burnin),((1+n1):(n1+n2))],1, .fun = identity)
  lambda1= alply( result$Avg_response_clinical[-(1:burnin),], 1,.fun=identity)
  lambda2= alply( result$Avg_response_control[-(1:burnin),], 1,.fun=identity)
  pi1=alply( result$pimat1[-(1:burnin),] , 1,.fun=identity) 
  pi2=alply( result$pimat2[-(1:burnin),] , 1,.fun=identity) 
  
  # hr1= parallel::mcmapply(function(timepoints,lambda, pi) sapply(timepoints, hazard.rate, lambda,pi),
  #                         lambda1,pi1,MoreArgs = list(timepoints=timepoints), mc.cores = parallel::detectCores()/2)
  # hr2= parallel::mcmapply(function(timepoints,lambda, pi) sapply(timepoints, hazard.rate, lambda,pi),
  #                         lambda2,pi1,MoreArgs = list(timepoints=timepoints), mc.cores = parallel::detectCores()/2)
  # hr1.mean=data.frame(HR= rowMeans(hr1), Time=timepoints,Arm="Treatment")
  # hr2.mean=data.frame(HR=rowMeans(hr2), Time=timepoints,Arm="Synthetic")
  # hr.mean.compare=rbind(hr1.mean,hr2.mean)
  
  dens1= parallel::mcmapply(function(timepoints,lambda, pi) sapply(timepoints, density.fn, lambda,pi),
                            lambda1,pi1,MoreArgs = list(timepoints=timepoints), mc.cores = parallel::detectCores()/2)
  surv1= parallel::mcmapply(function(timepoints,lambda, pi) sapply(timepoints, survival.fn, lambda,pi),
                            lambda1,pi1,MoreArgs = list(timepoints=timepoints), mc.cores = parallel::detectCores()/2)
  
  dens2= parallel::mcmapply(function(timepoints,lambda, pi) sapply(timepoints, density.fn, lambda,pi),
                            lambda2,pi1,MoreArgs = list(timepoints=timepoints), mc.cores = parallel::detectCores()/2)
  surv2= parallel::mcmapply(function(timepoints,lambda, pi) sapply(timepoints, survival.fn, lambda,pi),
                            lambda2,pi1,MoreArgs = list(timepoints=timepoints), mc.cores = parallel::detectCores()/2)
  
  
  
  dens1.mean=rowMeans(dens1)
  surv1.mean=rowMeans(surv1)
  
  dens2.mean=rowMeans(dens2)
  surv2.mean=rowMeans(surv2)
  
  hr1=dens1.mean/surv1.mean
  hr2=dens2.mean/surv2.mean
  
  # hazard.ratio=hr1/hr2
  
  
  hr1.dat=data.frame(HR= (hr1), Time=timepoints,Arm="Treatment")
  hr2.dat=data.frame(HR= (hr2), Time=timepoints,Arm="Synthetic")
  hr.mean.compare=rbind(hr1.dat,hr2.dat)
  
  dat.hrplot=melt( hr.mean.compare , id.vars = c("Time","Arm"))
  p.hazard.rate <- ggplot(dat.hrplot,  aes(x=Time,y=value, color=Arm)) +
    geom_smooth(method = "loess",fullrange=T,se=F)
  
  
  
  
  hazard.ratio=data.frame(HR= (hr1/hr2), Time=timepoints )
  p.hazard.ratio <- ggplot(hazard.ratio,  aes(x=Time,y=HR) ) +
    geom_smooth(method = "loess",fullrange=T)
  
  
  
  ###########for hazard-ratio with CI############
  hazard.ratio=data.frame(HR= ((dens1/surv1)/(dens2/surv2))[ ,round(seq(1,ncol(dens1), length.out = min(600,ncol(dens1))), 0) ], Time=timepoints )
  dat.hrplot=melt( hazard.ratio , id.vars = "Time")
  p.hazarad.ratio2 <- ggplot(dat.hrplot,  aes(Time,value)) +
    geom_smooth(method = "loess",fullrange=T,se=T)
  
  ###############################################
  
  list(p.hazard.rate=p.hazard.rate, p.hazarad.ratio=p.hazard.ratio, p.hazarad.ratio.withCI=p.hazarad.ratio2)
}

hazard_plots_lognormal=function(result, timepoints,parallel=F,  burnin=1e3){
  library(plyr); library(reshape2); library(ggplot2)
  # alloc_vars=result$`Allocation variables`+1 ###to change from C++ indexing to R
  # cluster1= alply( alloc_vars[-(1:burnin),1:n1],1, .fun = identity)
  # cluster2= alply( alloc_vars[-(1:burnin),((1+n1):(n1+n2))],1, .fun = identity)
  
  mu1=alply( result$Lognormal_params1[[1]][-(1:burnin),], 1,.fun=identity) 
  mu2=alply( result$Lognormal_params2[[1]][-(1:burnin),], 1,.fun=identity) 
  s1=alply( result$Lognormal_params1[[2]][-(1:burnin),], 1,.fun=purrr::compose(sqrt, identity))
  s2=alply( result$Lognormal_params2[[2]][-(1:burnin),], 1,.fun=purrr::compose(sqrt, identity)) 
  
  pi1=alply( result$pimat1[-(1:burnin),] , 1,.fun=identity) 
  pi2=alply( result$pimat2[-(1:burnin),] , 1,.fun=identity) 
  
  # hr1= parallel::mcmapply(function(timepoints,lambda, pi) sapply(timepoints, hazard.rate, lambda,pi),
  #                         lambda1,pi1,MoreArgs = list(timepoints=timepoints), mc.cores = parallel::detectCores()/2)
  # hr2= parallel::mcmapply(function(timepoints,lambda, pi) sapply(timepoints, hazard.rate, lambda,pi),
  #                         lambda2,pi1,MoreArgs = list(timepoints=timepoints), mc.cores = parallel::detectCores()/2)
  # hr1.mean=data.frame(HR= rowMeans(hr1), Time=timepoints,Arm="Treatment")
  # hr2.mean=data.frame(HR=rowMeans(hr2), Time=timepoints,Arm="Synthetic")
  # hr.mean.compare=rbind(hr1.mean,hr2.mean)
  
  # timepoints=log(seq(40,175,length.out=40) ) ##comment this
  timepoints=log(timepoints )
  ncores=ifelse(parallel,parallel::detectCores()/2,1)
  
  dens1= parallel::mcmapply(function(timepoints,mu,s, pi) sapply(timepoints, density.fn.lognorm, mu,s,pi),
                            mu1,s1,pi1,MoreArgs = list(timepoints=timepoints), mc.cores = ncores)
  surv1= parallel::mcmapply(function(timepoints,mu,s, pi) sapply(timepoints, survival.fn.lognorm, mu,s,pi),
                            mu1,s1,pi1,MoreArgs = list(timepoints=timepoints), mc.cores = ncores)
  
  dens2= parallel::mcmapply(function(timepoints,mu,s, pi) sapply(timepoints, density.fn.lognorm, mu,s,pi),
                            mu2,s2,pi1,MoreArgs = list(timepoints=timepoints), mc.cores = ncores)
  surv2= parallel::mcmapply(function(timepoints,mu,s, pi) sapply(timepoints, survival.fn.lognorm, mu,s,pi),
                            mu2,s2,pi1,MoreArgs = list(timepoints=timepoints), mc.cores = ncores)
  
  dens1.mean=rowMeans(dens1)
  surv1.mean=rowMeans(surv1)
  
  dens2.mean=rowMeans(dens2)
  surv2.mean=rowMeans(surv2)
  
  hr1=dens1.mean/surv1.mean
  hr2=dens2.mean/surv2.mean
  
  # hazard.ratio=hr1/hr2
  
  
  hr1.dat=data.frame(HR= (hr1), Time=exp(timepoints),Arm="Treatment")
  hr2.dat=data.frame(HR= (hr2), Time=exp(timepoints),Arm="Synthetic")
  hr.mean.compare=rbind(hr1.dat,hr2.dat)
  
  dat.hrplot=melt( hr.mean.compare , id.vars = c("Time","Arm"))
  p.hazard.rate <- ggplot(dat.hrplot,  aes(x=Time,y=value, color=Arm)) +
    geom_smooth(method = "loess",fullrange=T,se=F)
  
  hazard.ratio=data.frame(HR= (hr1/hr2), Time=exp(timepoints) )
  p.hazard.ratio <- ggplot(hazard.ratio,  aes(x=Time,y=HR) ) +
    geom_smooth(method = "loess",fullrange=T)#+ylim(c(0,2))
  
  
  
  ###########for hazard-ratio with CI############
  tmp.HR=((dens1/dens2)*(surv2/surv1))  #(dens1/surv1)/(dens2/surv2)
  
  HR.mean=rowMeans(tmp.HR)
  HR.qtiles=apply(tmp.HR,1,quantile,probs=c(.025,.975))
  
  # hazard.ratio=data.frame(HR= ( tmp.HR )[ ,round(seq(1,ncol(dens1), length.out = min(60,ncol(dens1))), 0) ], Time=exp(timepoints) )
  
  hazard.ratio=data.frame(HR= HR.mean, HR.band=t(HR.qtiles),  Time=(exp(timepoints)) )
  p.hazarad.ratio2=ggplot(data = hazard.ratio,  aes(Time,  HR)) + geom_line()+
    geom_ribbon(data=hazard.ratio,aes(ymin=HR.band.2.5.,ymax=HR.band.97.5.),alpha=0.3)+ylab("Hazard ratio")
  
  
  # dat.hrplot=melt( hazard.ratio , id.vars = "Time")
  # p.hazarad.ratio2 <- ggplot(dat.hrplot,  aes(Time,value)) +
  #   geom_smooth(method = "loess",fullrange=T,se=T)
  
  ###############################################
  
  list(p.hazard.rate=p.hazard.rate, p.hazarad.ratio=p.hazard.ratio, p.hazarad.ratio.withCI=p.hazarad.ratio2, Hazard.Ratio=tmp.HR)
}


gen_roc_for_covariates_bart=function(wts, eta.cat1, eta.cat2, test_prop=.15, q_seq, nmc=50, WR=T){
  #wts=weight of data points in synthetic arm
  #eta.cat1=clinical arm; eta.cat=synthetic arm
  #tst_prop= proportion of test observations to predict in random forest to make the roc
  #q_seq=x axis of the roc
  #nmc= number replicates the random forest is performed
  #WR=whether sampling synthetic data with replace
  nsamp1=nrow(eta.cat1)
  nsamp2=nrow(eta.cat2)
  if(length(wts)!= nsamp2)
    stop("length(wts)!=nrow(eta.cat2)")
  wts_random=rep(1,nsamp2)
  ntest=round(test_prop*nsamp1,0)
  
  library(foreach);library(doParallel)
  cl=makeCluster(parallel::detectCores()-1)
  registerDoParallel(cl)
  roc=foreach(q = q_seq,.combine = rbind,.packages = c("randomForest","caret","dbarts")) %dopar%{
    sampled.data=eta.cat2[ sample.int(nsamp2,size = nsamp1,prob=wts,replace =WR),]
    # sampled.data=eta.cat2[ subsample,]
    
    dat=rbind(eta.cat1,sampled.data)
    labs=factor(c(rep(1,nsamp1),rep(2,nsamp1)))
    
    # dat.imputed=round(na.roughfix(x=dat.factor,y=labs,iter=5)[,-1],0)
    dat.imputed=round(randomForest:: na.roughfix(dat),0)
    
    ind_list=replicate(nmc,c((1:nsamp1)[sample.int(nsamp1,size=ntest,replace = F) ],((nsamp1+1):(2*nsamp1))[sample.int(nsamp1,size=ntest,replace = F) ]),simplify = F)
    res=sapply(ind_list, function(inds,q, ntest,dat.imputed,labs){
      # inds=c((1:nsamp1)[sample.int(nsamp1,size=ntest,replace = F) ],((nsamp1+1):(2*nsamp1))[sample.int(nsamp1,size=ntest,replace = F) ])
      traindat=data.frame(x=dat.imputed[-inds,],y=labs[-inds])
      testdat=data.frame(x=dat.imputed[inds,])
      
      # rffit=randomForest(y~.,traindat,classwt=c(q,1-q))
      #########randomForest
      # rffit=randomForest(y~.,traindat)
      # classprob1=predict(rffit,testdat,type="prob")[,1]
      ##################
      
      #########BART
      rffit <- bart2(as.numeric(y)-1~. , data = traindat,test=testdat,combineChains=T,verbose = F)
      classprob1=1-colMeans(pnorm(rffit$yhat.test)); print(classprob1)
      ##################
      pred=as.factor(sapply(classprob1, function(pp,q) ifelse(pp>=q,1,2),q=q))
      # levels(pred)=c("1","2")
      trulab=as.factor(labs[inds])
      
      zz=numeric(2)
      zz[1]=caret::sensitivity(pred,trulab,positive=levels(trulab)[1])
      zz[2]=caret::specificity(pred,trulab,positive=levels(trulab)[1])
      names(zz)=c("Sensitivity","Specificity")
      zz=sapply(zz, function(z) ifelse(is.na(z),0,z))
      return(zz)
    },q=q,ntest=ntest,dat.imputed=dat.imputed, labs=labs)
    rowMeans(res)
  }
  roc=rbind(roc,diag(2))
  roc=data.frame(roc,Weight="IS")
  
  roc_random=foreach(q = q_seq,.combine = rbind,.packages = c("randomForest","caret","dbarts","ROCR")) %dopar%{
    sampled.data=eta.cat2[ sample.int(nsamp2,size = nsamp1,prob=wts_random,replace =WR),]
    
    dat=rbind(eta.cat1,sampled.data)
    labs=factor(c(rep(1,nsamp1),rep(2,nsamp1)))
    # dat.imputed=round(rfImpute(x=dat,y=labs,iter=5)[,-1],0)
    dat.imputed=round(na.roughfix(dat),0)
    
    ind_list=replicate(nmc,c((1:nsamp1)[sample.int(nsamp1,size=ntest,replace = F) ],((nsamp1+1):(2*nsamp1))[sample.int(nsamp1,size=ntest,replace = F) ]),simplify = F)
    res=sapply(ind_list, function(inds,q, ntest,dat.imputed,labs){
      # inds=c((1:nsamp1)[sample.int(nsamp1,size=ntest,replace = F) ],((nsamp1+1):(2*nsamp1))[sample.int(nsamp1,size=ntest,replace = F) ])
      traindat=data.frame(x=dat.imputed[-inds,],y=labs[-inds])
      testdat=data.frame(x=dat.imputed[inds,])
      
      # rffit=randomForest(y~.,traindat,classwt=c(q,1-q))
      #########randomForest
      # rffit=randomForest(y~.,traindat)
      # classprob1=predict(rffit,testdat,type="prob")[,1]
      ##################
      
      #########BART
      rffit <- bart2(as.numeric(y)-1~. , data = traindat,test=testdat,combineChains=T, verbose = F)
      classprob1=1-colMeans(pnorm(rffit$yhat.test)); print(classprob1)
      ##################
      pred=as.factor(sapply(classprob1, function(pp,q) ifelse(pp>=q,1,2),q=q))
      # levels(pred)=c("1","2")
      trulab=as.factor(labs[inds])
      
      zz=numeric(2)
      zz[1]=caret::sensitivity(pred,trulab,positive=levels(trulab)[1])
      zz[2]=caret::specificity(pred,trulab,positive=levels(trulab)[1])
      names(zz)=c("Sensitivity","Specificity")
      zz=sapply(zz, function(z) ifelse(is.na(z),0,z))
      return(zz)
    },q=q,ntest=ntest,dat.imputed=dat.imputed, labs=labs)
    rowMeans(res)
  }
  stopCluster(cl)
  roc_random=rbind(roc_random,diag(2))
  roc_random=data.frame(roc_random,Weight="Random")
  
  roc_merged<-rbind(roc,roc_random)
}

#' Calculate \emph{Average Treatment Effect (ATE)} from the output of \code{\link{cappmx_fit}}.
#'
#' @param result Output of \code{\link{cappmx_fit}}.
#' @param burnin Number of burnin samples to discard
#'
#' @return The estimated \emph{ATE} from the CAPPMx fit.
#' @export
average_trt_effect=function(result, burnin=200){
  mean(rowSums( ( (result$Lognormal_params1[[1]] -result$Lognormal_params2[[1]])
                  * result$pimat1)[-(1:burnin),] ) )
}

#' Fit CAPPMx Model
#' 
#' An implementation of the CAPPMx by \insertCite{chandra_GBM22;textual}{CAPPMx}. 
#' 
#' Fit the CAPPMx on treatment arm and RWD with survival endpoints.
#'  Currently the package only supports right-censored outcomes.
#'
#' @param cat_cov_trt The matrix of categorical variables in the Treatment arm.  
#' @param cont_cov_trt The matrix of continuous variables in the Treatment arm. 
#' @param response_trt The vector of \strong{\code{log}-transformed} survival times in the Treatment arm.
#' @param surv_ind_trt A logical vector of the same length as \code{response_trt} indicating whether the corresponding survival times is an observed failure or right-censored.
#' @param cat_cov_rwd The matrix of categorical variables in the RWD.  
#' @param cont_cov_rwd The matrix of continuous variables in the RWD.
#' @param response_rwd The vector of \strong{\code{log}-transformed} survival times in the RWD.
#' @param surv_ind_rwd A logical vector of the same length as \code{response_rwd} indicating whether the corresponding survival times is an observed failure or right-censored.
#' @param nmix Number of mixture components.
#' @param nrun Number of MCMC iterations.
#' @param burn Number of burn-in iterations.
#' @param thin The thinning interval.
#' @param del_range_response A vector of length 2 indicating the range of the Hamiltonian Monte Carlo (HMC) tuning parameter \eqn{\Delta t}  in the Leapgrog step for updating the hyperparameters of the response model. The lower bound MUST be positive.
#' @param nleapfrog_response A positive integer indicating the number of Leapfrog steps in the  HMC update step for updating the hyperparameters of the response model. 
#' @param del_range_alp1 A vector of length 2 indicating the range of the HMC tuning parameter \eqn{\Delta t}  for updating the Dirichlet concentration parameter \eqn{\alpha_1} in the mixture model of Treatment arm.
#' @param nleapfrog_alp1 A positive integer indicating the number of Leapfrog steps in the  HMC update step for updating \eqn{\alpha_1}. 
#' @param del_range_alp2 A vector of length 2 indicating the range of the HMC tuning parameter \eqn{\Delta t}  for updating the Dirichlet concentration parameter \eqn{\alpha_2} in the mixture model of Treatment arm.
#' @param nleapfrog_alp2 A positive integer indicating the number of Leapfrog steps in the  HMC update step for updating \eqn{\alpha_2}.
#'
#' @note
#' The following \strong{MUST BE} ensured 
#' \itemize{
#' \item{\code{nrow(cat_cov_trt)=nrow(cont_cov_trt)=length(response_trt)=length(surv_ind_trt)}}
#' \item{\code{nrow(cat_cov_rwd)=nrow(cont_cov_rwd)=length(response_rwd)=length(surv_ind_rwd)}}
#' \item{\code{ncol(cont_cov_trt)=ncol(cont_cov_rwd)}}
#' \item{\code{ncol(cat_cov_trt)=ncol(cat_cov_rwd)}}
#' \item{The category indicators of \code{cat_cov_trt} and \code{cat_cov_rwd} \strong{MUST BE} non-negative integers starting from 0. 
#' For example, a covariate with three categories must be indicated using 0,1 and 2.}
#' }
#' The package is under development. The user must format the inputs according to the above.
#' 
#' The HMC tuning parameters may need to be adjusted manually if the default settings  yield low acceptance rates.
#' 
#' @returns  A list with the following elements:
#' \describe{
#' \item{\code{pimat1}}{A matrix of MCMC samples of mixture weights \eqn{(\pi_{1,1},\dots,\pi_{1,\texttt{nmix}})} of the Treatment arm. Each row corresponds to a MCMC sample.}
#' \item{\code{pimat2}}{A matrix of MCMC samples of mixture weights \eqn{(\pi_{2,1},\dots,\pi_{2,\texttt{nmix}})} of the RWD. Each row corresponds to a MCMC sample.}
#' \item{\code{pimat2} }{A matrix of MCMC samples of mixture weights \eqn{(\pi_{2,1},\dots,\pi_{2,\texttt{nmix}})} of the RWD. Each row corresponds to a MCMC sample.}
#' \item{Weights2}{A matrix of MCMC samples of impotance resampling weights attached to the samples in the RWD. Each row corresponds to a MCMC sample. \code{ncol(Weights2)=} the  number of samples in the RWD.}
#' \item{Lognormal_params1}{A list of length 2-- the first element of the list is a \code{n_mc}\eqn{\times\texttt{nmix}}  matrix of MCMC samples containing the location parameters \eqn{(\mu_{1,1},\dots,\mu_{1,\texttt{nmix}})} of the location-scale mixture of normals  on the \code{log}-survival outcomes in the Treatment arm;
#'    the second element of the list is a \code{n_mc}\eqn{\times\texttt{nmix}} matrix of MCMC samples containing the scale parameters \eqn{(\sigma_{1,1}^{2},\dots,\sigma_{1,\texttt{nmix}}^{2} )} of the location-scale mixture of normals  on the \code{log}-survival outcomes in the Treatment arm. 
#'  We let \code{n_mc} denote the number of MCMC samples returned by the \code{cappmx_fit} function.}
#' \item{Lognormal_params2}{A list of length 2--organized in the same manner as \code{Lognormal_params1} but  corresponds to the mixture model on the survival outcomes of the RWD.}
#' \item{Unifs}{A \code{n_mc}\eqn{\times} \code{n} matrix for assessing model fit where \code{n=(length(response_trt)+length(response_rwd))}. 
#' In case the model fit is good, the distribution of the elements of each row of the \code{Unifs} should approximately be \eqn{Unif(0,1)}.}
#' \item{Lognormal_hyperparams}{A \code{n_mc}\eqn{\times}2 matrix of MCMC samples of the hyperparameters shared across the two arms in the response models.}
#' \item{Dirichlet_params}{A \code{n_mc}\eqn{\times}2 matrix of MCMC samples of \eqn{(\alpha_1,\alpha_2)}.}
#' \item{Acceptance rates}{A vector of length 3 indicating the acceptance rates of \eqn{\alpha_1}, \eqn{\alpha_2} and the \code{Lognormal_hyperparams}, respectively.}
#'}
#' @references
#' \insertAllCited{}
#' @export
cappmx_fit=function(cat_cov_trt=NULL,cont_cov_trt=NULL, response_trt, surv_ind_trt,
                    cat_cov_rwd=NULL,cont_cov_rwd=NULL, response_rwd, surv_ind_rwd,
                    nmix=15, nrun=5e3,burn=1e3,thin=5,
                    del_range_response=c(.005,.02)*15, nleapfrog_response=3,
                    del_range_alp1 = c(.1,.3)*5, nleapfrog_alp1 = 4,
                    del_range_alp2 = c(.1,.3)*3, nleapfrog_alp2 = 3){
  
  isnullcat1=is.null(cat_cov_trt); isnullcont1=is.null(cont_cov_trt); isnullcat2=is.null(cat_cov_rwd); isnullcont2=is.null(cont_cov_rwd)
  
  if(isnullcat1 & isnullcont1 & isnullcat2 & isnullcont2 )
    stop("All covariate matrices are null!")
  
  if(class(surv_ind_trt)!="logical" | class(surv_ind_rwd)!="logical" )
    stop("Survival indicators must be logical variables!")
  # if(ncol(cat_cov_trt)!=ncol(cat_cov_rwd)| ncol(cont_cov_trt)!=ncol(cont_cov_rwd)| 
  #    nrow(cat_cov_trt)!=nrow(cont_cov_trt) | )
  
  ######CHECKS FOR NULL COV MATRICES AND SETTING NULLS
  ######IN TRT ARM
  if(isnullcat1 ){
    if(isnullcont1 )
      stop("Both covariate matrices are null in the treatment arm!") else{
        cat_cov_trt=matrix(NA,nrow=nrow(cont_cov_trt),ncol=1)
      }
  }
  if(isnullcont1 ){
    if(isnullcat1)
      stop("Both covariate matrices are null in the treatment arm!") else{
        cont_cov_trt=matrix(NA,nrow=nrow(cat_cov_trt),ncol=1)
      }
  }
  ########################
  
  ######IN RWD######
  if(isnullcat2 ){
    if(isnullcont2 )
      stop("Both covariate matrices are null in the RWD!") else{
        cat_cov_rwd=matrix(NA,nrow=nrow(cont_cov_rwd),ncol=1)
      }
  }
  if(isnullcont2 ){
    if(isnullcat2)
      stop("Both covariate matrices are null in the RWD!") else{
        cont_cov_rwd=matrix(NA,nrow=nrow(cat_cov_rwd),ncol=1)
      }
  }
  ########################################################################
  
  eta.cont1=as.matrix(cont_cov_trt)
  eta.cont2=as.matrix(cont_cov_rwd)
  eta.cont=rbind(eta.cont1,eta.cont2)
  
  eta.cat1=as.matrix(cat_cov_trt); rownames(eta.cat1)=NULL ;nsamp1=nrow(eta.cat1)
  eta.cat2=as.matrix(cat_cov_rwd); rownames(eta.cat2)=NULL ;nsamp2=nrow(eta.cat2)
  eta.cat=rbind(eta.cat1,eta.cat2)
  ncats=apply(eta.cat,2, function(x) length(setdiff(unique(x),NA)))
  
  ################################################################
  #########GET THE INDICES OF THE NON-MISSING VALUES##############
  ################################################################
  
  ###############
  ###for the categorical covs 
  ###############
  non.na.inds=apply(eta.cat,1, function(eta) purrr::compose( which,"!",is.na)(eta)-1,simplify = F)
  
  ###############
  ###for the continuous covs 
  ###############
  non.na.inds_cont=apply(eta.cont,1, function(eta) purrr::compose( which,"!",is.na)(eta)-1,simplify = F)
  
  ##################################################################
  
  ###############Initial values for MCMC input###########
  ## For computing distance matrix
  
  hamm_dist=function(x,y){
    nonnax=which(!is.na(x))
    nonnay=which(!is.na(y))
    int=intersect(nonnax,nonnay)
    if(length(int)==0)
      return (length(x)) else return(e1071::hamming.distance(x[int],y[int]))
  }
  
  ds=matrix(0,nsamp2,nsamp2)
  if(!isnullcat2 ){
    for(i in 1:nsamp2){
      for(j in (i):nsamp2){
        if(j>nsamp2)
          print(j)
        ds[i,j]=hamm_dist(eta.cat2[i,],eta.cat2[j,])
        ds[j,i]=ds[i,j]
      }
    }
  }
  ds=as.dist(ds)
  
  if(!isnullcont2 ){
    x=eta.cont2; x[!is.finite(x)]=0
    ds=ds+ dist(x)
    rm(x)
  }
  
  fit <- hclust (ds , method = "ward.D2" )
  labels2 <- cutree ( fit , k =min(nmix,7)) -1
  
  #############NULL CHECKS###############
  if(isnullcat2 ){
    eta.imputed2=na.roughfix(eta.cont2)
    eta.imputed1=na.roughfix(eta.cont1)
  } else if(isnullcont2 ){
    eta.imputed2=na.roughfix(eta.cat2)
    eta.imputed1=na.roughfix(eta.cat1)
  } else if(!isnullcat1 & !isnullcont1 & !isnullcat2 & !isnullcont2 ){
    eta.imputed2=cbind(na.roughfix(eta.cat2),na.roughfix(eta.cont2))
    eta.imputed1=cbind(na.roughfix(eta.cat1),na.roughfix(eta.cont1))
  }  else stop("Check input covariates!")
  ############################################################
  
  rffit=randomForest::randomForest(x=eta.imputed2,y=as.factor(labels2), na.action=na.roughfix)
  labels1=as.numeric(predict(rffit,eta.imputed1))-1
  
  # cat_cov_trt,cont_cov_trt, response_trt, surv_ind_trt,
  # cat_cov_rwd,cont_cov_rwd, response_rwd, surv_ind_rwd
  
  # ind1=which(surv_ind_trt); ind2=which(surv_ind_rwd)
  log_observed_failure_times=(c( response_trt[surv_ind_trt], response_rwd[surv_ind_rwd]))
  
  ####log-normal hyperparameter tuning
  df0=1
  a0=10.1 ##must be greater than 1
  m= 2*(a0-1) #var(log_observed_failure_times)*(a0-1) ###mean of variance of outcomes
  s=5
  b_v=log1p(s/m^2)
  b_m=log(m)- b_v/2
  mu_m= mean(log_observed_failure_times)
  mu_v= 1
  
  common_atoms_cat_lognormal( nmix, ncat=ncats,
                              a0=a0, df0=df0, mu_m=mu_m, mu_v=mu_v, b_m=b_m, b_v=b_v,#normal and lognormal hyperparameters for the response
                              nrun=nrun, burn=burn, thin=thin,
                              eta.cat1  , eta.cat2,
                              eta.cont1, eta.cont2,
                              response_trt, surv_ind_trt,
                              response_rwd, surv_ind_rwd,
                              non.na.inds, non.na.inds_cont,
                              labels1, labels2, 
                              del_range_lognorm=del_range_response, nleapfrog_lognorm=nleapfrog_response,
                              alpha_hyper = c(1,10),del_range_alp1 = del_range_alp1, nleapfrog_alp1 = nleapfrog_alp1,
                              del_range_alp2 = del_range_alp2, nleapfrog_alp2 = nleapfrog_alp2)
}

#' Investigating Population Equivalence
#' 
#'  @details
#'  A randomly resampled set of covariates is selected from RWD with weights being the \emph{Importance-Resampling (IS)} ones
#'  derived from the output of \code{\link{cappmx_fit}}.
#'  The resampled subset and the covariate population in the treatment arm are merged 
#'  into a single dataset and then we try to classify patients in the merged sample as originally RWD or single-arm treatment cohort.
#'  For classification, we use Bayesian Additive Regression Tree (BART) \insertCite{bart2010}{CAPPMx} and report  the 
#'  area under the receiver operating characteristic curve (AUC) of the classification accuracy.
#'  For comparison, we also subsample randomly (instead of using the IS weights) from RWD
#'  to create another synthetic control arm and report the AUC.
#' 
#' @param cat_cov_trt The matrix of categorical variables in the Treatment arm.  
#' @param cont_cov_trt The matrix of continuous variables in the Treatment arm. 
#' @param cat_cov_rwd The matrix of categorical variables in the RWD.  
#' @param cont_cov_rwd The matrix of continuous variables in the RWD. 
#' @param result Output of \code{\link{cappmx_fit}} function.
#' @param burnin Number of burnin samples to discard from \code{result}.
#' @param cv_prop The value of \eqn{p} for a leave-\eqn{100\times p \%}-out 
#' cross-validation for computing out-of-sample classification accuracy.
#'
#' @return A vector of AUC values corresponding to the proposed IS and 
#' random resampling schemes, respectively.
#' @references
#' \insertAllCited{}
#' @export
get_auc=function(cat_cov_trt=NULL,cont_cov_trt=NULL,
                 cat_cov_rwd=NULL,cont_cov_rwd=NULL,
                 result, burnin=200,cv_prop=.15){
  isnullcat1=is.null(cat_cov_trt); isnullcont1=is.null(cont_cov_trt); isnullcat2=is.null(cat_cov_rwd); isnullcont2=is.null(cont_cov_rwd)
  if((isnullcat1&isnullcont1) | (isnullcat2&isnullcont2) )
    stop("Check input covariance matrices!")
  #############NULL CHECKS###############
  if(isnullcat2 ){
    X2=na.roughfix(cont_cov_rwd)
    X1=na.roughfix(cont_cov_trt)
  } else if(isnullcont2 ){
    X2=na.roughfix(cat_cov_rwd)
    X1=na.roughfix(cat_cov_trt)
  } else if(!isnullcat1 & !isnullcont1 & !isnullcat2 & !isnullcont2 ){
    X2=cbind(na.roughfix(cat_cov_rwd),na.roughfix(cont_cov_rwd))
    X1=cbind(na.roughfix(cat_cov_trt),na.roughfix(cont_cov_trt))
  }  else stop("Check input covariates!")
  ############################################################
  
  nsamp2=nrow(X2); n1= nrow(X1) #nsamp1
  WR=T
  ntest=round(cv_prop*n1,0)
  n.cv=round(1/cv_prop,0)
  
  ############FIND AUC FOR IMPORTANCE RESAMPLED DATA###########
  wts=colMeans(result$Weights2[-(1:burnin),])
  sampled.data=X2[ sample.int(nsamp2,size = nsamp1,prob=wts,replace =WR),]
  dat=rbind(X1,sampled.data)
  labs=factor(c(rep(1,nsamp1),rep(2,nsamp1)))
  auc_sampled=foreach(j=1:n.cv,.combine=c,.packages = c("dbarts","ROCR")) %do%{
    inds=c(sample.int(n1,size=ntest,replace = F) ,
           n1+sample.int(n1 ,size=ntest,replace = F) )
    traindat=data.frame(x=dat[-inds,],y=labs[-inds])
    testdat=data.frame(x=dat[inds,])
    #####################BART
    rf_model <- bart2(as.numeric(y)-1~. , data = traindat,test=testdat,combineChains=T)
    rf_prediction=colMeans(pnorm(rf_model$yhat.test))
    ########################
    
    ROC_rf=ROCR::prediction( rf_prediction, as.numeric(labs[inds])-1) ##for bart only
    auc.perf=performance(ROC_rf,measure = "auc")
    (auc.perf@y.values[[1]])
    ########################
  }
  
  AUC_adjusted=median(auc_sampled)
  rm(sampled.data,dat,labs,wts)
  #################################################################################
  
  ############FIND AUC FOR RANDOM RESAMPLED DATA###########
  wts=rep(1,nsamp2)
  sampled.data=X2[ sample.int(nsamp2,size = nsamp1,prob=wts,replace =WR),]
  dat=rbind(X1,sampled.data)
  labs=factor(c(rep(1,nsamp1),rep(2,nsamp1)))
  
  auc_sampled=foreach(j=1:n.cv,.combine=c,.packages = c("dbarts","ROCR")) %do%{
    inds=c(sample.int(n1,size=ntest,replace = F) ,
           n1+sample.int(n1 ,size=ntest,replace = F) )
    traindat=data.frame(x=dat[-inds,],y=labs[-inds])
    testdat=data.frame(x=dat[inds,])
    #####################BART
    rf_model <- bart2(as.numeric(y)-1~. , data = traindat,test=testdat,combineChains=T)
    rf_prediction=colMeans(pnorm(rf_model$yhat.test))
    ########################
    
    ROC_rf=ROCR::prediction( rf_prediction, as.numeric(labs[inds])-1) ##for bart only
    auc.perf=performance(ROC_rf,measure = "auc")
    (auc.perf@y.values[[1]])
    ########################
  }
  
  AUC_random=median(auc_sampled)
  rm(sampled.data,dat,labs,wts)
  #################################################################################
  aucs=c(AUC_adjusted,AUC_random)
  names(aucs)=c("IR_Adjusted", "Random")
  aucs
}

gen_non_lin_output=function(x1,bt, sig_sq=1){
  if(length(bt)<5) stop("length(bt)<5")
  # bt=c(runif(3,10,20), runif(2,-10,-5))
  p=ncol(x1);n=nrow(x1)
  y1=( abs(X1[,1]-2)<.5) *(X1[,2]>1.5) *bt[1] + ( (X1[,p]>=1)*(X1[,p-1]>=1)) *bt[2] +  
    ( (X1[,5])*(X1[,6])) *bt[3] + ( (abs(X1[,7]-1.5)>.1)*(abs(X1[,8]-1.5)>.1) ) *bt[4]+
    (X1[,3] >1.5)*(X1[,4] >1.5) *bt[5] +rnorm(nsamp1,mean = 0, sd=sqrt(sig_sq))
  # y1=bt[1] *(abs(x1[,1]-2)<.5)- bt[2]*(abs(x1[,2]-2)<.5)+  
  #   bt[3]*( (abs(x1[,3]-1)) *(abs(x1[,4]-1.5)^2 ) )+
  #   bt[4]* (abs(x1[,p-1]-1)<.1)*(abs(x1[,p]-1)<.1) +bt[5] *(abs(x1[,p-2]-1)<.1) +
  #   rnorm(n,mean=0,sd = sqrt(sig_sq))
  y1
}

gen_out=function(x1,bt, sig_sq=1){
  n=nrow(x1); p=ncol(x1)
  # y1=50* (abs(x1[,1]-2)< .9) *(abs(x1[,2]-2)< .9)- 50* (abs(x1[,3]-2)< .9)*(abs(x1[,4]-2)< .9) +
  #   250*(abs(x1[,5]-2)< .9)*(abs(x1[,6]-2)< .9) + rnorm(n,sd = sqrt(sig_sq))
  
  # bt1=runif(1,40,60); bt2=runif(1,40,60) ; bt3=runif(1,225,275) ; bt4=runif(1,-1,-5)
  bt[1]* (x1[,1]>1.25) * (x1[,2]>1.25)- bt[2]* (x1[,3]>1.25)*(x1[,4]>1.25) +
    bt[3]*(x1[,5]>1.25)*(x1[,6]>1.25) + bt[4]*(x1[,p]>=1)*(x1[,p-1]>=1)  + 
    rnorm(n,sd = sqrt(sig_sq))
}

