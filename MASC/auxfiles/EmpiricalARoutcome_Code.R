
estimate.tfe<-function(data,MAXT,nfactor=4){
  #Factor model estimation after removing time effects, following Bai (2009)
  #Also dividing by St. dev, so NORMALIZING:
  covdata<-data$outcome[1:MAXT,]
  W<-t(apply(covdata,MARGIN=1,function(x) (x - mean(x))/sd(x)))
  factors<-eigen(W%*%t(W)/(dim(covdata)[1]*dim(covdata)[2]))$vectors[,1:nfactor]*(dim(covdata)[1]^0.5)
  loadings<-t(t(factors)%*%W/dim(covdata)[1])
  fittedvals<-factors%*%t(loadings)
  fittedvals<- fittedvals*apply(covdata,MARGIN=1,FUN=function(x) sd(x)) +
    outer(apply(covdata,MARGIN=1,function(x) mean(x)),rep(1,dim(covdata)[2])) 
  residuals = covdata-fittedvals
  
  fittedvals.outcome<-fittedvals
  residuals.outcome<-residuals[grepl("gdpcap",row.names(residuals)),]
  
  
  rownames(factors)<-rownames(covdata)
  
  
  #Estimate AR coefficients for outcome path only
  #using the variance-covariance matrix without degree of freedom adjustment, as implied in Bai (2009) (pg 1250):
  arcoefs<-apply(residuals.outcome,MARGIN=2,FUN=function(x) ar(x,aic=FALSE,order.max=1,demean=FALSE)$ar)
  varcov<-diag(apply(residuals.outcome,MARGIN=2,FUN=function(x) ar(x,aic=FALSE,order.max=1,demean=FALSE)$var.pred))
  
  return(list(factors=factors,loadings=loadings,fittedvals=fittedvals.outcome,varcov=varcov,arcoef=arcoefs))
  
}



draw.ar<-function(r,MAXT,fittedval,arcoef,varcov.ar,
                  data.covariates,
                  draw=TRUE){
  
  errors <- mvrnorm(n=MAXT+burn_in,mu=rep(0,nrow(varcov.ar)),Sigma=varcov.ar)
  drawdata<-matrix(NA,nrow=nrow(fittedval),ncol=ncol(fittedval))
  Yres <-matrix(0,nrow=nrow(fittedval)+burn_in,ncol=ncol(fittedval))
  if(ncol(arcoef)==1){
    Yres[1,]<-errors[1,]
    for(t in 2:nrow(Yres)){
      Yres[t,]<-errors[t,]+arcoef*Yres[t-1,]
    }
  }
  if(ncol(arcoef)==2){
    Yres[1,]<-errors[1,]
    Yres[2,]<-errors[2,] + arcoef[,1]*Yres[1,]
    for(t in 3:nrow(Yres)){
      Yres[t,]<-errors[t,]+arcoef[,1]*Yres[t-1,]+arcoef[,2]*Yres[t-2,]
    }
  }
  for(t in 1:nrow(fittedval)){
    drawdata[t,]<- fittedval[t,]+Yres[burn_in+t,]
  }
  drawdata<-as.data.table(drawdata)
  temp<-as.data.table(t(drawdata))
  names(temp)<-rownames(fittedval)
  
  drawdata.covariates<-copy(data.covariates)
  setcolorder(drawdata.covariates,c("time","unit"))
  for(u in 1:length(unique(drawdata.covariates$unit))){
    unitnum<-unique(drawdata.covariates$unit)[u]
    drawdata.covariates[unit==unitnum,gdpcap:=drawdata[,u,with=FALSE]]
  }
  

  return(list(outcome=drawdata,
              covariates=drawdata.covariates))
}




#lambda=3.0: interactive FE model, no noise
#lambda=3.3: interactive FE model, white noise process (the one it implies)
###MAIN SIMULATION FUNCTION###
empirical_MC<-function(r=NA,data,philength=NA,Kvals=NA,firstforecast=NA,dump.info=FALSE,fold.errors=TRUE,includeoutput=FALSE,includedraw=FALSE){
  MAXT <- nrow(data$outcome)
  
  
  lambda<-data$lambda
  fitlambda<-floor(lambda)
  errorlambda<-round((lambda%%1)*10,2)
  if(fitlambda==3) {
    fittedval<-estimate.tfe(data=data,MAXT=MAXT)$fittedvals

  }

  if(errorlambda==0){
    arcoef<-matrix(0,nrow=ncol(data$outcome),ncol=2)
    varcov <- diag(rep(0,ncol(data$outcome)))
  }
  if(errorlambda==3){
    arcoef<-cbind(estimate.tfe(data=data,MAXT=MAXT)$arcoef,0)
    varcov<-    estimate.tfe(data=data,MAXT=MAXT)$varcov
  }
  if(includeoutput==TRUE) return(list(fittedval=fittedval,arcoefs=arcoef,varcov=varcov
                                      ))
  ##DRAWING FROM THE FITTED MODEL
  else{
    draw<-draw.ar(r,MAXT,fittedval=fittedval,arcoef=arcoef,varcov=varcov,
                  data.covariates=data$covariates)
    outtable<-data.table()
    for(iter in 1:ncol(draw$outcome)){
      drawdata<-copy(data)
      drawdata$treated<-as.matrix(draw$outcome[,iter,with=FALSE])
      drawdata$donors<-as.matrix(draw$outcome[,-iter,with=FALSE])
      drawdata$outcome<-as.matrix(cbind(drawdata$treated,drawdata$donors))
      drawdata$covariates<-copy(draw$covariates)
      drawdata$covariates[,unit:=unit+1]
      drawdata$covariates[unit==iter+1,unit:=1]
      
      allcovariates<-copy(drawdata$covariates)
      allcovariates<-allcovariates[time%in%data$ADH_fittimes,]
      allcovariates[,time:=time-min(time)+1]
      
      ADHdata<-copy(drawdata)
      ADHdata$outcome<-ADHdata$outcome[-(1:(min(ADHdata$ADH_fittimes)-1)),]
      ADHdata$donors<-ADHdata$donors[-(1:(min(ADHdata$ADH_fittimes)-1)),]
      ADHdata$treated<-as.matrix(ADHdata$treated[-(1:(min(ADHdata$ADH_fittimes)-1)),])
      ADHdata$treatment<-length(ADHdata$ADH_fittimes)+1
      ADHdata$covariates<-ADHdata$covariates[time%in%ADHdata$ADH_fittimes,]
      ADHdata$covariates[,time:=time-min(time)+1]
      ADHdata$firstforecast<-5
      
      ###SOLVING ESTIMATORS:
      output<-list()
      
      maxphi<-1e5
      maxphi.abadie<-finddist(estimator=solve.covreg,
                              data=ADHdata,
                              tune.pars.joint=lapply(maxphi, function(x) list(
                                K=1, phi=x, set.k=NA, min.preperiods=ADHdata$firstforecast,
                                Vfun=Cov.Vars,
                                maincovariates=allcovariates,
                                type="all"
                              )),
                              tune.pars.list=NULL,
                              model='regularize',
                              est.options=NULL,
                              only.synth=FALSE,only.match=FALSE,dump.info=dump.info)
      
      maxphi.abadie.outcome<-finddist(estimator=solve.locreg,
                              data=drawdata,
                              tune.pars.list=list(phi=maxphi,
                                                  K=1,
                                                  set.k=NA,
                                                  min.preperiods=firstforecast),
                              model='regularize',
                              only.synth=FALSE,only.match=FALSE,dump.info=FALSE,fold.errors=TRUE)
      
      output$reg<-solve.synth(data=ADHdata,
                              est.options=list(sigf.ipop=5,
                                               Margin.ipop=5e-04,
                                               Wbar=Wbar),
                              tune.pars=list(phi=0,K=1))
      output$reg$weights<-list(pars=output$reg$pars)
      
    


      
      outtable<-rbind(outtable,data.table(Placebo=data$names[iter],
                                          lambda=data$lambda,
                                          sim=r,
                                          Model="reg",
                                          MSPE4=mean(output$reg$pred.error[1:4]^2),
                                          MSPEAll=mean(output$reg$pred.error^2),
                                          fit=mean(output$reg$fit^2),
                                          cv.error=output$reg$cv.error,
                                          numweights=sum(output$reg$weights$pars>0.01),
                                          Phi=output$reg$tune.pars$phi,
                                          K=output$reg$tune.pars$K,
                                          weights=t(output$reg$weights$pars),
                                          pred.error=t(output$reg$pred.error)
      )
      ,fill=TRUE)
      
      
      
      phis<-  lapply(1:1,function(x) 0)
      output$reg.outcome<-cv.solver(estimator=solve.locreg,data=drawdata,model='regularize',
                            tune.pars.list=list(phi=phis,
                                                K=1,
                                                set.k=NA,
                                                min.preperiods=firstforecast), 
                            only.synth=TRUE,
                            only.match=FALSE,dump.info=TRUE)
      outtable<-rbind(outtable,data.table(Placebo=data$names[iter],
                                          lambda=data$lambda,
                                          sim=r,
                                          Model="reg.outcome",
                                          MSPE4=mean(output$reg.outcome$pred.error[1:4]^2),
                                          MSPEAll=mean(output$reg.outcome$pred.error^2),
                                          fit=mean(output$reg.outcome$fit^2),
                                          cv.error=output$reg.outcome$cv.error,
                                          numweights=sum(output$reg.outcome$weights$pars>0.01),
                                          Phi=output$reg.outcome$tune.pars$phi,
                                          K=output$reg.outcome$tune.pars$K,
                                          weights=t(output$reg.outcome$weights$pars),
                                          pred.error=t(output$reg.outcome$pred.error)
      )
      ,fill=TRUE)
      
      maxphi.abadiescweight<-finddist(estimator=solve.covreg,
                                      data=ADHdata,
                                      tune.pars.joint=lapply(maxphi, function(x) list(
                                        K=1, phi=x, set.k=NA, min.preperiods=ADHdata$firstforecast,
                                        scV=  unlist(output$reg$tune.pars$V),
                                        Vfun=Cov.Vars,
                                        maincovariates=allcovariates,
                                        type="all"
                                      )),
                                      tune.pars.list=NULL,
                                      model='regularize',
                                      est.options=list(sigf.ipop=5,
                                                       Margin.ipop=5e-04,
                                                       Wbar=Wbar),
                                      only.synth=FALSE,only.match=FALSE,dump.info=dump.info)
      
      
      output$match <- cv.solver(estimator=solve.locreg,
                                data=ADHdata,
                                model='analytic',
                                tune.pars.joint=lapply(Kvals,function(x) list(
                                  K=x, phi=1, set.k=NA, min.preperiods=ADHdata$firstforecast,
                                  Vfun=Cov.Vars,
                                  maincovariates=allcovariates,
                                  type="all"
                                )),
                                tune.pars.list=NULL,
                                est.options=list(Wbar=solve.covmatch),
                                only.synth=FALSE,only.match=TRUE,dump.info=dump.info)
      
      outtable<-rbind(outtable,data.table(Placebo=data$names[iter],
                                          lambda=data$lambda,
                                          sim=r,
                                          Model="match",
                                          MSPE4=mean(output$match$pred.error[1:4]^2),
                                          MSPEAll=mean(output$match$pred.error^2),
                                          fit=mean(output$match$fit^2),
                                          cv.error=output$match$cv.error,
                                          numweights=sum(output$match$weights$pars>0.01),
                                          Phi=output$match$tune.pars$phi,
                                          K=output$match$tune.pars$K,
                                          weights=t(output$match$weights$pars),
                                          pred.error=t(output$match$pred.error)
      ),fill=TRUE)
      
      
      phis<-  lapply(1:length(Kvals),function(x) 1)
      output$match.outcome <- cv.solver(estimator=solve.locreg,data=drawdata,
                                model='regularize',
                                tune.pars.list=list(phi=phis,K=Kvals,
                                                    set.k=NA,
                                                    min.preperiods=firstforecast),
                                only.synth=FALSE,only.match=TRUE,dump.info=dump.info)
      
      outtable<-rbind(outtable,data.table(Placebo=data$names[iter],
                                          lambda=data$lambda,
                                          sim=r,
                                          Model="match.outcome",
                                          MSPE4=mean(output$match.outcome$pred.error[1:4]^2),
                                          MSPEAll=mean(output$match.outcome$pred.error^2),
                                          fit=mean(output$match.outcome$fit^2),
                                          cv.error=output$match.outcome$cv.error,
                                          numweights=sum(output$match.outcome$weights$pars>0.01),
                                          Phi=output$match.outcome$tune.pars$phi,
                                          K=output$match.outcome$tune.pars$K,
                                          weights=t(output$match.outcome$weights$pars),
                                          pred.error=t(output$match.outcome$pred.error)
      ),fill=TRUE)
      
      
      
      
      phis<-seq(0,1,length.out=philength)
      output$SCAM.byphi<-cv.estimator.average(estimator=solve.synth,
                                              data=ADHdata,
                                              model='average',
                                              tune.pars=list(
                                                Kvals=Kvals, phis=phis, set.k=NA, min.preperiods=ADHdata$firstforecast,
                                                Vfun=Cov.Vars,
                                                maincovariates=allcovariates,
                                                type="all"
                                              ),
                                              est.options=list(Wbar=solve.covmatch,
                                                               sigf.ipop=5,
                                                               Margin.ipop=5e-04),
                                              only.synth=FALSE,only.match=FALSE,dump.info=dump.info,
                                              forecast.minlength=1,forecast.maxlength=1)
      
      outtable<-rbind(outtable,data.table(Placebo=data$names[iter],
                                          lambda=data$lambda,
                                          sim=r,
                                          Model="MASC",
                                          MSPE4=mean(output$SCAM.byphi$pred.error[1:4]^2),
                                          MSPEAll=mean(output$SCAM.byphi$pred.error^2),
                                          fit=mean(output$SCAM.byphi$fit^2),
                                          cv.error=min(output$SCAM.byphi$cv.error[[output$SCAM.byphi$tune.pars$K]]),
                                          numweights=sum(output$SCAM.byphi$weights$pars>0.01),
                                          Phi=output$SCAM.byphi$tune.pars$phi,
                                          K=output$SCAM.byphi$tune.pars$K,
                                          weights=t(output$SCAM.byphi$weights$pars),
                                          pred.error=t(output$SCAM.byphi$pred.error)
      ),fill=TRUE)
      outtable<-rbind(outtable,data.table(Placebo=data$names[iter],
                                          lambda=data$lambda,
                                          sim=r,
                                          Model="MASC.byphi",
                                          MSPE4=output$SCAM.byphi$all.results$pred.error.4,
                                          MSPEAll=output$SCAM.byphi$all.results$pred.error.all,
                                          fit=output$SCAM.byphi$all.results$fit,
                                          cv.error=output$SCAM.byphi$all.results$cv.error,
                                          Phi=rep(phis[[1]],length(Kvals)),
                                          K=rep(Kvals,each=length(phis[[1]]))
      ),fill=TRUE)
      
      
      output$SCAM.average.outcome<-cv.solver(estimator=solve.locreg,data=drawdata,
                                     tune.pars.list=list(phi=phis,K=Kvals,
                                                         set.k=NA,
                                                         min.preperiods=firstforecast),
                                     tune.pars.joint=NULL,
                                     only.synth=FALSE,only.match=FALSE,
                                     dump.info=dump.info,model='analytic',fold.errors=TRUE)
      outtable<-rbind(outtable,data.table(Placebo=data$names[iter],
                                          lambda=data$lambda,
                                          sim=r,
                                          Model="MASC.outcome",
                                          MSPE4=mean(output$SCAM.average.outcome$pred.error[1:4]^2),
                                          MSPEAll=mean(output$SCAM.average.outcome$pred.error^2),
                                          fit=mean(output$SCAM.average.outcome$fit^2),
                                          cv.error=output$SCAM.average.outcome$cv.error,
                                          numweights=sum(output$SCAM.average.outcome$weights$pars>0.01),
                                          Phi=output$SCAM.average.outcome$tune.pars$phi,
                                          K=output$SCAM.average.outcome$tune.pars$K,
                                          weights=t(output$SCAM.average.outcome$weights$pars),
                                          pred.error=t(output$SCAM.average.outcome$pred.error)
      ),fill=TRUE)
      
      phis<-  lapply(1:length(Kvals),function(x) seq(0,1,length.out=philength))
      
      output$SCAM.byphi.outcome<-cv.solver(estimator=solve.locreg,data=drawdata,
                                   tune.pars.list=list(phi=phis,K=Kvals,
                                                       set.k=NA,
                                                       min.preperiods=firstforecast),
                                   tune.pars.joint=NULL,
                                   only.synth=FALSE,only.match=FALSE,
                                   dump.info=dump.info,model='average',fold.errors=TRUE)
      
      outtable<-rbind(outtable,data.table(Placebo=data$names[iter],
                                          lambda=data$lambda,
                                          sim=r,
                                          Model="MASC.byphi.outcome",
                                          MSPE4=output$SCAM.byphi.outcome$all.results$pred.error.4,
                                          MSPEAll=output$SCAM.byphi.outcome$all.results$pred.error.all,
                                          fit=output$SCAM.byphi.outcome$all.results$fit,
                                          cv.error=output$SCAM.byphi.outcome$all.results$cv.error,
                                          Phi=rep(phis[[1]],length(Kvals)),
                                          K=rep(Kvals,each=length(phis[[1]]))
      ),fill=TRUE)
      
      
      phis<-  sapply(1:length(maxphi.abadie[[1]]),function(x) seq(0,maxphi.abadie[[x]],length.out=philength))
      output$SCAM.abadie<-cv.solver(estimator=solve.covreg,
                                    data=ADHdata,
                                    model='regularize',
                                    tune.pars.joint=lapply(phis,function(x) list(
                                      K=1, phi=x, set.k=NA, min.preperiods=ADHdata$firstforecast,
                                      Vfun=Cov.Vars,
                                      maincovariates=allcovariates,
                                      type="all"
                                    )),
                                    tune.pars.list=NULL,
                                    est.options=NULL,
                                    only.synth=FALSE,only.match=FALSE,dump.info=dump.info)
      
      outtable<-rbind(outtable,data.table(Placebo=data$names[iter],
                                          lambda=data$lambda,
                                          sim=r,
                                          Model="Penalized SC",
                                          MSPE4=mean(output$SCAM.abadie$pred.error[1:4]^2),
                                          MSPEAll=mean(output$SCAM.abadie$pred.error^2),
                                          fit=mean(output$SCAM.abadie$fit^2),
                                          cv.error=output$SCAM.abadie$cv.error,
                                          numweights=sum(output$SCAM.abadie$weights$pars>0.01),
                                          Phi=output$SCAM.abadie$tune.pars$phi,
                                          K=output$SCAM.abadie$tune.pars$K,
                                          weights=t(output$SCAM.abadie$weights$pars),
                                          pred.error=t(output$SCAM.abadie$pred.error)
      ),fill=TRUE)
      
      outtable<-rbind(outtable,data.table(Placebo=data$names[iter],
                                          lambda=data$lambda,
                                          sim=r,
                                          Model="Penalized SC.byphi",
                                          MSPE4=output$SCAM.abadie$all.results$pred.error.4,
                                          MSPEAll=output$SCAM.abadie$all.results$pred.error.all,
                                          fit=output$SCAM.abadie$all.results$fit,
                                          cv.error=output$SCAM.abadie$all.results$cv.error,
                                          Phi=seq(from=0,to=1,length.out=length(phis[[1]])),
                                          K=output$SCAM.abadie$tune.pars$K
      ),fill=TRUE)
      
      phis<-  lapply(1:length(maxphi.abadie.outcome[[1]]),function(x) seq(0,maxphi.abadie.outcome[[x]],length.out=philength))
      output$SCAM.abadie.outcome<-cv.solver(estimator=solve.locreg,data=drawdata,
                                    tune.pars.list=list(phi=phis,
                                                        K=1,
                                                        set.k=NA,
                                                        min.preperiods=firstforecast),
                                    tune.pars.joint=NULL,
                                    only.synth=FALSE,only.match=FALSE,
                                    dump.info=dump.info,model='regularize',fold.errors=TRUE)
      
      outtable<-rbind(outtable,data.table(Placebo=data$names[iter],
                                          lambda=data$lambda,
                                          sim=r,
                                          Model="Penalized SC.outcome",
                                          MSPE4=mean(output$SCAM.abadie.outcome$pred.error[1:4]^2),
                                          MSPEAll=mean(output$SCAM.abadie.outcome$pred.error^2),
                                          fit=mean(output$SCAM.abadie.outcome$fit^2),
                                          cv.error=output$SCAM.abadie.outcome$cv.error,
                                          numweights=sum(output$SCAM.abadie.outcome$weights$pars>0.01),
                                          Phi=output$SCAM.abadie.outcome$tune.pars$phi,
                                          K=output$SCAM.abadie.outcome$tune.pars$K,
                                          weights=t(output$SCAM.abadie.outcome$weights$pars),
                                          pred.error=t(output$SCAM.abadie.outcome$pred.error)
      ),fill=TRUE)
      
      outtable<-rbind(outtable,data.table(Placebo=data$names[iter],
                                          lambda=data$lambda,
                                          sim=r,
                                          Model="Penalized SC.byphi.outcome",
                                          MSPE4=output$SCAM.abadie.outcome$all.results$pred.error.4,
                                          MSPEAll=output$SCAM.abadie.outcome$all.results$pred.error.all,
                                          fit=output$SCAM.abadie.outcome$all.results$fit,
                                          cv.error=output$SCAM.abadie.outcome$all.results$cv.error,
                                          Phi=seq(from=0,to=1,length.out=length(phis[[1]])),
                                          K=output$SCAM.abadie.outcome$tune.pars$K
      ),fill=TRUE)
      
      phis<-  sapply(1:length(maxphi.abadiescweight[[1]]),function(x) seq(0,maxphi.abadiescweight[x],length.out=philength))
      output$SCAM.abadieSCweight<- cv.solver(estimator=solve.covreg,
                                             data=ADHdata,
                                             model='regularize',
                                             tune.pars.joint=lapply(phis,function(x) list(
                                               K=1, phi=x, set.k=NA, min.preperiods=ADHdata$firstforecast,
                                               scV=  unlist(output$reg$tune.pars$V),
                                               Vfun=Cov.Vars,
                                               maincovariates=allcovariates,
                                               type="all"
                                             )),
                                             tune.pars.list=NULL,
                                             # est.options=list(sigf.ipop=5,
                                             #                  Margin.ipop=5e-04,
                                             #                  Wbar=Wbar,
                                             #                  solver="ipop"),
                                             only.synth=FALSE,only.match=FALSE,dump.info=dump.info)
      
      
      outtable<-rbind(outtable,data.table(Placebo=data$names[iter],
                                          lambda=data$lambda,
                                          sim=r,
                                          Model="Penalized SC (SC weights)",
                                          MSPE4=mean(output$SCAM.abadieSCweight$pred.error[1:4]^2),
                                          MSPEAll=mean(output$SCAM.abadieSCweight$pred.error^2),
                                          fit=mean(output$SCAM.abadieSCweight$fit^2),
                                          cv.error=output$SCAM.abadieSCweight$cv.error,
                                          numweights=sum(output$SCAM.abadieSCweight$weights$pars>0.01),
                                          Phi=output$SCAM.abadieSCweight$tune.pars$phi,
                                          K=output$SCAM.abadieSCweight$tune.pars$K,
                                          weights=t(output$SCAM.abadieSCweight$weights$pars),
                                          pred.error=t(output$SCAM.abadieSCweight$pred.error)
      ),fill=TRUE)
      
      outtable<-rbind(outtable,data.table(Placebo=data$names[iter],
                                          lambda=data$lambda,
                                          sim=r,
                                          Model="Penalized SC (SC weights).byphi",
                                          MSPE4=output$SCAM.abadieSCweight$all.results$pred.error.4,
                                          MSPEAll=output$SCAM.abadieSCweight$all.results$pred.error.all,
                                          fit=output$SCAM.abadieSCweight$all.results$fit,
                                          cv.error=output$SCAM.abadieSCweight$all.results$cv.error,
                                          Phi=seq(from=0,to=1,length.out=length(phis[[1]])),
                                          K=output$SCAM.abadieSCweight$tune.pars$K
      ),fill=TRUE)
      
    }
    if(includedraw==TRUE){
      outtable<-list(outtable=outtable,drawdata=draw)
    }
    return(outtable)
  }
}

