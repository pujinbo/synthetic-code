
distcheck<-function(old,new){
  maindiff<-lapply(1:length(new$weights), function(x) sum((new$weights[[x]]$pars-old$weights[[x]]$pars)^2))
  cvdiffs<-unlist(lapply(1:length(new$weights), function(x) lapply(1:length(new$cvweights[[x]]), function(y) 
                sum((new$cvweights[[x]][[y]]$pars-old$cvweights[[x]][[y]]$pars)^2)
                                                            )))
  return(maindiff<=1e-4 & all(cvdiffs<=1e-4))
}
finddist<-function(estimator,data,tune.pars.list,tune.pars.joint=NULL,
                   model,only.synth,only.match,dump.info,fold.errors=FALSE,est.options=NULL){
  phivals<-NULL
  if(is.null(tune.pars.list)) tuneks <- tune.pars.joint
  else tuneks<-tune.pars.list$K
  
  for(k in 1:length(tuneks)){
    
    tunelistk<-tune.pars.list
    tunejointk<-tune.pars.joint
    
    if(is.null(tune.pars.list)) {
      tunejointk<-tunejointk[k]
      tunephi <- tunejointk[[1]]$phi
    }
    else {
      tunelistk$K<-tune.pars.list$K[k]
      tunephi <- tune.pars.list$phi
    } 
    
  weightsold<-cv.solver(estimator=estimator,data=data,model=model,
                        tune.pars.list=tunelistk,
                        tune.pars.joint=tunejointk,
                        only.synth=only.synth,
                        only.match=TRUE,dump.info=dump.info,fullweights=TRUE,
                        fold.errors=fold.errors,
                        est.options=est.options)$all.results[c('weights','cv.weights')]

  weightsnew<-cv.solver(estimator=estimator,data=data,model=model,
                        tune.pars.list=tunelistk,
                        tune.pars.joint=tunejointk,
                        only.synth=only.synth,
                        only.match=only.match,dump.info=dump.info,fullweights=TRUE,
                        fold.errors=fold.errors,
                        est.options=est.options)$all.results[c('weights','cv.weights')]

  while(!distcheck(weightsold,weightsnew)){
    print('WARNING, DISTANCE PENALTY WEIGHTS IMMEDIATELY CHANGE FROM INITIAL GUESS OF MAXIMAL PHI')
    print('INCREASING INITIAL GUESS')
    tunephi<-tunephi*10
    if(is.null(tune.pars.list)) tunejointk[[1]]$phi<-tunephi
    else tunelistk$phi<-tunephi
    
    
    weightsnew<-cv.solver(estimator=estimator,data=data,model=model,
                          tune.pars.list=tunelistk,
                          tune.pars.joint=tunejointk,
                          only.synth=only.synth,
                          only.match=only.match,dump.info=dump.info,fullweights=TRUE,
                          fold.errors=fold.errors,
                          est.options=est.options)$all.results[c('weights','cv.weights')]

  }
  while(distcheck(weightsnew,weightsold)&round(tunephi,5)>0){
    tunephi<-tunephi/5
    if(is.null(tune.pars.list)) tunejointk[[1]]$phi<-tunephi
    else tunelistk$phi<-tunephi
    weightsold<-weightsnew
    weightsnew<-cv.solver(estimator=estimator,data=data,model=model,
                          tune.pars.list=tunelistk,
                          tune.pars.joint=tunejointk,
                          only.synth=only.synth,
                          only.match=only.match,dump.info=dump.info,fullweights=TRUE,
                          fold.errors=fold.errors,
                          est.options=est.options)$all.results[c('weights','cv.weights')]
  }
  rm(weightsnew)
  rm(weightsold)
  phivals[k] <- tunephi*5

}
  return(phivals)
}

fold.forecast<-function(data,weights,fold.end,forecast.minlength,forecast.maxlength){
  forecast.end<- min((fold.end+forecast.maxlength),data$treatment-1)
  cv.interval<- (fold.end+forecast.minlength):(forecast.end)
  if(max(cv.interval)>=data$treatment) stop("One of the cross-validation folds is trying to forecast into treatment period")
  donors<-data$donors[cv.interval,]
  treated<-data$treated[cv.interval,]
  return(mean((treated-donors%*%weights)^2))
}

fold.forecast.average<-function(data, weights.SC, weights.match, phi, fold.end, forecast.minlength, forecast.maxlength){
  
  forecast.end<- min((fold.end+forecast.maxlength),data$treatment-1)
  cv.interval<- (fold.end+forecast.minlength):(forecast.end)
  if(max(cv.interval)>=data$treatment) stop("One of the cross-validation folds is trying to forecast into treatment period")
  donors<-data$donors[cv.interval,]
  treated<-data$treated[cv.interval,]
  return(mean((treated-donors%*%(weights.SC*(1-phi)+weights.match*phi))^2))
  
}


##This function takes as its argument data and a list of estimators, and returns cross-validation errors for those estimators.
cv.estimator<-function(estimator,data,tune.pars,only.match,only.synth,model,dump.info,fullweights=FALSE,
                       forecast.minlength, forecast.maxlength,foldtype,
                       est.options){
  min.preperiods=tune.pars$min.preperiods
  set.k=tune.pars$set.k
  
if((is.na(min.preperiods)&is.na(set.k))|(!is.na(min.preperiods)&!is.na(set.k)))  stop("You must specify precisely one of min.preperiods and set.k")
  if(!is.na(min.preperiods)) {
    if(min.preperiods+1+forecast.minlength > data$treatment) stop("No valid cross-validation folds: forecast.minlength is too high")
    set.k<-min.preperiods:(data$treatment-1-forecast.minlength)
  }
  
  if(model=="analytic"){
      Y_treated<-NULL
      Y_sc<-NULL
      Y_match<-NULL
      objweight<-NULL
      
      tune.pars$phi<-0
      SCfolds.loss.v<-NULL
      SCfolds.loss.w<-NULL
      for(k in set.k){
        if(foldtype != "leaveout"){
          treat<-k+1 
          keepobs<-1:nrow(data$donors)
        }
        if(foldtype=="leaveout"){
          treat<-data$treatment-1
          keepobs<-(1:nrow(data$donors))[-(k+1)]
        }
        folddata<-data
        folddata$treatment=treat
        folddata$donors<-data$donors[keepobs,]
        folddata$treated<-as.matrix(data$treated[keepobs,])
        treatinterval<-(k+forecast.minlength):min(k+forecast.maxlength,data$treatment-1)
      weight.match<-list(pars=est.options$Wbar(data=folddata,
                                               tune.pars=tune.pars))
      weight.SC<-estimator(data=folddata,
                           tune.pars=tune.pars,model="average",
                           only.match=FALSE,only.synth=TRUE, est.options=est.options)
      Y_treated<-c(Y_treated,data$treated[treatinterval,])
      Y_sc <- c(Y_sc,
                data$donors[treatinterval,]%*%weight.SC$pars)
      Y_match <- c(Y_match,
                data$donors[treatinterval,]%*%weight.match$pars)
      objweight<-c(objweight,rep(1/length(treatinterval),length(treatinterval)))
      SCfolds.loss.v<-c(SCfolds.loss.v, weight.SC$loss.v)
      SCfolds.loss.w<-c(SCfolds.loss.w, weight.SC$loss.w)
      }
      phi<-(objweight*(Y_treated-Y_sc))%*%(Y_match-Y_sc)/((objweight*(Y_match-Y_sc))%*%(Y_match-Y_sc))
      phi<-max(0,phi)
      phi<-min(1,phi)
      if(only.match==TRUE) phi<-1
      if(only.synth==TRUE) phi<-0
      
      cv.error<- sum(objweight*(Y_treated-phi*Y_match-(1-phi)*Y_sc)^2)/sum(objweight)
      finalweight.match<-list(pars=est.options$Wbar(data=data,
                              tune.pars=tune.pars))
      finalweight.SC<-estimator(data=data,
                           tune.pars=tune.pars,model="average",
                           only.match=FALSE,only.synth=TRUE,
                           est.options=est.options)

      output<-list(weights = list(pars=finalweight.match$pars*phi + finalweight.SC$pars*(1-phi)))
      output$pred.error<-data$treated[data$treatment:nrow(data$donors),] - data$donors[data$treatment:nrow(data$donors),]%*%output$weights$pars
      output$fit<-data$treated[1:(data$treatment-1)] - data$donors[1:(data$treatment-1),]%*%output$weights$pars
      output$tune.pars<-list(phi=phi,K=tune.pars$K)
      output$cv.error<-cv.error
      output$SC.loss.v<-weight.SC$loss.v
      output$SC.loss.w<-weight.SC$loss.w
      output$SCfolds.loss.v<-SCfolds.loss.v
      output$SCfolds.loss.w<-SCfolds.loss.w
      
      
  }
  else{
  treat<-data$treatment 
  keepobs<-1:nrow(data$donors)
  set.k<-set.k[which(data$treatment-set.k-1>=1)]
  output<-list(weights=estimator(data=data,model=model,
                                 tune.pars=tune.pars, only.match=only.match, only.synth=only.synth,
                                 est.options=est.options))
  output$pred.error<-data$treated[data$treatment:nrow(data$donors),] - data$donors[data$treatment:nrow(data$donors),]%*%output$weights$pars
  output$fit<-data$treated[1:(data$treatment-1)] - data$donors[1:(data$treatment-1),]%*%output$weights$pars
  output$tune.pars<-tune.pars
  output$cv.error<-NA
  output$fold.errors<-list(NA)
  output$fold.weights<-list(NA)
  position=1
  for(k in set.k){
    if(foldtype != "leaveout"){
     treat<-k+1 
     keepobs<-1:nrow(data$donors)
    }
    if(foldtype=="leaveout"){
      treat<-data$treatment-1
      keepobs<-(1:nrow(data$donors))[-(k+1)]
    }
    folddata<-copy(data)
    folddata$treatment<-treat
    folddata$donors<-data$donors[keepobs,]
    folddata$treated<-as.matrix(data$treated[keepobs,])
    output$fold.weights[[position]]<-estimator(data=folddata,
                                               tune.pars=tune.pars,model=model,
                                               only.match=only.match,only.synth=only.synth,
                                               est.options=est.options)
    output$fold.errors[[position]]<-fold.forecast(data=data,
                                                  weights=output$fold.weights[[position]]$pars,
                                                  fold.end=k,
                                                  forecast.minlength = forecast.minlength,forecast.maxlength = forecast.maxlength)
    output$SCfolds.loss.v<-c(output$SCfolds.loss.v, output$fold.weights[[position]]$loss.v)
    output$SCfolds.loss.w<-c(output$SCfolds.loss.w, output$fold.weights[[position]]$loss.w)
    
    
    position=position+1
  }
  output$SC.loss.v<-output$weights$loss.v
  output$SC.loss.w<-output$weights$loss.w
  
  
  names(output$fold.errors)<-paste0("fold ",set.k)
  names(output$fold.weights)<-paste0("fold ",set.k)
  output$cv.error <- mean(unlist(output$fold.errors)) 
  if(dump.info==FALSE) output$cv.specs<-list(set.k=set.k)
  if(dump.info==TRUE){
    output$fold.weights<-NULL
    output$tune.pars<-list(phi=output$tune.pars$phi,K=output$tune.pars$K)
  }
  }
  return(output)
  
}


cv.solver <- function(estimator,data,foldtype="rolling",tune.pars.list=list(phi=(0:100)/100,K=1:5),tune.pars.joint=NULL,
                      dump.info=FALSE,only.match=FALSE,only.synth=FALSE,model='regularize',fullweights=FALSE,fold.errors=FALSE,
                      forecast.minlength=1, forecast.maxlength = 1, est.options=list(Wbar=Wbar)){
  if((!is.null(tune.pars.list)&!is.null(tune.pars.joint))|(is.null(tune.pars.list)&is.null(tune.pars.joint))) stop("You must specify precisely one of min.preperiods and set.k")
  if(!is.null(tune.pars.list)) {
    tune.pars.joint<-list(NA)
    position=1
    if(model=="analytic") tune.pars.list$phi<-lapply(1:length(tune.pars.list$K), function(x) 0)
    for(d in 1:length(tune.pars.list$K)){
      phivals<-tune.pars.list$phi[[d]]
      for(a in 1:length(phivals)){
          for(g in 1:length(tune.pars.list$min.preperiods)){
            for(j in 1:length(tune.pars.list$set.k)){
              tune.pars.joint[[position]]<-list(phi=phivals[a],
                                                K=tune.pars.list$K[d],
                                                min.preperiods=tune.pars.list$min.preperiods[g],
                                                set.k=tune.pars.list$set.k[[j]],
                                                unit=tune.pars.list$unit
              )
              position=position+1
              
            }
          }
      }
    }
  }
  allresults<-lapply(tune.pars.joint,function(x) cv.estimator(estimator=estimator,data=data,foldtype=foldtype, tune.pars=x,
                                                              only.match=only.match,only.synth=only.synth,model=model,
                                                              dump.info=dump.info,fullweights=fullweights, 
                                                              forecast.minlength=forecast.minlength, forecast.maxlength=forecast.maxlength,
                                                              est.options=est.options))
  #allresults is a list of unnamed things, each with named list (weights, pred error, fold stuff, tune pars, cv errors)
  output<-allresults[[which.min(lapply(allresults,function(x) x$cv.error))]]
  output$all.results<-list()
  output$all.results$pred.error.4<-sapply(allresults,function(x) mean((x$pred.error[1:4])^2))
  output$all.results$pred.error.all<-sapply(allresults,function(x) mean((x$pred.error)^2))
  output$all.results$cv.error<-sapply(allresults,function(x) x$cv.error)
  output$all.results$fit<-sapply(allresults,function(x) mean(x$fit^2))
  output$all.results$pred.error.4.level<-sapply(allresults,function(x) mean((x$pred.error[1:4])))
  output$all.results$pred.error.4.abslevel<-sapply(allresults,function(x) mean(abs(x$pred.error[1:4])))
  output$all.results$SC.loss.v<-lapply(allresults, function(x) x$SC.loss.v)
  output$all.results$SC.loss.w<-lapply(allresults, function(x) x$SC.loss.w)
  output$all.results$SCfolds.loss.v<-lapply(allresults, function(x) x$SCfolds.loss.v)
  output$all.results$SCfolds.loss.w<-lapply(allresults, function(x) x$SCfolds.loss.w)

    if(fold.errors==TRUE) output$all.results$fold.errors<-sapply(allresults,function(x) x$fold.errors)

  if(fullweights==TRUE){
  output$all.results$cv.weights<-lapply(allresults,function(x) x$fold.weights)
  output$all.results$weights<-lapply(allresults,function(x) x$weights)
  }
  return(output)
}










##This function takes as its argument data and two estimators, returning the best convex combination of those two estimators over MSPE.
cv.estimator.average<-function(estimator,data,tune.pars,only.match,only.synth,model,dump.info,fullweights=FALSE,
                       forecast.minlength, forecast.maxlength,foldtype="rolling",
                       est.options){
  min.preperiods=tune.pars$min.preperiods
  set.k=tune.pars$set.k
  phis<-tune.pars$phis
  
  if((is.na(min.preperiods)&is.na(set.k))|(!is.na(min.preperiods)&!is.na(set.k)))  stop("You must specify precisely one of min.preperiods and set.k")
  if(!is.na(min.preperiods)) {
    if(min.preperiods+1+forecast.minlength > data$treatment) stop("No valid cross-validation folds: forecast.minlength is too high")
    set.k<-min.preperiods:(data$treatment-1-forecast.minlength)
  }
  
    treat<-data$treatment 
    keepobs<-1:nrow(data$donors)
    set.k<-set.k[which(data$treatment-set.k-1>=1)]
    output<-list(weights.SC=estimator(data=data,model=model,
                                   tune.pars=c(tune.pars,phi=0,K=1), only.match=only.match, only.synth=only.synth,
                                   est.options=est.options)$pars,
                 weights.match=list(NA))
    NNpos<-1
    for(NN in tune.pars$Kvals){
      tune.parsN<-tune.pars
      tune.parsN$K<-NN
    output$weights.match[[NNpos]] <-  est.options$Wbar(data=data,tune.pars=tune.parsN)
    NNpos<-NNpos+1
    }
    output$cv.error<-list(NA)
    output$fold.errors<-list(NA)
    output$fold.weights.SC<-list(NA)
    output$fold.weights.match<-list(NA)
    output$all.results<-list()
    output$all.results$pred.error.4<-NA
    for(NN in tune.pars$Kvals){
      output$fold.weights.match[[NN]]<-list(NA)
      output$fold.errors[[NN]]<-list(NA)
    }
    position=1
    output$phi.byk<-NULL
    for(k in set.k){
      if(foldtype != "leaveout"){
        treat<-k+1 
        keepobs<-1:nrow(data$donors)
      }
      if(foldtype=="leaveout"){
        treat<-data$treatment-1
        keepobs<-(1:nrow(data$donors))[-(k+1)]
      }
      folddata<-copy(data)
      folddata$treatment<-treat
      folddata$donors<-data$donors[keepobs,]
      folddata$treated<-as.matrix(data$treated[keepobs,])
      output$fold.weights.SC[[position]]<-estimator(data=folddata,
                                                 tune.pars=c(tune.pars,phi=0,K=1),
                                                 model=model,
                                                 only.match=only.match,only.synth=only.synth,
                                                 est.options=est.options)$modelweights$SC
      for(NN in tune.pars$Kvals){
        tune.parsN<-tune.pars
        tune.parsN$K<-NN
        output$fold.weights.match[[NN]][[position]]<-est.options$Wbar(data=folddata,
                                                                      tune.pars=tune.parsN)
      }

      position=position+1
    }

    allpos<-1
    NNpos<-1
    for(NN in tune.pars$Kvals){
    position=1
      for(phi in phis){
        output$fold.errors[[NNpos]][[position]]<-NA
        kpos=1
        for(k in set.k){
        output$fold.errors[[NNpos]][[position]]<-c(output$fold.errors[[NNpos]][[position]],fold.forecast.average(data=data,
                                                              weights.SC=output$fold.weights.SC[[kpos]],
                                                              weights.match=output$fold.weights.match[[NNpos]][[kpos]],
                                                              phi=phi,
                                                              fold.end=k,
                                                              forecast.minlength = forecast.minlength,forecast.maxlength = forecast.maxlength)
        )
        kpos=kpos+1
        }
        position<-position+1
        #Prediction errors for this value:
        weights<-output$weights.SC*(1-phi)+output$weights.match[[NNpos]]*phi
        output$all.results$pred.error.4[allpos]<-mean((data$treated[data$treatment:(data$treatment+3),] - data$donors[data$treatment:(data$treatment+3),]%*%weights)^2)
        output$all.results$pred.error.all[allpos]<-mean((data$treated[data$treatment:nrow(data$treated),] - data$donors[data$treatment:nrow(data$treated),]%*%weights)^2)
        output$all.results$fit[allpos]<-mean((data$treated[1:(data$treatment-1)] - data$donors[1:(data$treatment-1),]%*%weights)^2)
        output$all.results$pred.error.4.level[allpos]<-mean(data$treated[data$treatment:(data$treatment+3),] - data$donors[data$treatment:(data$treatment+3),]%*%weights)
        output$all.results$pred.error.4.abslevel[allpos]<-mean(abs(data$treated[data$treatment:(data$treatment+3),] - data$donors[data$treatment:(data$treatment+3),]%*%weights))
        
        allpos<-allpos+1
      }
    output$cv.error[[NNpos]]<-sapply(output$fold.errors[[NNpos]], function(x) mean(x,na.rm=TRUE)) 
    output$phi.byk<-c(output$phi.byk,phis[output$cv.error[[NNpos]]])
    
        
    NNpos<-NNpos+1
    }
    K<-tune.pars$Kvals[which.min(sapply(output$cv.error,function(x) min(x)))]
    phi<-phis[which.min(output$cv.error[[K]])]
    output$all.results$cv.error<-unlist(output$cv.error)
    output$weights<-list(pars=output$weights.SC*(1-phi)+output$weights.match[[K]]*phi)
    output$pred.error<-data$treated[data$treatment:nrow(data$donors),] - data$donors[data$treatment:nrow(data$donors),]%*%output$weights$pars
    output$fit<-data$treated[1:(data$treatment-1)] - data$donors[1:(data$treatment-1),]%*%output$weights$pars

    output$tune.pars<-list(phi=phi,K=K)
    
    if(dump.info==FALSE) output$cv.specs<-list(set.k=set.k)
    if(dump.info==TRUE){
      output$fold.weights<-NULL
    }
  
  return(output)
  
}

