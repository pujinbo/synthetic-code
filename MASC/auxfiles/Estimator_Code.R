source("gurobiSC.R")
library(Matrix)
##################################
#OPTIMIZATION PROBLEM
#################################

#Nearest Neighbor weight vector
Wbar<-function(data,tune.pars,...){
  dist <- apply((data$donors[1:(data$treatment-1),]-as.vector(data$treated)[1:(data$treatment-1)])^2,MARGIN=2,FUN=sum)
  W <-rep(0,length(dist))
  W[which(dist%in%sort(dist)[1:tune.pars$K])]<-1/length(which(dist%in%sort(dist)[1:tune.pars$K]))
  return(W)
}

#Function determining weights on covariates by their standard deviation
Cov.Vars<-function(input,covariates,...){

  if(input$type=="controls") output<-1/apply(covariates[,-1], 1, var)
  if(input$type=="all") output<- 1/apply(covariates, 1, var)
  if(!is.null(input$scV)) output <- output*input$scV
  return(output)
}
#Function determining weights on covariates, simply feeding in the weights from SC
Cov.Synth<-function(input,...){
  return(input$V)
}

#Function finding the K nearest neighbor(s), based on covariates
solve.covmatch<-function(data,tune.pars,...){
  treatment<-data$treatment
  covariates<-data$covariates[time < treatment,]
  covariates<-covariates[,lapply(.SD,mean,na.rm=TRUE),by=unit,.SDcols=names(covariates)[!names(covariates)%in%c("unit","time")]]
  setorder(covariates,unit)
  covariates[,unit:=NULL]
  covariates<-t(covariates)
  lowest  <- which(rownames(covariates)=="school.illit")
  highest <- which(rownames(covariates)=="school.high")
  covariates[lowest:highest,] <-  
    100 * scale(covariates[lowest:highest,],
                center=FALSE,
                scale=colSums(covariates[lowest:highest,])
    )

  dist <- apply(((covariates[,1]-covariates[,-1])^2),MARGIN=2,
                FUN=function(x) sum(x*(tune.pars$Vfun(c(tune.pars,treatment=treatment),covariates=covariates))))
  W <-rep(0,length(dist))
  W[which(dist%in%sort(dist)[1:tune.pars$K])]<-1/length(which(dist%in%sort(dist)[1:tune.pars$K]))
  return(pars=W)
}


#MASC estimator, using covariates in the same way as Abadie et al. (2010):
solve.synth <- function(data, est.options, tune.pars,...){
  treatment<-data$treatment
  covariates<-data$covariates[time < treatment,]
  covariates<-covariates[,lapply(.SD,mean,na.rm=TRUE),by=unit,.SDcols=names(covariates)[!names(covariates)%in%c("unit","time")]]
  setorder(covariates,unit)
  covariates[,unit:=NULL]
  covariates<-t(covariates)
  
  lowest  <- which(rownames(covariates)=="school.illit")
  highest <- which(rownames(covariates)=="school.high")
  covariates[lowest:highest,] <-  
    100 * scale(covariates[lowest:highest,],
                center=FALSE,
                scale=colSums(covariates[lowest:highest,])
    )
  

ADH<-gurobisynth(Z1=as.matrix(data$treated[1:(treatment-1),]),
                              Z0=data$donors[1:(treatment-1),],
                              X1=as.matrix(covariates[,1]),
                              X0=covariates[,-1],
                              sigf.ipop=est.options$sigf.ipop,
                              Margin.ipop=est.options$Margin.ipop,
                              custom.v=est.options$custom.v
                )

Match<-est.options$Wbar(data=data,tune.pars=tune.pars)                

                output<-list(pars=as.vector(ADH$solution.w)*(1-tune.pars$phi)+(Match)*tune.pars$phi)
                output$pred.error<-data$treated[data$treatment:nrow(data$donors),] - data$donors[data$treatment:nrow(data$donors),]%*%output$pars
                output$fit<-data$treated[1:(data$treatment-1)] - data$donors[1:(data$treatment-1),]%*%output$pars
                output$tune.pars<-list(phi=tune.pars$phi,K=tune.pars$K,V=ADH$solution.v)
                output$loss.v<-ADH$loss.v
                output$loss.w<-ADH$loss.w
                output$cv.error<-NA
                output$modelweights<-list(SC=as.vector(ADH$solution.w),Match=Match)
                
  return(output)
}

#Purely outcome-based MASC estimator
solve.covreg <- function(data, model, tune.pars, est.options, only.match, only.synth,...){
  treatment<-data$treatment
  covariates<-data$covariates[time < treatment,]
  covariates<-covariates[,lapply(.SD,mean,na.rm=TRUE),by=unit,.SDcols=names(covariates)[!names(covariates)%in%c("unit","time")]]
  setorder(covariates,unit)
  covariates[,unit:=NULL]
  covariates<-t(covariates)
  
  lowest  <- which(rownames(covariates)=="school.illit")
  highest <- which(rownames(covariates)=="school.high")
  covariates[lowest:highest,] <-  
    100 * scale(covariates[lowest:highest,],
                center=FALSE,
                scale=colSums(covariates[lowest:highest,])
    )
  covariates<-covariates*(tune.pars$Vfun(c(tune.pars,treatment=treatment),covariates=covariates)^0.5)
  
  
  output<-solve.locreg(data=list(treatment=nrow(covariates)+1,
                                 donors=covariates[,-1],
                                 treated=as.matrix(covariates[,1])
                                 ),
                       model=model,
                       tune.pars=tune.pars, only.match=only.match, only.synth=only.synth,
                       est.options=est.options,
                       checkhull=FALSE)
  output$loss.w<-output$objval

  
  output$tune.pars<-tune.pars

  return(output)
}

#Outcome-based Synthetic controls, local weights
solve.locreg <- function(data,tune.pars,only.synth=FALSE,only.match=FALSE,log=0,psdtol=1e-4,
                         BarConvTol=1e-8,BarIterLimit=1e5,model='regularize',
                         checkhull=TRUE,est.options=list(NULL),...){
  if(only.synth!=FALSE & only.match!=FALSE) stop("only.synth and only.match must take on logical values; only one may be TRUE.")
  ####Pulling objects out of lists###
  treatment<-data$treatment
  donors<-data$donors[1:(treatment-1),]
  treated<-as.matrix(data$treated[1:(treatment-1),])
  phi<-tune.pars$phi
  K <- tune.pars$K
  
  ###########################################      
  obj.synthetic<-function(donors,treated,treatment,only.match){
    if(only.match) return(0)
    else return(-2*t(treated)%*%donors)
  }
  objcon.synthetic<-function(donors,treated,treatment,only.match){
    if(only.match) return(0)
    else return(t(treated)%*%treated)
  }
  Q.synthetic<-function(donors,treated,treatment,only.match){
    if(only.match) return(0)
    else return(t(donors)%*%donors)
  }
  
  
  if(model=='regularize'){
    obj.matching<-function(donors,treated,treatment,only.synth){
      if(only.synth) return(0)
      else {
        return(t(apply((as.vector(treated)-donors)^2,MARGIN=2,FUN=sum)))
      }
    }
    objcon.matching<-function(donors,treated,treatment,only.synth){
      return(0)
    } 
    Q.matching<-function(donors,treated,treatment,only.synth){
      return(matrix(0,ncol(donors),ncol(donors)))
    } 
  }
  
  if(model=='average'){
    obj.matching<-function(data,tune.pars,only.synth){
      if(only.synth) return(rep(0,ncol(donors)))
      else {
        return((-2*Wbar(data=data,tune.pars=tune.pars)))
      }
    }
    objcon.matching<-function(data,tune.pars,only.synth){
      if(only.synth) return(0)
      else return(t(Wbar(data=data,tune.pars=tune.pars))%*%Wbar(data=data,tune.pars=tune.pars))
    } 
    Q.matching<-function(data,tune.pars,only.synth){
      if(only.synth) return(matrix(0,ncol(data$donors),ncol(data$donors)))
      else return(diag(ncol(data$donors)))
    } 
  }
  params<-list(pars=NA)
  
    if(model=='average'){ 
      modelreg<-list(
      A= matrix(c(rep(1,length(donors[1,])),rep(0,length(donors[1,])),rep(0,length(donors[1,])),rep(1,length(donors[1,]))),nrow=2,ncol=2*ncol(donors),byrow=TRUE),
      sense = c('=','='),
      rhs = c(1,1),
      lb = rep(0,2*length(donors[1,])),
      ub = rep(1,2*length(donors[1,])),
      obj = c(obj.synthetic(donors=donors,treated=treated,treatment=treatment,only.match=only.match) 
      ,obj.matching(data=data,tune.pars=tune.pars,only.synth=only.synth)),
      objcon = c(objcon.synthetic(donors=donors,treated=treated,treatment=treatment,only.match=only.match)
      ,objcon.matching(data=data,tune.pars=tune.pars,only.synth=only.synth)),
      Q = bdiag(Q.synthetic(donors=donors,treated=treated,treatment=treatment,only.match=only.match) 
      ,Q.matching(data=data,tune.pars=tune.pars,only.synth=only.synth)))
    }

      if(model=='regularize'){
          if(!is.null(tune.pars$unit)) synthmult<-(1-phi)
          else synthmult<-1
          weightlimit<-K
      modelreg<-list(
      A= t(as.matrix(rep(1,length(donors[1,])))),
      sense = '=',
      rhs = 1,
      lb = rep(0,length(donors[1,])),
      ub = rep(1/weightlimit,length(donors[1,])),
      obj = synthmult*obj.synthetic(donors=donors,treated=treated,treatment=treatment,only.match=only.match) 
      + phi*obj.matching(donors=donors,treated=treated,treatment=treatment,only.synth=only.synth),
      objcon = synthmult*objcon.synthetic(donors=donors,treated=treated,treatment=treatment,only.match=only.match)
      + phi*objcon.matching(donors=donors,treated=treated,treatment=treatment,only.synth=only.synth),
      Q = synthmult*Q.synthetic(donors=donors,treated=treated,treatment=treatment,only.match=only.match) 
      + phi*Q.matching(donors=donors,treated=treated,treatment=treatment,only.synth=only.synth))
      }
    continue<-0
    if(!is.null(est.options$solver)){
      if(est.options$solver=="ipop"){
      tparams<-ipop(
        c = modelreg$obj,
        H = 2*modelreg$Q,
        A = modelreg$A,
        b=1,
        l=modelreg$lb,
        u=modelreg$ub,
        r=0,
        sigf=est.options$sigf.ipop,
        margin=est.options$Margin.ipop,
        maxiter=1000
      )
      params$pars<- primal(tparams)
      params$objval<- modelreg$obj%*%params$pars + t(params$pars)%*%modelreg$Q%*%params$pars + modelreg$objcon
      }
    }
    else{
    tparams<-gurobi(modelreg,params=list(OutputFlag=log,PSDTol=psdtol,BarConvTol=BarConvTol,BarIterLimit=BarIterLimit))
    tparams<-rename(tparams,c('x'='pars'))
    if(tparams$status!='OPTIMAL'){
      print('WARNING: FAILED TO SOLVE PROBLEM')
    continue      
    }
    if(model=='regularize') params$pars<-tparams$pars
    if(model=='average')params$pars<-(1-phi)*tparams$pars[1:ncol(donors)]+phi*tparams$pars[(ncol(donors)+1):(2*ncol(donors))]
  if(log==0) suppressWarnings(file.remove("gurobi.log"))
  params$objval<-tparams$objval
    }
  return(params)
}





