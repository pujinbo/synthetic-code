
############PRELIMINARIES#############
library(reshape2)
library(data.table)
library(Synth)
library(plyr)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(forcats)
library(extrafont)
library(ggplot2)
library(tools)
library(textables)

setwd(codepath)
source("Estimator_Code.R")
source("Cross-validation_Code_byK.R")
source("plot_rules.R")

theme_update(strip.text.x=element_text(size=10))


#################SETTING UP DATASETS AND PARAMETERS##########################
#year-over-year GDP per capita changes in Spain:
#price index adjustment, 2010 to 1986 USD: (109.6/216.687)

data(basque)
basque<-as.data.table(basque)
basque<-basque[regionno!=1,]
basque[,regionname:= gsub(" (.*)","",regionname)]

data<-list()
data[[1]]<-list()
data[[1]]$outcome<-cbind(basque[regionno==17,gdpcap],t(reshape(basque[regionno!=17,.(regionno,year,gdpcap)], idvar='regionno', timevar='year',direction='wide')[,-"regionno",with=FALSE]))


setwd(basepath)
params<-readRDS('data_parameters.rds')

addtime<-c(1954)
setwd(outputpath)
philength<-100
Kvals<-1:10
forecast.minmaxmatch<-TRUE
forecast.minlength<-1
forecast.maxlength<- 1
position<-1
tune.pars.list<-list(phi=seq(0,1,length.out=philength),
                     bandwidth=1e20,
                     t.kernel=list(dnorm),K=Kvals,
                     Wbar=list(Wbar),
                     penalty='weight',
                     set.k=NA,
                     forecast.minmaxmatch=forecast.minmaxmatch,
                     forecast.minlength=forecast.minlength,
                     forecast.maxlength=forecast.maxlength
                     )


names <-list()
names[[1]]<- c(unique(basque[regionno==17,regionname]),unique(basque[regionno!=17,regionname]))
names[[1]][10]<-"Catalonia"
names[[1]][2]<-"Andalusia"
names[[1]][4]<-"Asturias"
names[[1]][5]<-"Balearic Isl"
names[[1]][6]<-"Canary Isl"
names[[1]][8]<-"Cast-Leon"
names[[1]][9]<-"Cast-Mancha"
names[[1]][11]<-"Valencia"
names[[1]][16]<-"Navarre"


nameappend <- c("Sp_")
graphnames<-c("Spain")
omits<-list()

treatment_times<-list()
treatment_times[[1]]<-16
control<-length(data)+1

rescale <- c(1000)

pregdp<-data[[1]]$outcome[treatment_times[[1]]-1,1]*rescale

secaxisey<-scale_y_continuous(breaks=pretty_breaks(n=7), sec.axis = sec_axis(~.*100/pregdp, name = "Average RMSPE \n (Pct of GDP per Capita, Basque 1969)", breaks=pretty_breaks(n=5)))
secaxisycv<-scale_y_continuous(breaks=pretty_breaks(n=7), sec.axis = sec_axis(~.*100/pregdp, name = "Difference in Average RMSPE \n (Pct of GDP per Capita, Basque 1969)", breaks=pretty_breaks(n=5)))
secaxisytreat<-scale_y_continuous(breaks=pretty_breaks(n=7), sec.axis = sec_axis(~.*100/pregdp, name = "Pct of GDP per Capita (Basque 1969)", breaks=pretty_breaks(n=5)))
secaxisyplacebo<-scale_y_continuous(breaks=pretty_breaks(n=4), sec.axis = sec_axis(~.*100/pregdp, name = "Pct of GDP per Capita (Basque 1969)", breaks=pretty_breaks(n=5)))

cfact_model_list<-list(c("ADH","SCAMADHAll","MatchADHAll","AbadieCov","AbadieCovSCweight"))
TE_models<-c("SCAMADHAll")

lines <- list(16)
ADH_fittimes<-list()
ADH_fittimes[[1]]<-6:15
firstforecast<-c(8)


pdf("application_results.pdf")

for(run in c(1)){
  omits[[run]]<-2:length(names[[run]])
  
  for(omit in omits[[run]]){
    control<-length(data)+1
    names[[control]]<-c(names[[run]][omit],names[[run]][-c(1,omit)])
    treatment_times[[control]]<-treatment_times[[run]]
    ADH_fittimes[[control]]<-ADH_fittimes[[run]]
    lines<-c(lines,lines[run])
    graphnames<-c(graphnames,paste0(nameappend[run],names[[run]][omit]))
    firstforecast[control]<-firstforecast[run]
    
    data[[control]]<-list()
    data[[control]]$outcome<-cbind(data[[run]]$outcome[,omit],data[[run]]$outcome[,-c(1,omit)])
  }
}
for(dt in c(1:length(data))){
  data[[dt]]$outcome<-data[[dt]]$outcome*rescale
}

############################################################
###########PLOTS:#################

ytitles<-c('Per Capita GDP (USD, 1986)')
treatvals<-c("Spain"=NA)
result<-data.table()
drop_placebos<-c("Balearic Isl","Extremadura","Madrid")
for(dt in c(1:length(params))){
  result1<- data.table(time=1:nrow(data[[dt]]$outcome)) 
  result1[,ADH := as.vector(data[[dt]]$outcome[,-1]%*%params[[dt]]$ADH$pars)]
  result1[,MatchADHAll := as.vector(data[[dt]]$outcome[,-1]%*%params[[dt]]$matchADHAll2$weights$pars)]
  result1[,SCAMADHAll := as.vector(data[[dt]]$outcome[,-1]%*%params[[dt]]$SCAM.CovAll2$weights$pars)]
  result1[,AbadieCov := as.vector(data[[dt]]$outcome[,-1]%*%params[[dt]]$SCAM.abadieCov$weights$pars)]
  result1[,AbadieCovSCweight := as.vector(data[[dt]]$outcome[,-1]%*%params[[dt]]$SCAM.abadieCovSCweight$weights$pars)]
  
  result1[,Treated := as.vector(data[[dt]]$outcome[,1])]
  result1<-melt(result1,id.vars="time")
  result1[,diff:=as.vector(data[[dt]]$outcome[,1])-value,by=variable]
  result1[variable%in%c("SCAM","Abadie","SC","Match")
         ,ddiff:=diff-result1[variable=="SCAM",diff],
         by=variable]
  result1[variable%in%c("SCAMADH","AbadieADH","AbadieCov","SCADH","MatchADH","SCAMADHCovs","ADH","SCAMADHV","MatchADHV","SCAMADHAll","MatchADHAll","AbadieCovSCweight")
         ,ddiff:=diff-result1[variable=="SCAMADHAll",diff],
         by=variable]
  
  result1[,post:=time>=lines[[1]]]
  result1[post==1,cum:=cumsum(diff),by=variable]
  result1[post==1,cumdiff:=cumsum(ddiff),by=variable]
  result1[,region:=names[[1]][dt]]
  result<-rbind(result,result1)
  }
  
  
  
  result[variable%in%c("ADH"),group:="SC with Covs"]
  result[variable%in%c("MatchADHAll"),group:="Matching with Covs"]
  result[variable%in%c("SCAMADHAll"),group:="MASC (both with Covs)"]
  result[variable%in%c("AbadieCov"),group:="Penalized SC with Covs"]
  result[variable%in%c("AbadieCovSCweight"),group:="Penalized SC with Covs (weighted)"]
  result[variable=="Treated",group:="Treated"]

  
  dt<-1


  outdata<-as.data.table(data[[dt]]$outcome)
  outdata[,t:=1:nrow(outdata)]
  outdata<-melt(outdata,id.vars='t')
  outdata[variable=='V1',Type:='Treated']
  outdata[variable!='V1',Type:='Control']

    #Plotting counterfactuals constructed by each estimator:
    result[region=="Basque",regiongroup:="Basque Country"]
    result[region!="Basque",regiongroup:="Control Regions"]
    result[,regiongroup:=factor(regiongroup,c("Control Regions","Basque Country"))]
    result[,region:=fct_relevel(region,"Basque",after=Inf)]
    for(mdl in cfact_model_list){
plot<-  (ggplot(data=result[variable%in%c(mdl,"Treated") & region=="Basque"
                            &!variable%in%c("SCAMADHAll","AbadieCovSCweight"),],aes(x=time+addtime[dt]))
        +  geom_line(aes(y=value,color=group,linetype=group),size=1.4)
        +geom_vline(xintercept=lines[[dt]]+addtime[dt],linetype='dashed',color='black')
        +graycolors
        +secaxisytreat
        +labs(x="Year",y=ytitles[dt]))

print(plot+   theme(legend.position="none"))



result[,labelvar:=variable]
result[variable=="ADH",labelvar:="Synthetic Control"]
result[variable=="SCAMADHAll",labelvar:="MASC"]
result[variable=="AbadieCovSCweight",labelvar:="PSC-AG"]
result[variable=="AbadieCov",labelvar:="PSC-S"]
result[variable=="MatchADHAll",labelvar:="Matching"]


plot1<- (ggplot(data=result[variable%in%mdl& ! region%in%c(drop_placebos) ,],aes(x=time+addtime[dt]))
         +  geom_line(aes(y=diff,color=interaction(regiongroup),
                          group=interaction(region)),size=1.4)
         +geom_vline(xintercept=lines[[dt]]+addtime[dt],linetype='dashed',color='black')
         +geom_hline(yintercept=0,linetype='dashed',color='black')
         +scale_color_manual(breaks=c(paste("Basque Country"),paste("Control Regions")),
                             values=c("black","gray80"),
                             labels=c("Basque Country","Placebos"))
         +secaxisyplacebo
         +labs(x="Year",y=ytitles[1])
         +facet_wrap(fct_relevel(labelvar,c("Synthetic Control","MASC","PSC-AG","PSC-S","Matching"),after=Inf)~.,
                     labeller= ,
                     scales="fixed",ncol=2)
         )

print(plot1+ 
        labs(color="",linetype="",shape="")+
      theme(legend.position=c(1,0), legend.justification=c(0.9,-0.3)))


  #exporting legend
  print(as_ggplot(get_legend(plot+guides(col=guide_legend(nrow=2))+ theme(legend.key.size=unit(3,"lines")))))
  #exporting legend
  print(as_ggplot(get_legend(plot1+labs(color="",linetype="",shape="")+theme(legend.direction="horizontal")+ theme(legend.key.size=unit(3,"lines")))))

  rm(plot1)
  }

    
######PLACEBO RESULTS######
#These placebos have relatively high prediction error
high_placebos<-c("Canary Isl","Cast-Mancha","Catalonia","Asturias","Valencia","Aragon","Cantabria","Balearic Isl","Extremadura","Madrid")

limpos<-0
for(drop_placebos in list(c("Balearic Isl","Extremadura","Madrid"),c("None"))){
  limpos<-limpos+1
  
  
plotdata<-cbind(graphnames[!graphnames%in%c("Spain")])
plotdata<-as.data.table(plotdata)
names(plotdata)<-"Placebo"
plotdata[grepl('Sp_',Placebo),dataset:='Spain']
plotdata[,Placebo:=gsub('Sp_','',Placebo)]
candnames<-NULL
mods<-c("ADH","SCAM.multiCV5","SCAM.multiCV3","SCAM.setCV5","SCAM.setCV3","matchADHAll2","SCAM.CovAll2","SCAM.byphiCovAll2","SCAM.abadieCov","SCAM.abadieCovSCweight")

placebonames<-unique(plotdata$Placebo)
##Making dataset of placebo results:
for(mdl in mods){
  if(mdl =="ADH")mdldata <- data.table(
    phi = sapply(params[!graphnames%in%c("Spain")]
                 ,function(x) x[[mdl]]$tune.pars$phi),
    K = sapply(params[!graphnames%in%c("Spain")]
               ,function(x) x[[mdl]]$tune.pars$K),
    minmspe4 = NA,
    minmspeall = NA,
    phi_mspe4 =NA,
    phi_mspeall =NA,
    K_mspe4 = NA,
    K_mspeall = NA,
    MSPEAll = sapply(params[!graphnames%in%c("Spain")]
                     ,function(x)mean(x[[mdl]]$pred.error^2)^0.5)*rescale,
    MSPE1 = sapply(params[!graphnames%in%c("Spain")]
                   ,function(x)mean(x[[mdl]]$pred.error[1]^2)^0.5)*rescale,
    MSPE4 = sapply(params[!graphnames%in%c("Spain")]
                   ,function(x)mean(x[[mdl]]$pred.error[1:4]^2)^0.5)*rescale,
    MSPE10 =   sapply(params[!graphnames%in%c("Spain")]
                      ,function(x)mean(x[[mdl]]$pred.error[1:10]^2)^0.5)*rescale,
    MSPE20 =   sapply(params[!graphnames%in%c("Spain")]
                      ,function(x)mean(x[[mdl]]$pred.error[1:20]^2)^0.5)*rescale,
    fit =   sapply(params[!graphnames%in%c("Spain")]
                   ,function(x)mean(x[[mdl]]$fit^2)^0.5)*rescale,
    numweights =   sapply(params[!graphnames%in%c("Spain")]
                          ,function(x)sum(x[[mdl]]$weights$pars>=0.01))
  )
  else mdldata<-data.table(
    phi = sapply(params[!graphnames%in%c("Spain")]
                 ,function(x) x[[mdl]]$tune.pars$phi),
    K = sapply(params[!graphnames%in%c("Spain")]
               ,function(x) x[[mdl]]$tune.pars$K),
    minmspe4 = sapply(params[!graphnames%in%c("Spain")]
                      ,function(x) min(x[[mdl]]$all.results$pred.error.4)^0.5)*rescale,
    minmspeall = sapply(params[!graphnames%in%c("Spain")]
                      ,function(x) min(x[[mdl]]$all.results$pred.error.all)^0.5)*rescale,
    phi_mspe4 =sapply(params[!graphnames%in%c("Spain")]
                      , function(x) seq(0,1,length.out=philength)[(which.min(x[[mdl]]$all.results$pred.error.4)-1)%%philength+1]),
    phi_mspeall =sapply(params[!graphnames%in%c("Spain")]
                      , function(x) seq(0,1,length.out=philength)[(which.min(x[[mdl]]$all.results$pred.error.all)-1)%%philength+1]),
    K_mspe4 = sapply(params[!graphnames%in%c("Spain")]
                     , function(x) Kvals[floor((which.min(x[[mdl]]$all.results$pred.error.4)-1)/philength)+1]),
    K_mspeall = sapply(params[!graphnames%in%c("Spain")]
                     , function(x) Kvals[floor((which.min(x[[mdl]]$all.results$pred.error.all)-1)/philength)+1]),
    MSPEAll = sapply(params[!graphnames%in%c("Spain")]
                     ,function(x)mean(x[[mdl]]$pred.error^2)^0.5)*rescale,
    MSPE1 = sapply(params[!graphnames%in%c("Spain")]
                   ,function(x)mean(x[[mdl]]$pred.error[1]^2)^0.5)*rescale,
    MSPE4 = sapply(params[!graphnames%in%c("Spain")]
                       ,function(x)mean(x[[mdl]]$pred.error[1:4]^2)^0.5)*rescale,
    MSPE10 =   sapply(params[!graphnames%in%c("Spain")]
                      ,function(x)mean(x[[mdl]]$pred.error[1:10]^2)^0.5)*rescale,
    MSPE20 =   sapply(params[!graphnames%in%c("Spain")]
                      ,function(x)mean(x[[mdl]]$pred.error[1:20]^2)^0.5)*rescale,
    fit =   sapply(params[!graphnames%in%c("Spain")]
                   ,function(x)mean(x[[mdl]]$fit^2)^0.5)*rescale,
    numweights =   sapply(params[!graphnames%in%c("Spain")]
                          ,function(x)sum(x[[mdl]]$weights$pars>=0.01))
  )
  mods<-c("ADH","SCAM.multiCV5","SCAM.multiCV3","SCAM.setCV5","SCAM.setCV3","matchADHAll2","SCAM.CovAll2","SCAM.byphiCovAll2","SCAM.abadieCov","SCAM.abadieCovSCweight")
  
  names(mdldata)<-paste(names(mdldata),c('ADH',"multiCV5","multiCV3",'setCV5','setCV3',"matchADHAll","averageCovAll","averageCovAllbyphi","abadieADHcovs","abadieADHcovsSCweight")[which(mods==mdl)],sep='.')
  plotdata<-cbind(plotdata,mdldata)
}

plotdata<-as.data.table(plotdata)

####Summary Table####
if(limpos == 1){
tabledata<-copy(plotdata)
vars<-c("Placebo","MSPEAll.ADH","MSPEAll.averageCovAll","MSPEAll.abadieADHcovs","MSPEAll.abadieADHcovsSCweight","MSPEAll.matchADHAll","fit.ADH","fit.averageCovAll",'fit.abadieADHcovs',"fit.abadieADHcovsSCweight",'fit.matchADHAll')
varsmin<-c("Placebo","MSPEAll.ADH","minmspeall.averageCovAllbyphi",'minmspeall.abadieADHcovs',"minmspeall.abadieADHcovsSCweight",'minmspeall.matchADHAll')
tabledatamin<-tabledata[!Placebo%in%drop_placebos,..varsmin]
tabledata<-tabledata[!Placebo%in%drop_placebos,..vars]

tab<-NULL
tab<-TR(c("","RMSPE","Pre-Period Fit"),cspan=c(1,5,5)) + 
  vspace(5) + 
  midrulep(list(c(2,6),c(7,11))) +
  TR(c("Placebo","SC","MASC","PSC-S", "PSC-AG","Matching","SC","MASC","PSC-S","PSC-AG","Matching"))+
  midrule()
for(i in 1:nrow(tabledata)){
  tab<-tab+TR(unlist(tabledata[i,1]))%:%TR(unlist(tabledata[i,-1]),dec=0) +
  TR("")%:%TR(unlist(tabledatamin[i,-1]), dec=0, se = TRUE)%:%TR(rep("",4))
}
#ADD AVERAGE
tab<- tab +
  midrule()+
  TR("\\textbf{Average}")%:%TR(apply(tabledata[,-1],MARGIN=2,FUN=mean),dec=0,surround="\\textbf{ %s}") +
  TR("")%:%TR(apply(tabledatamin[,-1],MARGIN=2,FUN=mean), dec=0, se = TRUE,surround="\\textbf{ %s}")%:%TR(rep("",5))

TS(tab, file=paste0("Placebo_Performances_",limpos),
   pretty_rules=T,
   header=c('r',rep('c',10)))

}

###Table of MSPEs####
if(limpos==2){
tabledata<-copy(plotdata)
vars<-c("Placebo",
        "MSPEAll.ADH","MSPEAll.averageCovAll","MSPEAll.abadieADHcovs","MSPEAll.abadieADHcovsSCweight","MSPEAll.matchADHAll",
        "MSPE4.ADH","MSPE4.averageCovAll","MSPE4.abadieADHcovs","MSPE4.abadieADHcovsSCweight","MSPE4.matchADHAll",
        "MSPE10.ADH","MSPE10.averageCovAll","MSPE10.abadieADHcovs","MSPE10.abadieADHcovsSCweight","MSPE10.matchADHAll",
        "fit.ADH","fit.averageCovAll","fit.abadieADHcovs","fit.abadieADHcovsSCweight","fit.matchADHAll"
        )
tabledata<-tabledata[!Placebo%in%drop_placebos,..vars]


tab<-NULL
tab<-TR(c("","Pre-Period Fit","4-Year RMSPE"),cspan=c(1,5,5)) + 
  vspace(2) + 
  midrulep(list(c(2,6),c(7,11))) +
  TR(c("Placebo","SC","MASC","PSC-S", "PSC-AG","Matching","SC","MASC","PSC-S","PSC-AG","Matching"))+
  midrule()

for(i in 1:nrow(tabledata)){
  tab<-tab+TR(unlist(tabledata[i,1]))%:%TR(unlist(tabledata[i,.(fit.ADH,fit.averageCovAll,fit.abadieADHcovs,fit.abadieADHcovsSCweight,fit.matchADHAll,
                                                                MSPE4.ADH,MSPE4.averageCovAll,MSPE4.abadieADHcovs,MSPE4.abadieADHcovsSCweight,MSPE4.matchADHAll)]),dec=0)
}

#ADD AVERAGE
tab<- tab +
  midrule()+
  TR("\\textbf{Average}")%:%TR(apply(tabledata[,.(fit.ADH,fit.averageCovAll,fit.abadieADHcovs,fit.abadieADHcovsSCweight,fit.matchADHAll,
                                                  MSPE4.ADH,MSPE4.averageCovAll,MSPE4.abadieADHcovs,MSPE4.abadieADHcovsSCweight,MSPE4.matchADHAll)],MARGIN=2,FUN=mean),dec=0,surround="\\textbf{ %s}")+ 
 midrule()
  


#20 and 28 year
tab<- tab+ 
  TR(c("","10-Year RMSPE","28-Year RMSPE"),cspan=c(1,5,5)) + 
  vspace(2) +
  midrulep(list(c(2,6),c(7,11))) +
  #TR(c("","Penalized SC","","","Penalized SC",""),cspan=c(3,2,1,2,2,1))+
  TR(c("Placebo","SC","MASC","PSC-S", "PSC-AG","Matching","SC","MASC","PSC-S","PSC-AG","Matching"))+
  midrule()

for(i in 1:nrow(tabledata)){
  tab<-tab+TR(unlist(tabledata[i,1]))%:%TR(unlist(tabledata[i,.(MSPE10.ADH,MSPE10.averageCovAll,MSPE10.abadieADHcovs,MSPE10.abadieADHcovsSCweight,MSPE10.matchADHAll,
                                                                MSPEAll.ADH,MSPEAll.averageCovAll,MSPEAll.abadieADHcovs,MSPEAll.abadieADHcovsSCweight,MSPEAll.matchADHAll)]),dec=0)
}

#ADD AVERAGE
tab<- tab +
  midrule()+
  TR("\\textbf{Average}")%:%TR(apply(tabledata[,.(MSPE10.ADH,MSPE10.averageCovAll,MSPE10.abadieADHcovs,MSPE10.abadieADHcovsSCweight,MSPE10.matchADHAll,
                                                   MSPEAll.ADH,MSPEAll.averageCovAll,MSPEAll.abadieADHcovs,MSPEAll.abadieADHcovsSCweight,MSPEAll.matchADHAll)],
                                     MARGIN=2,FUN=mean),dec=0,surround="\\textbf{ %s}") 


TS(tab, file=paste0("Placebo_RMSPEs_",limpos),
   pretty_rules=T,
   header=c('r',rep('c',10)))
}

}

###FIGURE-BASED RESULTS NOW:####
drop_placebos <- c("Balearic Isl","Extremadura","Madrid")

plotdata<-reshape(plotdata,direction='long',varying=names(plotdata)[-c(1:2)],sep='.')
temp<-plotdata[!Placebo%in%drop_placebos,lapply(.SD,mean),by=.(dataset,time),.SDcols=names(plotdata)[!names(plotdata)%in%c("dataset","time","Placebo")]]
temp[,Placebo:="Average"]
plotdata<-rbind(plotdata,temp)
rm(temp)
plotdata[,Placebo:=factor(Placebo,levels=c(placebonames,"Average"))]
plotdata[time=='ADH',time:='SC with Covs']
plotdata[time=='multiCV3',time:='multi-step ahead (3)']
plotdata[time=='multiCV5',time:='multi-step ahead (5)']
plotdata[time=='setCV5',time:='fixed window (5)']
plotdata[time=='setCV3',time:='fixed window (3)']


#restructuring dataset of placebo results for easier plotting:
meltdata<-melt(plotdata,id.vars=c('Placebo','dataset','time'))

meltcv<-copy(plotdata)
meltcv<-meltcv[,names(meltcv)%in%c("Placebo","dataset","time","MSPE1","MSPE4","MSPE10","MSPE20","MSPEAll"),with=FALSE]
setnames(meltcv,c("MSPEAll","time"),c("MSPE28","estimator"))
meltcv<-reshape(meltcv,varying = c("MSPE1","MSPE4","MSPE10","MSPE20","MSPE28"),direction="long",sep="SPE")
setnames(meltcv,c("M","time"),c("MSPE","period"))
nums<-unique(meltcv$period)
meltcv[,period:=as.character(period)]
meltcv[,period:=paste0(period," Years")]
meltcv[period=="1 Years",period:="1 Year"]
meltcv[,period:=factor(period,levels=c("1 Year",paste0(nums[-1]," Years")))]
meltcv<-meltcv[order(estimator,dataset,period,Placebo,),]
meltcv$MSPEdiff_MASC<-meltcv$MSPE-meltcv[estimator=="averageCovAll",MSPE]

dtnum<-0
for(dt in levels(as.factor(meltdata$dataset))){
  dtnum<-dtnum+1
   
     outcome<-"MSPE4" 

    meltdata<-meltdata[,Placebo:=factor(Placebo,levels=c(
      meltdata[variable==outcome&time=='averageCovAll'&dataset==dt&Placebo!="Average",unique(as.character(Placebo[order(value)]))],
                                        "Average"))]
    
    for(plotmodels in list(c("averageCovAll","SC with Covs","abadieADHcovsSCweight","abadieADHcovs","matchADHAll"))){

      outcomenum<-0

      ##AVERAGE RMSPE BY ESTIMATOR, OVER ALTERNATIVE PREDICTION WINDOWS:
      myplot<-(ggplot(data=meltcv[dataset==dt
                                    &Placebo == "Average"
                                    &estimator%in%plotmodels
                                  &period!="28 Years"
                                    ,]) + 
                 geom_col(aes(x=as.factor(period),y=as.numeric(MSPE), fill=factor(estimator,levels=plotmodels),  
                                      group = factor(estimator,levels=plotmodels)),
                                  color="black",
                                  position=position_dodge2(preserve="single")) + 
                 coord_cartesian(ylim=c(min(meltcv[dataset==dt
                                                   &Placebo == "Average"
                                                   &estimator%in%plotmodels
                                                   &period!="28 Years"
                                                   ,]$MSPE)*0.95,
                                        max(meltcv[dataset==dt
                                                   &Placebo == "Average"
                                                   &estimator%in%plotmodels
                                                   &period!="28 Years"
                                                   ,]$MSPE)))+
                 list(scale_fill_manual(breaks=c("averageCovAll","SC with Covs","abadieADHcovsSCweight","abadieADHcovs","matchADHAll"),
                                        values=c("abadieADHcovs"=colcols[4],
                                                 "abadieADHcovsSCweight"=colcols[3],
                                                 "SC with Covs"=colcols[2],"matchADHAll"="white",
                                                 "averageCovAll"=colcols[1]),
                                        name = element_blank(),
                                        labels=c("MASC","Synth. Control","PSC-AG","PSC-S","Matching")
                 ),
                 scale_color_manual(breaks=c("averageCovAll","SC with Covs","abadieADHcovsSCweight","abadieADHcovs","matchADHAll"),
                                    values=c("abadieADHcovs"=colcols[4],
                                             "abadieADHcovsSCweight"=colcols[3],
                                             "SC with Covs"=colcols[2],"matchADHAll"="white",
                                             "averageCovAll"=colcols[1]),
                                    name = element_blank(),
                                    labels=c("MASC","Synth. Control","PSC-AG","PSC-S","Matching")
                 )
                 )+
                 labs(x='Length of Post-Period Interval',y="Average RMSPE (GDP per capita)")+
                 secaxisey
               +facet_wrap(.~period,scales="free")
               +theme(
                 strip.background = element_blank(),
                 strip.text.x = element_blank()
               )
      )
      

      


      print(myplot+   theme(legend.position="none"))
      print(as_ggplot(get_legend(myplot+guides(fill=guide_legend(nrow=ceiling(length(plotmodels)/2)))+theme(legend.key.size=unit(2,"lines"),legend.direction="horizontal"))))
      
      
      
    }

    
     plotmodels<-c("MASC","Synth. Control","Penalized SC","Matching","SC with Covs")


   #Plotting diff in RMSPE between our procedure and alternative CV procedures for forecasting
   myplot<-(ggplot(data=meltcv[dataset==dt
                           & !Placebo%in%drop_placebos
                           &Placebo!="Average"
                           &estimator%in%c("multi-step ahead (5)","fixed window (5)",
                                           "multi-step ahead (3)","fixed window (3)")
                           ,]) +
           geom_col(data= meltcv[dataset==dt
                                 &estimator%in%c("multi-step ahead (5)","fixed window (5)",
                                                 "multi-step ahead (3)","fixed window (3)")
                                 & !Placebo%in%drop_placebos
                                &  Placebo=="Average"
                                ,],
                    aes(y=MSPEdiff_MASC,x=period,
                        group=factor(estimator,levels=c("multi-step ahead (3)","fixed window (3)",
                                                        "multi-step ahead (5)","fixed window (5)")),
                        fill=estimator),
                    position=position_dodge()
                    ) +
           geom_point(aes(y=MSPEdiff_MASC,x=period,
                          group=factor(estimator,levels=c("multi-step ahead (3)","fixed window (3)",
                                                          "multi-step ahead (5)","fixed window (5)"))),
                      position=position_dodge(width=1),
                      shape=1,size=2) +
           geom_hline(yintercept=0,linetype="dashed")+
           columncolors_cv+
           secaxisycv+
           textadjust+
           labs(y=paste0("Difference in Average RMSPE \n (GDP per capita)"),x="Length of Post-Period Interval")
         +theme(axis.title.y= element_text(size=15.5))
   )
   print(myplot+   theme(legend.position="none"))
   print(as_ggplot(get_legend(myplot+guides(fill=guide_legend(nrow=floor(length(plotmodels)/2)))+theme(legend.key.size=unit(2,"lines"),legend.direction="horizontal"))))
   

   
}
#fit and MSPE curves over values of phi, given matching estimator, averaging over placebos:
#excluding spain and any of the placebo regions we decided to drop:
fit_mspedata<-data.table()
for(regions in list(which(!names[[1]]%in%drop_placebos)[-1],9)){
  if(length(regions)!=1) regionname<-"Mean over Placebos"
  else regionname<-"Castile-La Mancha"
fit_mspedata<-rbind(fit_mspedata,
                    data.table(
                      Placebo=regionname,phi=seq(0,1,length.out=100),
      mspe.average = apply(sapply(params[regions], function(x) x$SCAM.byphiCovAll2$all.results$pred.error.all[(x$SCAM.byphiCovAll2$tune.pars$K-1)*philength+1:philength]^0.5)*rescale,
                                    MARGIN=1,FUN=mean),
      fit.average = apply(sapply(params[regions], function(x) x$SCAM.byphiCovAll2$all.results$fit[(x$SCAM.byphiCovAll2$tune.pars$K-1)*philength+1:philength]^0.5)*rescale,
      MARGIN=1,FUN=mean),
      mspe.abadie = apply(sapply(params[regions], function(x) x$SCAM.abadieCovSCweight$all.results$pred.error.all[(x$SCAM.abadieCovSCweight$tune.pars$K-1)*philength+1:philength]^0.5)*rescale,
                   MARGIN=1,FUN=mean),
      fit.abadie = apply(sapply(params[regions], function(x) x$SCAM.abadieCovSCweight$all.results$fit[(x$SCAM.abadieCovSCweight$tune.pars$K-1)*philength+1:philength]^0.5)*rescale,
                          MARGIN=1,FUN=mean)))
}


myplot<-(ggplot(data=fit_mspedata)
           +geom_line(aes(x=phi,y=mspe.average,group=Placebo,linetype=as.factor(Placebo)),size=0.8)
           +labs(x=expression(Tuning~Parameter~(phi)),y=paste('Prediction error (RMSPE)'),color='',linetype='')
           +expand_limits(y=c(0,0.285)*rescale)
            +scale_linetype_manual(values=c("Castile-La Mancha"="solid","Mean over Placebos"='dashed'))
)

print(myplot+   theme(legend.position="none"))

myplot<-(ggplot(data=fit_mspedata)
        +geom_line(aes(x=phi,y=mspe.abadie,group=Placebo,linetype=as.factor(Placebo)),size=0.8)
         +labs(x=expression(Tuning~Parameter~(pi)),y=paste('Prediction error (RMSPE)'),color='',linetype='')
        +expand_limits(y=c(0,0.285)*rescale)
        +scale_linetype_manual(values=c("Castile-La Mancha"="solid","Mean over Placebos"='dashed'))
)

print(myplot+   theme(legend.position="none"))

myplot<-(ggplot(data=fit_mspedata)
         +geom_line(aes(x=phi,y=fit.average,group=Placebo,linetype=Placebo),size=0.8)
         +labs(x=expression(Tuning~Parameter~(phi)),y=paste('Pre-period fit (RMSE)'),linetype='')
        +expand_limits(y=c(0,0.16)*rescale)
        +scale_linetype_manual(values=c("Castile-La Mancha"="solid","Mean over Placebos"='dashed'))
)
print(myplot+   theme(legend.position="none"))

myplot<-(ggplot(data=fit_mspedata)
         +geom_line(aes(x=phi,y=fit.abadie,group=Placebo,linetype=Placebo),size=0.8)
         +labs(x=expression(Tuning~Parameter~(pi)),y=paste('Pre-period fit (RMSE)'),linetype='')
        +expand_limits(y=c(0,0.16)*rescale)
        +scale_linetype_manual(values=c("Castile-La Mancha"="solid","Mean over Placebos"='dashed'))
)
print(myplot+   theme(legend.position="none"))

print(as_ggplot(get_legend(myplot+labs(shape="Estimator:")+guides(color=guide_legend(order=1))+theme(legend.key.size=unit(2,"lines"),legend.direction="horizontal"))))



dev.off()

