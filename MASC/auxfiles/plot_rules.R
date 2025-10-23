graycols<-gray.colors(n=2,start=0,end=0.6)
colcols<-c(gray.colors(n=4,start=0,end=0.8), hcl.colors(n=7,palette = "reds"))


graybreaks<-c("Treated","SC with Covs","Penalized SC with Covs","Matching with Covs")
graylabels<-c("Actual with Terrorism", "MASC/SC/PSC-AG", "PSC-S", "Matching")

colors<-list(scale_shape_manual(values=c("MASC"=15,"Penalized SC"=17,"Synth. Control"=18,"Matching"=1,"Controls"=NA,"Treated"=NA)),
             scale_linetype_manual(values=c("MASC"="dashed","Penalized SC"="dashed","Synth. Control"="dashed","Matching"="dashed","Controls"="solid","Treated"="solid")),
             guides(shape=guide_legend(nrow=3),linetype=guide_legend(nrow=3)))


graycolors<- list(scale_color_manual(breaks=graybreaks,
                                     values=c('Penalized SC with Covs'=graycols[2],
                                              "Treated"="black",
                                              "SC with Covs"=graycols[2],'Matching with Covs'=graycols[1]),
                                     labels=graylabels,
                                     name=element_blank()
                                     ),
                  scale_fill_manual(breaks=graybreaks,
                                    values=c('Penalized SC with Covs'=graycols[2],
                                             "Treated"="black",
                                             "SC with Covs"=graycols[2],'Matching with Covs'=graycols[1]),
                                    labels=graylabels,
                                    name=element_blank()
                                    ),
                  scale_linetype_manual(breaks=graybreaks,
                                        values=c('Penalized SC with Covs'="solid",
                                                 "Treated"="solid",
                                                 "SC with Covs"="dashed", "Matching with Covs"="dashed"),
                                        labels=graylabels,
                                        name=element_blank()
                                        ),
                  scale_shape_manual(breaks=graybreaks,
                                    values=c('Penalized SC with Covs'=NA,
                                             "Treated"=19,
                                             "SC with Covs"=NA, "Matching with Covs"=NA),
                                    labels=graylabels,
                                    name=element_blank())
)


columncolors_cv<-list(scale_fill_manual(breaks=c("multi-step ahead (3)","fixed window (3)",
                                                 "multi-step ahead (5)","fixed window (5)"),
                                     values=c("multi-step ahead (3)"=colcols[1], 'fixed window (3)'=colcols[2],
                                              "multi-step ahead (5)"=colcols[3], 'fixed window (5)'=colcols[4]),
                                     name = element_blank(),
                                     labels=c("multi-step ahead (3)","fixed window (3)",
                                              "multi-step ahead (5)","fixed window (5)")
)
)

theme_set(theme_bw())
theme_update(text= element_text(size=22,color='black'),axis.text.y = element_text(color='black',size=14),
             axis.text.x = element_text(color='black',angle = 0, hjust = 0.5,size=14),
             plot.title = element_text(size=22,color='black',hjust=0.5), panel.grid.major=element_blank(),
             legend.background=element_blank(),
             panel.grid.minor=element_blank())

addtoplot<-list(scale_shape_manual(breaks=c("Synthetic Control","MASC","Matching","Penalized SC"),values=c("MASC"=20,"Matching"=0,"Penalized SC"=5,"Synthetic Control"=2,"Controls"=NA,"Treated"=NA)))
textadjust<-theme(axis.text.x=element_text(angle=90,vjust=0.5))

