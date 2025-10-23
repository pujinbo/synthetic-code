rm(list=ls())

## please set the working directory
setwd("~/gitrepos/tscs/code/AJPS_tscs_replication/")

## check/create the folder to put outputs
if (file.exists(file.path(getwd(), "output"))){
    OUT_DIR <- file.path(getwd(), "output")
} else {
    dir.create(file.path(file.path(getwd(), "output")))
    OUT_DIR <- file.path(getwd(), "output")
}

##### Making Figure 1, the Treatment Variation Plot #####

library(PanelMatch)
library(ggplot2)
library(gridExtra)

load("Acemoglu.RData")
load("SS.RData")



###### Making Figure 1: the Treatment Variation Plot ######

ADis <- DisplayTreatment(unit.id = "wbcode2",
                         time.id = "year", 
                         xlab = "Year", ylab = "Countries", legend.position = "bottom",
                         legend.labels = c("Autocracy (Control)", "Democracy (Treatment)"),
                         title = "Democracy as the Treatment",
                         treatment = "dem", data = d2) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=0, size = 6.5, vjust=0.5)) + 
  scale_y_discrete(breaks = c(1960, 1970, 1980, 1990, 2000, 2010))



SDis <- DisplayTreatment(unit.id = "name", 
                         time.id = "year", legend.position = "bottom",
                         xlab = "Year", ylab = "Countries",
                         legend.labels = c("Peace (Control)", "War (Treatment)"),
                         y.size = 10,
                         title = "War as the Treatment",
                         treatment = "himobpopyear2p", data = d3) + 
  theme(axis.text.x = element_text(angle=0, size = 6.5, vjust=0.5), 
        axis.ticks.y = element_blank()) + 
  scale_y_discrete(breaks = c(1850, 1900, 1950, 2000))
g <- arrangeGrob(ADis, SDis, ncol = 2)

ggsave(file = "Figure1.pdf", height = 6, width = 10,
       g, path = OUT_DIR)



