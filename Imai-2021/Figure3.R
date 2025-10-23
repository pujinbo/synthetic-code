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


##### Making Figure 3, the Histogram #####
## For country name footnotes in the histogram section, see the end of the script ##

library(PanelMatch)
library(ggplot2)
library(gridExtra)

load("Acemoglu.RData")
load("SS.RData")


##############################################################
### Getting basic matched sets to produce histograms
## storing cov.formula
formula_Ace <- 
  ~ I(lag(y, 1:4)) + I(lag(Populationages014oftotal, 1:4)) + 
  I(lag(Populationages1564oftota, 1:4)) + I(lag(unrest, 1:4)) + 
  I(lag(tradewb, 1:4)) + I(lag(nfagdp, 1:4)) + I(lag(logpop, 1:4))

formula_SS <-
  ~ I(lag(topitaxrate2, 1:4)) + 
  I(lag(leftexec2, 1:4)) + I(lag(unisuffrage, 1:4)) + 
  I(lag(rgdppc, 1:4)) 

formula_Ace_l1 <- 
  ~ I(lag(y, 1)) + I(lag(Populationages014oftotal, 1)) + 
  I(lag(Populationages1564oftota, 1)) + I(lag(unrest, 1)) + 
  I(lag(tradewb, 1)) + I(lag(nfagdp, 1)) + I(lag(logpop, 1))

formula_SS_l1 <-
  ~ I(lag(topitaxrate2, 1)) + 
  I(lag(leftexec2, 1)) + I(lag(unisuffrage, 1)) + 
  I(lag(rgdppc, 1)) 


## getting matched sets for L4 Acemoglu
matches_L4_cbps <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
                              treatment = "dem", outcome = "y", 
                              refinement.method = "mahalanobis", 
                              qoi = "ate",
                              data = d2, 
                              match.missing = FALSE, 
                              covs.formula = formula_Ace,
                              size.match = 5)

## getting matched sets for L4 SS
matches_L4_cbps_SS <- PanelMatch(lag = 4, time.id = "year", unit.id = "ccode", 
                                 treatment = "himobpopyear2p", outcome = "topitaxrate2", 
                                 refinement.method = "mahalanobis", 
                                 qoi = "ate",
                                 data = d3, 
                                 match.missing = T, 
                                 covs.formula = formula_SS,
                                 size.match = 5)

## getting matched sets for L1 Acemoglu
matches_L1_cbps <- PanelMatch(lag = 1, time.id = "year", unit.id = "wbcode2", 
                              treatment = "dem", outcome = "y", 
                              refinement.method = "mahalanobis", 
                              qoi = "ate",
                              data = d2, 
                              match.missing = T, 
                              covs.formula = formula_Ace_l1,
                              size.match = 5)

## getting matched sets for L1 SS
matches_L1_cbps_SS <- PanelMatch(lag = 1, time.id = "year", unit.id = "ccode", 
                                 treatment = "himobpopyear2p", outcome = "topitaxrate2", 
                                 refinement.method = "mahalanobis", 
                                 qoi = "ate",
                                 data = d3, 
                                 match.missing = T, 
                                 covs.formula = formula_SS_l1,
                                 size.match = 5)

## extract objects to make histograms with
L4_cbps_att_sizes <- summary(matches_L4_cbps$att)$overview$matched.set.size
L4_cbps_atc_sizes <- summary(matches_L4_cbps$atc)$overview$matched.set.size


L4_cbps_att_sizes_SS <- summary(matches_L4_cbps_SS$att)$overview$matched.set.size
L4_cbps_atc_sizes_SS <- summary(matches_L4_cbps_SS$atc)$overview$matched.set.size


L1_cbps_att_sizes <- summary(matches_L1_cbps$att)$overview$matched.set.size
L1_cbps_atc_sizes <- summary(matches_L1_cbps$atc)$overview$matched.set.size


L1_cbps_att_sizes_SS <- summary(matches_L1_cbps_SS$att)$overview$matched.set.size
L1_cbps_atc_sizes_SS <- summary(matches_L1_cbps_SS$atc)$overview$matched.set.size

## Making the histograms with the objects extracted above
pdf(file = file.path(OUT_DIR, "Figure3.pdf"), width = 10, height = 6)
den.matrix <- matrix(1:4, nrow = 2, ncol = 2)
den.matrix[1,] <- 1:2
den.matrix[2,] <- 3:4
den.matrix <- 
  layout(
    den.matrix
  )
par(mar = c(2, 4, 2, 1), oma = c(2.5, 3, 1.5, 0))
hist(L4_cbps_att_sizes[L4_cbps_att_sizes>0], 
     col = "#ffc6c4", border = NA,
     freq = TRUE,
   
     breaks = c(1, 5, 10, 15, 20, 25,
                30, 35, 40, 45,
                50, 55, 60, 65, 70, 75, 80,
                85, 90, 95, 100, 105, 110, 115, 120, 125, 130),
     xlim = c(0, 130), ylim= c(0, 30), xlab = "",
     main = "Democratization",
     ylab = "")
lines(x = c(0,0), 
      y = c(0, length(L4_cbps_att_sizes[L4_cbps_att_sizes==0])), 
      lwd = 4,
      col = "#ffc6c4")

white.hist <- hist(L1_cbps_att_sizes[L1_cbps_att_sizes>0] , 
                   col = rgb(1,1,1,0),     
                   
                   breaks = c(1, 5, 10, 15, 20, 25,
                              30, 35, 40, 45,
                              50, 55, 60, 65, 70, 75, 80,
                              85, 90, 95, 100, 105, 110, 115, 120,
                              125, 130),
                   plot = FALSE, lty = 1)
lines(white.hist$breaks, c(0, white.hist$counts),
      type = 'S', lty = 1)
legend(x = 0, y = 30,  
       legend = c(paste0("Four Year Lags"), 
                  #                         " empty matched sets")
                  paste0("One Year Lag")),
      
       fill = c("#ffc6c4", "white"),
      
       xjust = 0,
      
       pt.cex = 1,
       bty = "n", ncol = 1, cex = 1, bg = "white")

mtext("Frequency", side = 2, line = 2.2, outer=FALSE,
      cex = 0.8)


hist(L4_cbps_atc_sizes[L4_cbps_atc_sizes>0], 
     col = "#ffc6c4", border = NA,
     freq = TRUE,
     # # main= "Average Treatment Effect \n for Treated",
     # breaks = seq(0, 130, by = 5),
     breaks = c(1, 5, 10, 15, 20, 25,
                30, 35, 40, 45,
                50, 55, 60, 65, 70, 75, 80,
                85, 90, 95, 100, 105, 110, 115, 120, 125, 130),
     xlim = c(0, 130), ylim= c(0, 30), xlab = "",
     main = "Authoritarian Reversal",
     ylab = "")
lines(x = c(0,0), 
      y = c(0, length(L4_cbps_atc_sizes[L4_cbps_atc_sizes==0])), 
      lwd = 4,
      col = "#ffc6c4")

white.hist <- plot(matches_L1_cbps$atc, 
                   col = rgb(1,1,1,0),      
                   # breaks = seq(0, 130, by = 5), 
                   breaks = c(1, 5, 10, 15, 20, 25,
                              30, 35, 40, 45,
                              50, 55, 60, 65, 70, 75, 80,
                              85, 90, 95, 100, 105, 110, 115, 120, 125, 130),
                   plot = FALSE, lty = 1)
lines(white.hist$breaks, c(0, white.hist$counts),
      type = 'S', lty = 1)

legend(x = 0, y = 30,  
       legend = c(paste0("Four Year Lags "), 
                  #                         " empty matched sets")
                  paste0("One Year Lag")),
       
       fill = c("#ffc6c4", "white"),
       
       xjust = 0,
       
       pt.cex = 1,
       bty = "n", ncol = 1, cex = 1, bg = "white")
mtext("Frequency", side = 2, line = 2.2, outer=FALSE,
      cex = 0.8)


hist(L4_cbps_att_sizes_SS[L4_cbps_att_sizes_SS > 0], 
     col = "#ffc6c4", border = NA,
     freq = TRUE,
     
     breaks = seq(1, 20, by = 1),
     
     xlim = c(0, 20), ylim= c(0, 30), xlab = "",
     main = "Starting War",
     ylab = "")
lines(x = c(0,0), 
      y = c(0, length(L4_cbps_att_sizes_SS[L4_cbps_att_sizes_SS == 0])), 
      lwd = 4,
      col = "#ffc6c4")

white.hist <- plot(matches_L1_cbps_SS$att, 
                   col = rgb(1,1,1,0),      
                   breaks = seq(1, 20, by = 1),
                   plot = FALSE, lty = 1)
lines(white.hist$breaks, c(0, white.hist$counts),
      type = 'S', lty = 1)

legend(x = 0, y = 30,  
       legend = c(paste0("Four Year Lags"), 
                  #                         " empty matched sets")
                  paste0("One Year Lag ")),
       
       fill = c("#ffc6c4", "white"),
       
       xjust = 0,
       
       pt.cex = 1,
       bty = "n", ncol = 1, cex = 1, bg = "white")
mtext("Frequency", side = 2, line = 2.2, outer=FALSE,
      cex = 0.8)



hist(L4_cbps_atc_sizes_SS[L4_cbps_atc_sizes_SS > 0], 
     col = "#ffc6c4", border = NA,
     freq = TRUE,
    
     breaks = seq(1, 20, by = 1),
    
     xlim = c(0, 20), ylim= c(0, 30), xlab = "",
     main = "Ending War",
     ylab = "")

white.hist <- hist(L1_cbps_att_sizes_SS[L1_cbps_att_sizes_SS>0], 
                   col = rgb(1,1,1,0), breaks = seq(1, 20, by = 1),
                   plot = FALSE, lty = 1)
lines(white.hist$breaks, c(0, white.hist$counts),
      type = 'S', lty = 1)
lines(x = c(-.1,-.1), 
      y = c(-.1, length(L4_cbps_atc_sizes_SS[L4_cbps_atc_sizes_SS == 0])), 
      lwd = 4,
      col = "#ffc6c4")

lines(x = c(0.3,0.3), 
      y = c(0.3, length(L1_cbps_atc_sizes_SS[L1_cbps_atc_sizes_SS == 0])), 
      lwd = 4,
      col = "black")

legend(x = 1.9, y = 30,  
       legend = c(paste0("Four Year Lags "), 
                  #                         " empty matched sets")
                  paste0("One Year Lag ")),
       # y.intersp = 2,
       fill = c("#ffc6c4", "white"),
       # x.intersp = 0.3,
       xjust = 0,
       # pch = c(22, 22), 
       pt.cex = 1,
       bty = "n", ncol = 1, cex = 1, bg = "white")
mtext("Frequency", side = 2, line = 2.2, outer=FALSE,
      cex = 0.8)
mtext("Acemoglu et al. (2018)", side = 2, line = 0, at = 0.75, outer = TRUE)
mtext("Scheve & Stasavage (2012)", side = 2, line = 0, at = 0.26, outer = TRUE)
mtext("Number of matched control units", side = 1, line = 0.5, at = 0.52, outer = TRUE)
dev.off()



