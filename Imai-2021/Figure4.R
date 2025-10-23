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


##### Making Figure 4, the Histogram #####
## For country name footnotes in the histogram section, see the end of the script ##

library(PanelMatch)
library(ggplot2)
library(gridExtra)

## This script takes the environment saved in Preparing_covariate_balance.R to make 
# scatter plots of covariate balance, allowing for treatment reversal

load("./PanelMatch_temp1.RData")

# setting it up 
den.matrix <- matrix(1:12, nrow = 3, ncol = 4)
den.matrix[1,] <- 1:4
den.matrix[2,] <- 5:8
den.matrix[3,] <- 9:12

xlim <- c(0, .8)
ylim <- c(0, .8)



pdf(file = file.path(OUT_DIR, "Figure4.pdf"),
    width = 9,
    height = 6)
par(oma = c(5, 10, 1.5, 0),
    mar = c(0.8, .9, 1.5, 0.45),
    mfrow = c(3,4),
    pty = "s")

for (i in 1:length(den.matrix)) {
  if (i == 1|i == 5|i == 9) {
    benchmark <- 
      all_nonrestricted$L1_Ace[ substr(names(all_nonrestricted$L1_Ace), 3, 3) == 0 & 
                                  substr(names(all_nonrestricted$L1_Ace), 5, 6) == "no" & 
                                  stringr::str_detect(names(all_nonrestricted$L1_Ace), "TRUE")]
    #all_nonrestricted$L1_Ace[[f]][grepl("Non", names(all_nonrestricted$L1_Ace[[f]]))]
    benchmark <- benchmark[[1]][1:(nrow(benchmark[[1]])-1 ),]
    legend_text <- c("Up to 5 matches",
                     "Up to 10 matches")
  } else if (i == 2|i == 6|i == 10) {
    benchmark <- 
      all_nonrestricted$L4_Ace[ substr(names(all_nonrestricted$L4_Ace), 3, 3) == 0 & 
                                  substr(names(all_nonrestricted$L1_Ace), 5, 6) == "no" & 
                                  stringr::str_detect(names(all_nonrestricted$L4_Ace), "TRUE")]
    
    benchmark <- benchmark[[1]][1:(nrow(benchmark[[1]])-1 ),]
  } else if (i == 3|i == 7| i == 11) {
    benchmark <- 
      all_nonrestricted$L1_SS[ substr(names(all_nonrestricted$L1_SS), 3, 3) == 0 & 
                                 substr(names(all_nonrestricted$L1_SS), 5, 6) == "no" & 
                                 stringr::str_detect(names(all_nonrestricted$L1_SS), "TRUE")]
    benchmark <- benchmark[[1]][1:(nrow(benchmark[[1]])-1 ),]
    legend_text <- c("Up to 1 matches",
                     "Up to 3 matches")
  } else if (i == 4|i == 8| i == 12) {
    benchmark <- 
      all_nonrestricted$L4_SS[ substr(names(all_nonrestricted$L4_SS), 3, 3) == 0 & 
                                 substr(names(all_nonrestricted$L1_SS), 5, 6) == "no" & 
                                 stringr::str_detect(names(all_nonrestricted$L4_SS), "TRUE")]
    benchmark <- benchmark[[1]][1:(nrow(benchmark[[1]])-1 ),]
  }
  if (i < 5) {
    #index <- i
    method <- "mahalanobis"
  } else if (i >= 5 & i <= 8) {
    #index <- i - 4
    method <- "CBPS.match"
  } else {
    #index <- i - 8
    method <- "CBPS.weight"
  }
  
  if (i == 1|i == 3|i == 5|i == 7|i == 9|i == 11) {
    specific_lag <- "L1"
  } else {
    specific_lag <- "L4"
  }
  
  if (i == 1|i == 2|i == 5|i == 6|i ==9|i == 10) {
    paper_name <- "Ace"
  } else {
    paper_name <- "SS"
  }
  
  if (i == 1|i == 3) {
    main_text <- "One Year Lag"
  } else if (i == 2|i == 4) {
    main_text <- "Four Year Lags"
  } else {
    main_text <- ""
  }
  
  sublist <- 
    all_nonrestricted[names(all_nonrestricted) == paste0(specific_lag, "_",
                                                         paper_name)][[1]]
  
  compared <- 
    sublist[stringr::str_detect(names(sublist), method) & 
              grepl("no ", names(sublist)) == FALSE &  
              substr(names(sublist), 3, 3) == 0]
  #all_nonrestricted[[index]][[f]][grepl(method,names(all_nonrestricted[[index]][[f]]))]
  compared <- sapply(compared, function(x) x <- x[1:(nrow(x)-1),])
  # compared <- compared[[1]][1:(nrow(compared[[1]])-1),]
  plot(abs(as.numeric(benchmark)), 
       abs(as.numeric(compared[,1])), pch = 1,
       xlab = "",
       ylab = "",
       xlim = xlim,
       ylim = ylim,
       main = main_text,
       font.main = 1)
  if (i < 9) {
    points(abs(as.numeric(benchmark)), 
           abs(as.numeric(compared[,2])),
           pch = "+")  
  }
  
  abline(h = 0, lty = "dashed")
  abline(0, 1, lty = 2, col = "red")
  if (i == 1|i == 3|i == 5|i == 7) {
    legend(x = 0, y = 0.8,  legend = legend_text,
           y.intersp = 0.65,
           x.intersp = 0.3,
           xjust = 0,
           pch = c(1, 3), pt.cex = 1,
           bty = "n", ncol = 1, cex = 1, bg = "white")
  }
  
  
}

mtext("Acemoglu et al. (2018)", line = 0, at = 0.28, outer = TRUE, cex = 1.2)
mtext("Scheve & Stasavage (2012)", line = 0, at = 0.77, outer = TRUE, cex = 1.2)
# mtext("CBPS Weighting", line = 0, at = 0.9, outer = TRUE, cex = 1.2)
# mtext("ATT", side = 2, line = 1.15, at = 0.71, outer = TRUE, cex = 1.2)
# mtext("ATC", side = 2, line = 1.15, at = 0.21, outer = TRUE, cex = 1.2)
mtext(1,
      text = "Standardized Mean Difference \n before Refinement",
      line = 3.5,
      at = 0.52, outer = TRUE, cex = 1)
mtext(2, text = "Standardized Mean Difference \n After Refinement",
      line = 4, outer = TRUE)

mtext(2, text = "Mahalanobis Distance \n Matching",
      line = .25, at = .82, outer = TRUE,
      cex = .8)

mtext(2, text = "Propensity Score \n Matching",
      line = .25, at = .5, outer = TRUE,
      cex = .8)

mtext(2, text = "Propensity Score \n Weighting",
      line = .25, at = .16, outer = TRUE,
      cex = .8)
dev.off()

