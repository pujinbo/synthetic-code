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


## This script takes the environment saved in Getting_PanelEstimates.R
load("./PanelMatch_temp2.RData")

## Scatter plots

# Balance scatterplot allowing for treatment reversal
# setting it up 
den.matrix <- matrix(1:12, nrow = 3, ncol = 4)
den.matrix[1,] <- 1:4
den.matrix[2,] <- 5:8
den.matrix[3,] <- 9:12

xlim <- c(0, .8)
ylim <- c(0, .8)

for (f in 1:5) {
  file = paste0(OUT_DIR, "/scatter_ATT_nonrestricted_listwise_nolagy_", 
                "F", f-1, ".pdf")
  pdf(file = file, 
      width = 9,
      height = 6)
  par(oma = c(5, 10, 1.5, 0),
      mar = c(0.8, .9, 1.5, 0.45),
      mfrow = c(3,4),
      pty = "s")
  
  for (i in 1:length(den.matrix)) {
    if (i == 1|i == 5|i == 9) {
      benchmark <- 
        all_nonrestricted$L1_Ace[ substr(names(all_nonrestricted$L1_Ace), 3, 3) == f-1 & 
                                    substr(names(all_nonrestricted$L1_Ace), 5, 6) == "no" & 
                                    stringr::str_detect(names(all_nonrestricted$L1_Ace), "TRUE")]
      #all_nonrestricted$L1_Ace[[f]][grepl("Non", names(all_nonrestricted$L1_Ace[[f]]))]
      benchmark <- benchmark[[1]][1:(nrow(benchmark[[1]])-1 ),]
      legend_text <- c("Up to 5 matches",
                       "Up to 10 matches")
    } else if (i == 2|i == 6|i == 10) {
      benchmark <- 
        all_nonrestricted$L4_Ace[ substr(names(all_nonrestricted$L4_Ace), 3, 3) == f-1 & 
                                    substr(names(all_nonrestricted$L1_Ace), 5, 6) == "no" & 
                                    stringr::str_detect(names(all_nonrestricted$L4_Ace), "TRUE")]
      
      benchmark <- benchmark[[1]][1:(nrow(benchmark[[1]])-1 ),]
    } else if (i == 3|i == 7| i == 11) {
      benchmark <- 
        all_nonrestricted$L1_SS[ substr(names(all_nonrestricted$L1_SS), 3, 3) == f-1 & 
                                   substr(names(all_nonrestricted$L1_SS), 5, 6) == "no" & 
                                   stringr::str_detect(names(all_nonrestricted$L1_SS), "TRUE")]
      benchmark <- benchmark[[1]][1:(nrow(benchmark[[1]])-1 ),]
      legend_text <- c("Up to 1 matches",
                       "Up to 3 matches")
    } else if (i == 4|i == 8| i == 12) {
      benchmark <- 
        all_nonrestricted$L4_SS[ substr(names(all_nonrestricted$L4_SS), 3, 3) == f-1 & 
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
                substr(names(sublist), 3, 3) == f-1]
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
}

## Balance scatterplot not allowing for treatment reversal
# setting it up 
den.matrix <- matrix(1:12, nrow = 3, ncol = 4)
den.matrix[1,] <- 1:4
den.matrix[2,] <- 5:8
den.matrix[3,] <- 9:12

xlim <- c(0, .8)
ylim <- c(0, .8)

for (f in 1:5) {
  file = paste0(OUT_DIR, "/scatter_ATT_restricted_listwise_nolagy_", 
                "F", f-1, ".pdf")
  pdf(file = file, 
      width = 9,
      height = 6)
  par(oma = c(5, 10, 1.5, 0),
      mar = c(0.8, .9, 1.5, 0.45),
      mfrow = c(3,4),
      pty = "s")
  
  for (i in 1:length(den.matrix)) {
    if (i == 1|i == 5|i == 9) {
      benchmark <- 
        all_restricted$L1_Ace[ substr(names(all_restricted$L1_Ace), 3, 3) == f-1 & 
                                 substr(names(all_restricted$L1_Ace), 5, 6) == "no" & 
                                 stringr::str_detect(names(all_restricted$L1_Ace), "TRUE")]
      #all_restricted$L1_Ace[[f]][grepl("Non", names(all_restricted$L1_Ace[[f]]))]
      benchmark <- benchmark[[1]][1:(nrow(benchmark[[1]])-1 ),]
      legend_text <- c("Up to 5 matches",
                       "Up to 10 matches")
    } else if (i == 2|i == 6|i == 10) {
      benchmark <- 
        all_restricted$L4_Ace[ substr(names(all_restricted$L4_Ace), 3, 3) == f-1 & 
                                 substr(names(all_restricted$L1_Ace), 5, 6) == "no" & 
                                 stringr::str_detect(names(all_restricted$L4_Ace), "TRUE")]
      
      benchmark <- benchmark[[1]][1:(nrow(benchmark[[1]])-1 ),]
    } else if (i == 3|i == 7| i == 11) {
      benchmark <- 
        all_restricted$L1_SS[ substr(names(all_restricted$L1_SS), 3, 3) == f-1 & 
                                substr(names(all_restricted$L1_SS), 5, 6) == "no" & 
                                stringr::str_detect(names(all_restricted$L1_SS), "TRUE")]
      benchmark <- benchmark[[1]][1:(nrow(benchmark[[1]])-1 ),]
      legend_text <- c("Up to 1 matches",
                       "Up to 3 matches")
    } else if (i == 4|i == 8| i == 12) {
      benchmark <- 
        all_restricted$L4_SS[ substr(names(all_restricted$L4_SS), 3, 3) == f-1 & 
                                substr(names(all_restricted$L1_SS), 5, 6) == "no" & 
                                stringr::str_detect(names(all_restricted$L4_SS), "TRUE")]
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
      method <- "CBPS.msm.weight"
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
      all_restricted[names(all_restricted) == paste0(specific_lag, "_",
                                                     paper_name)][[1]]
    
    compared <- 
      sublist[stringr::str_detect(names(sublist), method) & 
                grepl("no ", names(sublist)) == FALSE &  
                substr(names(sublist), 3, 3) == f-1]
    #all_restricted[[index]][[f]][grepl(method,names(all_restricted[[index]][[f]]))]
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
}

## Balance lineplot allowing for treatment reversal
for (f in 1:5) {
  max.lead <- f-1
  # some final preparation  
  all_lines_nonrestricted <- 
    list(all_nonrestricted$L4_Ace[stringr::str_detect(names(all_nonrestricted$L4_Ace), "no ") & 
                                    substr(names(all_nonrestricted$L4_Ace), 3, 3) == max.lead & 
                                    stringr::str_detect(names(all_nonrestricted$L4_Ace), "none")][[1]],
         all_nonrestricted$L4_Ace[stringr::str_detect(names(all_nonrestricted$L4_Ace), "no ") & 
                                    substr(names(all_nonrestricted$L4_Ace), 3, 3) == max.lead & 
                                    stringr::str_detect(names(all_nonrestricted$L4_Ace), "none") == FALSE][[1]],
         all_nonrestricted$L4_Ace[stringr::str_detect(names(all_nonrestricted$L4_Ace), "no ") == FALSE &
                                    stringr::str_detect(names(all_nonrestricted$L4_Ace), "mahalanobis 5") & 
                                    substr(names(all_nonrestricted$L4_Ace), 3, 3) == max.lead][[1]],
         all_nonrestricted$L4_Ace[stringr::str_detect(names(all_nonrestricted$L4_Ace), "no ") == FALSE &
                                    stringr::str_detect(names(all_nonrestricted$L4_Ace), "CBPS.match 5") & 
                                    substr(names(all_nonrestricted$L4_Ace), 3, 3) == max.lead][[1]],      
         all_nonrestricted$L4_Ace[stringr::str_detect(names(all_nonrestricted$L4_Ace), "no ") == FALSE &
                                    stringr::str_detect(names(all_nonrestricted$L4_Ace), "CBPS.weight 10") & 
                                    substr(names(all_nonrestricted$L4_Ace), 3, 3) == max.lead][[1]],
         
         all_nonrestricted_atc$L4_Ace[stringr::str_detect(names(all_nonrestricted$L4_Ace), "no ") & 
                                        substr(names(all_nonrestricted$L4_Ace), 3, 3) == max.lead & 
                                        stringr::str_detect(names(all_nonrestricted$L4_Ace), "none")][[1]],
         all_nonrestricted_atc$L4_Ace[stringr::str_detect(names(all_nonrestricted$L4_Ace), "no ") & 
                                        substr(names(all_nonrestricted$L4_Ace), 3, 3) == max.lead & 
                                        stringr::str_detect(names(all_nonrestricted$L4_Ace), "none") == FALSE][[1]],
         all_nonrestricted_atc$L4_Ace[stringr::str_detect(names(all_nonrestricted$L4_Ace), "no ") == FALSE &
                                        stringr::str_detect(names(all_nonrestricted$L4_Ace), "mahalanobis 5") & 
                                        substr(names(all_nonrestricted$L4_Ace), 3, 3) == max.lead][[1]],
         all_nonrestricted_atc$L4_Ace[stringr::str_detect(names(all_nonrestricted$L4_Ace), "no ") == FALSE &
                                        stringr::str_detect(names(all_nonrestricted$L4_Ace), "CBPS.match 5") & 
                                        substr(names(all_nonrestricted$L4_Ace), 3, 3) == max.lead][[1]],      
         all_nonrestricted_atc$L4_Ace[stringr::str_detect(names(all_nonrestricted$L4_Ace), "no ") == FALSE &
                                        stringr::str_detect(names(all_nonrestricted$L4_Ace), "CBPS.weight 10") & 
                                        substr(names(all_nonrestricted$L4_Ace), 3, 3) == max.lead][[1]],
         
         all_nonrestricted$L4_SS[stringr::str_detect(names(all_nonrestricted$L4_SS), "no ") & 
                                   substr(names(all_nonrestricted$L4_SS), 3, 3) == max.lead & 
                                   stringr::str_detect(names(all_nonrestricted$L4_SS), "none")][[1]],
         all_nonrestricted$L4_SS[stringr::str_detect(names(all_nonrestricted$L4_SS), "no ") & 
                                   substr(names(all_nonrestricted$L4_SS), 3, 3) == max.lead & 
                                   stringr::str_detect(names(all_nonrestricted$L4_SS), "none") == FALSE][[1]],
         all_nonrestricted$L4_SS[stringr::str_detect(names(all_nonrestricted$L4_SS), "no ") == FALSE &
                                   stringr::str_detect(names(all_nonrestricted$L4_SS), "mahalanobis 1") & 
                                   substr(names(all_nonrestricted$L4_SS), 3, 3) == max.lead][[1]],
         all_nonrestricted$L4_SS[stringr::str_detect(names(all_nonrestricted$L4_SS), "no ") == FALSE &
                                   stringr::str_detect(names(all_nonrestricted$L4_SS), "CBPS.match 1") & 
                                   substr(names(all_nonrestricted$L4_SS), 3, 3) == max.lead][[1]],      
         all_nonrestricted$L4_SS[stringr::str_detect(names(all_nonrestricted$L4_SS), "no ") == FALSE &
                                   stringr::str_detect(names(all_nonrestricted$L4_SS), "CBPS.weight 3") & 
                                   substr(names(all_nonrestricted$L4_SS), 3, 3) == max.lead][[1]]
         
    )
  
  
  ### setting it up ###
  plot_matrix <- matrix(1:length(all_lines_nonrestricted), nrow = 3)
  plot_matrix[1, ] <- 1:(length(all_lines_nonrestricted)/3)
  plot_matrix[2, ] <- (length(all_lines_nonrestricted)/3+1):(length(plot_matrix[1, ]) + length(plot_matrix[2, ]))
  plot_matrix[3, ] <- (length(plot_matrix[1, ]) + length(plot_matrix[2, ])+1):
    length(plot_matrix)
  
  
  file = paste0(OUT_DIR, "/Gaps_nonrestricted_listwise_nolagy_F", max.lead, ".pdf")
  pdf(file = file, width = 10, height = 6)
  layout(plot_matrix, heights = c(1,1,1.15))
  
  par(oma = c(2.5, 5, 2.5, 0))
  par(mar = c(1.5, 3.5, 2 , 1))
  for (i in 1:length(all_lines_nonrestricted)) {
    
    if (i ==11) {
      par(mar = c(1.5, 3.5, 4, 1))
    }
    plot(1, 1,
         type = "l", xlim = c(0.5, 4.5), ylim = c(-2, 2.5),
         xlab = "", ylab = "",
         xaxt = "n")
    abline(v = 4, lty = "dashed")
    abline(h = 0, lty = "dotted")
    axis(side = 1, 
         at = 1:4, 
         labels = c(-4, -3,
                    -2, -1))
    for (k in 1:ncol(all_lines_nonrestricted[[i]])) {
      if (i == 1) {
        # mtext("Democratization", side = 2, line = 4.0, outer=FALSE,
        #       cex = 0.8)
        mtext("Standardized Mean \n Differences for \n Democratization", 
              side = 2, line = 2.2, outer=FALSE,
              cex = 0.8)
      } else if (i == 2) {
        mtext(3, line = 0.3, text = "")
      } else if (i == 3) {
        mtext(3, line = 0.3, text = "",cex = 0.8)
      } else if (i == 4) {
        mtext(3, line = 0.3, text = "",cex = 0.8)
      } else if (i == 5) {
        mtext(3, line = 0.3, text = "")
      } else if (i == 6) {
        # mtext("Authoritarian Reversal", side = 2, line = 4.0, outer=FALSE,
        #       cex = 0.8)
        mtext("Standardized Mean \n Differences for \n Authoritarian Reversal", 
              side = 2, line = 2.2, outer=FALSE,
              cex = 0.8)
      } else if (i == 11) {
        # mtext("Starting War", side = 2, line = 4.0, outer=FALSE,
        #       cex = 0.8)
        mtext("Standardized Mean \n Differences for \n Starting War", 
              side = 2, line = 2.2, outer=FALSE,
              cex = 0.8)
      } else if (i == 13) {
        mtext(side = 3, line = .3, text = "",
              cex = 0.8)
      } else if (i == 14) {
        mtext(side = 3, line = .3, text = "",
              cex = 0.8)
      }
      if (k == 1){
        lines(1:4, all_lines_nonrestricted[[i]][,k][1:4],
              col = "black")
      } else {
        lines(1:4, all_lines_nonrestricted[[i]][,k][1:4],
              col = "grey")
      }
    }
  }
  
  mtext("Mahalanobis Distance \n Matching", line = -1.5, at = 0.52, outer = TRUE)
  mtext("Propensity Score \n Matching", line = -1.5, at = 0.72, outer = TRUE)
  mtext("Propensity Score \n Weighting", line = -1.5, at = 0.915, outer = TRUE)
  mtext("Before \n Matching", line = -1.5, at = 0.115, outer = TRUE)
  mtext("Before \n Refinement", line = -1.5, at = 0.317, outer = TRUE)
  mtext(2, text = "Acemoglu et al. (2018)", line = 2.6,
        at = 0.7, outer = TRUE, cex = 1)
  mtext(2, text = "Scheve & Stasavage (2012)", line = 2.6,
        at = 0.15, outer = TRUE, cex = 1)
  mtext(1, text = "Years relative to the administration of treatment", line = 1.2,
        at = 0.52, outer = TRUE, cex = 1)
  dev.off()
}

## Balance lineplot not allowing for treatment reversal
# restricted
for (f in 1:5) {
  max.lead <- f-1
  all_lines_restricted <- 
    list(all_restricted$L4_Ace[stringr::str_detect(names(all_restricted$L4_Ace), "no ") & 
                                 substr(names(all_restricted$L4_Ace), 3, 3) == max.lead & 
                                 stringr::str_detect(names(all_restricted$L4_Ace), "none")][[1]],
         all_restricted$L4_Ace[stringr::str_detect(names(all_restricted$L4_Ace), "no ") & 
                                 substr(names(all_restricted$L4_Ace), 3, 3) == max.lead & 
                                 stringr::str_detect(names(all_restricted$L4_Ace), "none") == FALSE][[1]],
         all_restricted$L4_Ace[stringr::str_detect(names(all_restricted$L4_Ace), "no ") == FALSE &
                                 stringr::str_detect(names(all_restricted$L4_Ace), "mahalanobis 5") & 
                                 substr(names(all_restricted$L4_Ace), 3, 3) == max.lead][[1]],
         all_restricted$L4_Ace[stringr::str_detect(names(all_restricted$L4_Ace), "no ") == FALSE &
                                 stringr::str_detect(names(all_restricted$L4_Ace), "CBPS.match 5") & 
                                 substr(names(all_restricted$L4_Ace), 3, 3) == max.lead][[1]],      
         all_restricted$L4_Ace[stringr::str_detect(names(all_restricted$L4_Ace), "no ") == FALSE &
                                 stringr::str_detect(names(all_restricted$L4_Ace), "CBPS.msm.weight 10") & 
                                 substr(names(all_restricted$L4_Ace), 3, 3) == max.lead][[1]],
         
         all_restricted_atc$L4_Ace[stringr::str_detect(names(all_restricted$L4_Ace), "no ") & 
                                     substr(names(all_restricted$L4_Ace), 3, 3) == max.lead & 
                                     stringr::str_detect(names(all_restricted$L4_Ace), "none")][[1]],
         all_restricted_atc$L4_Ace[stringr::str_detect(names(all_restricted$L4_Ace), "no ") & 
                                     substr(names(all_restricted$L4_Ace), 3, 3) == max.lead & 
                                     stringr::str_detect(names(all_restricted$L4_Ace), "none") == FALSE][[1]],
         all_restricted_atc$L4_Ace[stringr::str_detect(names(all_restricted$L4_Ace), "no ") == FALSE &
                                     stringr::str_detect(names(all_restricted$L4_Ace), "mahalanobis 5") & 
                                     substr(names(all_restricted$L4_Ace), 3, 3) == max.lead][[1]],
         all_restricted_atc$L4_Ace[stringr::str_detect(names(all_restricted$L4_Ace), "no ") == FALSE &
                                     stringr::str_detect(names(all_restricted$L4_Ace), "CBPS.match 5") & 
                                     substr(names(all_restricted$L4_Ace), 3, 3) == max.lead][[1]],      
         all_restricted_atc$L4_Ace[stringr::str_detect(names(all_restricted$L4_Ace), "no ") == FALSE &
                                     stringr::str_detect(names(all_restricted$L4_Ace), "CBPS.msm.weight 10") & 
                                     substr(names(all_restricted$L4_Ace), 3, 3) == max.lead][[1]],
         
         all_restricted$L4_SS[stringr::str_detect(names(all_restricted$L4_SS), "no ") & 
                                substr(names(all_restricted$L4_SS), 3, 3) == max.lead & 
                                stringr::str_detect(names(all_restricted$L4_SS), "none")][[1]],
         all_restricted$L4_SS[stringr::str_detect(names(all_restricted$L4_SS), "no ") & 
                                substr(names(all_restricted$L4_SS), 3, 3) == max.lead & 
                                stringr::str_detect(names(all_restricted$L4_SS), "none") == FALSE][[1]],
         all_restricted$L4_SS[stringr::str_detect(names(all_restricted$L4_SS), "no ") == FALSE &
                                stringr::str_detect(names(all_restricted$L4_SS), "mahalanobis 1") & 
                                substr(names(all_restricted$L4_SS), 3, 3) == max.lead][[1]],
         all_restricted$L4_SS[stringr::str_detect(names(all_restricted$L4_SS), "no ") == FALSE &
                                stringr::str_detect(names(all_restricted$L4_SS), "CBPS.match 1") & 
                                substr(names(all_restricted$L4_SS), 3, 3) == max.lead][[1]],      
         all_restricted$L4_SS[stringr::str_detect(names(all_restricted$L4_SS), "no ") == FALSE &
                                stringr::str_detect(names(all_restricted$L4_SS), "CBPS.msm.weight 3") & 
                                substr(names(all_restricted$L4_SS), 3, 3) == max.lead][[1]]
         
    )
  
  
  ### restricted all gaps ###
  plot_matrix <- matrix(1:length(all_lines_restricted), nrow = 3)
  plot_matrix[1, ] <- 1:(length(all_lines_restricted)/3)
  plot_matrix[2, ] <- (length(all_lines_restricted)/3+1):(length(plot_matrix[1, ]) + length(plot_matrix[2, ]))
  plot_matrix[3, ] <- (length(plot_matrix[1, ]) + length(plot_matrix[2, ])+1):
    length(plot_matrix)
  
  
  file = paste0(OUT_DIR, "/Gaps_restricted_listwise_nolagy_F", max.lead, ".pdf")
  pdf(file = file, width = 10, height = 6)
  layout(plot_matrix, heights = c(1,1,1.15))
  
  par(oma = c(2.5, 5, 2.5, 0))
  par(mar = c(1.5, 3.5, 2 , 1))
  for (i in 1:length(all_lines_restricted)) {
    
    if (i ==11) {
      par(mar = c(1.5, 3.5, 4, 1))
    }
    plot(1, 1,
         type = "l", xlim = c(0.5, 4.5), ylim = c(-2, 2.5),
         xlab = "", ylab = "",
         xaxt = "n")
    abline(v = 4, lty = "dashed")
    abline(h = 0, lty = "dotted")
    axis(side = 1, 
         at = 1:4, 
         labels = c(-4, -3,
                    -2, -1))
    for (k in 1:ncol(all_lines_restricted[[i]])) {
      if (i == 1) {
        # mtext("Democratization", side = 2, line = 4.0, outer=FALSE,
        #       cex = 0.8)
        mtext("Standardized Mean \n Differences for \n Democratization", 
              side = 2, line = 2.2, outer=FALSE,
              cex = 0.8)
      } else if (i == 2) {
        mtext(3, line = 0.3, text = "")
      } else if (i == 3) {
        mtext(3, line = 0.3, text = "",cex = 0.8)
      } else if (i == 4) {
        mtext(3, line = 0.3, text = "",cex = 0.8)
      } else if (i == 5) {
        mtext(3, line = 0.3, text = "")
      } else if (i == 6) {
        # mtext("Authoritarian Reversal", side = 2, line = 4.0, outer=FALSE,
        #       cex = 0.8)
        mtext("Standardized Mean \n Differences for \n Authoritarian Reversal", 
              side = 2, line = 2.2, outer=FALSE,
              cex = 0.8)
      } else if (i == 11) {
        # mtext("Starting War", side = 2, line = 4.0, outer=FALSE,
        #       cex = 0.8)
        mtext("Standardized Mean \n Differences for \n Starting War", 
              side = 2, line = 2.2, outer=FALSE,
              cex = 0.8)
      } else if (i == 13) {
        mtext(side = 3, line = .3, text = "",
              cex = 0.8)
      } else if (i == 14) {
        mtext(side = 3, line = .3, text = "",
              cex = 0.8)
      }
      if (k == 1){
        lines(1:4, all_lines_restricted[[i]][,k][1:4],
              col = "black")
      } else {
        lines(1:4, all_lines_restricted[[i]][,k][1:4],
              col = "grey")
      }
    }
  }
  
  mtext("Mahalanobis Distance \n Matching", line = -1.5, at = 0.52, outer = TRUE)
  mtext("Propensity Score \n Matching", line = -1.5, at = 0.72, outer = TRUE)
  mtext("Propensity Score \n Weighting", line = -1.5, at = 0.915, outer = TRUE)
  mtext("Before \n Matching", line = -1.5, at = 0.115, outer = TRUE)
  mtext("Before \n Refinement", line = -1.5, at = 0.317, outer = TRUE)
  mtext(2, text = "Acemoglu et al. (2018)", line = 2.6,
        at = 0.7, outer = TRUE, cex = 1)
  mtext(2, text = "Scheve & Stasavage (2012)", line = 2.6,
        at = 0.15, outer = TRUE, cex = 1)
  mtext(1, text = "Years relative to the administration of treatment", line = 1.2,
        at = 0.52, outer = TRUE, cex = 1)
  dev.off()
}

## Effect Acemoglu allowing for treatment reversal
## making PDFs
# L = 4 
number_total <- 10
nrow_plot <- 2
result.matrix <- matrix(1:number_total, nrow = nrow_plot, ncol = number_total/nrow_plot)
result.matrix[1,] <- 1:(number_total/nrow_plot)
result.matrix[2,] <- (number_total/nrow_plot + 1):(number_total/nrow_plot + number_total/nrow_plot)
post_periods <- 5
pdf(file = paste0(OUT_DIR, "/A_L4_con_all_listwise_nonrestricted_nolagy.pdf"), width = 10, height = 6)
layout(
  result.matrix
)
specific_lag <- "4"
par(mar = c(1.5, 2, 2 , 1), oma = c(4, 4,1.5, 0))
for (i in 1:10) {
  plot(1, xlab = "", ylab = "", 
       xaxt = "n", yaxt = "n",
       xlim = c(-.5,4.5),
       ylim = c(-.15, .05), pch = 16, cex = 1.4) 
  k <- 1
  if (i == 1|i == 6) {
    method <- "mahalanobis 5"
  } else if (i == 2|i == 7) {
    method <- "mahalanobis 10"
  } else if (i == 3|i == 8) {
    method <- "CBPS.match 5"
  } else if (i == 4|i == 9) {
    method <- "CBPS.match 10"
  } else {
    method <- "CBPS.weight 10"
  }
  # inside each plot
  #for (k in 1:length(all_L4)){
  qoi <- ifelse(i <= 5, "att", "atc")
  if (qoi == "att") {
    tmp_results <- results_Ace_att
  } else {
    tmp_results <- results_Ace_atc
  }
  # index <- ifelse(i<=5, i, i-5)
  one_row <- tmp_results[substr(names(tmp_results), 1,1) == specific_lag & 
                           stringr::str_detect(names(tmp_results), "no ") == FALSE & 
                           stringr::str_detect(names(tmp_results), method)]
  
  
  
  reverse <- ifelse(qoi == "atc", -1, 1)
  P <- sapply(one_row, function(x) x[[1]]) * reverse
  L <- sapply(one_row, function(x) x[[2]]) * reverse
  U <- sapply(one_row, function(x) x[[3]]) * reverse
  if (k == 1){
    pch <- 16
    x_c <- 0:(length(as.numeric(P))-1)
  } else {
    pch <- 25
    x_c <-seq(0.25, 4.25, 1) 
  }
  # pch <- ifelse(k == 1, 16, 25)
  # x_c <- ifelse(k == 1, 0:(length(as.numeric(P))-1), 
  #               seq(0.25, 4.25, 1))
  if (i == 3&k==2|i == 4&k==2|i == 8&k==2|i == 9&k==2) {
    P <- U <- L <- rep(NA, length(P))
  }
  points(x_c, P, pch = pch)
  for (j in x_c) {
    lines(c(j,j), c(L[j+1], U[j+1]))
  }
  abline(h = 0, lty = 3)
  
  if (i == 1) {
    mtext(3, text = "Up to 5 matches", line = 0.17)
    
    # legend(x = .6, y = .2,  
    #        legend = c("Treatment reversal \n allowed",
    #                   "Treatment reversal \n not allowed"),
    #        y.intersp = 2,
    #        # x.intersp = 0.3,
    #        xjust = 0,
    #        pch = c(16, 25), pt.cex = 1,
    #        bty = "n", ncol = 1, cex = 1, bg = "white")
    
  } else if (i == 2){
    mtext(3, text = "Up to 10 matches", line = 0.17)
  } else if (i == 3) {
    mtext(3, text = "Up to 5 matches", line = 0.17)
  } else if (i == 4) {
    mtext(3, text = "Up to 10 matches", line = 0.17)
  }
  
  if (i == 1 | i == 6){
    axis(side = 2, at = seq(-.3, .3, by = .05), 
         labels = c("-.3", "-.25", "-.2", "-.15", "-.1", "-.05", "0", 
                    ".05", ".1", ".15", ".2", ".25", ".3"),
         cex.axis = .8
         # seq(-.3, .3, by = .05))
    )
    # title(ylab = "Estimated Effects", line = 2.1)
    mtext("", side = 2, line = 2, outer=FALSE,
          cex = 0.8)
    # axis(2, at=i,labels=-10:15, col.axis="red", las=2)
    
  }
  axis(side = 1, at = 0:4, labels = 0:4)
  
  
}
mtext("Mahalanobis Matching", line = 0, at = 0.21, outer = TRUE, cex = 1.2)
mtext("Propensity Score Matching", line = 0, at = 0.61, outer = TRUE, cex = 1.2)
mtext("Propensity Score \n Weighting", line = -1.8, at = 0.9, outer = TRUE, cex = 1.2)
mtext("Estimated Effect of \n Democratization", side = 2, line = .3, at = 0.73, 
      outer = TRUE, cex = 1.2)
mtext("Estimated Effect of \n Authoritarian Reversal", side = 2, line = .3, at = 0.22, 
      outer = TRUE, cex = 1.2)
mtext(1, text = "Years relative to the administration of treatment", line = 1.2,
      at = 0.52, outer = TRUE, cex = 1.2)

dev.off()

# L = 1
number_total <- 10
nrow_plot <- 2
result.matrix <- matrix(1:number_total, nrow = nrow_plot, ncol = number_total/nrow_plot)
result.matrix[1,] <- 1:(number_total/nrow_plot)
result.matrix[2,] <- (number_total/nrow_plot + 1):(number_total/nrow_plot + number_total/nrow_plot)
post_periods <- 5
pdf(file = paste0(OUT_DIR, "/A_L1_con_all_listwise_nonrestricted_nolagy.pdf"), width = 10, height = 6)
layout(
  result.matrix
)
specific_lag <- "1"
par(mar = c(1.5, 2, 2 , 1), oma = c(4, 4,1.5, 0))
for (i in 1:10) {
  plot(1, xlab = "", ylab = "", 
       xaxt = "n", yaxt = "n",
       xlim = c(-.5,4.5),
       ylim = c(-.15, .05), pch = 16, cex = 1.4) 
  k <- 1
  if (i == 1|i == 6) {
    method <- "mahalanobis 5"
  } else if (i == 2|i == 7) {
    method <- "mahalanobis 10"
  } else if (i == 3|i == 8) {
    method <- "CBPS.match 5"
  } else if (i == 4|i == 9) {
    method <- "CBPS.match 10"
  } else {
    method <- "CBPS.weight 10"
  }
  # inside each plot
  #for (k in 1:length(all_L4)){
  qoi <- ifelse(i <= 5, "att", "atc")
  if (qoi == "att") {
    tmp_results <- results_Ace_att
  } else {
    tmp_results <- results_Ace_atc
  }
  # index <- ifelse(i<=5, i, i-5)
  one_row <- tmp_results[substr(names(tmp_results), 1,1) == specific_lag & 
                           stringr::str_detect(names(tmp_results), "no ") == FALSE & 
                           stringr::str_detect(names(tmp_results), method)]
  reverse <- ifelse(qoi == "atc", -1, 1)
  P <- sapply(one_row, function(x) x[[1]]) * reverse
  L <- sapply(one_row, function(x) x[[2]]) * reverse
  U <- sapply(one_row, function(x) x[[3]]) * reverse
  if (k == 1){
    pch <- 16
    x_c <- 0:(length(as.numeric(P))-1)
  } else {
    pch <- 25
    x_c <-seq(0.25, 4.25, 1) 
  }
  # pch <- ifelse(k == 1, 16, 25)
  # x_c <- ifelse(k == 1, 0:(length(as.numeric(P))-1), 
  #               seq(0.25, 4.25, 1))
  if (i == 3&k==2|i == 4&k==2|i == 8&k==2|i == 9&k==2) {
    P <- U <- L <- rep(NA, length(P))
  }
  points(x_c, P, pch = pch)
  for (j in x_c) {
    lines(c(j,j), c(L[j+1], U[j+1]))
  }
  abline(h = 0, lty = 3)
  
  if (i == 1) {
    mtext(3, text = "Up to 5 matches", line = 0.17)
    
    # legend(x = .6, y = .2,  
    #        legend = c("Treatment reversal \n allowed",
    #                   "Treatment reversal \n not allowed"),
    #        y.intersp = 2,
    #        # x.intersp = 0.3,
    #        xjust = 0,
    #        pch = c(16, 25), pt.cex = 1,
    #        bty = "n", ncol = 1, cex = 1, bg = "white")
    
  } else if (i == 2){
    mtext(3, text = "Up to 10 matches", line = 0.17)
  } else if (i == 3) {
    mtext(3, text = "Up to 5 matches", line = 0.17)
  } else if (i == 4) {
    mtext(3, text = "Up to 10 matches", line = 0.17)
  }
  
  if (i == 1 | i == 6){
    axis(side = 2, at = seq(-.3, .3, by = .05), 
         labels = c("-.3", "-.25", "-.2", "-.15", "-.1", "-.05", "0", 
                    ".05", ".1", ".15", ".2", ".25", ".3"),
         cex.axis = .8
         # seq(-.3, .3, by = .05))
    )
    # title(ylab = "Estimated Effects", line = 2.1)
    mtext("", side = 2, line = 2, outer=FALSE,
          cex = 0.8)
    # axis(2, at=i,labels=-10:15, col.axis="red", las=2)
    
  }
  axis(side = 1, at = 0:4, labels = 0:4)
  
  
}
mtext("Mahalanobis Matching", line = 0, at = 0.21, outer = TRUE, cex = 1.2)
mtext("Propensity Score Matching", line = 0, at = 0.61, outer = TRUE, cex = 1.2)
mtext("Propensity Score \n Weighting", line = -1.8, at = 0.9, outer = TRUE, cex = 1.2)
mtext("Estimated Effect of \n Democratization", side = 2, line = .3, at = 0.73, 
      outer = TRUE, cex = 1.2)
mtext("Estimated Effect of \n Authoritarian Reversal", side = 2, line = .3, at = 0.22, 
      outer = TRUE, cex = 1.2)
mtext(1, text = "Years relative to the administration of treatment", line = 1.2,
      at = 0.52, outer = TRUE, cex = 1.2)

dev.off()

## Effect Acemoglu not allowing for treatment reversal
## making the PDFs
# L = 4
number_total <- 10
nrow_plot <- 2
result.matrix <- matrix(1:number_total, nrow = nrow_plot, ncol = number_total/nrow_plot)
result.matrix[1,] <- 1:(number_total/nrow_plot)
result.matrix[2,] <- (number_total/nrow_plot + 1):(number_total/nrow_plot + number_total/nrow_plot)
post_periods <- 5
pdf(file = paste0(OUT_DIR, "/A_L4_con_all_listwise_restricted_nolagy.pdf"), width = 10, height = 6)
layout(
  result.matrix
)
specific_lag <- "4"
par(mar = c(1.5, 2, 2 , 1), oma = c(4, 4,1.5, 0))
for (i in 1:10) {
  plot(1, xlab = "", ylab = "", 
       xaxt = "n", yaxt = "n",
       xlim = c(-.5,4.5),
       ylim = c(-.25, .2), pch = 16, cex = 1.4) 
  k <- 2
  if (i == 1|i == 6) {
    method <- "mahalanobis 5"
  } else if (i == 2|i == 7) {
    method <- "mahalanobis 10"
  } else if (i == 3|i == 8) {
    method <- "CBPS.match 5"
  } else if (i == 4|i == 9) {
    method <- "CBPS.match 10"
  } else {
    method <- "CBPS.msm.weight 10"
  }
  # inside each plot
  #for (k in 1:length(all_L4)){
  qoi <- ifelse(i <= 5, "att", "atc")
  if (qoi == "att") {
    tmp_results <- results_Ace_att_restricted
  } else {
    tmp_results <- results_Ace_atc_restricted
  }
  # index <- ifelse(i<=5, i, i-5)
  one_row <- tmp_results[substr(names(tmp_results), 1,1) == specific_lag & 
                           stringr::str_detect(names(tmp_results), "no ") == FALSE & 
                           stringr::str_detect(names(tmp_results), method)]
  
  
  
  reverse <- ifelse(qoi == "atc", -1, 1)
  P <- sapply(one_row, function(x) x[[1]]) * reverse
  L <- sapply(one_row, function(x) x[[2]]) * reverse
  U <- sapply(one_row, function(x) x[[3]]) * reverse
  if (k == 1){
    pch <- 16
    x_c <- 0:(length(as.numeric(P))-1)
  } else {
    pch <- 25
    x_c <-seq(0.25, 4.25, 1) 
  }
  # pch <- ifelse(k == 1, 16, 25)
  # x_c <- ifelse(k == 1, 0:(length(as.numeric(P))-1), 
  #               seq(0.25, 4.25, 1))
  # if (i == 3&k==2|i == 4&k==2|i == 8&k==2|i == 9&k==2) {
  #   P <- U <- L <- rep(NA, length(P))
  # }
  points(x_c, P, pch = pch)
  for (j in x_c) {
    lines(c(j,j), c(L[j+1], U[j+1]))
  }
  abline(h = 0, lty = 3)
  
  if (i == 1) {
    mtext(3, text = "Up to 5 matches", line = 0.17)
    
    # legend(x = .6, y = .2,  
    #        legend = c("Treatment reversal \n allowed",
    #                   "Treatment reversal \n not allowed"),
    #        y.intersp = 2,
    #        # x.intersp = 0.3,
    #        xjust = 0,
    #        pch = c(16, 25), pt.cex = 1,
    #        bty = "n", ncol = 1, cex = 1, bg = "white")
    
  } else if (i == 2){
    mtext(3, text = "Up to 10 matches", line = 0.17)
  } else if (i == 3) {
    mtext(3, text = "Up to 5 matches", line = 0.17)
  } else if (i == 4) {
    mtext(3, text = "Up to 10 matches", line = 0.17)
  }
  
  if (i == 1 | i == 6){
    axis(side = 2, at = seq(-.3, .3, by = .05), 
         labels = c("-.3", "-.25", "-.2", "-.15", "-.1", "-.05", "0", 
                    ".05", ".1", ".15", ".2", ".25", ".3"),
         cex.axis = .8
         # seq(-.3, .3, by = .05))
    )
    # title(ylab = "Estimated Effects", line = 2.1)
    mtext("", side = 2, line = 2, outer=FALSE,
          cex = 0.8)
    # axis(2, at=i,labels=-10:15, col.axis="red", las=2)
    
  }
  axis(side = 1, at = 0:4, labels = 0:4)
  
  
}
mtext("Mahalanobis Matching", line = 0, at = 0.21, outer = TRUE, cex = 1.2)
mtext("Propensity Score Matching", line = 0, at = 0.61, outer = TRUE, cex = 1.2)
mtext("Propensity Score \n Weighting", line = -1.8, at = 0.9, outer = TRUE, cex = 1.2)
mtext("Estimated Effect of \n Democratization", side = 2, line = .3, at = 0.73, 
      outer = TRUE, cex = 1.2)
mtext("Estimated Effect of \n Authoritarian Reversal", side = 2, line = .3, at = 0.22, 
      outer = TRUE, cex = 1.2)
mtext(1, text = "Years relative to the administration of treatment", line = 1.2,
      at = 0.52, outer = TRUE, cex = 1.2)

dev.off()

# L = 1
number_total <- 10
nrow_plot <- 2
result.matrix <- matrix(1:number_total, nrow = nrow_plot, ncol = number_total/nrow_plot)
result.matrix[1,] <- 1:(number_total/nrow_plot)
result.matrix[2,] <- (number_total/nrow_plot + 1):(number_total/nrow_plot + number_total/nrow_plot)
post_periods <- 5
pdf(file = paste0(OUT_DIR, "/A_L1_con_all_listwise_restricted_nolagy.pdf"), width = 10, height = 6)
layout(
  result.matrix
)
specific_lag <- "1"
par(mar = c(1.5, 2, 2 , 1), oma = c(4, 4,1.5, 0))
for (i in 1:10) {
  plot(1, xlab = "", ylab = "", 
       xaxt = "n", yaxt = "n",
       xlim = c(-.5,4.5),
       ylim = c(-.2, .2), pch = 16, cex = 1.4) 
  k <- 2
  if (i == 1|i == 6) {
    method <- "mahalanobis 5"
  } else if (i == 2|i == 7) {
    method <- "mahalanobis 10"
  } else if (i == 3|i == 8) {
    method <- "CBPS.match 5"
  } else if (i == 4|i == 9) {
    method <- "CBPS.match 10"
  } else {
    method <- "CBPS.msm.weight 10"
  }
  # inside each plot
  #for (k in 1:length(all_L4)){
  qoi <- ifelse(i <= 5, "att", "atc")
  if (qoi == "att") {
    tmp_results <- results_Ace_att_restricted
  } else {
    tmp_results <- results_Ace_atc_restricted
  }
  # index <- ifelse(i<=5, i, i-5)
  one_row <- tmp_results[substr(names(tmp_results), 1,1) == specific_lag & 
                           stringr::str_detect(names(tmp_results), "no ") == FALSE & 
                           stringr::str_detect(names(tmp_results), method)]
  
  
  
  reverse <- ifelse(qoi == "atc", -1, 1)
  P <- sapply(one_row, function(x) x[[1]]) * reverse
  L <- sapply(one_row, function(x) x[[2]]) * reverse
  U <- sapply(one_row, function(x) x[[3]]) * reverse
  if (k == 1){
    pch <- 16
    x_c <- 0:(length(as.numeric(P))-1)
  } else {
    pch <- 25
    x_c <-seq(0.25, 4.25, 1) 
  }
  # pch <- ifelse(k == 1, 16, 25)
  # x_c <- ifelse(k == 1, 0:(length(as.numeric(P))-1), 
  #               seq(0.25, 4.25, 1))
  # if (i == 3&k==2|i == 4&k==2|i == 8&k==2|i == 9&k==2) {
  #   P <- U <- L <- rep(NA, length(P))
  # }
  points(x_c, P, pch = pch)
  for (j in x_c) {
    lines(c(j,j), c(L[j+1], U[j+1]))
  }
  abline(h = 0, lty = 3)
  
  if (i == 1) {
    mtext(3, text = "Up to 5 matches", line = 0.17)
    
    # legend(x = .6, y = .2,  
    #        legend = c("Treatment reversal \n allowed",
    #                   "Treatment reversal \n not allowed"),
    #        y.intersp = 2,
    #        # x.intersp = 0.3,
    #        xjust = 0,
    #        pch = c(16, 25), pt.cex = 1,
    #        bty = "n", ncol = 1, cex = 1, bg = "white")
    
  } else if (i == 2){
    mtext(3, text = "Up to 10 matches", line = 0.17)
  } else if (i == 3) {
    mtext(3, text = "Up to 5 matches", line = 0.17)
  } else if (i == 4) {
    mtext(3, text = "Up to 10 matches", line = 0.17)
  }
  
  if (i == 1 | i == 6){
    axis(side = 2, at = seq(-.3, .3, by = .05), 
         labels = c("-.3", "-.25", "-.2", "-.15", "-.1", "-.05", "0", 
                    ".05", ".1", ".15", ".2", ".25", ".3"),
         cex.axis = .8
         # seq(-.3, .3, by = .05))
    )
    mtext("", side = 2, line = 2, outer=FALSE,
          cex = 0.8)
    # axis(2, at=i,labels=-10:15, col.axis="red", las=2)
    
  }
  axis(side = 1, at = 0:4, labels = 0:4)
  
  
}
mtext("Mahalanobis Matching", line = 0, at = 0.21, outer = TRUE, cex = 1.2)
mtext("Propensity Score Matching", line = 0, at = 0.61, outer = TRUE, cex = 1.2)
mtext("Propensity Score \n Weighting", line = -1.8, at = 0.9, outer = TRUE, cex = 1.2)
mtext("Estimated Effect of \n Democratization", side = 2, line = .3, at = 0.73, 
      outer = TRUE, cex = 1.2)
mtext("Estimated Effect of \n Authoritarian Reversal", side = 2, line = .3, at = 0.22, 
      outer = TRUE, cex = 1.2)
mtext(1, text = "Years relative to the administration of treatment", line = 1.2,
      at = 0.52, outer = TRUE, cex = 1.2)

dev.off()

## Effect SS
## Allowing for treatment reversal
k <- 1
tmp_results <- results_SS_att
number_total <- 10
nrow_plot <- 2
result.matrix <- matrix(1:number_total, nrow = nrow_plot, ncol = number_total/nrow_plot)
result.matrix[1,] <- 1:(number_total/nrow_plot)
result.matrix[2,] <- (number_total/nrow_plot + 1):(number_total/nrow_plot + number_total/nrow_plot)
post_periods <- 5
pdf(file = paste0(OUT_DIR, "/SS_con_all_listwise_nonrestricted_nolagy.pdf"), width = 10, height = 6)
layout(
  result.matrix
)

par(mar = c(1.5, 2, 2 , 1), oma = c(4, 4,1.5, 0))
for (i in 1:10) {
  if(i <= 5) {
    specific_lag <- "1"
  } else {
    specific_lag <- "4"
  }
  
  if (i == 1|i == 6) {
    method <- "mahalanobis 1"
  } else if (i == 2|i == 7) {
    method <- "mahalanobis 3"
  } else if (i == 3|i == 8) {
    method <- "CBPS.match 1"
  } else if (i == 4|i == 9) {
    method <- "CBPS.match 3"
  } else {
    method <- "CBPS.weight 3"
  }
  plot(1, xlab = "", ylab = "", 
       xaxt = "n", yaxt = "n",
       xlim = c(-.5,4.5),
       pch = '',
       ylim = c(-4.5, 8.5), cex = 1.4) 
  
  one_row <- tmp_results[substr(names(tmp_results), 1,1) == specific_lag & 
                           stringr::str_detect(names(tmp_results), "no ") == FALSE & 
                           stringr::str_detect(names(tmp_results), method)]
  
  
  
  #reverse <- ifelse(qoi == "atc", -1, 1)
  P <- sapply(one_row, function(x) x[[1]]) 
  L <- sapply(one_row, function(x) x[[2]]) 
  U <- sapply(one_row, function(x) x[[3]]) 
  
  
  if (k == 1){
    pch <- 16
    x_c <- 0:(length(as.numeric(P))-1)
  } else {
    pch <- 25
    x_c <-seq(0.25, 4.25, 1) 
  }
  # pch <- ifelse(k == 1, 16, 25)
  # x_c <- ifelse(k == 1, 0:(length(as.numeric(P))-1), 
  #               seq(0.25, 4.25, 1))
  # if (i == 3&k==2|i == 4&k==2|i == 8&k==2|i == 9&k==2) {
  #   P <- U <- L <- rep(NA, length(P))
  # }
  points(x_c, P, pch = pch)
  for (j in x_c) {
    lines(c(j,j), c(L[j+1], U[j+1]))
  }
  abline(h = 0, lty = 3)
  #}
  if (i == 1) {
    mtext(3, text = "Up to 1 matches", line = 0.17)
    
    # legend(x = -.75, y = 17,  legend = c("Treatment reversal \n allowed",
    #                                      "Treatment reversal \n not allowed"),
    #        y.intersp = 1.5,
    #        # x.intersp = 0.3,
    #        xjust = 0,
    #        pch = c(16, 25), pt.cex = 1,
    #        bty = "n", ncol = 1, cex = 1, bg = "white")
    
  } else if (i == 2){
    mtext(3, text = "Up to 3 matches", line = 0.17)
  } else if (i == 3) {
    mtext(3, text = "Up to 1 matches", line = 0.17)
  } else if (i == 4) {
    mtext(3, text = "Up to 3 matches", line = 0.17)
  }
  
  if (i == 1 | i == 6){
    axis(side = 2, at = c(-8,
                          -4,
                          0,
                          4,
                          8,12,16), labels = c(-8,
                                               -4,
                                               0,
                                               4,
                                               8, 12,16))
    # title(ylab = "Estimated Effects", line = 2.1)
    mtext("", side = 2, line = 2, outer=FALSE,
          cex = 0.8)
    # axis(2, at=i,labels=-10:15, col.axis="red", las=2)
    
  }
  axis(side = 1, at = 0:4, labels = 0:4)
  
  
}
mtext("Mahalanobis Matching", line = 0, at = 0.22, outer = TRUE, cex = 1.2)
mtext("Propensity Score Matching", line = 0, at = 0.62, outer = TRUE, cex = 1.2)
mtext("Propensity Score \n Weighting", line = -1.8, at = 0.9, outer = TRUE, cex = 1.2)
# mtext("L = 1", side = 2, line = 1.15, at = 0.725, outer = TRUE, cex = 1.2)
# mtext("L = 4", side = 2, line = 1.15, at = 0.225, outer = TRUE, cex = 1.2)
mtext("Estimated effect of war", side = 2, line = .285, at = 0.73, outer = TRUE, cex = 1)
mtext("Estimated effect of war", side = 2, line = .285, at = 0.22, outer = TRUE, cex = 1)
mtext("One Year Lag", side = 2, line = 1.85, at = 0.73, outer = TRUE, cex = 1.2)
mtext("Four Year Lags", side = 2, line = 1.85, at = 0.22, outer = TRUE, cex = 1.2)
mtext(1, text = "Years relative to the administration of treatment", line = 1.2,
      at = 0.52, outer = TRUE, cex = 1)
dev.off()

## Not allowing for treatment reversal
k <- 2
tmp_results <- results_SS_att_restricted
number_total <- 10
nrow_plot <- 2
result.matrix <- matrix(1:number_total, nrow = nrow_plot, ncol = number_total/nrow_plot)
result.matrix[1,] <- 1:(number_total/nrow_plot)
result.matrix[2,] <- (number_total/nrow_plot + 1):(number_total/nrow_plot + number_total/nrow_plot)
post_periods <- 5
pdf(file = paste0(OUT_DIR, "/SS_con_all_listwise_restricted_nolagy.pdf"), width = 10, height = 6)
layout(
  result.matrix
)

par(mar = c(1.5, 2, 2 , 1), oma = c(4, 4,1.5, 0))
for (i in 1:10) {
  if(i <= 5) {
    specific_lag <- "1"
  } else {
    specific_lag <- "4"
  }
  
  if (i == 1|i == 6) {
    method <- "mahalanobis 1"
  } else if (i == 2|i == 7) {
    method <- "mahalanobis 3"
  } else if (i == 3|i == 8) {
    method <- "CBPS.match 1"
  } else if (i == 4|i == 9) {
    method <- "CBPS.match 3"
  } else {
    method <- "CBPS.msm.weight 3"
  }
  plot(1, xlab = "", ylab = "", 
       xaxt = "n", yaxt = "n",
       xlim = c(-.5,4.5),
       pch = '',
       ylim = c(-6, 16), cex = 1.4) 
  
  one_row <- tmp_results[substr(names(tmp_results), 1,1) == specific_lag & 
                           stringr::str_detect(names(tmp_results), "no ") == FALSE & 
                           stringr::str_detect(names(tmp_results), method)]
  
  
  
  #reverse <- ifelse(qoi == "atc", -1, 1)
  P <- sapply(one_row, function(x) x[[1]]) 
  L <- sapply(one_row, function(x) x[[2]]) 
  U <- sapply(one_row, function(x) x[[3]]) 
  
  
  if (k == 1){
    pch <- 16
    x_c <- 0:(length(as.numeric(P))-1)
  } else {
    pch <- 25
    x_c <-seq(0.25, 4.25, 1) 
  }
  # pch <- ifelse(k == 1, 16, 25)
  # x_c <- ifelse(k == 1, 0:(length(as.numeric(P))-1), 
  #               seq(0.25, 4.25, 1))
  # if (i == 3&k==2|i == 4&k==2|i == 8&k==2|i == 9&k==2) {
  #   P <- U <- L <- rep(NA, length(P))
  # }
  points(x_c, P, pch = pch)
  for (j in x_c) {
    lines(c(j,j), c(L[j+1], U[j+1]))
  }
  abline(h = 0, lty = 3)
  #}
  if (i == 1) {
    mtext(3, text = "Up to 1 matches", line = 0.17)
    
    # legend(x = -.75, y = 17,  legend = c("Treatment reversal \n allowed",
    #                                      "Treatment reversal \n not allowed"),
    #        y.intersp = 1.5,
    #        # x.intersp = 0.3,
    #        xjust = 0,
    #        pch = c(16, 25), pt.cex = 1,
    #        bty = "n", ncol = 1, cex = 1, bg = "white")
    
  } else if (i == 2){
    mtext(3, text = "Up to 3 matches", line = 0.17)
  } else if (i == 3) {
    mtext(3, text = "Up to 1 matches", line = 0.17)
  } else if (i == 4) {
    mtext(3, text = "Up to 3 matches", line = 0.17)
  }
  
  if (i == 1 | i == 6){
    axis(side = 2, at = c(-8,
                          -4,
                          0,
                          4,
                          8,12,16), labels = c(-8,
                                               -4,
                                               0,
                                               4,
                                               8, 12,16))
    # title(ylab = "Estimated Effects", line = 2.1)
    mtext("", side = 2, line = 2, outer=FALSE,
          cex = 0.8)
    # axis(2, at=i,labels=-10:15, col.axis="red", las=2)
    
  }
  axis(side = 1, at = 0:4, labels = 0:4)
  
  
}
mtext("Mahalanobis Matching", line = 0, at = 0.22, outer = TRUE, cex = 1.2)
mtext("Propensity Score Matching", line = 0, at = 0.62, outer = TRUE, cex = 1.2)
mtext("Propensity Score \n Weighting", line = -1.8, at = 0.9, outer = TRUE, cex = 1.2)
# mtext("L = 1", side = 2, line = 1.15, at = 0.725, outer = TRUE, cex = 1.2)
# mtext("L = 4", side = 2, line = 1.15, at = 0.225, outer = TRUE, cex = 1.2)
mtext("Estimated effect of war", side = 2, line = .285, at = 0.73, outer = TRUE, cex = 1)
mtext("Estimated effect of war", side = 2, line = .285, at = 0.22, outer = TRUE, cex = 1)
mtext("One Year Lag", side = 2, line = 1.85, at = 0.73, outer = TRUE, cex = 1.2)
mtext("Four Year Lags", side = 2, line = 1.85, at = 0.22, outer = TRUE, cex = 1.2)
mtext(1, text = "Years relative to the administration of treatment", line = 1.2,
      at = 0.52, outer = TRUE, cex = 1)
dev.off()


## clean up and rename output file names
system("mv output/scatter_ATT_restricted_listwise_nolagy_F1.pdf output/Section_S4_Figure_S1.pdf")
system("mv output/scatter_ATT_restricted_listwise_nolagy_F2.pdf output/Section_S4_Figure_S2.pdf")
system("mv output/scatter_ATT_restricted_listwise_nolagy_F3.pdf output/Section_S4_Figure_S3.pdf")
system("mv output/scatter_ATT_restricted_listwise_nolagy_F4.pdf output/Section_S4_Figure_S4.pdf")
system("mv output/Gaps_restricted_listwise_nolagy_F1.pdf output/Section_S4_Figure_S5.pdf")
system("mv output/Gaps_restricted_listwise_nolagy_F2.pdf output/Section_S4_Figure_S6.pdf")
system("mv output/Gaps_restricted_listwise_nolagy_F3.pdf output/Section_S4_Figure_S7.pdf")
system("mv output/Gaps_restricted_listwise_nolagy_F4.pdf output/Section_S4_Figure_S8.pdf")
system("mv output/A_L1_con_all_listwise_nonrestricted_nolagy.pdf output/Section_S5_Figure_S1.pdf")
system("rm output/Gaps*")
system("rm output/scatter*")
system("rm output/SS_*")
system("rm output/A_*")


