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

## This script takes the environment saved in Preparing_covariate_balance.R to make 
# scatter plots of covariate balance, allowing for treatment reversal

load("./PanelMatch_temp1.RData")


max.lead <- 0
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


file = file.path(OUT_DIR, "Figure5.pdf")
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

