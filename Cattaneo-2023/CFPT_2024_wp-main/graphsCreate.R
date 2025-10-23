###############################################################################
# Load stuff 
pacman::p_load(ggplot2, Qtools, Hmisc, abind, Matrix, nloptr, reshape2, latex2exp,
               dplyr, tidyr, doSNOW, doRNG, parallel, CVXR, haven, fastDummies, stringr)

path.code <- "/Users/fpalomba/Dropbox/projects/scpi/packages/R/scpi/R"
path <- "/Users/fpalomba/Dropbox/projects/scpi/CFPT_2022_application/"


path.data <- paste0(path,"data/")
path.fig <- paste0(path,"fig/")
path.tab <- paste0(path,"tab/")

dpi <- 300

theme_set(theme_bw())

source(paste0(path.code, "/supporting_functions.R"))
source(paste0(path.code, "/scest.R"))
source(paste0(path.code, "/scpi.R"))
source(paste0(path.code, "/scdata.R"))
source(paste0(path.code, "/scdataMulti.R"))
source(paste0(path.code, "/scplot.R"))
source(paste0(path.code, "/scplotMulti.R"))
source(paste0(path.code, "/scest_methods.R"))
source(paste0(path.code, "/scpi_methods.R"))

###############################################################################
# Case-study - one feature
###############################################################################

load(paste0(path,"fig/case_study/workspace_casestudy.Rdata"))

path.code <- "/Users/fpalomba/Dropbox/projects/scpi/packages/R/scpi/R"
path <- "/Users/fpalomba/Dropbox/projects/scpi/CFPT_2022_application/"

path.data <- paste0(path,"data/")
path.fig <- paste0(path,"fig/case_study/")
path.tab <- paste0(path,"tab/")

theme_set(theme_bw())

############
# L1-L2
############

p <- scplotMulti(res.avg.l1l2, type = "series", joint = TRUE)
pp <- p$plot_out_gau + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.justification = c(0, 1), legend.position = c(0, 1),
        legend.background = element_rect(fill='transparent', colour = NA),
        legend.box.background = element_rect(fill='transparent', colour = NA)) + 
  ylab("(log) GDP per capita (thousand US dollars)") + xlab("year") + ggtitle("")

ggsave(filename = paste0(path.fig, "ave_ridge",".png"), plot = pp,
       width = 10, height = 6, dpi = dpi)

p <- scplotMulti(res.avg.l1l2, type = "treatment", joint = TRUE)
pp <- p$plot_out_gau + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  ylab(TeX("$\\widehat{\\tau}_{\\cdot t, 1991}$")) + xlab("year") + ggtitle("")

ggsave(filename = paste0(path.fig, "ave_ridge_treat",".png"), plot = pp,
       width = 10, height = 6, dpi = dpi)

w <- res.avg.l1l2$est.results$w
aux <- data.frame(w=w,
                  donor=unlist(lapply(strsplit(names(w), "\\."), "[[", 2)),
                  treated=unlist(lapply(strsplit(names(w), "\\."), "[[", 1)))

p <- ggplot() + 
  geom_point(data=aux, aes(x=donor, y=w, size=abs(w))) + xlab("") + ylab("Weight") +
  geom_hline(yintercept=0, linetype="dotted") +
  theme(legend.position = "none",
        axis.text.x=element_text(angle = 45, hjust=1, size=14),
        axis.text.y =element_text(size=13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = paste0(path.fig, "donors_ridge",".png"), plot = p,
       width = 10, height = 6, dpi = dpi)


############
# SIMPLEX
############

p <- scplotMulti(res.avg, type = "series", joint = TRUE)
pp <- p$plot_out_gau + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.justification = c(0, 1), legend.position = c(0, 1),
        legend.background = element_rect(fill='transparent', colour = NA),
        legend.box.background = element_rect(fill='transparent', colour = NA)) + 
  ylab("(log) GDP per capita (thousand US dollars)") + xlab("year") + ggtitle("")

ggsave(filename = paste0(path.fig, "ave_simplex",".png"), plot = pp,
       width = 10, height = 6, dpi = dpi)

p <- scplotMulti(res.avg, type = "treatment", joint = TRUE)
pp <- p$plot_out_gau + 
  scale_size_continuous(range=c(3,4), breaks = c(9,10)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text=element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  ylab(TeX("$\\widehat{\\tau}_{\\cdot t, 1991}$")) + xlab("year") + ggtitle("")

ggsave(filename = paste0(path.fig, "ave_simplex_treat",".png"), plot = pp,
       width = 10, height = 6, dpi = dpi)

w <- res.avg$est.results$w
aux <- data.frame(w=w,
                  donor=unlist(lapply(strsplit(names(w), "\\."), "[[", 2)),
                  treated=unlist(lapply(strsplit(names(w), "\\."), "[[", 1)))

p <- ggplot() + 
  geom_point(data=aux, aes(x=donor, y=w, size=w)) + xlab("") + ylab("Weight") +
  geom_hline(yintercept=0, linetype="dotted") +
  theme(legend.position = "none",
        axis.text.x=element_text(angle = 45, hjust=1, size=14),
        axis.text.y =element_text(size=13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = paste0(path.fig, "donors_simplex",".png"), plot = p,
       width = 10, height = 6, dpi = dpi)


###############################################################################
# Case-study - two features
###############################################################################

load(paste0(path,"fig/case_study/workspace_casestudy_excluded.Rdata"))

path.code <- "/Users/fpalomba/Dropbox/projects/scpi/packages/R/scpi/R"
path <- "/Users/fpalomba/Dropbox/projects/scpi/CFPT_2022_application/"

path.data <- paste0(path,"data/")
path.fig <- paste0(path,"fig/case_study/")
path.tab <- paste0(path,"tab/")

theme_set(theme_bw())

############
# L1-L2
############

p <- scplotMulti(res.avg.l1l2, type = "series", joint = TRUE)
pp <- p$plot_out_gau + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.justification = c(0, 1), legend.position = c(0, 1),
        legend.background = element_rect(fill='transparent', colour = NA),
        legend.box.background = element_rect(fill='transparent', colour = NA)) + 
  ylab("(log) GDP per capita (thousand US dollars)") + xlab("year") + ggtitle("")

ggsave(filename = paste0(path.fig, "excluded_ave_ridge",".png"), plot = pp,
       width = 10, height = 6, dpi = dpi)

p <- scplotMulti(res.avg.l1l2, type = "treatment", joint = TRUE)
pp <- p$plot_out_gau + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  ylab(TeX("$\\widehat{\\tau}_{\\cdot t, 1991}$")) + xlab("year") + ggtitle("")

ggsave(filename = paste0(path.fig, "excluded_ave_ridge_treat",".png"), plot = pp,
       width = 10, height = 6, dpi = dpi)

w <- res.avg.l1l2$est.results$w
aux <- data.frame(w=w,
                  donor=unlist(lapply(strsplit(names(w), "\\."), "[[", 2)),
                  treated=unlist(lapply(strsplit(names(w), "\\."), "[[", 1)))

p <- ggplot() + 
  geom_point(data=aux, aes(x=donor, y=w, size=abs(w))) + xlab("") + ylab("Weight") +
  geom_hline(yintercept=0, linetype="dotted") +
  theme(legend.position = "none",
        axis.text.x=element_text(angle = 45, hjust=1, size=14),
        axis.text.y =element_text(size=13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = paste0(path.fig, "excluded_donors_ridge",".png"), plot = p,
       width = 10, height = 6, dpi = dpi)


############
# SIMPLEX
############

p <- scplotMulti(res.avg, type = "series", joint = TRUE)
pp <- p$plot_out_gau + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.justification = c(0, 1), legend.position = c(0, 1),
        legend.background = element_rect(fill='transparent', colour = NA),
        legend.box.background = element_rect(fill='transparent', colour = NA)) + 
  ylab("(log) GDP per capita (thousand US dollars)") + xlab("year") + ggtitle("")

ggsave(filename = paste0(path.fig, "excluded_ave_simplex",".png"), plot = pp,
       width = 10, height = 6, dpi = dpi)

p <- scplotMulti(res.avg, type = "treatment", joint = TRUE)
pp <- p$plot_out_gau + 
  scale_size_continuous(range=c(3,4), breaks = c(9,10)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text=element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  ylab(TeX("$\\widehat{\\tau}_{\\cdot t, 1991}$")) + xlab("year") + ggtitle("")

ggsave(filename = paste0(path.fig, "excluded_ave_simplex_treat",".png"), plot = pp,
       width = 10, height = 6, dpi = dpi)

w <- res.avg$est.results$w
aux <- data.frame(w=w,
                  donor=unlist(lapply(strsplit(names(w), "\\."), "[[", 2)),
                  treated=unlist(lapply(strsplit(names(w), "\\."), "[[", 1)))

p <- ggplot() + 
  geom_point(data=aux, aes(x=donor, y=w, size=w)) + xlab("") + ylab("Weight") +
  geom_hline(yintercept=0, linetype="dotted") +
  theme(legend.position = "none",
        axis.text.x=element_text(angle = 45, hjust=1, size=14),
        axis.text.y =element_text(size=13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = paste0(path.fig, "excluded_donors_simplex",".png"), plot = p,
       width = 10, height = 6, dpi = dpi)



###############################################################################
# Each continent
###############################################################################

ncols.list <- list('Africa' = 2, 'Europe' = 2, 'Asia' = 2, 'South America' = 2, 'North America' = 2)
axis.text.x.list <- list('Africa' = 10, 'Europe' = 13, 'Asia' = 10, 'South America' = 10, 'North America' = 13)
leg.just.list.y <- list('Africa' = 0.1, 'Europe' = 1, 'Asia' = 0.3, 'South America' = 0.15, 'North America' = 1)
leg.just.list.x <- list('Africa' = 0, 'Europe' = 0, 'Asia' = 0, 'South America' = 0, 'North America' = 0)
leg.text.list <- list('Africa' = 8, 'Europe' = 10, 'Asia' = 10, 'South America' = 10, 'North America' = 10)

for (model in c("simplex", "ridge")) {

  load(paste0(path,"fig/",model,"/workspace_", model,".Rdata"))
  
  path.code <- "/Users/fpalomba/Dropbox/projects/scpi/packages/R/scpi/R"
  path <- "/Users/fpalomba/Dropbox/projects/scpi/CFPT_2022_application/"
  
  cores <- 1

  path.data <- paste0(path,"data/")
  path.fig <- paste0(path,"fig/",model,"/")
  path.tab <- paste0(path,"tab/")
  
  theme_set(theme_bw())
  source(paste0(path.code, "/scplotMulti.R"))
  
  for (cont in c("Africa", "Europe", "South America", "North America", "Asia")) {
    
    #####################################################
    # UNIT-TIME PREDICTAND
    #####################################################
    
    p <- scplotMulti(res[[cont]], type = "series", joint = TRUE, scales = "free_y", ncols = ncols.list[[cont]])
    pp <- p$plot_out_gau + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            axis.text=element_text(size=11),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13),
            legend.key = element_rect(colour = NA, fill = NA),
            legend.text = element_text(size=leg.text.list[[cont]]),
            legend.justification = c(leg.just.list.x[[cont]], leg.just.list.y[[cont]]),
            legend.position = c(leg.just.list.x[[cont]], leg.just.list.y[[cont]]),
            legend.background = element_rect(fill='transparent', colour = NA),
            legend.box.background = element_rect(fill='transparent', colour = NA)) + 
      ylab("(log) GDP per capita (thousand US dollars)") + xlab("year") + ggtitle("")
    
    ggsave(filename = paste0(path.fig, str_replace(cont, " ", "_"), "_series_gau", ".png"), plot = pp,
           width = 8, height = 8, dpi = dpi)
    
    p <- scplotMulti(res[[cont]], type = "treatment", joint = TRUE,  ncols = ncols.list[[cont]])
    pp <- p$plot_out_gau + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.background = element_blank(),
            axis.text=element_text(size=11),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13)) +
      ylab(TeX("$\\widehat{\\tau}_{it}$")) + xlab("year") + ggtitle("")
    
    ggsave(filename = paste0(path.fig, str_replace(cont, " ", "_"), "_treat_gau", ".png"), plot = pp,
           width = 8, height = 8, dpi = dpi)
    
    
    w <- res[[cont]]$est.results$w
    aux <- data.frame(w=w,
                      donor=unlist(lapply(strsplit(names(w), "\\."), "[[", 2)),
                      treated=unlist(lapply(strsplit(names(w), "\\."), "[[", 1)))
    
    p <- ggplot() + 
      geom_point(data=aux, aes(x=donor, y=w, size=abs(w))) + xlab("") + ylab("Weight") +
      geom_hline(yintercept=0, linetype="dotted") +
      facet_wrap(~treated, ncol=ncols.list[[cont]]) +
      theme(legend.position = "none",
            strip.background = element_blank(),
            axis.text.x=element_text(angle = 90, hjust=1, vjust=0.2, size=axis.text.x.list[[cont]]),
            axis.text.y=element_text(size=13),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13),
            panel.grid.minor = element_blank())
    
    ggsave(filename = paste0(path.fig, str_replace(cont, " ", "_"), "_donors",".png"), plot = p,
           width = 10, height = 8, dpi = dpi)
    

    #####################################################
    # UNIT AVERAGE OVER TIME PREDICTAND
    #####################################################
        
    p <- scplotMulti(res.unit[[cont]], type = "series", scales = "free_y", joint = TRUE, ncols = ncols.list[[cont]])
    pp <- p$plot_out_gau + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text=element_text(size=11),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13),
            legend.key = element_rect(colour = NA, fill = NA),
            legend.text = element_text(size=leg.text.list[[cont]]),
            legend.justification = c(leg.just.list.x[[cont]], leg.just.list.y[[cont]]),
            legend.position = c(leg.just.list.x[[cont]], leg.just.list.y[[cont]]),
            legend.background = element_rect(fill='transparent', colour = NA),
            legend.box.background = element_rect(fill='transparent', colour = NA)) + 
      ylab("(log) GDP per capita (thousand US dollars)") + xlab("year") + ggtitle("")
    
    ggsave(filename = paste0(path.fig, str_replace(cont, " ", "_"), "_ate_series_gau",".png"), plot = pp,
           width = 8, height = 8, dpi = dpi)
    
    p <- scplotMulti(res.unit[[cont]], type = "treatment", joint = TRUE, ncols = ncols.list[[cont]])
    pp <- p$plot_out_gau + ggtitle("") + ylab(TeX("$\\widehat{\\tau}_{i\\cdot}$")) +
      scale_size_continuous(range=c(3,4), breaks = c(9,10)) + xlab("") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.y = element_blank(),
            strip.background = element_blank(),
            axis.text=element_text(size=11),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13),
            legend.position = "none")
    
    
    ggsave(filename = paste0(path.fig, str_replace(cont, " ", "_"), "_ate_treat_gau", ".png"), plot = pp,
           width = 8, height = 8, dpi = dpi)
    
    w <- res.unit[[cont]]$est.results$w
    aux <- data.frame(w=w,
                      donor=unlist(lapply(strsplit(names(w), "\\."), "[[", 2)),
                      treated=unlist(lapply(strsplit(names(w), "\\."), "[[", 1)))
    
    p <- ggplot() + 
      geom_point(data=aux, aes(x=donor, y=w, size=abs(w))) + xlab("") + ylab("Weight") +
      geom_hline(yintercept=0, linetype="dotted") +
      facet_wrap(~treated, ncol=ncols.list[[cont]]) +
      theme(legend.position = "none",
            strip.background = element_blank(),
            axis.text.x=element_text(angle = 90, hjust=1, vjust=0.2, size=axis.text.x.list[[cont]]),
            axis.text.y=element_text(size=13),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13),
            panel.grid.minor = element_blank())
    
    ggsave(filename = paste0(path.fig, str_replace(cont, " ", "_"), "_ate_donors",".png"), plot = p,
           width = 10, height = 8, dpi = dpi)
    
    #####################################################
    # AVERAGE OVER UNITS PREDICTAND
    #####################################################
    
    p <- scplotMulti(res.time[[cont]], type = "series", scales = "free_y", joint = TRUE, ncols = ncols.list[[cont]])
    pp <- p$plot_out_gau + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text=element_text(size=11),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13),
            legend.key = element_rect(colour = NA, fill = NA),
            legend.text = element_text(size=leg.text.list[[cont]]),
            legend.justification = c(0, 1),
            legend.position = c(0, 1),
            legend.background = element_rect(fill='transparent', colour = NA),
            legend.box.background = element_rect(fill='transparent', colour = NA)) + 
      ylab("(log) GDP per capita (thousand US dollars)") + xlab("year") + ggtitle("")
    
    ggsave(filename = paste0(path.fig, str_replace(cont, " ", "_"), "_att_series_gau",".png"), plot = pp,
           width = 10, height = 6, dpi = dpi)
    
    p <- scplotMulti(res.time[[cont]], type = "treatment", joint = TRUE, ncols = ncols.list[[cont]])
    pp <- p$plot_out_gau + ggtitle("") + ylab(TeX("$\\widehat{\\tau}_{\\cdot t}$")) +
      scale_size_continuous(range=c(3,4), breaks = c(9,10)) + xlab("") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            strip.background = element_blank(),
            axis.text=element_text(size=11),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13))
    
    
    ggsave(filename = paste0(path.fig, str_replace(cont, " ", "_"), "_att_treat_gau", ".png"), plot = pp,
           width = 10, height = 6, dpi = dpi)
    
    w <- res.time[[cont]]$est.results$w
    aux <- data.frame(w=w,
                      donor=unlist(lapply(strsplit(names(w), "\\."), "[[", 2)),
                      treated=unlist(lapply(strsplit(names(w), "\\."), "[[", 1)))
    
    p <- ggplot() + 
      geom_point(data=aux, aes(x=donor, y=w, size=abs(w))) + xlab("") + ylab("Weight") +
      geom_hline(yintercept=0, linetype="dotted") +
      facet_wrap(~treated, ncol=ncols.list[[cont]]) +
      theme(legend.position = "none",
            strip.background = element_blank(),
            axis.text.x=element_text(angle = 90, hjust=1, vjust=0.2, size=axis.text.x.list[[cont]]),
            axis.text.y=element_text(size=13),
            axis.title.x = element_text(size = 13),
            axis.title.y = element_text(size = 13),
            panel.grid.minor = element_blank())
    
    ggsave(filename = paste0(path.fig, str_replace(cont, " ", "_"), "_att_donors",".png"), plot = p,
           width = 10, height = 8, dpi = dpi)
    
  }
  
}
