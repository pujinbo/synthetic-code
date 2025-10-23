###############################################################################
## Replication File for Cattaneo, Feng, Palomba, and Titiunik (2022)
###############################################################################


######################################################################################
# Three predictands: unit-time, time, and unit
######################################################################################

# Load stuff 

pacman::p_load(scpi, haven)
path <- 'YOUR_PATH_HERE'

sims <- 200

path.data <- paste0(path,"data/")
path.fig <- paste0(path,"fig/")

set.seed(8894)

data <- haven::read_dta(paste0(path.data, "final_data.dta"))
data <- subset(data, restricted2 == 1)
data$lgdp <- log(data$rgdppp)
post.est <- 10

units <- df <- est <- res <- list()
df.unit <- est.unit <- res.unit <- list()
df.time <- est.time <- res.time <- list()

features <- list(c("lgdp","lsc"))
covs.adj <- list("lgdp" = c("constant","trend"), "lsc" = c("constant","trend"))

for (cont in c("Africa", "Asia", "Europe", "South America", "North America")) {
  print(paste0("Currently working on ", cont))
  
  print("Estimating whole treatment effect time series")
  
  units[[cont]] <- unique(subset(data, continent == cont & treated == 1 & trDate <= 1992)$countryname)
  if (cont == "Europe") {
    units[[cont]] <- units[[cont]][!(units[[cont]] %in% c("Slovenia"))]
  }
  
  df[[cont]] <- scdataMulti(data, id.var = "countryname", outcome.var = "lgdp", 
                            treatment.var = "liberalization", time.var = "year", constant = TRUE, 
                            cointegrated.data = TRUE, post.est = post.est, units.est = units[[cont]],
                            features=features, cov.adj = covs.adj, anticipation = 1)
  
  res[[cont]] <- scpi(df[[cont]], sims = sims, cores = cores, w.constr = list("name" = "L1-L2"),
                      u.order = 1, u.lags = 1)
  
  # Analysis on average treatment
  print("Estimating average unit treatment effect")
  
  df.unit[[cont]] <- scdataMulti(data, id.var = "countryname", outcome.var = "lgdp", effect = "unit",
                                 treatment.var = "liberalization", time.var = "year", constant = TRUE, 
                                 cointegrated.data = TRUE, post.est = post.est, units.est = units[[cont]],
                                 features = features, cov.adj = covs.adj, anticipation = 1)
  
  res.unit[[cont]] <- scpi(df.unit[[cont]], sims = sims, cores = cores, w.constr = list("name" = "L1-L2"),
                           u.order = 1, u.lags = 1)
  
  
  # Analysis on average treatment on the treated
  print("Estimating average treatment effect")
  # if (cont == "Europe") {
  #   units[[cont]] <- unique(subset(data, continent == cont & treated == 1 & trDate <= 1992)$countryname)
  # }
  df.time[[cont]] <- scdataMulti(data, id.var = "countryname", outcome.var = "lgdp", effect = "time",
                                 treatment.var = "liberalization", time.var = "year", constant = T, 
                                 cointegrated.data = T, post.est = post.est, units.est = units[[cont]],
                                 features=features, cov.adj = covs.adj, anticipation = 1)
  
  res.time[[cont]] <- scpi(df.time[[cont]], sims = sims, cores = cores, w.constr = list("name" = "L1-L2"),
                           u.order = 1, u.lags = 1)
  
}

save.image(file = paste0(path.fig, "workspace_ridge.Rdata"))


######################################################################################
# Fourth predictand - Case study
######################################################################################

data <- haven::read_dta(paste0(path.data, "final_data.dta"))
data <- subset(data, restricted2 == 1)
data$lgdp <- log(data$rgdppp)
post.est <- 10

# One Feature
eu.tr.91 <- unique(subset(data, trDate %in% c(1991) & continent == "Europe")$countryname)

data.co <- subset(data, !(countryname %in% eu.tr.91))
data.tr <- subset(data, countryname %in% eu.tr.91)

data.tr.agg <- aggregate(data.tr[c("lgdp", "liberalization", "lsc")], 
                         by=list(data.tr$year), FUN=mean, na.rm = TRUE)
names(data.tr.agg) <- c("year","lgdp", "liberalization", "lsc")
data.tr.agg$countryname <- "Average Unit"
data.tr.agg$liberalization <- 1*(data.tr.agg$liberalization > 0)

data.agg <- rbind(data.tr.agg, data.co[names(data.tr.agg)])

features <- list(c("lgdp"))
covs.adj <- list(c("constant","trend"))

df.agg <- scdataMulti(data.agg, id.var = "countryname", outcome.var = "lgdp",
                      treatment.var = "liberalization", time.var = "year", constant = FALSE, 
                      cointegrated.data = T, post.est = 10, units.est = "Average Unit", 
                      features=features, cov.adj=covs.adj) 

res.avg.l1l2 <- scpi(df.agg, sims = sims, cores = 1, w.constr = list("name" = "L1-L2"),
                     u.order = 1, u.lags = 1)

res.avg <- scpi(df.agg, sims = sims, cores = 1, w.constr = list("name" = "simplex"),
                u.order = 1, u.lags = 1)

save.image(file = paste0(path.fig, "workspace_casestudy_excluded.Rdata"))


# two features
features <- list(c("lgdp","lsc"))
covs.adj <- list(c("constant","trend"))

df.agg <- scdataMulti(data.agg, id.var = "countryname", outcome.var = "lgdp",
                      treatment.var = "liberalization", time.var = "year", constant = TRUE, 
                      cointegrated.data = T, post.est = 10, units.est = "Average Unit", 
                      features=features, cov.adj=covs.adj) 

res.avg.l1l2 <- scpi(df.agg, sims = sims, cores = 1, w.constr = list("name" = "L1-L2"),
                     u.order = 1, u.lags = 1)

res.avg <- scpi(df.agg, sims = sims, cores = 1, w.constr = list("name" = "simplex"),
                u.order = 1, u.lags = 1)

save.image(file = paste0(path.fig, "workspace_casestudy.Rdata"))


