#Per Capita cigarrette sales (in packs)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
rm(list = ls())
library(matrixStats)
load("ca_comp.rdata")
ca_dat <- mixtape::smoking

years <- unique(ca_dat$year)
t_0   <- min(years) + t_0
t     <- max(years)

# Creating new dataframe

graph_dat <- tibble(
  "year" = years,
  "model" = "BL-TVP",
  "estimate" = rowMedians(y.pred),
  "lower_95" = rowQuantiles(y.pred, probs = .025),
  "upper_95" = rowQuantiles(y.pred, probs = .975),
  "actual" = y
)

graph_dat <- rbind(
  graph_dat,
  tibble(
    "year" = years,
    "model" = "CI",
    "estimate" = cim_pred[,1],
    "lower_95" = cim_pred[,2],
    "upper_95" = cim_pred[,3],
    "actual" = y
  ),
  tibble(
    "year" = years,
    "model" = "CI-TVP",
    "estimate" = cim_pred_tvp[,1],
    "lower_95" = cim_pred_tvp[,2],
    "upper_95" = cim_pred_tvp[,3],
    "actual" = y
  ),
  tibble(
    "year" = years,
    "model" = "SC",
    "estimate" = sc_pred,
    "lower_95" = sc_pred,
    "upper_95" = sc_pred,
    "actual" = y
  ),
  # tibble(
  #   "year" = years,
  #   "model" = "SC-Cattaneo",
  #   "estimate" = c(result$est.results$Y.pre.fit,result$est.results$Y.post.fit),
  #   "lower_95" = c(rep(NA,19),result$inference.results$CI.all.qreg[,1]),
  #   "upper_95" = c(rep(NA,19),result$inference.results$CI.all.qreg[,2]),
  #   "actual" = y
  # ),
  tibble(
    "year"=years,
    "model"="DM-LFM",
    "estimate"=dmlfm$estimated_counterfactual,
    "lower_95"=dmlfm$counterfactual_ci_l,
    "upper_95"=dmlfm$counterfactual_ci_u,
    "actual"=y
  ),
  tibble(
    "year"=years,
    "model"="ArCo",
    "estimate"=c(arco_fit$fitted.values,arco_fit$cf),
    "lower_95"=c(rep(NA,18),lb),
    "upper_95"=c(rep(NA,18),ub),
    "actual"=y
  ),
  tibble(
    "year"=years,
    "model"="BSCM-Horseshoe",
    "estimate"=bayesreg_pred_dat$pred,
    "lower_95"=bayesreg_pred_dat$lower_95,
    "upper_95"=bayesreg_pred_dat$upper_95,
    "actual"=y
  )
)

graph_dat %>%
  mutate(
    estimate = actual - estimate,
    lower_95 = actual - lower_95,
    upper_95 = actual - upper_95,
    actual = actual - actual
  ) %>%
  ggplot(aes(x = year, y = estimate)) +
  geom_line() +
  geom_line(aes(x = year, y = actual), linetype = "dashed") +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = .2) +
  ylab("Per Capita Cigarrette Sales (in Packs)") +
  xlab("Year") +
  geom_vline(xintercept = t_0 - 1, linetype = "dotted") +
  coord_cartesian(ylim = c(-80, 40)) +
  facet_wrap( ~ model,ncol=2) +
  theme_minimal(base_size=14)
  ggsave("../pictures_tables/ca_g2.pdf",
         width = 5, height = 8, dpi = 150, units = "in", device='png')
