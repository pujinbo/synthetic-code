# Set working directory to source file location
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(magrittr)
library(tidysynth)
library(ggpubr)
library(readr)
notvp_placebo<-read_csv("notvp_placebo_1000.csv")

ca_dat <- mixtape::smoking
ca_data<-tidysynth::smoking
years <- unique(ca_dat$year)
t_0<-1989
state<-  "tennessee"
state_fancy<- "Tennessee"
sim_num<-1

notvp_placebo %<>% 
  mutate(period=period+1969,
         model=ifelse(model=="BLTVP","BL-TVP",model)) 
# General Placebo test here
(ex<-notvp_placebo %>% 
  filter(unit==state & sim==sim_num) %>% 
  group_by(model,period) %>% 
  summarise(sim_y=mean(sim_y),
            pred_y=mean(prediction),
            lower_95=mean(lower_95, probs=.025),
            upper_95=mean(upper_95, probs=.975)
            ) %>% 
  ggplot(aes(x=period,y=pred_y))+
  geom_line(lty="dashed",col="black")+
  geom_ribbon(aes(ymin=lower_95,ymax=upper_95), alpha=.2)+
  geom_line(aes(x=period,y=sim_y),col="black")+
  geom_vline(xintercept=1988, lty="dotted")+
  labs(x="Year",
       y="Cigarettes per Capita")+
  theme_minimal(base_size=24)+
    coord_cartesian(ylim = c(40, 150)) +
  facet_wrap(~model) 
)
ggsave(paste0("graphs/",state, "_example.pdf"),
         width = 12,
         height = 6)

# new graph idea ----------------------------------------------------------

bltvp<-notvp_placebo %>% 
  filter(model=="BL-TVP") %>% 
  dplyr::select(period,sim,unit,"bltvp"=prediction,
                "bltvp_upper"=upper_95,
                "bltvp_lower"=lower_95)


# cred ratio to msfe ration
notvp_placebo %>% 
  left_join(bltvp) %>% 
  filter(
    model!="BL-TVP") %>% 
  filter(period > t_0) %>%
  group_by(model,unit) %>%
  summarise(
    cred_spread_bltvp = mean(bltvp_upper - bltvp_lower, na.rm=T),
    cred_spread = mean(upper_95 - lower_95, na.rm=T),
    msfe_bltvp=mean((bltvp-sim_y)^2,na.rm=T),
    msfe=mean((prediction-sim_y)^2,na.rm=T)
  ) %>% 
  mutate(cred_ratio=log(cred_spread/cred_spread_bltvp),
         msfe_ratio=log(msfe/msfe_bltvp),
         state_graph=ifelse(state=="Tennessee","Tennessee","NA")
  ) %>% 
  mutate(cred_ratio=ifelse(model!="SC",cred_ratio,0)) %>% 
  ggplot(aes(x=cred_ratio, y=msfe_ratio))+
  geom_point(size=3)+
  ylim(-2.5,2.5)+
  xlim(-2.5,2.5)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  facet_wrap(~model,ncol=2)+
  theme_minimal(base_size=24)+
  theme(legend.position = "none")+
  scale_color_manual(values=c("black","red"))+
  labs(y="log(MSFE Alternative Model)-log(MSFE BL-TVP)",
       x="log(CI Spread Alternative Model)-log(CI Spread BL-TVP)")
ggsave(paste0("graphs/", "log_ratio.pdf"),
       width = 12,
       height = 12)

# Tables ------------------------------------------------------------------


msfe_table<-notvp_placebo %>% 
  left_join(bltvp) %>% 
  filter(period > t_0) %>%
  group_by(model,unit) %>%
  summarise(
    msfe = mean((sim_y - prediction) ^ 2,na.rm=T)
  ) %>% 
  pivot_wider(names_from = model, values_from = msfe) %>% 
  janitor::clean_names()
  

cred_spread_table<-notvp_placebo %>% 
  left_join(bltvp) %>% 
  filter(period > t_0) %>%
  group_by(model,unit) %>%
  summarise(
    cred_spread = mean(upper_95 - lower_95, na.rm=T)
  ) %>% 
  pivot_wider(names_from = model, values_from = cred_spread) %>% 
  janitor::clean_names() %>% 
  mutate(sc=NA)

means<-cred_spread_table %>% 
  rename_at(vars(-unit), ~ paste0("cred_", .x)) %>% 
  left_join(msfe_table %>% 
              rename_at(vars(-unit), ~ paste0("msfe_", .x))) %>% 
  dplyr::select(-unit) %>% 
  summarise(across(everything(),mean,na.rm=T)
  )

medians<-cred_spread_table %>% 
  rename_at(vars(-unit), ~ paste0("cred_", .x)) %>% 
  left_join(msfe_table %>% 
              rename_at(vars(-unit), ~ paste0("msfe_", .x))) %>% 
  dplyr::select(-unit) %>% 
  summarise(across(everything(),median,na.rm=T)
  )


table<-cred_spread_table %>%  
  rename_at(vars(-unit), ~ paste0("cred_", .x)) %>% 
  left_join(msfe_table %>% 
              rename_at(vars(-unit), ~ paste0("msfe_", .x))) %>%
  left_join(ca_data %>%
              dplyr::select(state) %>%
              mutate(unit=snakecase::to_snake_case(state))
  ) %>%
  dplyr::select(-unit) %>% 
  relocate(state,.before=1) %>%
  rbind(.,c("state"="Average",means)) %>%
  rbind(.,c("state"="Median",medians))%>%
  distinct()

View(table)

table %>% 
  write_csv("sim_table.csv")
