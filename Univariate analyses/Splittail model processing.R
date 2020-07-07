require(brms)
require(dplyr)
require(tidybayes)
require(tibble)
require(purrr)
require(ggplot2)
source("Univariate analyses/Survey assessment functions.R")

# Full model -------------------------------------------------------------
load("Univariate analyses/Splittail models/Splittail full model.Rds")

Full_eval<-model_diagnose(model)

Full_change<-Post_processor(model, max_year=2018, model_name="Full")

rm(model)


# Missing seasons from month cuts -----------------------------------------
load("Univariate analyses/Split data.Rds")

Months<-Data_split%>%
  filter(Year<=2018)%>%
  group_by(Month_num, Year, Season, .drop=F)%>%
  summarise(N=n(), .groups="drop")%>%
  filter(N<5)

# Reduced models ----------------------------------------------------------

N_station<-c(10,5,3,2)
N_month<-c(2,1)

Reduced_models<-tibble(Replicate=sequence(N_station), N_station=rep(N_station, times=N_station))%>%
  mutate(File=paste0("Splittail ", 1/N_station, " station cut ", Replicate, " of ", N_station, ".Rds"))%>%
  bind_rows(tibble(N_month=rep(N_month, each=3), Replicate=rep(1:3, 2))%>%
              mutate(File=paste0("Splittail ", 3-N_month, " month cut ", Replicate, " of ", 3, ".Rds"))
  )

Reduced_model_processor<- function(file){
  
  load(file.path("Univariate analyses", "Splittail models", file))
  Reduced_eval<-model_diagnose(model)
  
  Reduced_change<-Post_processor(model, max_year=2018, model_name=file, Intervals=Full_change)
  rm(model)
  message(paste("Finished:", file))
  return(list(Reduced_prob=Reduced_change, Errors=Reduced_eval))
}

Reduced_probs<-map(Reduced_models$File, Reduced_model_processor)

#save(Reduced_probs, Full_eval, file="Univariate analyses/Splittail model processing.Rds")

Reduced_probs_extracted<-map_dfr(1:nrow(Reduced_models), ~Reduced_probs[[.x]]$Reduced_prob%>%
                                   mutate(N_station=Reduced_models$N_station[.x], 
                                          Replicate=Reduced_models$Replicate[.x],
                                          N_month=Reduced_models$N_month[.x]))%>%
  mutate(Cut_type=if_else(is.na(N_station), "Month", "Station"),
         Cut=if_else(Cut_type=="Month", (3-N_month)/3, 1/N_station))%>%
  select(-N_station, -N_month, -N)
# No Bulk_ESS, Tail_ESS, or Rhat issues from any models

# How much do the replicates differ?
p_rep<-ggplot(Reduced_probs_extracted, aes(x=Year, y=Prob_local, color=as.factor(Replicate), group=Replicate))+
  geom_line()+
  facet_grid(Cut_type*round(Cut,2)~Season)+
  scale_color_brewer(palette="Paired", aesthetics = c("fill", "color"))+
  coord_cartesian(expand=0, ylim=c(0,1))+
  ylab("Proportional overlap")+
  theme_bw()+
  theme(strip.background=element_blank(), legend.position="none", text=element_text(size=8), panel.spacing.y = unit(0.5, "lines"))

ggsave(p_rep, file="Univariate analyses/Figures/Splittail reduced model replicates.png", device="png", units="in", width=6, height=6)
# Summarise results for each year and season ------------------------------

Reduced_probs_year<-Reduced_probs_extracted%>%
  group_by(Cut_type, Cut, Year, Season)%>%
  summarise(across(c(Prob_global, Prob_local), list(mean=mean, sd=sd)), .groups="drop")

p_ribbon<-ggplot(Reduced_probs_year, aes(x=Year, y=Prob_local_mean, fill=Cut, 
                                         shape=Cut_type, group=interaction(Cut,Cut_type)))+
  geom_ribbon(aes(ymin=Prob_local_mean-Prob_local_sd, ymax=Prob_local_mean+Prob_local_sd), alpha=0.4)+
  geom_line(aes(color=Cut))+
  facet_grid(Cut_type~Season)+
  scale_color_viridis_c(aesthetics = c("fill", "color"), name="Proportion removed", direction = -1,
                        guide=guide_colorbar(barheight=0.5, title.position="top", title.hjust=0.5, direction="horizontal"))+
  coord_cartesian(expand=0, ylim=c(0,1))+
  ylab("Proportional overlap")+
  theme_bw()+
  theme(strip.background=element_blank(), text=element_text(size=8), panel.spacing.y = unit(0.5, "lines"), legend.position = c(0.5, 1.12),
        legend.margin=margin(0,0,0,0), legend.spacing=unit(0, "lines"), plot.margin = margin(30,0,0,0))

ggsave(p_ribbon, file="Univariate analyses/Figures/Splittail reduced model ribbon.png", device="png", units="in", width=6, height=4)
# Summarise results for each season ---------------------------------------

Reduced_probs_season<-Reduced_probs_extracted%>%
  group_by(Cut_type, Cut, Season)%>%
  summarise(across(c(Prob_global, Prob_local), list(mean=mean, sd=sd)), .groups="drop")

p_points<-ggplot(Reduced_probs_season, aes(x=Season, y=Prob_local_mean, ymin=Prob_local_mean-Prob_local_sd, 
                                 ymax=Prob_local_mean+Prob_local_sd, color=Cut, fill=Cut, 
                                 shape=Cut_type, group=interaction(Cut,Cut_type)))+
  geom_pointrange(position=position_dodge(width=0.3), shape=21, color="black", size=0.3, stroke=0.3)+
  facet_grid(Cut_type~.)+
  scale_color_viridis_c(aesthetics = c("fill", "color"), name="Proportion removed", direction = -1,
                        guide=guide_colorbar(barheight=0.5, title.position="top", title.hjust=0.5, direction="horizontal"))+
  coord_cartesian(expand=0, ylim=c(0,1))+
  ylab("Proportional overlap")+
  theme_bw()+
  theme(strip.background=element_blank(), text=element_text(size=8), panel.spacing.y = unit(0.5, "lines"), legend.position = c(0.5, 1.08),
        legend.margin=margin(0,0,0,0), legend.spacing=unit(0, "lines"), plot.margin = margin(30,0,0,0))

ggsave(p_points, file="Univariate analyses/Figures/Splittail reduced model summarized.png", device="png", units="in", width=4, height=4)

# TODO
# 1) For month cuts, remove season x year combos not present in reducd datasets