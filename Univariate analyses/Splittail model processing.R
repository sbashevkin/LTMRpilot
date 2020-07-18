require(brms)
require(dplyr)
require(tidybayes)
require(tibble)
require(purrr)
require(ggplot2)
require(patchwork)
require(stringr)
require(tidyr)
source("Univariate analyses/Survey assessment functions.R")

# Full model -------------------------------------------------------------
load("Univariate analyses/Splittail models/Splittail full model.Rds")

Full_eval<-model_diagnose(model)

Full_change<-Post_processor(model, max_year=2018, model_name="Full")

rm(model)

save(Full_change, file="Univariate analyses/Full model local trend.Rds")

# Missing seasons from month cuts -----------------------------------------
load("Univariate analyses/Split data.Rds")

Months<-Data_split%>%
  filter(Year<=2018)%>%
  group_by(Month_num, Year, Season, .drop=FALSE)%>%
  summarise(N=n(), .groups="drop")%>%
  complete(Month_num=1:3, Year=1985:2018, Season)%>%
  filter(is.na(N))

# Months not represented are:

Season_removals<-Months%>%
  mutate(Cut=2/3)%>%
  rename(Replicate=Month_num)%>%
  select(-N)%>%
  bind_rows(Months%>%
              group_by(Year, Season)%>%
              mutate(N=n())%>%
              ungroup()%>%
              filter(N>1)%>%
              group_by(Year, Season)%>%
              summarise(Replicate=(1:3)[which(!1:3%in%unique(Month_num))], .groups="drop")%>%
              mutate(Cut=1/3))%>%
  mutate(Cut_type="Month",
         Remove=T)

# Reduced models ----------------------------------------------------------

N_station<-c(10,5,3,2, 1.5)
station_reps<-c(10,5,3,2, 3)
N_month<-c(2,1)

Reduced_models<-tibble(Replicate=sequence(station_reps), N_station=rep(N_station, times=station_reps))%>%
  mutate(File=paste0("Splittail ", 1/N_station, " station cut ", Replicate, " of ", if_else(N_station==1.5, 3, N_station), ".Rds"))%>%
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

load("Univariate analyses/Full model local trend.Rds")
Reduced_probs<-map(Reduced_models$File, Reduced_model_processor)
# No Bulk_ESS, Tail_ESS, or Rhat issues from any models

#save(Reduced_probs, file="Univariate analyses/Reduced model proportions.Rds")

Reduced_probs_extracted<-map_dfr(1:nrow(Reduced_models), ~Reduced_probs[[.x]]$Reduced_prob%>%
                                   mutate(N_station=Reduced_models$N_station[.x], 
                                          Replicate=Reduced_models$Replicate[.x],
                                          N_month=Reduced_models$N_month[.x]))%>%
  mutate(Cut_type=if_else(is.na(N_station), "Month", "Station"),
         Cut=if_else(Cut_type=="Month", (3-N_month)/3, 1/N_station))%>%
  select(-N_station, -N_month, -N)%>%
  left_join(Season_removals, by=c("Year", "Season", "Replicate", "Cut", "Cut_type"))%>%
  mutate(Remove=replace_na(Remove, FALSE))%>%
  mutate(across(c(Prob_global, Prob_local), ~if_else(Remove, NA_real_, .x)))%>%
  select(-Remove)%>%
  mutate(Season=factor(Season, levels=c("Winter", "Spring", "Summer", "Fall")))

# How much do the replicates differ?
p_rep<-ggplot(Reduced_probs_extracted, aes(x=Year, y=Prob_local, color=as.factor(Replicate), group=Replicate))+
  geom_line()+
  facet_grid(Cut_type*round(Cut,2)~Season)+
  scale_color_brewer(palette="Paired", aesthetics = c("fill", "color"))+
  coord_cartesian(expand=0, ylim=c(0,1))+
  ylab("Proportional overlap with full model")+
  theme_bw()+
  theme(strip.background=element_blank(), legend.position="none", text=element_text(size=8), panel.spacing.y = unit(0.5, "lines"))

ggsave(p_rep, file="Univariate analyses/Figures/Splittail reduced model replicates.png", device="png", units="in", width=6, height=6)
# Summarise results for each year and season ------------------------------

Reduced_probs_year<-Reduced_probs_extracted%>%
  group_by(Cut_type, Cut, Year, Season)%>%
  summarise(across(c(Prob_global, Prob_local), list(mean=mean, sd=sd), na.rm=T), .groups="drop")

p_ribbon<-ggplot(Reduced_probs_year, aes(x=Year, y=Prob_local_mean, fill=Cut, 
                                         shape=Cut_type, group=interaction(Cut,Cut_type)))+
  geom_ribbon(aes(ymin=Prob_local_mean-Prob_local_sd, ymax=Prob_local_mean+Prob_local_sd), alpha=0.4)+
  geom_line(aes(color=Cut))+
  facet_grid(Cut_type~Season)+
  scale_color_viridis_c(aesthetics = c("fill", "color"), name="Proportion removed", direction = -1,
                        guide=guide_colorbar(barheight=0.5, title.position="top", title.hjust=0.5, direction="horizontal"))+
  coord_cartesian(expand=0, ylim=c(0,1))+
  ylab("Proportional overlap with full model")+
  theme_bw()+
  theme(strip.background=element_blank(), text=element_text(size=8), panel.spacing.y = unit(0.5, "lines"), legend.position = c(0.5, 1.12),
        legend.margin=margin(0,0,0,0), legend.spacing=unit(0, "lines"), plot.margin = margin(30,0,0,0))

ggsave(p_ribbon, file="Univariate analyses/Figures/Splittail reduced model ribbon.png", device="png", units="in", width=6, height=4)
# Summarise results for each season ---------------------------------------

Reduced_probs_season<-Reduced_probs_extracted%>%
  group_by(Cut_type, Cut, Season)%>%
  summarise(across(c(Prob_global, Prob_local), list(mean=mean, sd=sd), na.rm=T), .groups="drop")

p_points<-ggplot(Reduced_probs_season, aes(x=Cut, y=Prob_local_mean, ymin=Prob_local_mean-Prob_local_sd, 
                                           ymax=Prob_local_mean+Prob_local_sd, color=Cut, fill=Cut, 
                                           shape=Cut_type, group=interaction(Cut,Cut_type)))+
  geom_pointrange(color="black", size=0.3, stroke=0.3, position=position_dodge(width=0.05))+
  facet_grid(~Season)+
  geom_hline(yintercept=0.95, linetype=2)+
  scale_color_viridis_c(aesthetics = c("fill", "color"), name="Proportion removed", direction = -1,
                        guide="none")+
  scale_y_continuous(expand=expansion(0,0), limits=c(0,1))+
  scale_shape_manual(values=c(21, 24), name="Cut type", guide=guide_legend(override.aes = list(stroke=0.5, linetype=0, fill="black")))+
  ylab("Proportional overlap with full model")+
  xlab("Proportional reduction in sampling effort")+
  theme_bw()+
  theme(strip.background=element_blank(), text=element_text(size=8), legend.position=c(0.9, 0.2), legend.background=element_rect(color="black"))

ggsave(p_points, file="Univariate analyses/Figures/Splittail reduced model summarized.png", device="png", units="in", width=5, height=3)

# Same plot with global prob
p_points_global<-ggplot(Reduced_probs_season, aes(x=Season, y=Prob_global_mean, ymin=Prob_global_mean-Prob_global_sd, 
                                 ymax=Prob_global_mean+Prob_global_sd, color=Cut, fill=Cut, 
                                 shape=Cut_type, group=interaction(Cut,Cut_type)))+
  geom_pointrange(position=position_dodge(width=0.5), shape=21, color="black", size=0.3, stroke=0.3)+
  facet_grid(Cut_type~.)+
  scale_color_viridis_c(aesthetics = c("fill", "color"), name="Proportion removed", direction = -1,
                        guide=guide_colorbar(barheight=0.5, title.position="top", title.hjust=0.5, direction="horizontal"))+
  coord_cartesian(expand=0, ylim=c(0,1))+
  ylab("Proportional overlap")+
  theme_bw()+
  theme(strip.background=element_blank(), text=element_text(size=8), panel.spacing.y = unit(0.5, "lines"), legend.position = c(0.5, 1.08),
        legend.margin=margin(0,0,0,0), legend.spacing=unit(0, "lines"), plot.margin = margin(40,0,0,0))

ggsave(p_points_global, file="Univariate analyses/Figures/Splittail reduced model summarized global.png", device="png", units="in", width=3, height=4)

# Distribution plots ------------------------------------------------------

Distributional_model_plotter<-function(file, Full_post){
  file2<-str_remove(file, fixed(".Rds"))
  load(file.path("Univariate analyses", "Splittail models", file))
  Data<-model_predictor(model)
  rm(model)
  p<-Distribution_plotter(Full_post, Data, "Change_local")+
    ylab("Local trend")
  ggsave(p, file=paste0("Univariate analyses/Figures/Distribution plots/", file2, " local trend distribution.png"), device="png", units="in", width=4, height=4)
  rm(p)
  
  p<-Distribution_plotter(Full_post, Data, ".value")+
    ylab("Predicted count")
  ggsave(p, file=paste0("Univariate analyses/Figures/Distribution plots/", file2, " predicted count distribution.png"), device="png", units="in", width=4, height=4)
  rm(p)
  
  p<-Distribution_plotter(Full_post, Data, "Change_global")+
    ylab("Global standardized trend")
  ggsave(p, file=paste0("Univariate analyses/Figures/Distribution plots/", file2, " global standardized trend distribution.png"), device="png", units="in", width=4, height=4)
  rm(p)
  
  return(message(paste("Finished:", file)))
}

load("~/LTMRpilot/Univariate analyses/Splittail models/Splittail full model.Rds")
Full_post<-model_predictor(model)

Distributional_model_plotter(Reduced_models$File[1], Full_post)

map(Reduced_models$File, ~Distributional_model_plotter(.x, Full_post))

# Ribbon plots
load(file.path("Univariate analyses", "Splittail models", Reduced_models%>%filter(N_station==2 & Replicate==1)%>%pull(File)))
Data<-model_predictor(model)
rm(model)
p1<-Ribbon_plotter(Full_post, Data, ".value")+ylab("Predicted count")+theme(legend.position=c(0.92, 0.85), legend.background = element_rect(color="black"))
p2<-Ribbon_plotter(Full_post, Data, "Change_local")+ylab("Local trend")+theme(legend.position = "none")

p<-p1/p2+plot_annotation(tag_levels = "A", tag_suffix = ")")
ggsave(p, file=paste0("Univariate analyses/Figures/Distribution plots/", str_remove(Reduced_models%>%filter(N_station==2 & Replicate==1)%>%pull(File), fixed(".Rds")), 
                      " ribbon example.png"), device="png", units="in", width=8, height=8)


# Correlates of proportion overlap ----------------------------------------

Full_sum<-Full_post%>%
  filter(Year!=min(Year))%>%
  group_by(Year, Season)%>%
  summarise(Count_full=mean(.value), Change_local_full=mean(Change_local), .groups="drop")

Data_sum<-Data_split%>%
  group_by(Year, Season)%>%
  summarise(N=n(), .groups="drop")

water_year<-waterYearType::water_year_indices%>%
  select(location, WY, Yr_type, Index)%>%
  pivot_wider(names_from = location, values_from = c(Yr_type, Index))%>%
  rename(Sac_WY_type=`Yr_type_Sacramento Valley`, SJ_WY_type=`Yr_type_San Joaquin Valley`,
         Sac_Index=`Index_Sacramento Valley`, SJ_Index=`Index_San Joaquin Valley`)

Reduced_correlates<-left_join(Reduced_probs_extracted, Full_sum)%>%
  left_join(Data_sum)%>%
  mutate(WY=if_else(Season=="Fall", Year, Year-1))%>%
  left_join(water_year)

ggplot(Reduced_correlates, aes(x=Count_full, y=Prob_local, color=factor(Replicate)))+
  geom_point()+
  facet_grid(Cut_type~Cut)+
  theme_bw()

# No relationships of proportional overlaps with sample size, splittail abundance, local trend, or water year