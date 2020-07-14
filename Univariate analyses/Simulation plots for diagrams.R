Season_means<-Data%>%
  group_by(Season)%>%
  summarise(mean=mean(Count), .groups="drop")

Data_sim<-expand.grid(Season=unique(Data$Season), Year=unique(Data$Year))%>%
  as_tibble()%>%
  left_join(Season_means, by="Season")%>%
  mutate(Sim_1=rpois(n(), lambda=mean),
         Sim_2=rpois(n(), lambda=mean),
         Sim_3=rpois(n(), lambda=mean),
         Sim_4=rpois(n(), lambda=mean),
         Sim_5=rpois(n(), lambda=mean),
         Sim_6=rpois(n(), lambda=mean),
         Sim_7=rpois(n(), lambda=mean))%>%
  select(-mean)%>%
  mutate(across(starts_with("Sim"), list(sd=sqrt)))%>%
  rename_with(.cols=starts_with("Sim") & !ends_with("_sd"), ~paste0(.x, "_mean"))%>%
  pivot_longer(cols=starts_with("Sim"), names_to=c("Simulation", ".value"), names_sep=6, values_ptypes=list(numeric(), numeric()))%>%
  mutate(Month=case_when(
    Season=="Winter" ~ 1,
    Season=="Spring" ~ 4,
    Season=="Summer" ~ 7,
    Season=="Fall" ~ 10
  ),
  Date=parse_date_time(paste(Month, "1", Year, sep="/"), "%m/%d/%Y"))%>%
  filter(Year>=2010)

p_reduced_sim<-ggplot(filter(Data_sim, Simulation!="Sim_7_"), aes(x=Date, y=mean))+
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), alpha=0.3, fill="darkorange1")+
  geom_line(color="darkorange1")+
  facet_wrap(~Simulation, ncol=2)+
  ylab("Catch")+
  theme_bw()+
  theme(panel.grid=element_blank(), strip.text = element_blank(), strip.background = element_blank(), axis.text=element_blank(), axis.ticks = element_blank())

p_full_sim<-ggplot(filter(Data_sim, Simulation=="Sim_7_"), aes(x=Date, y=mean))+
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), alpha=0.3, fill="dodgerblue3")+
  geom_line(color="dodgerblue3")+
  ylab("Catch")+
  theme_bw()+
  theme(panel.grid=element_blank(), strip.text = element_blank(), strip.background = element_blank(), axis.text=element_blank(), axis.ticks = element_blank())

ggsave(p_reduced_sim, file="Univariate analyses/Figures/Reduced_sim.png", device="png", height=3, width=2, units="in")

ggsave(p_full_sim, file="Univariate analyses/Figures/Full_sim.png", device="png", height=2, width=2, units="in")

Full_overlay<-map_dfr(1:6, ~filter(Data_sim, Simulation=="Sim_7_")%>%mutate(Simulation=paste0("Sim_", .x, "_")))

p_overlay<-ggplot(filter(Data_sim, Simulation!="Sim_7_"), aes(x=Date, y=mean))+
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), alpha=0.3, fill="darkorange1")+
  geom_line(color="darkorange1")+
  geom_ribbon(data=Full_overlay, aes(x=Date, y=mean, ymin=mean-sd, ymax=mean+sd), fill="dodgerblue3", alpha=0.3)+
  geom_line(data=Full_overlay, aes(x=Date, y=mean), color="dodgerblue3")+
  facet_wrap(~Simulation, ncol=2)+
  ylab("Catch")+
  theme_bw()+
  theme(panel.grid=element_blank(), strip.text = element_blank(), strip.background = element_blank(), axis.text=element_blank(), axis.ticks = element_blank())

ggsave(p_overlay, file="Univariate analyses/Figures/Reduced_sim_overlay.png", device="png", height=3, width=2, units="in")
