require(LTMRdata)
require(dplyr)
require(mgcv)
require(lubridate)
require(tidyr)
require(brms)
require(gamm4)
require(ggplot2)
require(purrr)
require(future)
source("Univariate analyses/Survey assessment functions.R")

Data <- LTMRpilot(convert_lengths=TRUE, remove_unconverted_lengths = TRUE, size_cutoff=40)%>%
  zero_fill(species="Leptocottus armatus", remove_unknown_lengths = TRUE, univariate = TRUE)%>%
  filter((Method=="Midwater trawl" | Source=="Suisun") & year(Date)>=1980 & year(Date)<=2019)%>%
  mutate(Tow_volume=if_else(Source=="Suisun", Tow_area*1.5, Tow_volume))%>%
  drop_na(Sal_surf, Temp_surf, Count, Latitude, Longitude, Date, Tow_volume)%>%
  group_by(Source, Station, Latitude, Longitude, Date, Tow_volume, SampleID, Sal_surf, Temp_surf)%>%
  summarise(Count=sum(Count), Length=mean(Length, na.rm=T), .groups="drop")%>%
  group_by(Source, Station, Latitude, Longitude, Date)%>%
  summarise(Count=sum(Count), Tow_volume=sum(Tow_volume), Sal_surf=mean(Sal_surf), 
            Temp_surf=mean(Temp_surf), Length=mean(Length, na.rm=T), .groups="drop")%>%
  mutate(Julian_day=yday(Date),
         Year=year(Date),
         Month=month(Date),
         CPUE=Count/Tow_volume,
         Station_fac=factor(Station))%>%
  mutate(Year=if_else(Month==12, Year+1, Year))%>%
  mutate(Year_fac=ordered(Year),
         Month_fac=ordered(Month),
         ID=1:nrow(.),
         Season=case_when(
           Month%in%c(12,1,2) ~ "Winter",
           Month%in%c(3,4,5) ~ "Spring",
           Month%in%c(6,7,8) ~ "Summer",
           Month%in%c(9,10,11) ~ "Fall"
         ))%>%
  filter(Season=="Fall")%>%
  mutate_at(vars(Length, Longitude, Latitude, Year, Julian_day, Sal_surf, Temp_surf, Tow_volume), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T)))%>% # Create centered and standardized versions of covariates
  arrange(Date)%>%
  filter(Year>=1985)%>% # Remove early years before Suisun survey was regular
  droplevels()


Station_splits<-Data%>%
  select(Station)%>%
  distinct()%>%
  mutate(Group_10=random_groups(12, nrow(.), 10),
         Group_5=random_groups(123, nrow(.), 5),
         Group_3=random_groups(1234, nrow(.), 3),
         Group_2=random_groups(12345, nrow(.), 2),
         Group_2.3=random_groups(123456, nrow(.), 3)) # For 2/3 station cuts

Data_split<-Data%>%
  left_join(Station_splits, by="Station")

save(Data_split, file="Univariate analyses/Sculpin split data.Rds")

# Full model --------------------------------------------------------------

iterations<-5e3
warmup<-iterations/4
plan(multiprocess)

model<-brm(as.integer(round(Count)) ~ Tow_volume_s + Year_fac + (1|Station_fac) + (1|ID),
           family=poisson, data=Data_split,
           prior=prior(normal(0,5), class="Intercept")+
             prior(normal(0,5), class="b")+
             prior(cauchy(0,5), class="sd"),
           chains=3, cores=3,
           iter = iterations, warmup = warmup)
model<-add_criterion(model, c("waic", "loo"))
model<-add_criterion(model, c("kfold"), chains=1)

model2<-brm(as.integer(round(Count)) ~ Tow_volume_s + Year_fac + (1|Station_fac),
           family=poisson, data=Data_split,
           prior=prior(normal(0,5), class="Intercept")+
             prior(normal(0,5), class="b")+
             prior(cauchy(0,5), class="sd"),
           chains=3, cores=3,
           iter = iterations, warmup = warmup)
model2<-add_criterion(model2, c("waic", "loo"))
model2<-add_criterion(model2, c("kfold"), chains=1)

model3<-brm(as.integer(round(Count)) ~ Tow_volume_s + Year_fac + (1|Station_fac),
            family=negbinomial, data=Data_split,
            prior=prior(normal(0,5), class="Intercept")+
              prior(normal(0,5), class="b")+
              prior(cauchy(0,5), class="sd"),
            chains=3, cores=3,
            iter = iterations, warmup = warmup)
model3<-add_criterion(model3, c("waic", "loo"))
model3<-add_criterion(model3, c("kfold"), chains=1)
