require(LTMRdata)
require(dplyr)
require(mgcv)
require(lubridate)
require(tidyr)
require(brms)
require(gamm4)
require(ggplot2)
require(purrr)
source("Univariate analyses/Survey assessment functions.R")

Data <- LTMRpilot(convert_lengths=TRUE, remove_unconverted_lengths = TRUE, size_cutoff=40)%>%
  zero_fill(species="Pogonichthys macrolepidotus", remove_unknown_lengths = TRUE, univariate = TRUE)%>%
  filter(Method=="Otter trawl" & year(Date)>=1980 & year(Date)<=2019)%>%
  drop_na(Sal_surf, Temp_surf, Count, Latitude, Longitude, Date, Tow_area)%>%
  group_by(Source, Station, Latitude, Longitude, Date, Tow_area, SampleID, Sal_surf, Temp_surf)%>%
  summarise(Count=sum(Count), Length=mean(Length, na.rm=T), .groups="drop")%>%
  group_by(Source, Station, Latitude, Longitude, Date)%>%
  summarise(Count=sum(Count), Tow_area=sum(Tow_area), Sal_surf=mean(Sal_surf), 
            Temp_surf=mean(Temp_surf), Length=mean(Length, na.rm=T), .groups="drop")%>%
  mutate(Julian_day=yday(Date),
         Year=year(Date),
         Month=month(Date),
         CPUE=Count/Tow_area,
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
  mutate_at(vars(Length, Longitude, Latitude, Year, Julian_day, Sal_surf, Temp_surf, Tow_area), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T)))%>% # Create centered and standardized versions of covariates
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
  left_join(Station_splits, by="Station")%>%
  mutate(Month_num=case_when(
    Month%in%c(12,3,6,9) ~ 1,
    Month%in%c(1,4,7,10) ~ 2,
    Month%in%c(2,5,8,11) ~ 3,
  ))

save(Data_split, file="Univariate analyses/Split data.Rds")