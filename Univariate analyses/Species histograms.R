require(LTMRdata)
require(ggridges)
require(ggplot2)
require(dplyr)
require(tidyr)
require(dtplyr)
require(lubridate)

Data <- LTMRpilot(convert_lengths=TRUE, remove_unconverted_lengths = TRUE, size_cutoff=40)%>%
  zero_fill(species=NULL, remove_unknown_lengths = TRUE, univariate = TRUE)%>%
  filter(Method=="Otter trawl")%>%
  drop_na(Sal_surf, Temp_surf, Count, Latitude, Longitude, Date)%>%
  lazy_dt()%>%
  group_by(Taxa, Source, Station, Latitude, Longitude, Date, Tow_area, SampleID, Tow_area, Sal_surf, Temp_surf)%>%
  summarise(Count=sum(Count), Length=mean(Length, na.rm=TRUE))%>%
  ungroup()%>%
  group_by(Taxa, Source, Station, Latitude, Longitude, Date)%>%
  summarise(Count=sum(Count), Tow_area=sum(Tow_area), Sal_surf=mean(Sal_surf), Temp_surf=mean(Temp_surf), Length=mean(Length, na.rm=TRUE))%>%
  ungroup()%>%
  as_tibble()%>%
  mutate(CPUE=Count/Tow_area,
         Year=year(Date))%>%
  drop_na(CPUE)

Stations<-Data%>%
  group_by(Taxa, Source, Station)%>%
  summarize(Detection_prop=length(which(CPUE>0))/n())%>%
  filter(Detection_prop>1/12)

Data2<-left_join(Data, Stations, by=c("Taxa", "Source", "Station"))%>%
  filter(!is.na(Detection_prop))

Species_data_quality<-Data2%>%
  group_by(Taxa)%>%
  summarise(N_0=length(which(CPUE==0)), N_stations=length(unique(Station)), CPUE_total_sum=sum(CPUE))%>%
  ungroup()%>%
  left_join(LTMRdata::Species%>%select(ScientificName, CommonName)%>%distinct(),by=c("Taxa"="ScientificName"))%>%
  select(Taxa, CommonName, N_0, N_stations, CPUE_total_sum)%>%arrange(CPUE_total_sum, N_stations, N_0)

Data2%>%
  group_by(Taxa)%>%
  summarise(N_0=length(which(CPUE==0)), N_stations=length(unique(Station)), CPUE_total_sum=sum(CPUE))
  
ggplot(Stations, aes(x=Avg_CPUE, y=Taxa))+
    geom_density_ridges()
  
ggplot(Data2, aes(x=CPUE, y=Taxa))+
  geom_density_ridges(rel_min_height = 0.005)+
  scale_x_continuous(expand = c(0.01, 0), limits=c(0,1))
