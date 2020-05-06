require(LTMRdata)
require(dplyr)
require(mgcv)
require(lubridate)

Data <- LTMRpilot(convert_lengths=TRUE, remove_unconverted_lengths = TRUE, size_cutoff=40)%>%
  zero_fill(species="Pogonichthys macrolepidotus", remove_unknown_lengths = TRUE, univariate = TRUE)%>%
  group_by(Source, Station, Latitude, Longitude, Date, Datetime, SampleID, Tow_area, Tow_volume, Sal_surf, Temp_surf, Method)%>%
  summarise(Count=sum(Count))%>%
  ungroup()%>%
  mutate(Julian_day=yday(Date),
         Year=year(Date))

m<-gam(log(Count+1) ~ te(Latitude, Longitude, Year, d=c(2,1)) + s(Julian_day, bs="cc"), family=tw(), data=filter(Data, Method=="Otter trawl"))
