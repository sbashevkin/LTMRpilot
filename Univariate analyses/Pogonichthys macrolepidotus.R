require(LTMRdata)
require(dplyr)
require(mgcv)
require(lubridate)
require(tidyr)
require(brms)

Data <- LTMRpilot(convert_lengths=TRUE, remove_unconverted_lengths = TRUE, size_cutoff=40)%>%
  zero_fill(species="Pogonichthys macrolepidotus", remove_unknown_lengths = TRUE, univariate = TRUE)%>%
  group_by(Source, Station, Latitude, Longitude, Date, Datetime, SampleID, Tow_area, Tow_volume, Sal_surf, Temp_surf, Method)%>%
  summarise(Count=sum(Count))%>%
  ungroup()%>%
  drop_na(Sal_surf, Temp_surf, Count, Latitude, Longitude, Date)%>%
  mutate(Julian_day=yday(Date),
         Year=year(Date),
         Month=month(Date))%>%
  mutate_at(vars(Longitude, Latitude, Year, Julian_day, Sal_surf, Temp_surf), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T))) # Create centered and standardized versions of covariates

mtw<-gam(log(Count+1) ~ te(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(20, 10)) + s(Julian_day_s, bs="cc") + Source, 
         family=tw(), data=filter(Data, Method=="Otter trawl"))

mzip<-gam(as.integer(round(Count)) ~ te(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(20, 10)) + s(Julian_day_s, bs="cc") + Source, 
          family=ziP(), data=filter(Data, Method=="Otter trawl"))

mzip2<-gam(list(as.integer(round(Count)) ~ te(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(20, 10), bs="ds") + s(Julian_day_s, bs="cc", m=1) + Source, ~s(Sal_surf_s, m=1)), 
          family=ziplss(), data=filter(Data, Method=="Otter trawl"))

mtw2<-gam(log(Count+1) ~ te(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(20, 10)) + s(Julian_day_s, bs="cc") + Source + s(Sal_surf_s) + s(Temp_surf_s), 
        family=tw(), data=filter(Data, Method=="Otter trawl"))

mtw3<-gam(log(Count+1) ~ te(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(20, 10)) + s(Julian_day_s, bs="cc") + Source + poly(Sal_surf_s, 2) + s(Temp_surf_s), 
          family=tw(), data=filter(Data, Method=="Otter trawl"))

mtw4<-gam(log(Count+1) ~ te(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(20, 10)) + s(Month, bs="cc") + Source + poly(Sal_surf_s, 2) + s(Temp_surf_s), 
          family=tw(), data=filter(Data, Method=="Otter trawl"))

iterations <- 5e3
warmup <- iterations/4
mbayes<- brm(Count ~ t2(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(15, 10)) + s(Julian_day_s, bs="cc") + Source + s(Sal_surf_s) + s(Temp_surf_s),
             family=hurdle_gamma(), data=filter(Data, Method=="Otter trawl"),
             prior=prior(normal(0,10), class="Intercept")+
               prior(normal(0,5), class="b"),
             chains=1, cores=1,
             iter = iterations, warmup = warmup)

mbayes2<- brm(Count ~ t2(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(15, 10)) + s(Julian_day_s, bs="cc") + Source + s(Sal_surf_s) + s(Temp_surf_s),
             family=hurdle_lognormal(), data=filter(Data, Method=="Otter trawl"),
             prior=prior(normal(0,10), class="Intercept")+
               prior(normal(0,5), class="b")+
               prior(cauchy(0,5), class="sigma"),
             chains=1, cores=1,
             iter = iterations, warmup = warmup)

mbayes3<- brm(Count ~ t2(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(15, 10)) + s(Julian_day_s, bs="cc") + Source + s(Sal_surf_s) + s(Temp_surf_s),
              family=hurdle_poisson(), data=filter(Data, Method=="Otter trawl"),
              prior=prior(normal(0,10), class="Intercept")+
                prior(normal(0,5), class="b")+
                prior(cauchy(0,5), class="sigma"),
              chains=1, cores=1,
              iter = iterations, warmup = warmup)

Data_length <- LTMRpilot(convert_lengths=TRUE, remove_unconverted_lengths = TRUE, size_cutoff=40)%>%
  zero_fill(species="Pogonichthys macrolepidotus", remove_unknown_lengths = TRUE, univariate = TRUE)%>%
  drop_na(Sal_surf, Temp_surf, Latitude, Longitude, Date)%>%
  mutate(Julian_day=yday(Date),
         Year=year(Date),
         Month=month(Date))%>%
  mutate_at(vars(Longitude, Latitude, Year, Julian_day, Sal_surf, Temp_surf, Length), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T))) # Create centered and standardized versions of covariates


mbayes3<- brm(bf(Count ~ t2(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(15, 10)) + s(Julian_day_s, bs="cc") + Source + s(Sal_surf_s) + s(Temp_surf_s) + Length_s,
                 hu~s(Sal_surf_s)),
              family=hurdle_poisson(), data=filter(Data, Method=="Otter trawl"),
              prior=prior(normal(0,10), class="Intercept")+
                prior(normal(0,5), class="b")+
                prior(cauchy(0,5), class="sigma"),
              chains=1, cores=1,
              iter = iterations, warmup = warmup)
