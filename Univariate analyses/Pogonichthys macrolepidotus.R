require(LTMRdata)
require(dplyr)
require(mgcv)
require(lubridate)
require(tidyr)
require(brms)
require(gamm4)
require(ggplot2)
require(purrr)

Data <- LTMRpilot(convert_lengths=TRUE, remove_unconverted_lengths = TRUE, size_cutoff=40)%>%
  zero_fill(species="Pogonichthys macrolepidotus", remove_unknown_lengths = TRUE, univariate = TRUE)%>%
  filter(Method=="Otter trawl")%>%
  drop_na(Sal_surf, Temp_surf, Count, Datetime, Latitude, Longitude, Date)%>%
  group_by(Source, Station, Latitude, Longitude, Date, Tow_area, SampleID, Tow_area, Sal_surf, Temp_surf)%>%
  summarise(Count=sum(Count), Length=mean(Length, na.rm=T))%>%
  ungroup()%>%
  group_by(Source, Station, Latitude, Longitude, Date)%>%
  summarise(Count=sum(Count), Tow_area=sum(Tow_area), Sal_surf=mean(Sal_surf), Temp_surf=mean(Temp_surf), Length=mean(Length, na.rm=T))%>%
  ungroup()%>%
  mutate(Julian_day=yday(Date),
         Year=year(Date),
         Month=month(Date),
         CPUE=Count/Tow_area,
         Station_fac=factor(Station),
         Year_fac=ordered(Year),
         Month_fac=ordered(Month))%>%
  mutate_at(vars(Length, Longitude, Latitude, Year, Julian_day, Sal_surf, Temp_surf, Tow_area), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T)))%>% # Create centered and standardized versions of covariates
  arrange(Date)

Stations<-Data%>%
  group_by(Source, Year, Station)%>%
  summarize(Avg_CPUE=mean(CPUE))%>%
  ungroup()%>%
  group_by(Source, Station)%>%
  summarize(Avg_CPUE=mean(Avg_CPUE))%>%
  ungroup()%>%
  filter(Avg_CPUE>0.005)

Data2<-Data%>%
  filter(Station%in%unique(Stations$Station))%>%
  droplevels()


# Trying to get a model to fit well ---------------------------------------

# log(CPUE+1), numeric year, tw family
m2 <- bam(log(CPUE+1) ~ s(Year_s) + s(Julian_day_s) + s(Station_fac, bs="re"), family=tw, data=Data2, method="fREML", discrete=T, nthreads=4)
m2b <- bam(log(CPUE+1) ~ s(Year_s, k=30) + s(Julian_day_s) + s(Station_fac, bs="re"), family=tw, data=Data2, method="fREML", discrete=T, nthreads=4)

# CPUE, numeric year, tw family
m22 <- bam(CPUE ~ s(Year) + s(Julian_day) + s(Station_fac, bs="re"), family=tw, data=Data2, method="fREML", discrete=T, nthreads=4)

# Count, Tow_area, numeric year, tw family
m23 <- bam(Count ~ Tow_area_s + s(Year_s) + s(Julian_day_s) + s(Station_fac, bs="re"), family=tw, data=Data2, method="fREML", discrete=T, nthreads=4)

# log(Count+1), Tow_area, numeric year, tw family
m24 <- bam(log(Count+1) ~ Tow_area + s(Year) + s(Julian_day) + s(Station_fac, bs="re"), family=tw, data=Data2, method="fREML", discrete=T, nthreads=4)

# round(Count), Tow_area, numeric year, negative binomial family
m25 <- bam(as.integer(round(Count)) ~ Tow_area_s + s(Year_s) + s(Julian_day_s) + s(Station_fac, bs="re"), family=nb, data=Data2, method="fREML", discrete=T, nthreads=4)

# round(Count), Tow_area, numeric year, poisson family
m26 <- bam(as.integer(round(Count)) ~ Tow_area_s + s(Year_s) + s(Julian_day_s) + s(Station_fac, bs="re"), family=poisson, data=Data2, method="fREML", discrete=T, nthreads=4)

# round(Count), Tow_area, numeric year, quasipoisson family
m27 <- bam(as.integer(round(Count)) ~ Tow_area_s + s(Year_s) + s(Julian_day_s) + s(Station_fac, bs="re"), family=quasipoisson, data=Data2, method="fREML", discrete=T, nthreads=4)

# round(Count), Tow_area, numeric year, zero inflated poisson family
m28 <- bam(as.integer(round(Count)) ~ Tow_area_s + s(Year_s) + s(Julian_day_s) + s(Station_fac, bs="re"), family=ziP, data=Data2, method="fREML", discrete=T, nthreads=4)

# log(CPUE+1), factor year, tw family
m29 <- bam(log(CPUE+1) ~ Year_fac + s(Julian_day_s) + s(Station_fac, bs="re"), family=tw, data=Data2, method="fREML", discrete=T, nthreads=4)

# CPUE, factor year, tw family
m210 <- bam(CPUE ~ Year_fac + s(Julian_day_s) + s(Station_fac, bs="re"), family=tw, data=Data2, method="fREML", discrete=T, nthreads=4)

# log(Count+1), Tow_area, factor year, tw family
m211 <- bam(log(Count+1) ~ Tow_area_s + Year_fac + s(Julian_day_s) + s(Station_fac, bs="re"), family=tw, data=Data2, method="fREML", discrete=T, nthreads=4)

# Count, Tow_area, factor year, tw family
m212 <- bam(Count ~ Tow_area_s + Year_fac + s(Julian_day_s) + s(Station_fac, bs="re"), family=tw, data=Data2, method="fREML", discrete=T, nthreads=4)

# round(Count), Tow_area, factor year, negative binomial family
m213 <- bam(as.integer(round(Count)) ~ Tow_area_s + Year_fac + s(Julian_day_s) + s(Station_fac, bs="re"), family=nb, data=Data2, method="fREML", discrete=T, nthreads=4)

# round(Count), Tow_area, factor year, poisson family
m214 <- bam(as.integer(round(Count)) ~ Tow_area_s + Year_fac + s(Julian_day_s) + s(Station_fac, bs="re"), family=poisson, data=Data2, method="fREML", discrete=T, nthreads=4)

# round(Count), Tow_area, factor year, quasipoisson family
m215 <- bam(as.integer(round(Count)) ~ Tow_area_s + Year_fac + s(Julian_day_s) + s(Station_fac, bs="re"), family=quasipoisson, data=Data2, method="fREML", discrete=T, nthreads=4)

# round(Count), Tow_area, factor year, zero inflated poisson family
m216 <- bam(as.integer(round(Count)) ~ Tow_area_s + Year_fac + s(Julian_day_s) + s(Station_fac, bs="re"), family=ziP, data=Data2, method="fREML", discrete=T, nthreads=4)

# Count, Tow_area, Spatiotemporal smoothers, tw family
m217 <- bam(Count ~ Tow_area_s + te(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(20, 10)) + s(Julian_day_s, bs="cc"), family=tw, data=Data2, method="fREML", discrete=T, nthreads=4)

# log(CPUE+1), numeric year, tw family, temporal autocorrelation
m218 <- gamm(log(CPUE+1) ~ s(Year_s) + s(Julian_day_s), random =list(Station_fac=~1), family=tw, correlation = corCAR1(form = ~ Date|Station), data=Data2)

# lmer
m219 <- lmer(CPUE^(1/3) ~ Year_fac + Month_fac + (1|Station), data=Data2)

# CPUE^(1/3), numeric year, tw family
m220<- bam(CPUE^(1/3) ~ s(Year_s) + s(Julian_day_s) + s(Station_fac, bs="re"), family=tw, data=Data2, method="fREML", discrete=F, nthreads=4)

# round(Count), Tow_area, numeric year, average length, zero-inflated poisson with salinity predicting presence
m221<- gam(list(as.integer(round(Count)) ~ Tow_area_s + s(Year_s) + s(Julian_day_s) + s(Length_s) + s(Station_fac, bs="re"), ~s(Sal_surf_s)), family=ziplss, data=mutate(Data2, Length_s=replace_na(Length_s, 0)), method="REML")

# CPUE^(1/3), numeric year, year and julian_day by station, tw family
m222<- bam(CPUE^(1/3) ~ s(Year_s, by=Station_fac) + s(Julian_day_s, by=Station_fac) + s(Station_fac, bs="re"), family=tw, data=Data2, method="fREML", discrete=T, nthreads=4)

# CPUE^(1/3), factor year and station with interaction, julian_day by station, tw family
m223<- bam(CPUE^(1/3) ~ Year_fac*Station_fac + s(Julian_day_s, by=Station_fac), family=tw, data=Data2, method="fREML", discrete=T, nthreads=4)

# CPUE^(1/3), 3D smooth with julian day, factor year, and factor station, tw family
m224 <- bam(CPUE^(1/3) ~ te(Julian_day_s, Year_fac, Station_fac, bs=c("cc", "fs", "fs")), family=tw, data=Data2, method="fREML", discrete=T, nthreads=4)

# CPUE^(1/3), 3D smooth with julian day, factor year, and factor station, tw family
m225 <- bam(CPUE^(1/3) ~ te(Julian_day_s, Year_fac, Station_fac, bs=c("cc", "re", "re")), family=tw, data=Data2, method="fREML", discrete=T, nthreads=4)

fit_plot<-function(model, X, title, subtitle, type="fit"){
  if(type=="fit"){
    fit<-fitted(model, type="response")
  }
  if(type=="predict"){
    fit<-predict(model, type="response")
  }
  p<-ggplot(data=Data2)+
    geom_point(aes(x=X, y=fit), alpha=0.2)+
    geom_abline(intercept = 0, slope=1, color="red", size=1)+
    ggtitle(label=title, subtitle = subtitle)+
    xlab("Response")+
    ylab("Model predicted fitted value")+
    coord_fixed()+
    theme_bw()
  
  return(p)
}



p<-pmap(list(model=list(m2, m22, m23, m24, m25, m26, m27, m28, m29, m210, m211, m212, m213, m214, m215, m216, m217, m218$gam, m219, m220, m221, m222, m223, m224),
             X=list(log(Data2$CPUE+1), 
                    Data2$CPUE, 
                    Data2$Count, 
                    log(Data2$Count+1), 
                    as.integer(round(Data2$Count)), 
                    as.integer(round(Data2$Count)), 
                    as.integer(round(Data2$Count)), 
                    as.integer(round(Data2$Count)), 
                    log(Data2$CPUE+1), 
                    Data2$CPUE, 
                    log(Data2$Count+1), 
                    Data2$Count, 
                    as.integer(round(Data2$Count)), 
                    as.integer(round(Data2$Count)), 
                    as.integer(round(Data2$Count)), 
                    as.integer(round(Data2$Count)),
                    Data2$Count,
                    log(Data2$CPUE+1),
                    Data2$CPUE^(1/3),
                    Data2$CPUE^(1/3),
                    as.integer(round(Data2$Count)),
                    Data2$CPUE^(1/3),
                    Data2$CPUE^(1/3),
                    Data2$CPUE^(1/3)),
             title=list('log(CPUE+1), numeric year, tw family', 
                        'CPUE, numeric year, tw family', 
                        'Count, Tow_area, numeric year, tw family',
                        'log(Count+1), Tow_area, numeric year, tw family', 
                        'round(Count), Tow_area, numeric year, negative binomial family',
                        'round(Count), Tow_area, numeric year, poisson family', 
                        'round(Count), Tow_area, numeric year, quasipoisson family',
                        'round(Count), Tow_area, numeric year, zero inflated poisson family', 
                        'log(CPUE+1), factor year, tw family',
                        'CPUE, factor year, tw family', 
                        'log(Count+1), Tow_area, factor year, tw family', 
                        'Count, Tow_area, factor year, tw family',
                        'round(Count), Tow_area, factor year, negative binomial family', 
                        'round(Count), Tow_area, factor year, poisson family',
                        'round(Count), Tow_area, factor year, quasipoisson family', 
                        'round(Count), Tow_area, factor year, zero inflated poisson family',
                        'Count, Tow_area, Spatiotemporal smoothers, tw family',
                        'log(CPUE+1), numeric year, tw family, temporal autocorrelation',
                        'Simple linear mixed effects model',
                        'CPUE^(1/3), numeric year, tw family',
                        'round(Count), Tow_area, numeric year, average length, zero-inflated poisson with salinity predicting presence',
                        'CPUE^(1/3), numeric year, year and julian_day by station, tw family',
                        'CPUE^(1/3), factor year and station with interaction, julian_day by station, tw family',
                        'CPUE^(1/3), 3D smooth with julian day, factor year, and factor station, tw family'),
             subtitle=c('bam(log(CPUE+1) ~ s(Year_s) + s(Julian_day_s) + s(Station_fac, bs="re"), family=tw)',
                        'bam(CPUE ~ s(Year) + s(Julian_day) + s(Station_fac, bs="re"), family=tw)',
                        'bam(Count ~ Tow_area_s + s(Year_s) + s(Julian_day_s) + s(Station_fac, bs="re"), family=tw)',
                        'bam(log(Count+1) ~ Tow_area + s(Year) + s(Julian_day) + s(Station_fac, bs="re"), family=tw)',
                        'bam(as.integer(round(Count)) ~ Tow_area_s + s(Year_s) + s(Julian_day_s) + s(Station_fac, bs="re"), family=nb)',
                        'bam(as.integer(round(Count)) ~ Tow_area_s + s(Year_s) + s(Julian_day_s) + s(Station_fac, bs="re"), family=poisson)',
                        'bam(as.integer(round(Count)) ~ Tow_area_s + s(Year_s) + s(Julian_day_s) + s(Station_fac, bs="re"), family=quasipoisson)',
                        'bam(as.integer(round(Count)) ~ Tow_area_s + s(Year_s) + s(Julian_day_s) + s(Station_fac, bs="re"), family=ziP)',
                        'bam(log(CPUE+1) ~ Year_fac + s(Julian_day_s) + s(Station_fac, bs="re"), family=tw)',
                        'bam(CPUE ~ Year_fac + s(Julian_day_s) + s(Station_fac, bs="re"), family=tw)',
                        'bam(log(Count+1) ~ Tow_area_s + Year_fac + s(Julian_day_s) + s(Station_fac, bs="re"), family=tw)',
                        'bam(Count ~ Tow_area_s + Year_fac + s(Julian_day_s) + s(Station_fac, bs="re"), family=tw)',
                        'bam(as.integer(round(Count)) ~ Tow_area_s + Year_fac + s(Julian_day_s) + s(Station_fac, bs="re"), family=nb)',
                        'bam(as.integer(round(Count)) ~ Tow_area_s + Year_fac + s(Julian_day_s) + s(Station_fac, bs="re"), family=poisson)',
                        'bam(as.integer(round(Count)) ~ Tow_area_s + Year_fac + s(Julian_day_s) + s(Station_fac, bs="re"), family=quasipoisson)',
                        'bam(as.integer(round(Count)) ~ Tow_area_s + Year_fac + s(Julian_day_s) + s(Station_fac, bs="re"), family=ziP)',
                        'bam(Count ~ Tow_area_s + te(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(20, 10)) + s(Julian_day_s, bs="cc"), family=tw',
                        'gamm(log(CPUE+1) ~ s(Year_s) + s(Julian_day_s), random =list(Station_fac=~1), family=tw, correlation = corCAR1(form = ~ Date|Station))',
                        'lmer(CPUE^(1/3) ~ Year_fac + Month_fac + (1|Station))',
                        'bam(CPUE^(1/3) ~ s(Year_s) + s(Julian_day_s) + s(Station_fac, bs="re"), family=tw)',
                        'gam(list(as.integer(round(Count)) ~ Tow_area_s + s(Year_s) + s(Julian_day_s) + s(Length_s) + s(Station_fac, bs="re"), ~s(Sal_surf_s)), family=ziplss)',
                        'bam(CPUE^(1/3) ~ s(Year_s, by=Station_fac) + s(Julian_day_s, by=Station_fac) + s(Station_fac, bs="re"), family=tw)',
                        'bam(CPUE^(1/3) ~ Year_fac*Station_fac + s(Julian_day_s, by=Station_fac), family=tw)',
                        'bam(CPUE^(1/3) ~ te(Julian_day_s, Year_fac, Station_fac, bs=c("cc", "fs", "fs")), family=tw)'),
             type=c(rep("fit", 20), "predict", rep("fit", 3))),
        fit_plot)


walk2(p, paste("model", 1:24), ~ggsave(filename = paste0("Univariate analyses/Figures/", .y, ".png"), plot = .x, device = "png", width=15, height=6, units="in"))





# Including length --------------------------------------------------------

Data_length <- LTMRpilot(convert_lengths=TRUE, remove_unconverted_lengths = TRUE, size_cutoff=40)%>%
  zero_fill(species="Pogonichthys macrolepidotus", remove_unknown_lengths = TRUE, univariate = TRUE)%>%
  filter(Method=="Otter trawl")%>%
  drop_na(Sal_surf, Temp_surf, Latitude, Longitude, Date, Tow_area)%>%
  mutate(Julian_day=yday(Date),
         Year=year(Date),
         Month=month(Date),
         CPUE=Count/Tow_area,
         Station_fac=factor(Station),
         Year_fac=ordered(Year),
         Month_fac=ordered(Month))%>%
  mutate_at(vars(Length, Longitude, Latitude, Year, Julian_day, Sal_surf, Temp_surf, Tow_area), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T)))%>% # Create centered and standardized versions of covariates
  mutate(Length_s=replace_na(Length_s, 0))%>%
  filter(Station%in%Stations$Station)


# Length
m221<- gam(list(as.integer(round(Count)) ~ s(Year_s) + s(Julian_day_s, bs="cc") + s(Length_s, k=100) + s(Station_fac, bs="re"), ~s(Sal_surf_s)), family=ziplss, data=Data_length, method="REML")




# Old models --------------------------------------------------------------


mtw<-gam(log(Count+1) ~ Tow_area_s + te(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(20, 10)) + s(Julian_day_s, bs="cc") + Source, 
         family=tw(), data=filter(Data, Method=="Otter trawl"))

mzip<-gam(as.integer(round(Count)) ~ Tow_area_s + te(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(20, 10)) + s(Julian_day_s, bs="cc") + Source, 
          family=ziP(), data=filter(Data, Method=="Otter trawl"))

mzip2<-gam(list(as.integer(round(Count)) ~ Tow_area_s + te(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(20, 10), bs="ds") + s(Julian_day_s, bs="cc", m=1) + Source, ~s(Sal_surf_s, m=1)), 
           family=ziplss(), data=filter(Data, Method=="Otter trawl"))

mnb<-bam(as.integer(round(Count)) ~ Tow_area_s + te(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(20, 10), bs="ds") + s(Julian_day_s, bs="cc", m=1) + Source, 
         family=nb(), data=filter(Data, Method=="Otter trawl"), method="fREML", discrete=T, nthreads=4)

mnb2<-bam(as.integer(round(Count)) ~ Tow_area_s + te(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(20, 10), bs="ds") + s(Julian_day_s, bs="cc", m=1) + Source, 
          family=nb(), data=filter(Data, Method=="Otter trawl" & Count>0), method="fREML", discrete=T, nthreads=4)

mtw2<-gam(log(Count+1) ~ Tow_area_s + te(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(20, 10)) + s(Julian_day_s, bs="cc") + Source + s(Sal_surf_s) + s(Temp_surf_s), 
          family=tw(), data=filter(Data, Method=="Otter trawl"))

mtw3<-gam(log(Count+1) ~ Tow_area_s + te(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(20, 10)) + s(Julian_day_s, bs="cc") + Source + poly(Sal_surf_s, 2) + s(Temp_surf_s), 
          family=tw(), data=filter(Data, Method=="Otter trawl"))

mtw4<-gam(log(Count+1) ~ Tow_area_s + te(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(20, 10)) + s(Month, bs="cc") + Source + s(Sal_surf_s) + s(Temp_surf_s), 
          family=tw(), data=filter(Data, Method=="Otter trawl"))

# Remove stations where the species has never been seen
# Make Year a factor
# Try other species that are less "boom or bust"
# Try simpler models like gamm(CPUE ~ Year + (1|Station))

# Take splittail, eliminate stations where they're caught less than 15 times per year, normalize data, then build simple gamm(CPUE ~ Year + Source+ 1|Station)

iterations <- 5e3
warmup <- iterations/4
mbhg<- brm(Count ~ Tow_area_s + t2(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(15, 10)) + s(Julian_day_s, bs="cc") + Source + s(Sal_surf_s) + s(Temp_surf_s),
           family=hurdle_gamma(), data=filter(Data, Method=="Otter trawl"),
           prior=prior(normal(0,10), class="Intercept")+
             prior(normal(0,5), class="b"),
           chains=1, cores=1, control=list(adapt_delta=0.9),
           iter = iterations, warmup = warmup)

mbhg2<- brm(CPUE ~ t2(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(15, 10)) + s(Julian_day_s, bs="cc") + Source,
            family=hurdle_gamma(), data=filter(Data, Method=="Otter trawl"),
            prior=prior(normal(0,10), class="Intercept")+
              prior(normal(0,5), class="b"),
            chains=1, cores=1, control=list(adapt_delta=0.9),
            iter = iterations, warmup = warmup)

mbhln<- brm(Count ~ Tow_area_s + t2(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(15, 10)) + s(Julian_day_s, bs="cc") + Source + s(Sal_surf_s) + s(Temp_surf_s),
            family=hurdle_lognormal(), data=filter(Data, Method=="Otter trawl"),
            prior=prior(normal(0,10), class="Intercept")+
              prior(normal(0,5), class="b")+
              prior(cauchy(0,5), class="sigma"),
            chains=1, cores=1,
            iter = iterations, warmup = warmup)

mbnb<- brm(as.integer(round(Count)) ~ Tow_area_s + t2(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(15, 10)) + s(Julian_day_s, bs="cc") + Source + s(Sal_surf_s) + s(Temp_surf_s),
           family=hurdle_negbinomial(), data=filter(Data, Method=="Otter trawl"),
           prior=prior(normal(0,10), class="Intercept")+
             prior(normal(0,5), class="b"),
           chains=1, cores=1,
           iter = iterations, warmup = warmup)

## Bayesian length model

mbayes4<- brm(bf(as.integer(round(Count)) ~ Tow_area_s + t2(Latitude_s, Longitude_s, Year_s, d=c(2,1), k=c(15, 10)) + s(Julian_day_s, bs="cc") + Source + poly(Sal_surf_s, 2) + s(Temp_surf_s) + Length_s,
                 hu~poly(Sal_surf_s, 2)),
              family=hurdle_poisson(), data=filter(Data_length, Method=="Otter trawl"),
              prior=prior(normal(0,10), class="Intercept")+
                prior(normal(0,5), class="b"),
              chains=1, cores=1,
              iter = iterations, warmup = warmup)


# Aggregating data --------------------------------------------------------

Data_ag<-Data%>%
  group_by(Station, Station_fac, Year, Latitude, Longitude, Source)%>%
  summarize(Count=sum(Count), Tow_area=sum(Tow_area))%>%
  ungroup()%>%
  mutate(CPUE=Count/Tow_area)%>%
  mutate_at(vars(Longitude, Latitude, Year, Tow_area), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T)))

Data_ag2<-Data_ag%>%
  filter(Station%in%unique(Stations$Station))%>%
  droplevels()

m3<-gam(CPUE ~ s(Year) + s(Station_fac, bs="re"), data=Data_ag2, family=gaussian)

m3.2<-gam(CPUE ~ s(Year) + s(Station_fac, bs="re"), data=Data_ag2, family=tw)

m3.3<-gam(CPUE ~ s(Year) + s(Station_fac, bs="re"), data=Data_ag2, family=scat)
