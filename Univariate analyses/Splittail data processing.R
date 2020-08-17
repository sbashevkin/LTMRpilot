require(LTMRdata)
require(dplyr)
require(lubridate)
require(tidyr)
source("Univariate analyses/Survey assessment functions.R") # for random_groups function

Data <- LTMRpilot(convert_lengths=TRUE, remove_unconverted_lengths = TRUE, size_cutoff=40)%>% # Generate data
  zero_fill(species="Pogonichthys macrolepidotus", remove_unknown_lengths = TRUE, univariate = TRUE)%>% # Fill in 0s
  filter(Method=="Otter trawl" & year(Date)>=1980 & year(Date)<=2019)%>% # Select otter trawl data within date window
  drop_na(Sal_surf, Temp_surf, Count, Latitude, Longitude, Date, Tow_area)%>% # Remove rows with NAs in key variables
  group_by(Source, Station, Latitude, Longitude, Date, Tow_area, SampleID, Sal_surf, Temp_surf)%>%
  summarise(Count=sum(Count), Length=mean(Length, na.rm=T), .groups="drop")%>% # Sum catches of each length class within a sample
  group_by(Source, Station, Latitude, Longitude, Date)%>%
  summarise(Count=sum(Count), Tow_area=sum(Tow_area), Sal_surf=mean(Sal_surf), 
            Temp_surf=mean(Temp_surf), Length=mean(Length, na.rm=T), .groups="drop")%>% # Take averages when multiple samples were collected in a day
  mutate(Julian_day=yday(Date), 
         Year=year(Date), 
         Month=month(Date),
         CPUE=Count/Tow_area, # Calculate CPUE
         Station_fac=factor(Station))%>% # Create factor version of Station
  mutate(Year=if_else(Month==12, Year+1, Year))%>% # Move Decembers to the following year so winters are all in the same year
  mutate(Year_fac=ordered(Year), # Create ordered factor version of year
         Month_fac=ordered(Month), # Create ordered factor version of year
         ID=1:nrow(.), # Create observation level identifier (ID)
         Season=case_when( # Define seasons
           Month%in%c(12,1,2) ~ "Winter",
           Month%in%c(3,4,5) ~ "Spring",
           Month%in%c(6,7,8) ~ "Summer",
           Month%in%c(9,10,11) ~ "Fall"
         ))%>%
  mutate_at(vars(Length, Longitude, Latitude, Year, Julian_day, Sal_surf, Temp_surf, Tow_area), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T)))%>% # Create centered and standardized versions of covariates
  arrange(Date)%>% # sort data by date
  filter(Year>=1985)%>% # Remove early years before Suisun survey was regular
  droplevels() # Remove empty factor levels


# Randomly divide data for spatial data reduction scenarios. Randomly divide sampling stations into 10, 5, 3, 2, and 3 groups for the simulations
# of 1/10, 1/5, 1/3, 1/2, and 2/3 reductions in spatial sampling effort
Station_splits<-Data%>%
  select(Station)%>%
  distinct()%>%
  mutate(Group_10=random_groups(12, nrow(.), 10),
         Group_5=random_groups(123, nrow(.), 5),
         Group_3=random_groups(1234, nrow(.), 3),
         Group_2=random_groups(12345, nrow(.), 2),
         Group_2.3=random_groups(123456, nrow(.), 3)) # For 2/3 station cuts

# Add spatial sampling groups to dataset and also divide data by month for temporal sampling effort reduction scenarios
Data_split<-Data%>%
  left_join(Station_splits, by="Station")%>%
  mutate(Month_num=case_when(
    Month%in%c(12,3,6,9) ~ 1,
    Month%in%c(1,4,7,10) ~ 2,
    Month%in%c(2,5,8,11) ~ 3,
  ))

save(Data_split, file="Univariate analyses/Split data.Rds")