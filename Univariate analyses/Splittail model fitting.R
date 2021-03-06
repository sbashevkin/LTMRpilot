# Model fitting script
# Due to the long computational demand of this script, it is designed to be divided among multiple computers
# After running lines 7-15, each chunk can be run separately
# To run the whole script at once, the code should be modified to move all reduced model model fitting within the
# for loops as updates to the Full model. 

require(brms)

# Enter the folder where you've saved the "Split data.Rds" file and where you would like the models saved
save_folder<-"Univariate analyses"

load(file.path(save_folder, "Split data.Rds"))

iterations<-5e3
warmup<-iterations/4

# Full model --------------------------------------------------------------

model<-brm(as.integer(round(Count)) ~ Tow_area_s + Year_fac*Season + (1|Station_fac) + (1|ID),
            family=poisson, data=Data_split,
            prior=prior(normal(0,5), class="Intercept")+
              prior(normal(0,5), class="b")+
              prior(cauchy(0,5), class="sd"),
            chains=3, cores=3, control=list(max_treedepth=15),
            iter = iterations, warmup = warmup)
save(model, file=file.path(save_folder, paste0("Splittail full model.Rds")), compress="xz")

# 10% station cut ---------------------------------------------------------

N<-10

# Fit first model so we only need to compile once
model1<-brm(as.integer(round(Count)) ~ Tow_area_s + Year_fac*Season + (1|Station_fac) + (1|ID),
            family=poisson, data=subset(Data_split, Group_10!=1),
            prior=prior(normal(0,5), class="Intercept")+
              prior(normal(0,5), class="b")+
              prior(cauchy(0,5), class="sd"),
            chains=3, cores=3, control=list(max_treedepth=15),
            iter = iterations, warmup = warmup)
model<-model1 # Ensure all saved models have same name and model1 is not overwritten
save(model, file=file.path(save_folder, paste0("Splittail ", 1/N, " station cut ", "1", " of ", N, ".Rds")), compress="xz")
rm(model)
gc()
# Fit remaining models
for(i in 2:N){
  model<-update(model1, newdata=subset(Data_split, Group_10!=i),
                chains=3, cores=3, control=list(max_treedepth=15),
                iter = iterations, warmup = warmup)
  save(model, file=file.path(save_folder, paste0("Splittail ", 1/N, " station cut ", i, " of ", N, ".Rds")), compress="xz")
  rm(model)
  gc()
}  
warn<-warnings()
if(!is.null(warn)){
  save(warn, file=file.path(save_folder, paste0("Splittail ", 1/N, " station cut warnings.Rds")), compress="xz")
}

rm(model1)

# 20% station cut ---------------------------------------------------------

N<-5

# Fit first model so we only need to compile once
model1<-brm(as.integer(round(Count)) ~ Tow_area_s + Year_fac*Season + (1|Station_fac) + (1|ID),
            family=poisson, data=subset(Data_split, Group_5!=1),
            prior=prior(normal(0,5), class="Intercept")+
              prior(normal(0,5), class="b")+
              prior(cauchy(0,5), class="sd"),
            chains=3, cores=3, control=list(max_treedepth=15),
            iter = iterations, warmup = warmup)
model<-model1 # Ensure all saved models have same name and model1 is not overwritten
save(model, file=file.path(save_folder, paste0("Splittail ", 1/N, " station cut ", "1", " of ", N, ".Rds")), compress="xz")
rm(model)
gc()
# Fit remaining models
for(i in 2:N){
  model<-update(model1, newdata=subset(Data_split, Group_5!=i),
                chains=3, cores=3, control=list(max_treedepth=15),
                iter = iterations, warmup = warmup)
  save(model, file=file.path(save_folder, paste0("Splittail ", 1/N, " station cut ", i, " of ", N, ".Rds")), compress="xz")
  rm(model)
  gc()
}  
warn<-warnings()
if(!is.null(warn)){
  save(warn, file=file.path(save_folder, paste0("Splittail ", 1/N, " station cut warnings.Rds")), compress="xz")
}
rm(model1)

# 33% & 50% station cuts ---------------------------------------------------------

N<-3

# Fit first model so we only need to compile once
model1<-brm(as.integer(round(Count)) ~ Tow_area_s + Year_fac*Season + (1|Station_fac) + (1|ID),
            family=poisson, data=subset(Data_split, Group_3!=1),
            prior=prior(normal(0,5), class="Intercept")+
              prior(normal(0,5), class="b")+
              prior(cauchy(0,5), class="sd"),
            chains=3, cores=3, control=list(max_treedepth=15),
            iter = iterations, warmup = warmup)
model<-model1 # Ensure all saved models have same name and model1 is not overwritten
save(model, file=file.path(save_folder, paste0("Splittail ", 1/N, " station cut ", "1", " of ", N, ".Rds")), compress="xz")
rm(model)
gc()

# Fit remaining models
for(i in 2:N){
  model<-update(model1, newdata=subset(Data_split, Group_3!=i),
                chains=3, cores=3, control=list(max_treedepth=15),
                iter = iterations, warmup = warmup)
  save(model, file=file.path(save_folder, paste0("Splittail ", 1/N, " station cut ", i, " of ", N, ".Rds")), compress="xz")
  rm(model)
  gc()
}  
warn<-warnings()
if(!is.null(warn)){
  save(warn, file=file.path(save_folder, paste0("Splittail ", 1/N, " station cut warnings.Rds")), compress="xz")
}

N<-2
# Fit 50% cut models
for(i in 1:N){
  model<-update(model1, newdata=subset(Data_split, Group_2!=i),
                chains=3, cores=3, control=list(max_treedepth=15),
                iter = iterations, warmup = warmup)
  save(model, file=file.path(save_folder, paste0("Splittail ", 1/N, " station cut ", i, " of ", N, ".Rds")), compress="xz")
  rm(model)
  gc()
}  
warn<-warnings()
if(!is.null(warn)){
  save(warn, file=file.path(save_folder, paste0("Splittail ", 1/N, " station cut warnings.Rds")), compress="xz")
}
rm(model1)


# 67% (2/3) station cuts --------------------------------------------------------

N<-3

# Fit first model so we only need to compile once
model1<-brm(as.integer(round(Count)) ~ Tow_area_s + Year_fac*Season + (1|Station_fac) + (1|ID),
            family=poisson, data=subset(Data_split, Group_2.3==1),
            prior=prior(normal(0,5), class="Intercept")+
              prior(normal(0,5), class="b")+
              prior(cauchy(0,5), class="sd"),
            chains=3, cores=3, control=list(max_treedepth=15),
            iter = iterations, warmup = warmup)
model<-model1 # Ensure all saved models have same name and model1 is not overwritten
save(model, file=file.path(save_folder, paste0("Splittail ", 2/N, " station cut ", "1", " of ", N, ".Rds")), compress="xz")
rm(model)
gc()

# Fit remaining models
for(i in 2:N){
  model<-update(model1, newdata=subset(Data_split, Group_2.3==i),
                chains=3, cores=3, control=list(max_treedepth=15),
                iter = iterations, warmup = warmup)
  save(model, file=file.path(save_folder, paste0("Splittail ", 2/N, " station cut ", i, " of ", N, ".Rds")), compress="xz")
  rm(model)
  gc()
}  
warn<-warnings()
if(!is.null(warn)){
  save(warn, file=file.path(save_folder, paste0("Splittail ", 2/N, " station cut warnings.Rds")), compress="xz")
}

# 1 & 2 month cuts ------------------------------------------------------------

# Fit first model so we only need to compile once
model1<-brm(as.integer(round(Count)) ~ Tow_area_s + Year_fac*Season + (1|Station_fac) + (1|ID),
            family=poisson, data=subset(Data_split, Month_num!=1),
            prior=prior(normal(0,5), class="Intercept")+
              prior(normal(0,5), class="b")+
              prior(cauchy(0,5), class="sd"),
            chains=3, cores=3, control=list(max_treedepth=15),
            iter = iterations, warmup = warmup)
model<-model1 # Ensure all saved models have same name and model1 is not overwritten
save(model, file=file.path(save_folder, paste0("Splittail 1 month cut ", "1", " of ", 3, ".Rds")), compress="xz")
rm(model)
gc()

# Fit remaining models
for(i in 2:3){
  model<-update(model1, newdata=subset(Data_split, Month_num!=i),
                chains=3, cores=3, control=list(max_treedepth=15),
                iter = iterations, warmup = warmup)
  save(model, file=file.path(save_folder, paste0("Splittail 1 month cut ", i, " of ", 3, ".Rds")), compress="xz")
  rm(model)
  gc()
}  
warn<-warnings()
if(!is.null(warn)){
  save(warn, file=file.path(save_folder, paste0("Splittail 1 month cut warnings.Rds")), compress="xz")
}

# Fit 2-month cut models
for(i in 1:3){
  model<-update(model1, newdata=subset(Data_split, Month_num==i),
                chains=3, cores=3, control=list(max_treedepth=15),
                iter = iterations, warmup = warmup)
  save(model, file=file.path(save_folder, paste0("Splittail 2 month cut ", i, " of ", 3, ".Rds")), compress="xz")
  rm(model)
  gc()
}  
warn<-warnings()
if(!is.null(warn)){
  save(warn, file=file.path(save_folder, paste0("Splittail 2 month cut warnings.Rds")), compress="xz")
}
rm(model1)
