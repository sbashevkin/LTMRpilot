require(brms)
require(dplyr)
require(tidybayes)
require(tibble)
require(purrr)
require(ggplot2)
require(patchwork)
require(stringr)
require(tidyr)
require(furrr)
source("Univariate analyses/Survey assessment functions.R")

Data_split<-readRDS("Univariate analyses/Split data.Rds")

iterations<-5e3
warmup<-iterations/4

model_offset<-brm(as.integer(round(Count)) ~ offset(log(Tow_area)) + Year_fac*Season + (1|Station_fac) + (1|ID),
                  family=poisson, data=Data_split,
                  prior=prior(normal(0,5), class="Intercept")+
                    prior(normal(0,5), class="b")+
                    prior(cauchy(0,5), class="sd"),
                  chains=3, cores=3, control=list(max_treedepth=15),
                  iter = iterations, warmup = warmup,
                  backend = "cmdstanr", threads = threading(5))
#saveRDS(model_offset, file.path("Univariate analyses", "Splittail models", "Splittail full model offset.Rds"))
model_offset<-readRDS(file.path("Univariate analyses", "Splittail models", "Splittail full model offset.Rds"))
model_offset<-add_criterion(model_offset, c("waic", "loo"), cores=10)

model<-readRDS(file.path("Univariate analyses", "Splittail models", "Splittail full model.Rds"))
model<-add_criterion(model, c("waic", "loo"), cores=10)

# Compare models with loo and waic
loo_compare(model, model_offset, criterion="waic")
loo_compare(model, model_offset, criterion="loo")
# the model without an offset (estimating a linear coefficient for Tow_area_s) is better by both waic and loo, but the se of the 
# difference is greater than the difference in both cases, and there were some issues
# with the calculation of both metrics, so moving to the very computationally expensive kfold


model_offset<-add_criterion(model_offset, criterion="kfold", chains = 1, threads=threading(15))
model<-add_criterion(model, criterion="kfold", chains = 1, threads=threading(15))
# took 2 days
# Load the model results from above through here from this file: https://deltacouncil.box.com/shared/static/sncj3f0hrs23llov94vlqoz3xllu1qga.rdata
# load("offset comparison full model.RData")

loo_compare(model, model_offset, criterion="kfold")
# Now the model with the offset is better, but the se of the difference is still much greater than the difference

# Do the credible intervals of the predictions from the 2 models differ?

Full_change<-readRDS("Univariate analyses/Full model local trend.Rds") # From the "Splittail model processing.R" script
Full_change_offset<-Post_processor(model_offset, max_year=2018, model_name="Full")

# Compare upper and lower limits between the 2 models
plot(Full_change$Change_local.lower, Full_change_offset$Change_local.lower)
plot(Full_change$Change_local.upper, Full_change_offset$Change_local.upper)

# Credible intervals of the local trend are extremely similar between the 2 models

# Let's use a simple regression to see how similar they are

summary(lm(Full_change$Change_local.lower~Full_change_offset$Change_local.lower))
# R2 = 0.9992, slope=0.9996521, median residual = 0.000527

summary(lm(Full_change$Change_local.upper~Full_change_offset$Change_local.upper))
# R2 = 0.9991, slope=0.998238, median residual = 0.000360