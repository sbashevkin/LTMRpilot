require(brms)
require(brms)
require(tibble)
require(tidyr)
require(purrr)
require(dplyr)
require(ggplot2)

# Fit models --------------------------------------------------------------

# Load data

Data_split<-readRDS("Univariate analyses/Split data.Rds")

iterations<-5e3
warmup<-iterations/4

# Fit 2 types of models to decompose variance. 

# Decompose variance into year, season, station, and ID
model_var<-brm(as.integer(round(Count)) ~ Tow_area_s + (1|Year_fac) + (1|Season) + (1|Station_fac) + (1|ID),
               family=poisson, data=Data_split,
               prior=prior(normal(0,5), class="Intercept")+
                 prior(normal(0,5), class="b")+
                 prior(cauchy(0,5), class="sd"),
               chains=3, cores=3, control=list(adapt_delta=0.99),
               iter = iterations, warmup = warmup)

# Decompose variance into year, month, station, and ID
# This will be used in write-up since it aligns more closely to other methods (months were removed in sampling reduction scenarios, not seasons)
model_var2<-brm(as.integer(round(Count)) ~ Tow_area_s + (1|Year_fac) + (1|Month) + (1|Station_fac) + (1|ID),
                family=poisson, data=Data_split,
                prior=prior(normal(0,5), class="Intercept")+
                  prior(normal(0,5), class="b")+
                  prior(cauchy(0,5), class="sd"),
                chains=3, cores=3, control=list(adapt_delta=0.9),
                iter = iterations, warmup = warmup,
                backend = "cmdstanr", threads = threading(2))
model_var2<-add_criterion(model_var2, c("waic", "loo"), cores=6)

# Try with an offset for tow area
model_var3<-brm(as.integer(round(Count)) ~ offset(log(Tow_area)) + (1|Year_fac) + (1|Month) + (1|Station_fac) + (1|ID),
                family=poisson, data=Data_split,
                prior=prior(normal(0,5), class="Intercept")+
                  prior(cauchy(0,5), class="sd"),
                chains=3, cores=3, control=list(adapt_delta=0.9),
                iter = iterations, warmup = warmup,
                backend = "cmdstanr", threads = threading(2))
model_var3<-add_criterion(model_var3, c("waic", "loo"), cores=6)

p<-pp_check(model_var2)
p
p+scale_x_continuous(trans="log1p")
prop_zero <- function(x) mean(x == 0)
pp_check(model_var2, type="stat", stat=prop_zero, resp=resp)

# Save both models
saveRDS(model_var2, file=file.path("Univariate analyses", "Splittail models", "variance model.Rds"), compress="xz")


# Create plots ------------------------------------------------------------

# Load models
model_var2<-readRDS(file.path("Univariate analyses", "Splittail models", "variance model.Rds"))

# Summarise model parameters and convert to data frame
sum<-summary(model_var2)

sum2<-enframe(sum$random)%>%
  unnest(cols=value)%>%
  mutate(name=recode(name, Year_fac="Year", Station_fac="Station", ID="Sample"))%>%
  mutate(name=factor(name, levels=c("Sample", "Station", "Month", "Year")))

p_var<-ggplot(sum2, aes(y=name, x=Estimate, xmin=`l-95% CI`, xmax=`u-95% CI`))+
  geom_pointrange()+
  geom_vline(xintercept=0)+
  ylab("Variable")+
  xlab("Variance parameter estimate (± 95% CI)")+
  theme_bw()
ggsave(p_var, file="Univariate analyses/Figures/Variance plot.png", device="png", height=5, width=5, units="in")
