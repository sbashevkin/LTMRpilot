require(brms)
require(dplyr)
require(ggplot2)
require(tidyr)
require(tidybayes)

newdata<-expand.grid(Tow_area_s=0,  Year_fac=factor(1980:2019), Season=c("Winter", "Spring", "Summer", "Fall"))

pred_full<-fitted(mbrm7, newdata = newdata, re_formula=NA)
pred_fullb<-fitted(mbrm7, newdata = newdata, re_formula=NA, summary=F)

mbrm7_test<-brm(as.integer(round(Count)) ~ Tow_area_s + Year_fac*Season + (1|Station_fac) + (1|ID),
                family=poisson, data=filter(Data, Station_fac%in%sample(unique(Station_fac), 60)),
                prior=prior(normal(0,10), class="Intercept")+
                  prior(normal(0,5), class="b")+
                  prior(cauchy(0,5), class="sd"),
                chains=1, cores=1,
                iter = iterations, warmup = warmup)

pred_test<-fitted(mbrm7_test, newdata = newdata, re_formula=NA)

mbrm7_test2<-brm(as.integer(round(Count)) ~ Tow_area_s + Year_fac*Season + (1|Station_fac) + (1|ID),
                 family=poisson, data=filter(Data, Station_fac%in%sample(unique(Station_fac), 40)),
                 prior=prior(normal(0,10), class="Intercept")+
                   prior(normal(0,5), class="b")+
                   prior(cauchy(0,5), class="sd"),
                 chains=1, cores=1,
                 iter = iterations, warmup = warmup)

pred_test2<-fitted(mbrm7_test2, newdata = newdata, re_formula=NA)
pred_test2b<-fitted(mbrm7_test2, newdata = newdata, re_formula=NA, summary=F)

newdata2<-newdata%>%
  mutate(Full_pred=pred_full[,1],
         Full_SD=pred_full[,2],
         Full_l95=pred_full[,3],
         Full_u95=pred_full[,4],
         Full_SD2=Full_u95-Full_l95,
         Test_pred=pred_test[,1],
         Test_SD=pred_test[,2],
         Test_l95=pred_test[,3],
         Test_u95=pred_test[,4],
         Test_SD2=Test_u95-Test_l95,
         Test2_pred=pred_test2[,1],
         Test2_SD=pred_test2[,2],
         Test2_l95=pred_test2[,3],
         Test2_u95=pred_test2[,4],
         Test2_SD2=Test2_u95-Test2_l95,
         Year=as.numeric(as.character(Year_fac)))

newdata_plot<-newdata2%>%
  pivot_longer(cols=c(Full_pred, Full_SD, Full_l95, Full_u95, Test_pred, Test_SD, Test_l95, Test_u95, Test2_pred, Test2_SD, Test2_l95, Test2_u95),
               names_to=c("Model", ".value"), names_sep="_")

ggplot(filter(newdata_plot, Model%in%c("Full", "Test2")), aes(x=Year, y=pred, ymin=l95, ymax=u95, color=Model, fill=Model))+
  geom_ribbon(alpha=0.2)+
  geom_line()+
  scale_y_continuous(expand=expansion(0,0))+
  facet_wrap(~Season)+
  scale_fill_discrete(labels=c("Full", "Reduced"), aesthetics = c("colour", "fill"))+
  ylab("Predicted value")+
  theme_bw()+
  theme(panel.grid=element_blank(), text=element_text(size=18), legend.position=c(0.9, 0.9), legend.background = element_rect(color="black"))


ggplot(newdata_plot, aes(x=Year, y=pred, ymin=l95, ymax=u95, color=Model, fill=Model))+
  geom_pointrange(position=position_dodge(width=1))+
  scale_y_continuous(expand=expansion(0,0))+
  facet_wrap(~Season)+
  theme_bw()+
  theme(panel.grid=element_blank())

newdata_diff<-newdata2%>%
  mutate(across(c(Test_pred, Test2_pred), list(diff=~(.x-Full_pred)/Full_pred)))%>%
  mutate(across(c(Test_SD, Test2_SD), list(diff=~(.x-Full_SD)/Full_SD)))%>%
  mutate(across(c(Test_SD2, Test2_SD2), list(diff=~(.x-Full_SD2)/Full_SD2)))

ggplot(newdata_diff)+
  geom_tile(aes(x=Year, y=Season, fill=Test2_pred_diff))+
  scale_fill_viridis_c(labels=scales::percent, name="Change in predicted value")+
  coord_cartesian(expand=0)+
  theme(text=element_text(size=18))

ggplot(newdata_diff)+
  geom_tile(aes(x=Year, y=Season, fill=Test2_SD_diff))+
  scale_fill_viridis_c(labels=scales::percent, name="Change in SD")+
  coord_cartesian(expand=0)+
  theme(text=element_text(size=18))


# Using change in abundance -----------------------------------------------
Posterior_differ<-function(model, newdata, name){
  pred<-fitted(model, newdata = newdata, re_formula=NA, summary=F)
  col_names<-set_names(c("Pred", "Pred_lag", "Pred_change"), paste0(name, c("", "_lag", "_change")))
  out<-newdata%>%
    mutate(Row=row_number(),
           Year=as.numeric(as.character(Year_fac)))%>%
    rowwise()%>%
    mutate(Pred=list(pred[,Row]))%>%
    ungroup()%>%
    group_by(Season)%>%
    mutate(Pred_lag=lag(Pred, order_by=Year))%>%
    filter(Year!=min(Year))%>%
    ungroup()%>%
    rowwise()%>%
    mutate(Pred_change=list(Pred_lag-Pred))%>%
    rename(!!!col_names)%>%
    mutate(across(all_of(c(names(col_names))), list(mean=mean, sd=sd)))%>%
    ungroup()%>%
    select(-all_of(names(col_names)), -Row, -Tow_area_s, -Year_fac)
}

Change<-Posterior_differ(mbrm7, newdata, "Full")%>%
  left_join(Posterior_differ(mbrm7_test, newdata, "Reduced1"), by=c("Season", "Year"))%>%
  left_join(Posterior_differ(mbrm7_test2, newdata, "Reduced2"), by=c("Season", "Year"))%>%
  mutate(across(c(Reduced1_change_mean, Reduced2_change_mean), list(diff=~.x-Full_change_mean)))%>%
  mutate(across(c(Reduced1_change_sd, Reduced2_change_sd), list(diff=~.x-Full_change_sd)))

ggplot(filter(Change, Year>1985))+
  geom_tile(aes(x=Year, y=Season, fill=Reduced2_change_mean_diff))+
  scale_fill_gradient2(name="Change in \npredicted value")+
  coord_cartesian(expand=0)+
  theme(text=element_text(size=18))

Change_plot<-Change%>%
  select(Year, Season, Reduced1_change_mean, Reduced2_change_mean, Full_change_mean, Reduced1_change_sd, Reduced2_change_sd, Full_change_sd)%>%
  pivot_longer(cols=c(Reduced1_change_mean, Reduced2_change_mean, Full_change_mean, Reduced1_change_sd, Reduced2_change_sd, Full_change_sd),
               names_to=c("Model","change", ".value"), names_sep="_")%>%
  select(-change)

ggplot(filter(Change_plot, Model%in%c("Full", "Reduced2") & Year>1982), aes(x=Year, y=mean, ymin=mean-sd, ymax=mean+sd, color=Model, fill=Model))+
  geom_ribbon(alpha=0.2)+
  geom_line()+
  geom_hline(yintercept=0)+
  scale_y_continuous(expand=expansion(0,0))+
  facet_wrap(~Season)+
  scale_fill_discrete(labels=c("Full", "Reduced"), aesthetics = c("colour", "fill"))+
  ylab("'Population' change from prior year")+
  theme_bw()+
  theme(panel.grid=element_blank(), text=element_text(size=18), legend.position=c(0.6, 0.9), legend.background = element_rect(color="black"))

# Using tidybayes ---------------------------------------------------------


Effort<-Data%>%group_by(Season, Year)%>%summarise(N=n(), .groups="drop")

Post_data<-newdata%>%
  mutate(Year=as.numeric(as.character(Year_fac)))%>%
  filter(Year>=1985)%>%
  add_fitted_draws(mbrm7, re_formula=NA, scale="response")%>%
  ungroup()%>%
  mutate(Model="Full",
         Mean=mean(.value))%>%
  group_by(Season, .draw)%>%
  mutate(Lag=lag(.value, order_by=Year))%>%
  mutate(Change_global=(.value-Lag)/(Mean),
         Change_local=(.value-Lag)/(.value+Lag))%>%
  ungroup()%>%
  filter(Year!=min(Year))%>%
  bind_rows(newdata%>%
              mutate(Year=as.numeric(as.character(Year_fac)))%>%
              filter(Year>=1985)%>%
              add_fitted_draws(mbrm7_test, re_formula=NA, scale="response")%>%
              ungroup()%>%
              mutate(Model="Reduced",
                     Mean=mean(.value))%>%
              group_by(Season, .draw)%>%
              mutate(Lag=lag(.value, order_by=Year))%>%
              mutate(Change_global=(.value-Lag)/(Mean),
                     Change_local=(.value-Lag)/(.value+Lag))%>%
              ungroup()%>%
              filter(Year!=min(Year))
            
  )
  

Intervals<-Post_data%>%
  filter(Model=="Full")%>%
  group_by(Season, Year, Year_fac)%>%
  median_qi(Change_local, Change_global, .width = 0.95)%>%
  select(Season, Year, Year_fac, Change_local.lower, Change_local.upper, Change_global.lower, Change_global.upper)%>%
  ungroup()

Post_data_probs<-Post_data%>%
  filter(Model=="Reduced")%>%
  left_join(Intervals, by=c("Year_fac", "Year", "Season"))%>%
  mutate(IN_global=if_else(Change_global > Change_global.lower & Change_global <= Change_global.upper, 1, 0),
         IN_local=if_else(Change_local > Change_local.lower & Change_local <= Change_local.upper, 1, 0))%>%
  group_by(Year, Year_fac, Season)%>%
  summarise(N=n(), Prob_global=sum(IN_global)/N, Prob_local=sum(IN_local)/N, .groups="drop")

ggplot(Post_data)+
  stat_slab(aes(x=Year_fac, y=Change_local, fill = Model), alpha=0.5)+
  geom_point(data=Post_data_probs2, aes(x=Year_fac, y=Prob))+
  geom_hline(yintercept = c(-.1, .1), linetype = "dashed") +
  facet_wrap(~Season)+
  ylab("Change in abundance (Standardized by local magnitude)")+
  xlab("Year")+
  #coord_cartesian(ylim=c(-5,5))+
  scale_x_discrete(breaks=unique(Post_data$Year), labels = if_else(unique(Post_data$Year)%% 2 == 0, as.character(unique(Post_data$Year)), ""))+
  scale_fill_manual(values = c("dodgerblue3", "firebrick1"), aesthetics = c("fill", "color"))+
  theme_bw()+
  theme(panel.grid=element_blank(), text=element_text(size=18), axis.text.x=element_text(angle=45, hjust=1))
  
ggplot(Post_data, aes(x=Year_fac, y=Change_local, fill = Model))+
  stat_slab(alpha=0.5)+
  ylab("Change in abundance (F(t)-F(t-1))/(F(t)+F(t-1))")+
  geom_hline(yintercept = c(-.1, .1), linetype = "dashed") +
  facet_wrap(~Season)+
  scale_x_discrete(breaks=unique(Post_data$Year), labels = if_else(unique(Post_data$Year)%% 2 == 0, as.character(unique(Post_data$Year)), ""))+
  scale_fill_manual(values = c("gray80", "skyblue"), aesthetics = c("fill", "color"))+
  theme_bw()+
  theme(panel.grid=element_blank(), axis.text.x=element_text(angle=45, hjust=1))

ggplot(Post_data_probs, aes(x=Year, y=Prob))+
  geom_point()+
  facet_wrap(~Season)+
  scale_x_discrete(breaks=unique(Post_data_probs$Year), labels = if_else(unique(Post_data_probs$Year)%% 2 == 0, as.character(unique(Post_data_probs$Year)), ""))+
  theme_bw()+
  theme(panel.grid=element_blank(), axis.text.x=element_text(angle=45, hjust=1))
