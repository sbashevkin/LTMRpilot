require(brms)
require(dplyr)
require(ggplot2)
require(tidyr)

newdata<-expand.grid(Tow_area_s=0,  Year_fac=factor(1980:2019), Season=c("Winter", "Spring", "Summer", "Fall"), scale='response')

pred_full<-fitted(mbrm7, newdata = newdata, re_formula=NA)
pred_full2<-fitted(mbrm7, newdata = newdata, re_formula=NA, summary=F)

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
