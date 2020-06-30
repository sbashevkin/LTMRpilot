Post_processor<-function(model, min_year=1985, model_name=NULL, Intervals=NULL){
  out<-expand.grid(Tow_area_s=0,  Year_fac=factor(1980:2019), Season=c("Winter", "Spring", "Summer", "Fall"))%>%
    mutate(Year=as.numeric(as.character(Year_fac)))%>%
    filter(Year>=min_year)%>%
    add_fitted_draws(model, re_formula=NA, scale="response")%>%
    ungroup()%>%
    mutate(Model=model_name,
           Mean=mean(.value))%>%
    group_by(Season, .draw)%>%
    mutate(Lag=lag(.value, order_by=Year))%>%
    mutate(Change_global=(.value-Lag)/(Mean),
           Change_local=(.value-Lag)/(.value+Lag))%>%
    ungroup()%>%
    filter(Year!=min(Year))%>%
    {if(model_name=="Full"){
      group_by(., Season, Year, Year_fac)%>%
        median_qi(Change_local, Change_global, .width = 0.95)%>%
        select(Season, Year, Year_fac, Change_local.lower, Change_local.upper, Change_global.lower, Change_global.upper)%>%
        ungroup()
    } else{
      left_join(., Intervals, by=c("Year_fac", "Year", "Season"))%>%
        mutate(IN_global=if_else(Change_global > Change_global.lower & Change_global <= Change_global.upper, 1, 0),
               IN_local=if_else(Change_local > Change_local.lower & Change_local <= Change_local.upper, 1, 0))%>%
        group_by(Year, Year_fac, Season)%>%
        summarise(N=n(), Prob_global=sum(IN_global)/N, Prob_local=sum(IN_local)/N, .groups="drop")
    }}
 return(out) 
}

random_groups<-function(seed, N_total, N_groups){
  set.seed(seed)
  Groups<-rep(1:N_groups, floor(N_total/N_groups))
  if(length(Groups)!=N_total){
    Groups<-c(Groups, sample(1:N_groups, size=(N_total-length(Groups)), replace=F))
  }
  out<-sample(Groups, size=N_total, replace=F)
  set.seed(NULL)
  return(out)
}

Station_splits<-Data%>%
  select(Station)%>%
  distinct()%>%
  mutate(Group_10=random_groups(12, nrow(.), 10),
         Group_5=random_groups(123, nrow(.), 5),
         Group_3=random_groups(1234, nrow(.), 3),
         Group_2=random_groups(12345, nrow(.), 2))

Data_split<-Data%>%
  left_join(Station_splits, by="Station")%>%
  mutate(Month_num=case_when(
    Month%in%c(12,3,6,9) ~ 1,
    Month%in%c(1,4,7,10) ~ 2,
    Month%in%c(2,5,8,11) ~ 3,
  ))

save(Data_split, file="Split data.Rds")

# Check reduced models
model_diagnose<-function(model){
  Bulk_ESS<-map_dbl(parnames(model)[-grep("^r_", parnames(model))], ~rstan::ess_bulk(as.array(model)[,,.x]))
  Tail_ESS<-map_dbl(parnames(model)[-grep("^r_", parnames(model))], ~rstan::ess_tail(as.array(model)[,,.x]))
  Rhat<-rhat(model, pars=parnames(model)[-grep("^r_", parnames(model))])
  out<-tibble(Bulk_ESS, Tail_ESS, Rhat)
  return(out)
}

# Should be at least 100 ESS per chain, Rhat <1.05