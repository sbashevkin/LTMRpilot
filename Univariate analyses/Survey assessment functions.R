# Generate random groups in a repeatable manner 
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

# Create posterior predictions and calculate change from one year to next within each season
model_predictor<-function(model, max_year=2018){
  out<-expand.grid(Tow_area_s=0,  Year_fac=factor(1985:max_year), Season=c("Winter", "Spring", "Summer", "Fall"))%>%
    mutate(Year=as.numeric(as.character(Year_fac)))%>%
    add_fitted_draws(model, re_formula=NA, scale="response")%>%
    ungroup()%>%
    mutate(Mean=mean(.value))%>%
    group_by(Season, .draw)%>%
    mutate(Lag=lag(.value, order_by=Year))%>%
    mutate(Change_global=(.value-Lag)/(Mean),
           Change_local=(.value-Lag)/(.value+Lag))%>%
    ungroup()
  
  return(out)
}

# For the full model, calculates the 95% intervals
# For the reduced models, calculates the proportion of posterior samples within the 95% intervals of the full model
# model_name = "Full", or "Reduced"
# Intervals (i.e., output of this function when model_name="Full") should be supplied when model_name="Reduced"
Post_processor<-function(model, max_year=2018, model_name=NULL, Intervals=NULL){
  out<-model_predictor(model, max_year=max_year)%>%
    filter(Year!=min(Year))%>%
    mutate(Model=model_name)%>%
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

# Check models for any diagnostic issues (Rhat and effective sample size)
# ESS should be at least 100 per chain (300 total), Rhat should be <1.05.
model_diagnose<-function(model){
  sum<-summary(model)
  out<-tibble(par=c(names(sum$fixed[,1]), names(sum$random)), 
              Bulk_ESS=c(sum$fixed[,"Bulk_ESS"], sum$random$ID[,"Bulk_ESS"], sum$random$Station_fac[,"Bulk_ESS"]),
              Tail_ESS=c(sum$fixed[,"Tail_ESS"], sum$random$ID[,"Tail_ESS"], sum$random$Station_fac[,"Tail_ESS"]),
              Rhat=c(sum$fixed[,"Rhat"], sum$random$ID[,"Rhat"], sum$random$Station_fac[,"Rhat"]))
  Bulk_ESS<-filter(out, Bulk_ESS<300)%>%
    pull(par)
  Tail_ESS<-filter(out, Tail_ESS<300)%>%
    pull(par)
  Rhat<-filter(out, Rhat>1.05)%>%
    pull(par)
  
  if(length(Bulk_ESS)>0){
    message(paste("WARNING: Bulk_ESS is too low in the parameters:", paste(Bulk_ESS, collapse=", ")))
  }
  
  if(length(Tail_ESS)>0){
    message(paste("WARNING: Tail_ESS is too low in the parameters:", paste(Tail_ESS, collapse=", ")))
  }
  
  if(length(Rhat)>0){
    message(paste("WARNING: Rhat is too high in the parameters:", paste(Rhat, collapse=", ")))
  }
  
  if(length(c(Bulk_ESS, Tail_ESS, Rhat)>0)){
    return(out)
  }else{
    print("All good, no warnings!")
  }
}


# Function to visualize posterior distributions and predictions
# plot_type can be either "distributions" or "time series"
# when plot_type="distributions", y should be supplied as either "Change_local" or "Change_global"

Posterior_plotter<-function(model_full, model_reduced, plot_type, y=NULL, max_year=2018){
  if(is.character(model_reduced)){
    load(file.path("Univariate analyses", "Splittail models", model_reduced))
    model_reduced<-model
    rm(model)
  }
  
  Data<-model_predictor(model_full, max_year=max_year)%>%
    mutate(Model="Full")%>%
    bind_rows(model_predictor(model_reduced, max_year=max_year)%>%
    mutate(Model="Reduced"))
  
  if(plot_type=="distributions"){
    p<-ggplot(Data)+
      stat_slab(aes(x=Year_fac, y=.data[[y]], fill = Model), alpha=0.5)+
      facet_wrap(~Season)+
      ylab(paste0("Change in abundance ", if_else(y=="Change_local", "(Standardized by local magnitude)", "(Standardized by global mean)")))+
      xlab("Year")+
      {if(y=="Change_global"){
        coord_cartesian(ylim=c(-5,5))
      }}+
      scale_x_discrete(breaks=unique(Data$Year), labels = if_else(unique(Data$Year)%% 2 == 0, as.character(unique(Data$Year)), ""))+
      scale_fill_manual(values = c("dodgerblue3", "firebrick1"), aesthetics = c("fill", "color"))+
      theme_bw()+
      theme(panel.grid=element_blank(), text=element_text(size=18), axis.text.x=element_text(angle=45, hjust=1))
  }
  
  if(plot_type=="time series"){
    Data<-Data%>%
      group_by(Season, Year, Year_fac, Model)%>%
      mean_qi(.value)%>%
      mutate(Model=factor(Model, levels=c("Full", "Reduced")))
    
    p<-ggplot(Data, aes(x=Year, y=.value, ymin=.lower, ymax=.upper, color=Model, fill=Model))+
      geom_ribbon(alpha=0.2)+
      geom_line()+
      scale_y_continuous(expand=expansion(0,0))+
      facet_wrap(~Season)+
      scale_fill_manual(aesthetics = c("colour", "fill"), values = c("dodgerblue3", "firebrick1"))+
      ylab("Predicted value")+
      theme_bw()+
      theme(panel.grid=element_blank(), text=element_text(size=18), legend.position=c(0.1, 0.9), legend.background = element_rect(color="black"))
  }
  
  return(p) 
}

Distribution_plotter<-function(Full_post, Reduced_post, Y){
  Data<-bind_rows(Full_post%>%
                    mutate(Model="Full"),
                  Reduced_post%>%
                    mutate(Model="Reduced"))
  
  if(Y%in%c("Change_local", "Change_global")){
    Data<-Data%>%
      filter(Year!=min(Year))
  }
  
  if(Y=="Change_local"){
    ylims<-c(-1,1)
  } else{
    ylims<-Data%>%
      group_by(Model, Season, Year)%>%
      summarise(L95=quantile(.data[[Y]], probs=0.025), U95=quantile(.data[[Y]], probs=0.975), .groups="drop")
    
    ylims=c(min(ylims$L95), max(ylims$U95))
  }
  
  p<-ggplot(Data, aes(x=Year_fac, y=.data[[Y]], fill = Model))+
    stat_slab(alpha=0.5)+
    {if(Y%in%c("Change_local", "Change_global")){
      geom_hline(yintercept = 0, linetype = "dashed")
    }}+
    coord_cartesian(ylim=ylims)+
    facet_wrap(~Season)+
    scale_x_discrete(breaks=seq(1985, 2020, by=5))+
    scale_fill_manual(values = c("dodgerblue3", "firebrick1"), aesthetics = c("fill", "color"))+
    xlab("Year")+
    theme_bw()+
    theme(panel.grid=element_blank(), text=element_text(size=8), axis.text.x=element_text(angle=45, hjust=1), 
          strip.background=element_blank(), legend.position="none")
  
  return(p)
}

Ribbon_plotter<-function(Full_post, Reduced_post, Y){
  Data<-bind_rows(Full_post%>%
                    mutate(Model="Full"),
                  Reduced_post%>%
                    mutate(Model="Reduced"))
  
  if(Y%in%c("Change_local", "Change_global")){
    Data<-Data%>%
      filter(Year!=min(Year))
  }
  
  Data<-Data%>%
    group_by(Model, Year, Year_fac, Season)%>%
    summarise(Count=mean(.data[[Y]]), L95=quantile(.data[[Y]], probs=0.025), U95=quantile(.data[[Y]], probs=0.975), .groups="drop")
  
  p<-ggplot(Data, aes(x=Year, y=Count, ymin=L95, ymax=U95, fill = Model, group=Model))+
    geom_ribbon(alpha=0.2)+
    geom_line(aes(color=Model))+
    {if(Y%in%c("Change_local", "Change_global")){
      geom_hline(yintercept = 0, linetype = "dashed")
    }}+
    coord_cartesian(expand=0)+
    facet_wrap(~Season, nrow=1)+
    scale_x_continuous(breaks=seq(1985, 2020, by=5), limits=c(1985, 2018))+
    scale_fill_manual(values = c("dodgerblue3", "firebrick1"), aesthetics = c("fill", "color"))+
    xlab("Year")+
    theme_bw()+
    theme(text=element_text(size=8), axis.text.x=element_text(angle=45, hjust=1), strip.background=element_blank())
  
  return(p)
}
