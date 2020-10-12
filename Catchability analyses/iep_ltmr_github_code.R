#<><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
#<> Citation <>#
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
#<> IEP Long-term Survey Review Team. 2020. Interagency Ecological Program long-term survey designs to 
  #   provide data relevant to community status and trends in the San Francisco Estuary: a pilot review 
  #   of the Bay Study, UC Davis Suisun Marsh, and Fall Midwater Trawl Surveys, 2020. IEP Technical Report 

#<> Chapter 4: Catchability <>#

#<><><><><><><><><><><><>#
#<> Set Working Directory and Load Libraries <>#
#<><><><><><><><><><><><>#
setwd("~/IEP Survey Review/catchability github")

library(lubridate)
library(dplyr)
library(tidyr)
library(jagsUI)
library(ggplot2)
library(plotly)
library(grid)
library(mgcv)
library(gridExtra)

#<><><><><><><><><><><><>#
#<> Functions <>#
#<><><><><><><><><><><><>#
logit <- function(x){
  log(x / (1 - x))
}
ci05 <- function(x) {
  quantile(x, prob = 0.025)
}
ci25 <- function(x) {
  quantile(x, prob = 0.25)
}
ci75 <- function(x) {
  quantile(x, prob = 0.75)
}
ci95 <- function(x) {
  quantile(x, prob = 0.975)
}
firstNA <- function(x) {
  max(which(!is.na(x)))
}
MyNorm <- function(x) { 
  (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
} # Covariate standardizing function
MyNorm_integ <- function(x, y, z){ #Normalization Function for integrated data#
  (x - y) / z # x = observed data, y = total mean, z = total sd
}

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
#<> Workflow 4.1: Bias Caused by Ignoring the Observation Process in Survey Data <>#
  #<> R and JAGS code to simulate and analyze a standard Binomial N-mixture model based on Royle (2004).

  #<> Royle, J.A. 2004. N-Mixture Models for Estimating Population Size from Spatially Replicated Counts. 
    #<>   Biometrics 60(1): 108-115.
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
#<><><><><><><><><><><><>#
#<> Data Simulation <>#
#<><><><><><><><><><><><>#
set.seed(10)
nSites <- 100                                                               # Number of sampling stations
nRep <- 3                                                                   # Number of replicate samples collected within a site (repeat samples during closure)

p.intercept <- logit(0.6)                                                   # Intercept of detection efficiency set at 0.6 for the mean
p.covariate <- 2                                                            # Covariate effect on the detection (availability) process
N.intercept <- log(20)                                                      # Intercept of abundance mean set at 20
covariate <- rnorm(nSites, mean = 0, sd = 1)                                # Simulate covariate already standardized for analysis

N <- rpois(nSites, (exp(N.intercept)))                                      # Latent state abundance simulation

counts <- matrix(NA, nrow = nSites, ncol = nRep)                            # Placeholder
for (r in 1:nRep){
  counts[, r] <- rbinom(nSites,                                             # Binomial N-mixture simulation
                        N,
                        plogis(p.intercept + p.covariate * covariate))
} #r

#<><><><><><><><><><><><>#
#<> JAGS Models <>#
#<><><><><><><><><><><><>#
#<> Model 1: GLM of Hierarchical Data - Equation 3 in Report
sink("CH_4.1_model1.txt")
cat("
model{
  # Priors #
  b.0 ~ dnorm(0, 0.001)
  b.cov ~ dnorm(0, 0.001)

  # Likelihood and Constraints #
  for (i in 1:nSites){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- b.0 + b.cov * covariate[i]
    for (r in 1:nRep){
      count[i, r] ~ dpois(lambda[i])
    }
  }
} # # End model
", fill = TRUE)
sink()

#<> Model 2: Binomial N-mixture with covariate on N not p - Equation 4 in Report
sink("CH_4.1_model2.txt")
cat("
model{
  # Priors #
  b.0 ~ dnorm(0, 0.001)
  b.cov ~ dnorm(0, 0.001)
  a.0 ~ dnorm(0, 0.1)
  
  # Likelihood and Constraints #
  for (i in 1:nSites){
    logit(p[i]) <- a.0
  
    log(lambda[i]) <- b.0 + b.cov * covariate[i]
    N[i] ~ dpois(lambda[i])
    
    for (r in 1:nRep){
      count[i, r] ~ dbin(p[i], N[i])
    }
  }
} # # End model
", fill = TRUE)
sink()

#<> Model 3: Binomial N-mixture with covariate on p not N (TRUTH) - Equation 4 in Report
sink("CH_4.1_model3.txt")
cat("
model{
  # Priors #
  b.0 ~ dnorm(0, 0.001)
  a.cov ~ dnorm(0, 0.001)
  a.0 ~ dnorm(0, 0.1)
  
  # Likelihood and Constraints #
  for (i in 1:nSites){
    logit(p[i]) <- a.0 + a.cov * covariate[i]
  
    log(lambda[i]) <- b.0
    N[i] ~ dpois(lambda[i])
    
    for (r in 1:nRep){
      count[i, r] ~ dbin(p[i], N[i])
    }
  }
} # # End model
", fill = TRUE)
sink()

#<><><><><><><><><><><><>#
#<> Analyze Models <>#
#<><><><><><><><><><><><>#
# Load Data Into List for JAGS #
dat <- list(count = counts,
            nSites = nSites,
            nRep = nRep,
            covariate = covariate)

# Initial Values for JAGS #
inits <- function(){
  list(N = N)
}

#<> Parameters To Monitor <>#
param <- c("a.0", "a.cov", 
           "b.0", "b.cov",
           "N")

#<> Jags simulation setup <>#
nb <- 10000
ni <- 100000 + nb
nt <- 50
nc <- 3

#<> Vector of models to run and save <>#
model_ref <- c("CH_4.1_model1.txt",
               "CH_4.1_model2.txt",
               "CH_4.1_model3.txt")

model_name <- c("CH_4.1_model1",
               "CH_4.1_model2",
               "CH_4.1_model3")

#<> Pass Models Through JAGS and Save <>#
start <- Sys.time()
for (mod in 1:length(model_ref)){
  out <- jagsUI(dat, inits, param, model.file = model_ref[mod],
                n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
  
  save(out, file = paste(model_name[mod], "RData", sep = "."))
} #mod
end <- Sys.time()
end - start # 5.1 Minutes

#<><><><><><><><><><><><>#
#<> Load Results and Make Figure 4.3 <>#
#<><><><><><><><><><><><>#
cov_pred <- seq(from = min(covariate), to = max(covariate), length = 20)

PredLatentAVG <- array(NA, dim = c((((ni-nb) / nt) * nc), length(cov_pred), length(model_name)),      # Abundance Prediction Placeholder
                       dimnames = list(paste("Iter", 1:(((ni-nb) / nt) * nc), sep = "_"),
                                       paste("Cov", cov_pred, sep = "_"),
                                       model_name))

for (z in 1:length(model_name)){   
  load(file = paste(model_name[z], ".RData", sep = ""))              # Load models                                
  
  for (n in 1:length(cov_pred)){
    PredLatentAVG[, n, z] <- rpois(length(out$sims.list$b.0), exp(out$sims.list$b.0 +       # Predict abundance for each model
                                    ifelse(is.null(out$sims.list$b.cov), 0, out$sims.list$b.cov) * cov_pred[n]))
  }
}

PredLatentAVG <- as.data.frame.table(PredLatentAVG) %>%            # Put predictions into a dataframe for plotting
  group_by(Var2, Var3) %>%
  summarize(Mean = mean(Freq),
            Lower = ci05(Freq),
            Upper = ci95(Freq)) %>%
  ungroup() %>%
  separate(Var2, c("Remove", "Covariate"), sep = "_") %>%
  mutate(Model = ifelse(Var3 == "CH_4.1_model1", "Model 1",
                        ifelse(Var3 == "CH_4.1_model2", "Model 2", "Model 3")),
         Covariate = round(as.numeric(Covariate), digits = 1))
  
abundance_plot <- data.frame(Covariate = rep(covariate, 2),      # observed vs truth catch data frame for plotting
                             Data = rep(c("Observed", "Truth"), each = length(covariate)),
                             Mean = c(counts[, 1], N))          # just take one of the replicate passes for simplicity

figure_4.3 <- ggplot(data = PredLatentAVG,
       aes(x = Covariate, y = Mean)) +
#  ylim(0, 85) +
  labs(x = "Covariate", y = "Abundance") +
  geom_point(data = abundance_plot,
             aes(x = Covariate, y = Mean, shape = Data), size = 4, stroke = 1.5) +
  scale_shape_manual(values = c(1, 16)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Model), alpha = 0.3) +
  scale_fill_manual(values = c("#CC79A7", "#0072B2", "#009E73"),  #c("red", "blue", "green"),
                    labels = c("Model 1", "Model 2", "Model 3")) +
  theme(legend.position = c(0.15, 0.85), 
        axis.text = element_text(size = 12),
        legend.key = element_blank(),
        axis.title.y = element_text(vjust = 3),
        axis.title.x = element_text(vjust = -2.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        legend.box.background = element_blank(),
        text = element_text(size = 16),
        plot.title = element_text(vjust = -18, hjust = 0.95)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))


#<><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
#<> Workflow 4.2: Decomposition of Catchability and Integrating Common Fisheries Data Sets <>#
  #<> Code is based on Hostetter et al. (2019) which was adapted from Amundson et al. 2014 and 
  #<> Kery and Royle (2016) <>#

  #<> Hostetter, N.J., Gardner, B., Sillett, T.S., Pollock, K.H., and Simons, T.R. 2019. 
    #<> An integrated model decomposing the components of detection probability and abundance in
    #<> unmarked populations. Ecosphere 10(3): e02586. doi:10.1002/ecs2.2586.
  #<> Amundson, C.L., Royle, J.A., and Handel, C.M. 2014. A hierarchical model combining distance sampling 
    #<> and time removal to estimate detection probability during avian point counts. Auk 131(4): 476-494. 
    #<> doi:10.1642/auk-14-11.1.
  #<> Kéry, M., and Royle, J.A. 2016. Applied Hierarchical Modeling in Ecology: Analysis of distribution, 
    #<> abundance and species richness in R and BUGS: Volume 1: Prelude and Static Models. Academic Press.
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
#<><><><><><><><><><><><>#
#<> JAGS Models <>#
#<><><><><><><><><><><><>#
#<> Codend Only Analysis <>#
cat("
  model {
    ##################################
    ### Priors ###
    ##################################
    gamma0 ~ dnorm(0, 0.01)
    gamma3 ~ dnorm(0, 0.01)
    gamma4 ~ dnorm(0, 0.01)
    
    ##################################
    ### Covered Cod End ###
    ##################################
    for (x in 1:nCover){
      logit(pDet[x]) <- gamma0 +
        gamma3 * covariate3_codend[x] + 
        gamma4 * covariate4_codend[x]
      
      counts_trawl[x] ~ dbin(pDet[x], counts_total[x])
    } #x
  } #end model
", file = "CH_4.2_Codend.txt")   


#<> Binomial N-mixture Only Analysis - Appropriate <>#
cat("
  model {
    ##################################
    ### Priors ###
    ##################################
    theta0 ~ dnorm(0, 0.01)
    theta2 ~ dnorm(0, 0.01)
    theta3 ~ dnorm(0, 0.01)
    theta4 ~ dnorm(0, 0.01)
    
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)

    ##################################
    ### Available and Retention Probabilities for Full Data ###
    ##################################
    for (i in 1:nSites){
      for (z in 1:nPrimary){
         for (r in 1:nSecondary){
          logit(p[i, z, r]) <- theta0 + 
                                  theta2 * covariate2_replicate[i, z, r] +
                                  theta3 * covariate3_replicate[i, z, r] +
                                  theta4 * covariate4_replicate[i, z, r]
         } #r
      } #z
    } #i
        
    ##################################
    ### Standard N-mixture ###
    ##################################
    for (i in 1:nSites){
      for (z in 1:nPrimary){
        for (r in 1:nSecondary){
          y[i, z, r] ~ dbin(p[i, z, r], N.close[i, z])
        } #r
    N.close[i, z] ~ dpois(lambda[i, z])
    log(lambda[i, z]) <- beta0 + 
                         beta1 * covariate1_replicate[i, z]
      } #z
    } #i
  } #end model
", file = "CH_4.2_Binomial_Appropriate.txt")   


#<> Binomial N-mixture Only Analysis - Ambiguous (Months as Closed Replicates) <>#
cat("
  model {
    ##################################
    ### Priors ###
    ##################################
    theta0 ~ dnorm(0, 0.01)
    theta2 ~ dnorm(0, 0.01)
    theta3 ~ dnorm(0, 0.01)
    theta4 ~ dnorm(0, 0.01)
    
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)

    ##################################
    ### Available and Retention Probabilities for Full Data ###
    ##################################
    for (i in 1:nSites){
      for (z in 1:nPrimary){
          logit(p[i, z]) <- theta0 +
                                  theta2 * covariate2_replicate[i, z, 1] +
                                  theta3 * covariate3_replicate[i, z, 1] +
                                  theta4 * covariate4_replicate[i, z, 1]
      } #z
    } #i
        
    ##################################
    ### Standard N-mixture ###
    ##################################
    for (i in 1:nSites){
      for (z in 1:nPrimary){
          y[i, z, 1] ~ dbin(p[i, z], N.single[i])
      } #z
      N.single[i] ~ dpois(lambda[i])
      log(lambda[i]) <- beta0 + 
                        beta1 * covariate1_replicate[i, 1]
    } #i
  } #end model
", file = "CH_4.2_Binomial_Ambiguous.txt")


#<> Integrated N-mixture Analysis - Appropriate <>#
cat("
  model {
    ##################################
    ### Priors ###
    ##################################
    gamma0 ~ dnorm(0, 0.01)
    gamma3 ~ dnorm(0, 0.01)
    gamma4 ~ dnorm(0, 0.01)
    
    alpha0 ~ dnorm(0, 0.01)
    alpha2 ~ dnorm(0, 0.01)

    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)

    ##################################
    ### Covered Cod End ###
    ##################################
    for (x in 1:nCover){
      logit(pDet[x]) <- gamma0 +
        gamma3 * covariate3_codend[x] + # shared with binomial
        gamma4 * covariate4_codend[x] # shared with binomial
      
      counts_trawl[x] ~ dbin(pDet[x], counts_total[x])
    } #x
    
    ##################################
    ### Available and Retention Probabilities for Full Data ###
    ##################################
    for (i in 1:nSites){
      for (z in 1:nPrimary){
        for (r in 1:nSecondary){
          p[i, z, r] <- pAvail[i, z, r] * pRet[i, z, r]
          
          logit(pAvail[i, z, r]) <- alpha0 + 
                                    alpha2 * covariate2_replicate[i, z, r]
          
          logit(pRet[i, z, r]) <- gamma0 +
                                  gamma3 * covariate3_replicate[i, z, r] +
                                  gamma4 * covariate4_replicate[i, z, r]
        } #r
      } #z
    } #i
        
    ##################################
    ### Standard N-mixture ###
    ##################################
    for (i in 1:nSites){
      for (z in 1:nPrimary){
        for (r in 1:nSecondary){
          y[i, z, r] ~ dbin(p[i, z, r], N.close[i, z])
        } #r
    N.close[i, z] ~ dpois(lambda[i, z])
    log(lambda[i, z]) <- beta0 + 
                         beta1 * covariate1_replicate[i, z]
      } #z
    } #i
  } #end model
", file = "CH_4.2_Integrated_Appropriate.txt")


#<> Integrated N-mixture Analysis - Ambiguous (Months as Closed Replicates) <>#
cat("
  model {
    ##################################
    ### Priors ###
    ##################################
    gamma0 ~ dnorm(0, 0.01)
    gamma3 ~ dnorm(0, 0.01)
    gamma4 ~ dnorm(0, 0.01)
    
    alpha0 ~ dnorm(0, 0.01)
    alpha2 ~ dnorm(0, 0.01)

    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)

    ##################################
    ### Covered Cod End ###
    ##################################
    for (x in 1:nCover){
      logit(pDet[x]) <- gamma0 +
        gamma3 * covariate3_codend[x] + # shared with binomial
        gamma4 * covariate4_codend[x] # shared with binomial
      
      counts_trawl[x] ~ dbin(pDet[x], counts_total[x])
    } #x
    
    ##################################
    ### Available and Retention Probabilities for Full Data ###
    ##################################
    for (i in 1:nSites){
      for (z in 1:nPrimary){
          p[i, z] <- pAvail[i, z] * pRet[i, z]
          
          logit(pAvail[i, z]) <- alpha0 + 
                                    alpha2 * covariate2_replicate[i, z, 1]
          
          logit(pRet[i, z]) <- gamma0 +
                                  gamma3 * covariate3_replicate[i, z, 1] +
                                  gamma4 * covariate4_replicate[i, z, 1]
      } #z
    } #i
        
    ##################################
    ### Standard N-mixture ###
    ##################################
    for (i in 1:nSites){
      for (z in 1:nPrimary){
          y[i, z, 1] ~ dbin(p[i, z], N.single[i])
      } #z
    N.single[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 + 
                         beta1 * covariate1_replicate[i, 1]
    } #i
  } #end model
", file = "CH_4.2_Integrated_Ambiguous.txt") 

#<> Model References <> #
model_ref <- paste("CH_4.2_", c("Codend", "Binomial_Appropriate", "Binomial_Ambiguous",
                                "Integrated_Appropriate", "Integrated_Ambiguous"), sep = "")

#<><><><><><><><><><><><>#
#<> Data Simulation <>#
#<><><><><><><><><><><><>#
nb <- 1000                                                                  # Parameters for JAGS 
ni <- 10000 + nb
nt <- 5
nc <- 3

nIter <- 100                                                                # Number of simulations
nSites <- 200 
nPrimary <- 3                                                               # Number of months within a season: Current assumption is that these are also closed replicates in the FMWT analysis
nSecondary <- 3                                                             # Number of true closed replicates

#<><> Observation and Latent State Processes <><>#
p.det <- 0.7                                                                # probability a fish is retained once in contact with the gear: Mean intercept not on logit scale
gamma.cov3 <- 0.8                                                           # Covariate 3 effect - Fish Length effect on the retention
gamma.cov4 <- 1.5                                                           # Covariate 4 effect - Biomass effect on retention

p.avail <- 0.4                                                              # probability a fish is present and available to the sampling gear: Mean intercept value not on logit scale
b.cov2 <- 0.5                                                               # Covariate 2 effect - Turbidity effect on availability

lambda <- 20                                                                # Super population abundance rate: Equivalent to an entire region within the delta
b.cov1 <- -1                                                                # Covariate 1 effect - Salt effect on abundance/survival

#<> Closed Cod end Parameters for N-mixture portion <>#
nClosed <- 100
covariate3_codend <- rbinom(nClosed, 1, prob = 0.5)                         # Covariate 3 for codend data - length dummy variable 
covariate4_codend <- rbinom(nClosed, 1, prob = 0.5)                         # Covariate 4 for codend data - Biomass dummy variable

pDet <- plogis(qlogis(p.det) + 
                 gamma.cov3 * covariate3_codend +                               # retention of just closed codend data
                 gamma.cov4 * covariate4_codend)

Navail <- rbinom(n = nClosed, size = lambda, prob = p.avail)                # Number fish available to the gear                            
counts_trawl <- rbinom(n = nClosed, size = Navail, prob = pDet)             # Independent closed cod end data generation based on a conceptual availability of fish
counts_cover <- Navail - counts_trawl                                       # Observed counts in the cod end cover

#<> Closed Cod end Parameters for N-mixture portion <>#
param <-                                                                    # Parameters to Monitor
  c("gamma0", "gamma3", "gamma4",                                           # p retained parameters                                  
    "alpha0", "alpha2",                                                     # p available parameters
    "theta0", "theta2", "theta3", "theta4",                                 # specific to the binomial model 
    "beta0", "beta1")                                                       # N parameters

codend_results <- array(NA, dim = c(nIter, 3, 4),                           # Placeholders for Results
                        dimnames = list(paste("Iter", 1:nIter, sep = "_"),
                                        param[1:3],
                                        c("Median", "Lower", "Upper", "Truth")))

binomial_appropriate_results <- array(NA, dim = c(nIter, 6, 4),   # Placeholder for binomial N-mixture only results
                          dimnames = list(paste("Iter", 1:nIter, sep = "_"),
                                          param[6:11],
                                          c("Median", "Lower", "Upper", "Truth")))

binomial_ambiguous_results <- array(NA, dim = c(nIter, 6, 4),   # Placeholder for binomial N-mixture only results
                                      dimnames = list(paste("Iter", 1:nIter, sep = "_"),
                                                      param[6:11],
                                                      c("Median", "Lower", "Upper", "Truth")))

integrate_appropriate_results <- array(NA, dim = c(nIter, 7, 4),     # Placeholder for Integrated model results
                        dimnames = list(paste("Iter", 1:nIter, sep = "_"),
                                        param[c(1:5, 10:11)],
                                        c("Median", "Lower", "Upper", "Truth")))

integrate_ambiguous_results <- array(NA, dim = c(nIter, 7, 4),     # Placeholder for Integrated model results
                                       dimnames = list(paste("Iter", 1:nIter, sep = "_"),
                                                       param[c(1:5, 10:11)],
                                                       c("Median", "Lower", "Upper", "Truth")))

#<> Simulate data for Replicated Data set <>#
start <- Sys.time()
for (iter in 1:nIter){
  print(iter)                                                                 # Keep track of iteration currently being simulated
  
  #<> Binomial N-mixture portion for availability <>#
  covariate2_rep <- rnorm(nSites, 0, 1)                                # Simulate covariate 2 (turbidity) for latent state variable
  covariate1_rep <- rnorm(nSites, 0, 1)                                # Simulate covariate 1 (salt) for latent state variable
  
  covariate2_replicate <- array(NA, dim = c(length(covariate2_rep), nPrimary, nSecondary))
  covariate2_replicate[, , 1] <- covariate2_rep
  
  covariate1_replicate <-  array(NA, dim = c(length(covariate1_rep), nPrimary))
  covariate1_replicate[, 1] <- covariate1_rep
  
  for (z in 2:nPrimary){
    covariate2_replicate[, , z] <- rnorm(nSites, covariate2_replicate[, , z - 1], 1) # make it dependent on previous years observation
    covariate1_replicate[, z] <- rnorm(nSites, covariate1_replicate[, z - 1], 1) # make it dependent on previous years observation
  } #z
  
  covariate3_replicate <- array(NA, dim = c(nSites, nPrimary, nSecondary))
  covariate3_replicate[, , 1] <- matrix(rbinom((nSites * nPrimary), 1, 0.5), ncol = nPrimary)
  covariate3_replicate[, , 2] <- covariate3_replicate[, , 1]
  covariate3_replicate[, , 3] <- covariate3_replicate[, , 1]
  
  covariate4_replicate <- array(NA, dim = c(nSites, nPrimary, nSecondary))
  covariate4_replicate[, , 1] <- matrix(rbinom((nSites * nPrimary), 1, 0.5), ncol = nPrimary)
  covariate4_replicate[, , 2] <- covariate4_replicate[, , 1]
  covariate4_replicate[, , 3] <- covariate4_replicate[, , 1]
  
  N <- matrix(rpois(nSites * nPrimary, exp(log(lambda) + 
                                            b.cov1 * covariate1_replicate)), ncol = nPrimary)                  # covariate 1 (Salt) Effect on abundance
  
  p_avail <- plogis(qlogis(p.avail) + b.cov2 * covariate2_replicate)                     # P availability ~ covariates
  p_det <- plogis(qlogis(p.det) + gamma.cov3 * covariate3_replicate +  
                    gamma.cov4 * covariate4_replicate)                                  # P detect ~ covariates via standard replicat survey
  p_all <- p_avail * p_det                                                    # Total P
  
  y <- p_all                                                                  # Observation Place Holder
  
    for (rep in 1:nSecondary) {
      for (primary in 1:nPrimary){
        y[, primary, rep] <-                                                  # observed counts for all subsequent primary periods
          rbinom(n = nSites, size = N[, primary], prob = p_all[, primary, rep])
      } #primary
    }
 
    ### Load Data Into JAGS List ###
    dat <- list(
      # For Spatial-Temporal Replicates #
      y = y,
      nSites = dim(y)[1],
      nPrimary = dim(y)[2],
      nSecondary = dim(y)[3], #either single secondary period or 3
      covariate2_replicate = covariate2_replicate,
      covariate3_replicate = covariate3_replicate,
      covariate4_replicate = covariate4_replicate,
      covariate1_replicate = covariate1_replicate,
      
      # For Covered Cod End #
      counts_trawl = counts_trawl,
      counts_total = (counts_cover + counts_trawl), 
      nCover = length(counts_trawl),
      covariate3_codend = covariate3_codend,
      covariate4_codend = covariate4_codend
    )
    
    ### Initial Values ###
    init <- function(){
      list(
        gamma0 = qlogis(p.det),
        gamma3 = gamma.cov3,
        gamma4 = gamma.cov4,
        
        alpha0 = qlogis(p.avail),
        alpha2 = b.cov2,
        
        N.close = apply(y, c(1, 2), max),
        N.single = apply(y, c(1), max),
        beta0 = log(lambda),
        beta1 = b.cov1
      )
    }
    
    truth <- c(qlogis(p.det), gamma.cov3, gamma.cov4, 
               qlogis(p.avail), b.cov2, 
               qlogis(p.det * p.avail), b.cov2, gamma.cov3, gamma.cov4, 
               log(lambda), b.cov1)
    
    ### Run JAGS ###
    #<> Closed Codend Results <>#
    out_codend <- jagsUI(dat, init, param, paste(model_ref[1], "txt", sep = "."), 
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
    
    codend_results[iter, , ] <- if(max(unlist(out_codend$Rhat)) > 1.1){    # only keep analyses where all parameters are less than or equal to 1.1
    NA }
    else {
      cbind(c(out_codend$q50$gamma0, out_codend$q50$gamma3,
              out_codend$q50$gamma4),
            c(out_codend$q2.5$gamma0, out_codend$q2.5$gamma3,
              out_codend$q2.5$gamma4),
            c(out_codend$q97.5$gamma0, out_codend$q97.5$gamma3,
              out_codend$q97.5$gamma4),
            truth[1:3])
    }
    
    #<> Only Binomial N-Mixture Model with appropriate design <>#
    out_nmix_appropriate <- jagsUI(dat, init, param, paste(model_ref[2], "txt", sep = "."),  
                                 n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
    
    binomial_appropriate_results[iter, ,] <- if(max(unlist(out_nmix_appropriate$Rhat)) > 1.1){    # only keep analyses where all parameters are less than or equal to 1.1
    NA }
    else {
      cbind(c(out_nmix_appropriate$q50$theta0, out_nmix_appropriate$q50$theta2, out_nmix_appropriate$q50$theta3,
              out_nmix_appropriate$q50$theta4, out_nmix_appropriate$q50$beta0, out_nmix_appropriate$q50$beta1),
            c(out_nmix_appropriate$q2.5$theta0, out_nmix_appropriate$q2.5$theta2, out_nmix_appropriate$q2.5$theta3,
              out_nmix_appropriate$q2.5$theta4, out_nmix_appropriate$q2.5$beta0, out_nmix_appropriate$q2.5$beta1),
            c(out_nmix_appropriate$q97.5$theta0, out_nmix_appropriate$q97.5$theta2, out_nmix_appropriate$q97.5$theta3,
              out_nmix_appropriate$q97.5$theta4, out_nmix_appropriate$q97.5$beta0, out_nmix_appropriate$q97.5$beta1), 
            truth[c(6:11)])
    }

    #<> Only Binomial N-Mixture Model with ambiguous design <>#
    out_nmix_ambiguous <- jagsUI(dat, init, param, paste(model_ref[3], "txt", sep = "."),  
                                   n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
    
    binomial_ambiguous_results[iter, ,] <- if(max(unlist(out_nmix_ambiguous$Rhat)) > 1.1){    # only keep analyses where all parameters are less than or equal to 1.1
    NA }
    else {
      cbind(c(out_nmix_ambiguous$q50$theta0, out_nmix_ambiguous$q50$theta2, out_nmix_ambiguous$q50$theta3,
              out_nmix_ambiguous$q50$theta4, out_nmix_ambiguous$q50$beta0, out_nmix_ambiguous$q50$beta1),
            c(out_nmix_ambiguous$q2.5$theta0, out_nmix_ambiguous$q2.5$theta2, out_nmix_ambiguous$q2.5$theta3,
              out_nmix_ambiguous$q2.5$theta4, out_nmix_ambiguous$q2.5$beta0, out_nmix_ambiguous$q2.5$beta1),
            c(out_nmix_ambiguous$q97.5$theta0, out_nmix_ambiguous$q97.5$theta2, out_nmix_ambiguous$q97.5$theta3,
              out_nmix_ambiguous$q97.5$theta4, out_nmix_ambiguous$q97.5$beta0, out_nmix_ambiguous$q97.5$beta1), 
            truth[c(6:11)])
    }   
    
    #<> Integrated N-Mixture Model with appropriate design <>#
    out_int_appropriate <- jagsUI(dat, init, param, paste(model_ref[4], "txt", sep = "."),  
                                 n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
    
    integrate_appropriate_results[iter, ,] <- if(max(unlist(out_int_appropriate$Rhat)) > 1.1){    # only keep analyses where all parameters are less than or equal to 1.1
    NA }
    else {
      cbind(c(out_int_appropriate$q50$gamma0, out_int_appropriate$q50$gamma3, out_int_appropriate$q50$gamma4,
              out_int_appropriate$q50$alpha0, out_int_appropriate$q50$alpha2, out_int_appropriate$q50$beta0,
              out_int_appropriate$q50$beta1),
            c(out_int_appropriate$q2.5$gamma0, out_int_appropriate$q2.5$gamma3, out_int_appropriate$q2.5$gamma4,
              out_int_appropriate$q2.5$alpha0, out_int_appropriate$q2.5$alpha2, out_int_appropriate$q2.5$beta0,
              out_int_appropriate$q2.5$beta1),
            c(out_int_appropriate$q97.5$gamma0, out_int_appropriate$q97.5$gamma3, out_int_appropriate$q97.5$gamma4,
              out_int_appropriate$q97.5$alpha0, out_int_appropriate$q97.5$alpha2, out_int_appropriate$q97.5$beta0,
              out_int_appropriate$q97.5$beta1),
            truth[c(1:5, 10:11)])
    }      
    
    #<> Integrated N-Mixture Model with ambiguous design <>#
    out_int_ambiguous <- jagsUI(dat, init, param, paste(model_ref[5], "txt", sep = "."),  
                                  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
    
    integrate_ambiguous_results[iter, ,] <- if(max(unlist(out_int_ambiguous$Rhat)) > 1.1){    # only keep analyses where all parameters are less than or equal to 1.1
    NA }
    else {
      cbind(c(out_int_ambiguous$q50$gamma0, out_int_ambiguous$q50$gamma3, out_int_ambiguous$q50$gamma4,
              out_int_ambiguous$q50$alpha0, out_int_ambiguous$q50$alpha2, out_int_ambiguous$q50$beta0,
              out_int_ambiguous$q50$beta1),
            c(out_int_ambiguous$q2.5$gamma0, out_int_ambiguous$q2.5$gamma3, out_int_ambiguous$q2.5$gamma4,
              out_int_ambiguous$q2.5$alpha0, out_int_ambiguous$q2.5$alpha2, out_int_ambiguous$q2.5$beta0,
              out_int_ambiguous$q2.5$beta1),
            c(out_int_ambiguous$q97.5$gamma0, out_int_ambiguous$q97.5$gamma3, out_int_ambiguous$q97.5$gamma4,
              out_int_ambiguous$q97.5$alpha0, out_int_ambiguous$q97.5$alpha2, out_int_ambiguous$q97.5$beta0,
              out_int_ambiguous$q97.5$beta1),
            truth[c(1:5, 10:11)])
    }
} #iter
end <- Sys.time()
end - start # 10.3 hours for 100 simulations

#save(codend_results, file = "codend_results_1_100.RData")
#save(binomial_appropriate_results, file = "binomial_appropriate_results_101_200.RData")
#save(binomial_ambiguous_results, file = "binomial_ambiguous_results_201_210.RData")
#save(integrate_appropriate_results, file = "integrate_appropriate_results_101_250.RData")
#save(integrate_ambiguous_results, file = "integrate_ambiguous_results_101_250.RData")
    
#<><><><><><><><><><><><>#
#<> Load Results and Make Figure 4.5 and Table 4.1 <>#
#<><><><><><><><><><><><>#
#<><><> Covered Codend Results <><><>#
load("codend_results_1_100.RData")
codend_results <- as.data.frame.table(codend_results) %>%
  spread(., Var3, Freq) %>%
  na.omit(.) %>%
  mutate(Model = "Cod end",
         Secondary = "Appropriate",
         Difference = Median - Truth,
         RBias = (Median - Truth) / Truth) %>%
  select(Model, Iteration = Var1, Parameter = Var2, Secondary, Median:Truth) %>%
  rbind.data.frame(., .) %>%
  transform(Secondary = rep(c("Appropriate", "Ambiguous"), each = (nrow(.))/ 2))
  
#<><><> Binomial N-Mixture Only Results <><><>#
load("binomial_appropriate_results_1_320.RData")
binomial_approp <- as.data.frame.table(binomial_appropriate_results) %>%
  na.omit(.) %>%
  spread(., Var3, Freq) %>%
  arrange(Var1, Var2) %>%
  group_by(Var2) %>%
  mutate(ID = as.numeric(as.factor(as.numeric(Var1)))) %>% # Makes iterations consecutive to filter the first 100
  ungroup() %>%
  filter(ID < 101) %>%
  mutate(Model = "Binomial",
         Secondary = "Appropriate",
         Difference = Median - Truth,
         RBias = (Median - Truth) / Truth) %>%
  select(Model, Iteration = Var1, Parameter = Var2, Secondary, Median:Truth)

load("binomial_ambiguous_results_1_320.RData")
binomial_amb <- as.data.frame.table(binomial_ambiguous_results) %>%
  na.omit(.) %>%
  spread(., Var3, Freq) %>%
  arrange(Var1, Var2) %>%
  group_by(Var2) %>%
  mutate(ID = as.numeric(as.factor(as.numeric(Var1)))) %>% # Makes iterations consecutive to filter the first 100
  ungroup() %>%
  filter(ID < 101) %>%
  mutate(Model = "Binomial",
         Secondary = "Ambiguous",
         Difference = Median - Truth,
         RBias = (Median - Truth) / Truth) %>%
  select(Model, Iteration = Var1, Parameter = Var2, Secondary, Median:Truth)


#<><><> Integrated Model Results <><><>#
load("integrate_appropriate_results_1_100.RData")
int_approp1 <- as.data.frame.table(integrate_appropriate_results)

load("integrate_appropriate_results_101_250.RData")
integrate_approp <- as.data.frame.table(integrate_appropriate_results) %>%
  rbind.data.frame(int_approp1, .) %>%
  na.omit(.) %>%
  spread(., Var3, Freq) %>%
  arrange(Var1, Var2) %>%
  group_by(Var2) %>%
  mutate(ID = as.numeric(as.factor(as.numeric(Var1)))) %>% # Makes iterations consecutive to filter the first 100
  ungroup() %>%
  filter(ID < 101) %>%
  mutate(Model = "Integrated",
         Secondary = "Appropriate",
         Difference = Median - Truth,
         RBias = (Median - Truth) / Truth) %>%
  select(Model, Iteration = Var1, Parameter = Var2, Secondary, Median:Truth)

load("integrate_ambiguous_results_1_100.RData")
int_amb1 <- as.data.frame.table(integrate_ambiguous_results)

load("integrate_ambiguous_results_101_250.RData")
integrate_amb <- as.data.frame.table(integrate_ambiguous_results) %>%
  rbind.data.frame(int_amb1, .) %>%
  na.omit(.) %>%
  spread(., Var3, Freq) %>%
  arrange(Var1, Var2) %>%
  group_by(Var2) %>%
  mutate(ID = as.numeric(as.factor(as.numeric(Var1)))) %>% # Makes iterations consecutive to filter the first 100
  ungroup() %>%
  filter(ID < 101) %>%
  mutate(Model = "Integrated",
         Secondary = "Ambiguous",
         Difference = Median - Truth,
         RBias = (Median - Truth) / Truth) %>%
  select(Model, Iteration = Var1, Parameter = Var2, Secondary, Median:Truth)


#<><><> All Model Results <><><>#
all_results <- rbind.data.frame(codend_results, binomial_approp, binomial_amb,
                                integrate_approp, integrate_amb)

sample_size <- aggregate(Median ~ Model + Parameter + Secondary, all_results, length) # check sample size

#<> RMSE <>#
rmse <- mutate(all_results, 
               Difference = Median - Truth,
               RMSE = Difference^2) %>%
  group_by(Model, Parameter, Secondary) %>%
  summarize(RMSE = sqrt(sum(RMSE) / 100)) %>%
  ungroup() %>%
  data.frame(.)

rmse_table <- spread(rmse, Secondary, RMSE)

#<> 95% CI Coverage (Amundson et al. 2014) <>#
coverage <- mutate(all_results,
                   Coverage = ifelse(Truth < Upper & Truth > Lower, 1, 0)) %>%
  group_by(Model, Parameter, Secondary) %>%
  summarize(Coverage = (sum(Coverage) / 100)) %>%
  ungroup() %>%
  data.frame(.)

#<> Table 4.1 <>#
table_4.1 <- merge(rmse, coverage, all = T) %>%
  gather(., Response, Value, RMSE:Coverage) %>%
  transform(Secondary = factor(Secondary, levels = c("Appropriate", "Ambiguous")),
            Response = factor(Response, levels = c("Coverage", "RMSE")),
            Model = factor(Model, levels = c("Cod end", "Binomial", "Integrated"))) %>%
  spread(., Response, Value) %>%
  arrange(Parameter, Model, Secondary)

#<> Figure 4.5 <>#
full_results_fig <- all_results %>%
  mutate(Parameter = recode(Parameter, "theta0" = "theta[0]", "theta2" = "theta[Covariate2]",
                            "theta3" = "theta[Covariate3]", "theta4" = "theta[Covariate4]",
                            "gamma0" = "gamma[0]", "gamma3" = "gamma[Covariate3]", 
                            "gamma4" = "gamma[Covariate4]",
                            "alpha0" = "alpha[0]", "alpha2" = "alpha[Covariate2]",
                            "beta0" = "beta[0]", "beta1" = "beta[Covariate1]")) %>%
  transform(Parameter = factor(Parameter, levels = c("gamma[0]", "gamma[Covariate3]", "gamma[Covariate4]",
                                                     "alpha[0]", "alpha[Covariate2]", " ",
                                                     "beta[0]", "beta[Covariate1]", "  ")),
            Secondary = factor(Secondary, levels = c("Appropriate", "Ambiguous")))

figure_4.5 <- ggplot(data = filter(full_results_fig, Model == "Integrated"), 
                     aes(x = Secondary, y = Median)) + #, fill = Secondary
  geom_hline(aes(yintercept = Truth), lwd = 2, lty = 3) +
  geom_boxplot(fill = "light gray") +
  labs(x = "Secondary sampling period") +
  facet_wrap(. ~ Parameter, ncol = 3, scale = "free_y", drop = FALSE,
             labeller = label_parsed) +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.y = element_text(vjust = 3),
        axis.title.x = element_text(vjust = -2.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        legend.position = "none", #c(0.08, 0.8), 
        legend.key = element_blank(),
        legend.spacing.x = unit(1.0, 'mm'),
        legend.spacing.y = unit(1.0, "mm"), # adjusts position between legend themes (title and rest)
        legend.key.size = unit(0.7, "lines"), # Adjusts spacing between legend symbols
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.box.background = element_rect(colour = "black", size = 1),
        legend.box = "horizontal",
        plot.title = element_text(vjust = -18, hjust = 0.95)) +
  theme(plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))

figure_4.5 <- ggplotGrob(figure_4.5)
rm_grobs <- figure_4.5$layout$name %in% c("panel-2-3", "panel-3-3", "strip-t-3-2", "strip-t-3-3")
figure_4.5$grobs[rm_grobs] <- NULL
figure_4.5$layout <- figure_4.5$layout[!rm_grobs, ]
figure_4.5$layout[figure_4.5$layout$name == "axis-b-3-3", c("t", "b")] <- c(11, 11)
grid.newpage()
grid.draw(figure_4.5)


#<><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
#<> Workflow 4.3: Application of the Integrated Abundance Model to the FMWT <>#
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
#<><><><><><><><><><><><><><><><><><>#
#<> Longterm Data Cleaning: FMWT - Length is FL <>#
#<><><><><><><><><><><><><><><><><><>#
require(LTMRdata)

ltmr <- data.frame(LTMRdata::FMWT) %>%
  transform(Taxa = ifelse(is.na(Taxa) & Length_NA_flag == "No fish caught", "None", Taxa)) %>%
  mutate(Counts = ifelse(Taxa == "None", 0, Count),
         Taxa = gsub(" ", "_", Taxa), # replace space with _
         Year = year(Date),
         Month = month(Date),
         Season = ifelse(Month > 2 & Month < 6, "Spring",
                         ifelse(Month > 5 & Month < 9, "Summer",
                                ifelse(Month > 8 & Month < 12, "Fall", "Winter"))),
         WaterYear = ifelse(Month > 11, Year + 1, Year), # california water year runs from Oct1 to Sept 30. However only use December because it is the start of the next season (Sept is still fall, which would force fall to be split between water years)
         Size = ifelse(Length < 50.1, "Small", "Large"),
         Region = as.numeric(Station),
         RegionR = ifelse(Region < 200 & Region > 99, "South Bay",
                          ifelse(Region > 199 & Region < 300, "Central Bay",
                                 ifelse(Region > 299 & Region < 400, "San Pablo",
                                        ifelse(Region < 600 & Region > 499, "Honker",
                                               ifelse(Region < 712 & Region > 699, "Lower Sacramento",
                                                      ifelse(Region < 800 & Region > 711, "Upper Sacramento",
                                                             ifelse(Region < 900 & Region > 799, "Lower San Joaquin",
                                                                    ifelse(Region > 899, "Upper San Joaquin", 
                                                                           ifelse(Region < 100, "Upper Sacramento", "Suisun")))))))))) %>% 
  filter(!is.na(Tow_volume), # remove samples without tow volume information
         WaterYear < 2019,
         WaterYear > 1989) %>% # remove 2019 and after because full sample has yet to be collected and before 1990 because that is what polansky et al. (2019) used
  select(Source:Date, Year = WaterYear, Month:Season, Region = RegionR, Tow_volume, Tow_direction, SampleID:Secchi, Taxa, Size, Length, Counts)

#<> LMTR Counts <>#
cluster <- c("Morone_saxatilis", "Platichthys_stellatus", "Alosa_sapidissima", #striped bass, starry flounder, american shad, sac sucker, topsmelt, pacific sardine
             "Catostomus_occidentalis", "Atherinops_affinis", "Sardinops_sagax") # Fall only cluster if using more than one species

ltmr_count <- select(ltmr, Station:Station, Date:Region, Taxa:Size, Counts) %>%
  group_by(Station, Date, Year, Month, Season, Region, Taxa, Size) %>%
  summarize(Counts = sum(Counts)) %>%
  ungroup() %>%
  spread(., Size, Counts) %>%
  select(-`<NA>`) %>%
  transform(Large = replace_na(Large, 0),
            Small = replace_na(Small, 0)) %>%
  gather(., Size, Counts, Large:Small) %>%
  spread(., Taxa, Counts) %>%
  select(Station:Size, cluster) %>%    # select attributes of the data frame and species identified by the cluster object above
  transform(Morone_saxatilis = replace_na(Morone_saxatilis, 0),
            Platichthys_stellatus = replace_na(Platichthys_stellatus, 0),
            Alosa_sapidissima = replace_na(Alosa_sapidissima, 0),
            Catostomus_occidentalis = replace_na(Catostomus_occidentalis, 0),
            Atherinops_affinis = replace_na(Atherinops_affinis, 0),
            Sardinops_sagax = replace_na(Sardinops_sagax, 0),
            Size = factor(Size, levels = c("Small", "Large"))) %>%
  arrange(Station, Date, Region, Size)

#<> Stations To Remove <>#
station_check <- distinct(select(ltmr_count, Station, Date)) %>%
  mutate(count = 1) %>%
  group_by(Station) %>%
  summarize(Count = sum(count)) %>%
  ungroup()

figure_4.7 <- ggplot(data = station_check) +
  geom_histogram(aes(x = Count), fill = "light blue", colour = "black") +
  labs(x = "Number of samples collected", y = "Frequency") +
  theme(legend.position = c(0.85, 0.15), 
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(vjust = 3),
        axis.title.x = element_text(vjust = -2.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        legend.key = element_blank(),
        legend.text = element_text(size = 14),
        legend.box.background = element_rect(colour = "black", size = 1),
        text = element_text(size = 16),
        plot.title = element_text(vjust = -18, hjust = 0.95),
        plot.margin = unit(c(0, 0, 0.05, 0.03), 'null'))

#<> Remove Stations <>#
station_check <- filter(station_check, Count > 99)
ltmr_count <- filter(ltmr_count, Station %in% station_check$Station) 

#<><><> LMTR Covariates <><><>#
# Replace all missing values with means based on an order of grouping
# (1) Station, Year, and Month; (2) Region, Year, and Month; (3) Station, Year, and Season; 
# (4) Region, Year, and Season; (5) Station and Year; (6) Region and Year
ltmr_cov <- select(ltmr, Station:Tow_direction, Tide:Secchi, Length) %>%
  mutate(Length = replace_na(Length, 0)) %>%
  group_by(Station, Latitude, Longitude, Date, Year, Season, Month, Region) %>% # Replace Covariates that are missing with means by station, and month
  summarize(Tow_volume = mean(Tow_volume, na.rm = T),
            Tow_direction = unique(Tow_direction),
            Tide = unique(Tide),
            Salinity = mean(Sal_surf, na.rm = T),
            Temperature = mean(Temp_surf, na.rm = T),
            Secchi = mean(Secchi, na.rm = T),
            Length = sum(Length, na.rm = T)) %>% # Proxy of biomass in the net based on summed lengths
  ungroup() %>%
  data.frame(.) %>%
  group_by(Region, Year, Season, Month) %>% # Replace NA's with Region, year, season and Month specific means (No Station unlike the previous grouping)
  mutate(Salt_group = mean(Salinity, na.rm = T),
         Temp_group = mean(Temperature, na.rm = T),
         Secchi_group = mean(Secchi, na.rm = T)) %>%
  ungroup() %>%
  transform(Salinity = ifelse(is.na(Salinity), Salt_group, Salinity),
            Temperature = ifelse(is.na(Temperature), Temp_group, Temperature),
            Secchi = ifelse(is.na(Secchi), Secchi_group, Secchi)) %>%
  select(Station:Length) %>%
  group_by(Station, Year, Season) %>% # Replace NA's that weren't replaced by region: Use Station, Year, Season
  mutate(Salt_group = mean(Salinity, na.rm = T),
         Temp_group = mean(Temperature, na.rm = T),
         Secchi_group = mean(Secchi, na.rm = T)) %>%
  ungroup() %>%
  transform(Salinity = ifelse(is.na(Salinity), Salt_group, Salinity),
            Temperature = ifelse(is.na(Temperature), Temp_group, Temperature),
            Secchi = ifelse(is.na(Secchi), Secchi_group, Secchi)) %>%
  select(Station:Length) %>%
  group_by(Region, Year, Season) %>% # Replace NA's that weren't replaced: Use Region, Year, Season: All upper and lower sac in 2009 for salinity
  mutate(Salt_group = mean(Salinity, na.rm = T),
         Temp_group = mean(Temperature, na.rm = T),
         Secchi_group = mean(Secchi, na.rm = T)) %>%
  ungroup() %>%
  transform(Salinity = ifelse(is.na(Salinity), Salt_group, Salinity),
            Temperature = ifelse(is.na(Temperature), Temp_group, Temperature),
            Secchi = ifelse(is.na(Secchi), Secchi_group, Secchi)) %>%
  select(Station:Length) %>%
  group_by(Station, Year) %>% # Replace NA's that weren't replaced: Use Station and Year
  mutate(Salt_group = mean(Salinity, na.rm = T),
         Temp_group = mean(Temperature, na.rm = T),
         Secchi_group = mean(Secchi, na.rm = T)) %>%
  ungroup() %>%
  group_by(Region, Year) %>% # Replace NA's that weren't replaced: Use Region and Year: All upper and lower sac in 2009 for Turbidity (winter)
  mutate(Salt_group = mean(Salinity, na.rm = T),
         Temp_group = mean(Temperature, na.rm = T),
         Secchi_group = mean(Secchi, na.rm = T)) %>%
  ungroup() %>%
  transform(Salinity = ifelse(is.na(Salinity), Salt_group, Salinity),
            Temperature = ifelse(is.na(Temperature), Temp_group, Temperature),
            Secchi = ifelse(is.na(Secchi), Secchi_group, Secchi)) %>%
  select(Station:Length)

#<><><> LMTR Brought Back to One Data Frame <><><>#
ltmr_full <- merge(ltmr_count, ltmr_cov, all.x = T) %>%
  group_by(Station, Year, Season, Size) %>%
  mutate(Rep = row_number()) %>%
  ungroup() %>%
  select(Station:Region, Latitude:Longitude, Rep, Tow_volume:Secchi, TotalFL = Length,
         Size:Sardinops_sagax) %>%
  transform(Station = as.numeric(Station)) %>%
  data.frame(.)

#<><><><><><><><><><><><><><><><><><>#
#<> Covered Cod End Data <>#
#<><><><><><><><><><><><><><><><><><>#
codend <- read.csv("mitchel_covered_codend.csv", header = T) %>%
  filter(GearTypeCode == "FMWTC2B", # 2 types of trawls, this is the more common (2671 vs 304)
         Volume > 0) %>% # remove samples where volume was not recorded
  mutate(Size = ifelse(ForkLength < 50.1 & ForkLength > 0, "Small",
                       ifelse(ForkLength > 50, "Large", "None"))) %>%
  filter(Size != "None")

#<><><> Density Covariate (Sum of Fork Lengths) <><><>#
density <- select(codend, SampleTimeStart, TowDuration:Volume, CodendType, 
                  OrganismCode:Catch, TotalFL = ForkLength) %>%
  group_by(SampleTimeStart, CodendType) %>%
  summarize(TotalFL = sum(TotalFL)) %>%
  ungroup() %>%
  filter(CodendType == "Inside")

#<><><> Volume <><><>#
volume_cover <- select(codend, SampleTimeStart, Volume, CodendType) %>%
  group_by(SampleTimeStart, CodendType) %>%
  summarize(Volume = mean(Volume)) %>%
  ungroup() %>%
  filter(CodendType == "Inside")

#<><><> All Fish General <><><>#
all_fish <- select(codend, SampleTimeStart, CodendType, Size) %>%
  mutate(Count = 1) %>%
  group_by(SampleTimeStart, CodendType, Size) %>%
  summarize(Count = sum(Count)) %>%
  spread(., CodendType, Count) %>%
  mutate(Species = "All",
         Outside = replace_na(Outside, 0),
         Inside = replace_na(Inside, 0),
         Total = Inside + Outside,
         Proportion = Inside / Total) %>%
  data.frame(.) %>%
  merge(., select(density, -CodendType), all.x = T) %>%
  merge(., select(volume_cover, -CodendType), all.x = T)

#<> Covered Codend Total Fish Selectivity Curve <>#
mitchell_selectivity_plot <- ggplot(data = all_fish, aes(x = TotalFL, y = Proportion, colour = Size)) +
  geom_point(size = 10, alpha = 0.5) +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),
              lty = 1, lwd = 2,
              se = F) +
  labs(x = "Cumulative fish length (mm)", y = "Proportion of fish captured by the trawl") +
  theme(legend.position = c(0.85, 0.15), 
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(vjust = 3),
        axis.title.x = element_text(vjust = -2.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        legend.key = element_blank(),
        legend.text = element_text(size = 14),
        legend.box.background = element_rect(colour = "black", size = 1),
        text = element_text(size = 16),
        plot.title = element_text(vjust = -18, hjust = 0.95),
        plot.margin = unit(c(0, 0, 0.05, 0.03), 'null'))

#<><><><><><><><><><><><><><><><><><>#
#<> Create American Shad Capture History From LTMR Data <>#
#<><><><><><><><><><><><><><><><><><>#
#<><><> American Shad <><><>#
Ashad_ltmr <- select(ltmr_full, Station, Year, Season:Region, Rep, Size, Alosa_sapidissima) %>%
  spread(., Rep, Alosa_sapidissima)

#<> Capture Histories and Grouping Factors <>#
Ashad_ch <- select(Ashad_ltmr, `1`:`5`) %>%
  data.frame(.)
ch <- array(NA, dim = c(nrow(Ashad_ch), ncol(Ashad_ch))) #Won't work unless data is formatted this way
ch <- as.matrix(Ashad_ch)
ch[, 1] <- as.integer(Ashad_ch[, 1])
ch[, 2] <- as.integer(Ashad_ch[, 2])
ch[, 3] <- as.integer(Ashad_ch[, 3])
ch[, 4] <- as.integer(Ashad_ch[, 4])
ch[, 5] <- as.integer(Ashad_ch[, 5])

na_ref <- apply(Ashad_ch, 1, firstNA) # Reference to Last observation of the replicate used for Nmixture
Ashad_station <- Ashad_ltmr$Station
Ashad_year <- Ashad_ltmr$Year
Ashad_season <- Ashad_ltmr$Season
Ashad_region <- Ashad_ltmr$Region
Ashad_size <- Ashad_ltmr$Size

#<> Covariates by Replicates <>#
Ashad_volume <- select(ltmr_full, Station, Year, Season:Region, Rep, Size, Tow_volume) %>%
  spread(., Rep, Tow_volume)

Ashad_turb <- select(ltmr_full, Station, Year, Season:Region, Rep, Size, Secchi) %>%
  spread(., Rep, Secchi)

Ashad_tide <- select(ltmr_full, Station, Year, Season:Region, Rep, Size, Tide) %>%
  spread(., Rep, Tide)

Ashad_tow <- select(ltmr_full, Station, Year, Season:Region, Rep, Size, Tow_direction) %>%
  spread(., Rep, Tow_direction)

Ashad_density <- select(ltmr_full, Station, Year, Season:Region, Rep, Size, TotalFL) %>%
  spread(., Rep, TotalFL)

Ashad_salt <- select(ltmr_full, Station, Year, Season:Region, Rep, Size, Salinity) %>%
  spread(., Rep, Salinity)

#<><><><><><><><><><><><><><><><><><>#
#<> JAGS Models <>#
#<><><><><><><><><><><><><><><><><><>#
#<> Model 1 (Equation 9): GLM Without Integration <>#
cat("
model {
  # Priors #
  beta0 ~ dnorm(0, 0.001)
  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0, 0.001)
  beta3 ~ dnorm(0, 0.001)
  beta4 ~ dnorm(0, 0.001)
  beta5 ~ dnorm(0, 0.001)

  tauS <- pow(sigmaS, -2)
  sigmaS ~ dunif(0, 1000)
  for (z in 1:nStation){
    epsS[z] ~ dnorm(0, tauS)
  }
  
  tauY <- pow(sigmaY, -2)
  sigmaY ~ dunif(0, 1000)
  for (t in 1:nYear){
    epsY[t] ~ dnorm(0, tauY)
  }
  
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 1000)

  # Likelihood and Constraints #
  for (i in 1:nobs){
    eps[i] ~ dnorm(0, tau)
    
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 +
        beta1 * sizeL[i] +
        beta2 * densityL[i, 1] +
        beta3 * effortL[i, 1] +
        beta4 * turbidityL[i, 1] +
        beta5 * saltL[i] +
        epsS[station[i]] +      # random station effect
        epsY[year[i]] +         # random year effect
        eps[i]                  # random observation effect

    for (z in 1:nona[i]){ # used since replicate numbers per row is unequal
      count[i, z] ~ dpois(lambda[i])
      
  #######################################
  ### BAYESIAN P-VALUE ###
  #######################################
      y_new[i, z] ~ dpois(lambda[i])
      resid_actual[i, z] <- pow((count[i, z] - N[i]), 2)
      resid_new[i, z] <- pow((y_new[i, z] - N[i]), 2)
    } #z
    sum_act[i] <- sum(resid_actual[i, 1:nona[i]])
    sum_sim[i] <- sum(resid_new[i, 1:nona[i]])
  } #i
  fit <- sum(sum_act)
  fit_new <- sum(sum_sim)
  c_hat <- fit_new + 0.1 / fit + 0.1 # Need to add a small constant to both to make sure fit is not 0
  bpv <- step(fit_new - fit)
} #end model
", file = "CH_4.3_model1.txt")


#<> Model 2 (Equation 10): Integrated model with Turbidity on P_available <>#
cat("
model {
##################################
### Shared Parameters ###
##################################
  gamma0 ~ dnorm(0, 0.001)
  gamma1 ~ dnorm(0, 0.001)
  gamma2 ~ dnorm(0, 0.001)

##################################
### Covered Cod End ###
##################################
for (x in 1:nCover){
   logit(pDet[x]) <- gamma0 +
                     gamma1 * sizeT[x] +
                     gamma2 * densityT[x]

   cover_trawl[x] ~ dbin(pDet[x], cover_total[x])
} #x

##################################
### FMWT Binomial N-Mixture ###
##################################
  # Priors #
  alpha0 ~ dnorm(0, 0.001)
  alpha1 ~ dnorm(0, 0.001)
  alpha2 ~ dnorm(0, 0.001)
  
  beta0 ~ dnorm(0, 0.001)
  beta1 ~ dnorm(0, 0.001)
  
  tauS <- pow(sigmaS, -2)
  sigmaS ~ dunif(0, 1000)
  for (z in 1:nStation){
    epsS[z] ~ dnorm(0, tauS)
  }
  
  tauY <- pow(sigmaY, -2)
  sigmaY ~ dunif(0, 1000)
  for (t in 1:nYear){
    epsY[t] ~ dnorm(0, tauY)
  }
  
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 1000)
  tauP <- pow(sigmaP, -2)
  sigmaP ~ dunif(0, 1000)

  # Likelihood and Constraints #
  for (i in 1:nobs){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 +
      beta1 * saltL[i] +
      epsS[station[i]] +      # random station effect
      epsY[year[i]] +         # random year effect
      eps[i]                  # random observation effect
    
    eps[i] ~ dnorm(0, tau)

    for (z in 1:nona[i]){ # used since replicate numbers per row is unequal
      logit(pRet[i, z]) <- gamma0 + # must predict what retention efficiency would be for LTMR data: same as pDet
                    gamma1 * sizeL[i] +
                    gamma2 * densityL[i, z]

      p[i, z] <- pAvail[i, z] * pRet[i, z] # full detection efficiency is the product of retention efficency and availabiltiy

      logit(pAvail[i, z]) <- alpha0 + # Constrain availability by additional covariates we believe will effect just availability
        alpha1 * effortL[i, z] +
        alpha2 * turbidityL[i, z] +
        epsP[i, z]
        
      epsP[i, z] ~ dnorm(0, tauP)

      count[i, z] ~ dbin(p[i, z], N[i])
      
  #######################################
  ### BAYESIAN P-VALUE ###
  #######################################
      y_new[i, z] ~ dbin(p[i, z], N[i])
      e_count[i, z] <- p[i, z] * N[i]
      chi_actual[i, z] <- pow((count[i, z] - e_count[i, z]), 2) / (e_count[i, z] + 0.1)
      chi_sim[i, z] <- pow((y_new[i, z] - e_count[i, z]), 2) / (e_count[i, z] + 0.1)
    } #z
    sum_act[i] <- sum(chi_actual[i, 1:nona[i]])
    sum_sim[i] <- sum(chi_sim[i, 1:nona[i]])
  } #i
  fit <- sum(sum_act)
  fit_new <- sum(sum_sim)
  c_hat <- fit_new / fit
  bpv <- step(fit_new - fit)
} #end model
", file = "CH_4.3_model2.txt")


#<> Model 3 (Equation 11): Integrated model with Turbidity on Abundance <>#
cat("
model {
##################################
### Shared Parameters ###
##################################
  gamma0 ~ dnorm(0, 0.001)
  gamma1 ~ dnorm(0, 0.001)
  gamma2 ~ dnorm(0, 0.001)

##################################
### Covered Cod End ###
##################################
for (x in 1:nCover){
   logit(pDet[x]) <- gamma0 +
                     gamma1 * sizeT[x] +
                     gamma2 * densityT[x]

   cover_trawl[x] ~ dbin(pDet[x], cover_total[x])
} #x

##################################
### FMWT Binomial N-Mixture ###
##################################
  # Priors #
  alpha0 ~ dnorm(0, 0.001)
  alpha1 ~ dnorm(0, 0.001)

  beta0 ~ dnorm(0, 0.001)
  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0, 0.001)
  
  tauS <- pow(sigmaS, -2)
  sigmaS ~ dunif(0, 1000)
  for (z in 1:nStation){
    epsS[z] ~ dnorm(0, tauS)
  }
  
  tauY <- pow(sigmaY, -2)
  sigmaY ~ dunif(0, 1000)
  for (t in 1:nYear){
    epsY[t] ~ dnorm(0, tauY)
  }
  
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 1000)
  tauP <- pow(sigmaP, -2)
  sigmaP ~ dunif(0, 1000)

  # Likelihood and Constraints #
  for (i in 1:nobs){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 +
      beta1 * saltL[i] +
      beta2 * turbidityL[i, 1] +
      epsS[station[i]] +      # random station effect
      epsY[year[i]] +         # random year effect
      eps[i]                  # random observation effect
    
    eps[i] ~ dnorm(0, tau)

    for (z in 1:nona[i]){ # used since replicate numbers per row is unequal
      logit(pRet[i, z]) <- gamma0 + # must predict what retention efficiency would be for LTMR data: same as pDet
                    gamma1 * sizeL[i] +
                    gamma2 * densityL[i, z]

      p[i, z] <- pAvail[i, z] * pRet[i, z] # full detection efficiency is the product of retention efficency and availabiltiy

      logit(pAvail[i, z]) <- alpha0 + # Constrain availability by additional covariates we believe will effect just availability
        alpha1 * effortL[i, z] +
        epsP[i, z]
        
      epsP[i, z] ~ dnorm(0, tauP)

      count[i, z] ~ dbin(p[i, z], N[i])
      
  #######################################
  ### BAYESIAN P-VALUE ###
  #######################################
      y_new[i, z] ~ dbin(p[i, z], N[i])
      e_count[i, z] <- p[i, z] * N[i]
      chi_actual[i, z] <- pow((count[i, z] - e_count[i, z]), 2) / (e_count[i, z] + 0.1)
      chi_sim[i, z] <- pow((y_new[i, z] - e_count[i, z]), 2) / (e_count[i, z] + 0.1)
    } #z
    sum_act[i] <- sum(chi_actual[i, 1:nona[i]])
    sum_sim[i] <- sum(chi_sim[i, 1:nona[i]])
  } #i
  fit <- sum(sum_act)
  fit_new <- sum(sum_sim)
  c_hat <- fit_new / fit
  bpv <- step(fit_new - fit)
} #end model
", file = "CH_4.3_model3.txt")

#<><><><><><><><><><><><><><><><><><>#
#<> Load Data And Run Jags <>#
#<><><><><><><><><><><><><><><><><><>#
#<> Allow for standardizing density for both data sets: cod end and LTMR <>#
density_all <- gather(Ashad_density, Rep, Density, `1`:`5`) %>%
  select(Density) %>%
  data.frame(.) %>%
  c(., as.vector(data.frame(all_fish$TotalFL))) %>%
  unlist(., use.names = F)

#<> Read In Jags Data <>#
dat <- list(count = as.matrix(ch), 
            nona = as.vector(na_ref), 
            nobs = nrow(Ashad_ch),
            year = as.numeric(as.factor(Ashad_year)),
            nYear = max(as.numeric(as.factor(Ashad_year))),
            station = as.numeric(as.factor(Ashad_station)),
            nStation = max(as.numeric(as.factor(Ashad_station))),
            region = as.numeric(as.factor(Ashad_region)),
            nRegion = max(as.numeric(as.factor(Ashad_region))),
            saltL = as.vector(MyNorm(rowMeans(select(Ashad_salt, `1`:`5`), na.rm = T))),
            sizeL = model.matrix(~ Ashad_size)[, 2], # Large fish reference
            densityL = MyNorm_integ(as.matrix(select(Ashad_density, `1`:`5`)),
                                    mean(density_all, na.rm = T), 
                                    sd(density_all, na.rm = T)),
            effortL = MyNorm(as.matrix(select(Ashad_volume, `1`:`5`))),
            turbidityL = MyNorm(as.matrix(select(Ashad_turb, `1`:`5`))),
            
            cover_trawl = all_fish$Inside,
            cover_total = all_fish$Inside + all_fish$Outside,
            nCover = nrow(all_fish),
            sizeT = model.matrix(~ all_fish$Size -1)[, 1],
            densityT = MyNorm_integ(all_fish$TotalFL,
                                    mean(density_all, na.rm = T), 
                                    sd(density_all, na.rm = T)))

#<> Inits <>#
inits <- function(){
  list(N = round(apply(Ashad_ch, 1, max, na.rm = T), digits = 0) + 2)
}

#<> Parameters to monitor <>#
param <- c("alpha0", "alpha1", "alpha2",
           "beta0", "beta1", "beta2", "beta3", "beta4", "beta5",
           "gamma0", "gamma1", "gamma2",
           "sigmaS", "sigmaY", "sigmaR", "sigma", "sigmaP",
           "c_hat", "bpv", "fit", "fit_new"
           #           "chi_actual",
           #           "chi_sim"
)

#<> Jags simulation setup <>#
nb <- 10000
ni <- 60000 + nb
nt <- 30
nc <- 3

#<> Analysis <>#
start <- Sys.time()
out_mod1 <- jagsUI(dat, inits, param, model.file = "CH_4.3_model1.txt",
                   n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
out_mod2 <- jagsUI(dat, inits, param, model.file = "CH_4.3_model2.txt",
                   n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
out_mod3 <- jagsUI(dat, inits, param, model.file = "CH_4.3_model3.txt",
                   n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
end <- Sys.time()
end - start

print(out_mod1, digits = 3)
print(out_mod2, digits = 3)
print(out_mod3, digits = 3)

#save(out_mod1, file = "CH_4.3_model1.RData")
#save(out_mod2, file = "CH_4.3_model2.RData")
#save(out_mod3, file = "CH_4.3_model3.RData")

#<><><><><><><><><><><><><><><><><><>#
#<> Process Model Results and Figures 4.8 - 4.10 <>#
#<><><><><><><><><><><><><><><><><><>#
dat_pred <- data.frame(turb = seq(from = min(dat$turbidityL, na.rm = T), 
                                  to = max(dat$turbidityL, na.rm = T), length = 40),
                       salt = seq(from = min(dat$saltL, na.rm = T), 
                                  to = max(dat$saltL, na.rm = T), length = 40),
                       effort = seq(from = min(dat$effortL, na.rm = T), 
                                    to = max(dat$effortL, na.rm = T), length = 40),
                       density = seq(from = min(dat$densityL, na.rm = T), 
                                     to = max(dat$densityL, na.rm = T), length = 40))

pred_N <- array(NA, dim = c(nrow(dat_pred), 5, 5),
                dimnames = list(paste("Cov", 1:40, sep = "_"),
                                c("Mean", "Low", "Up", "Turbidity", "Salinity"),
                                c("Model_1_Turbidity", "Model_1_Salinity", "Model_2_Salinity", 
                                  "Model_3_Salinity", "Model_3_Turbidity")))
pred_p <- array(NA, dim = c(nrow(dat_pred), 5, 3), dimnames = list(paste("Cov", 1:40, sep = "_"),
                                                                   c("Mean", "Low", "Up", "Effort", "Turbidity"),
                                                                   c("Model_2_Effort", "Model_2_Turbidity", "Model_3_Effort")))
pred_ret <- array(NA, dim = c(nrow(dat_pred), 4, 4), dimnames = list(paste("Cov", 1:40, sep = "_"),
                                                                     c("Mean", "Low", "Up", "Biomass"),
                                                                     c("Model_2_sizesmall", "Model_2_sizelarge", 
                                                                       "Model_3_sizesmall", "Model_3_sizelarge")))


load(file = "CH_4.3_model1.RData") # Basic Count Model
for (x in 1:nrow(dat_pred)){
  pred_N[x, 1, 1] <- mean(exp(out_mod1$sims.list$beta0 + out_mod1$sims.list$beta1 + out_mod1$sims.list$beta2 + 
                                out_mod1$sims.list$beta3 * mean(dat_pred$effort) + 
                                out_mod1$sims.list$beta4 * dat_pred$turb[x] + out_mod1$sims.list$beta5 * mean(dat_pred$salt)))
  pred_N[x, 1, 2] <- mean(exp(out_mod1$sims.list$beta0 + out_mod1$sims.list$beta1 + out_mod1$sims.list$beta2 + 
                                out_mod1$sims.list$beta3 * mean(dat_pred$effort) + 
                                out_mod1$sims.list$beta4 * mean(dat_pred$turb) + out_mod1$sims.list$beta5 * dat_pred$salt[x]))
  
  pred_N[x, 2, 1] <- ci05(exp(out_mod1$sims.list$beta0 + out_mod1$sims.list$beta1 + out_mod1$sims.list$beta2 + 
                                out_mod1$sims.list$beta3 * mean(dat_pred$effort) + 
                                out_mod1$sims.list$beta4 * dat_pred$turb[x] + out_mod1$sims.list$beta5 * mean(dat_pred$salt)))
  pred_N[x, 2, 2] <- ci05(exp(out_mod1$sims.list$beta0 + out_mod1$sims.list$beta1 + out_mod1$sims.list$beta2 + 
                                out_mod1$sims.list$beta3 * mean(dat_pred$effort) + 
                                out_mod1$sims.list$beta4 * mean(dat_pred$turb) + out_mod1$sims.list$beta5 * dat_pred$salt[x]))
  
  pred_N[x, 3, 1] <- ci95(exp(out_mod1$sims.list$beta0 + out_mod1$sims.list$beta1 + out_mod1$sims.list$beta2 + 
                                out_mod1$sims.list$beta3 * mean(dat_pred$effort) + 
                                out_mod1$sims.list$beta4 * dat_pred$turb[x] + out_mod1$sims.list$beta5 * mean(dat_pred$salt)))
  pred_N[x, 3, 2] <- ci95(exp(out_mod1$sims.list$beta0 + out_mod1$sims.list$beta1 + out_mod1$sims.list$beta2 + 
                                out_mod1$sims.list$beta3 * mean(dat_pred$effort) + 
                                out_mod1$sims.list$beta4 * mean(dat_pred$turb) + out_mod1$sims.list$beta5 * dat_pred$salt[x]))
}

load(file = "CH_4.3_model2.RData") # N-mixture with turb on detection
for (x in 1:nrow(dat_pred)){
  pred_N[x, 1, 3] <- mean(exp(out_mod2$sims.list$beta0 + out_mod2$sims.list$beta1 * dat_pred$salt[x]))
  pred_N[x, 2, 3] <- ci05(exp(out_mod2$sims.list$beta0 + out_mod2$sims.list$beta1 * dat_pred$salt[x]))
  pred_N[x, 3, 3] <- ci95(exp(out_mod2$sims.list$beta0 + out_mod2$sims.list$beta1 * dat_pred$salt[x]))
  
  pred_p[x, 1, 1] <- mean(plogis(out_mod2$sims.list$alpha0 + out_mod2$sims.list$alpha1 * dat_pred$effort[x] + 
                                   out_mod2$sims.list$alpha2 * mean(dat_pred$turb)))
  pred_p[x, 1, 2] <- mean(plogis(out_mod2$sims.list$alpha0 + out_mod2$sims.list$alpha1 * mean(dat_pred$effort) + 
                                   out_mod2$sims.list$alpha2 * dat_pred$turb[x]))
  pred_p[x, 2, 1] <- ci05(plogis(out_mod2$sims.list$alpha0 + out_mod2$sims.list$alpha1 * dat_pred$effort[x] + 
                                   out_mod2$sims.list$alpha2 * mean(dat_pred$turb)))
  pred_p[x, 2, 2] <- ci05(plogis(out_mod2$sims.list$alpha0 + out_mod2$sims.list$alpha1 * mean(dat_pred$effort) + 
                                   out_mod2$sims.list$alpha2 * dat_pred$turb[x]))
  pred_p[x, 3, 1] <- ci95(plogis(out_mod2$sims.list$alpha0 + out_mod2$sims.list$alpha1 * dat_pred$effort[x] + 
                                   out_mod2$sims.list$alpha2 * mean(dat_pred$turb)))
  pred_p[x, 3, 2] <- ci95(plogis(out_mod2$sims.list$alpha0 + out_mod2$sims.list$alpha1 * mean(dat_pred$effort) + 
                                   out_mod2$sims.list$alpha2 * dat_pred$turb[x]))
  
  pred_ret[x, 1, 1] <- mean(plogis(out_mod2$sims.list$gamma0 + out_mod2$sims.list$gamma2 * dat_pred$density[x]))  
  pred_ret[x, 2, 1] <- ci05(plogis(out_mod2$sims.list$gamma0 + out_mod2$sims.list$gamma2 * dat_pred$density[x]))  
  pred_ret[x, 3, 1] <- ci95(plogis(out_mod2$sims.list$gamma0 + out_mod2$sims.list$gamma2 * dat_pred$density[x]))  
  pred_ret[x, 1, 2] <- mean(plogis(out_mod2$sims.list$gamma0 + out_mod2$sims.list$gamma1 +
                                     out_mod2$sims.list$gamma2 * dat_pred$density[x]))  
  pred_ret[x, 2, 2] <- ci05(plogis(out_mod2$sims.list$gamma0 + out_mod2$sims.list$gamma1 +
                                     out_mod2$sims.list$gamma2 * dat_pred$density[x]))  
  pred_ret[x, 3, 2] <- ci95(plogis(out_mod2$sims.list$gamma0 + out_mod2$sims.list$gamma1 +
                                     out_mod2$sims.list$gamma2 * dat_pred$density[x])) 
}

load(file = "CH_4.3_model3.RData") # N-mixture with turb on abundance
for (x in 1:nrow(dat_pred)){
  pred_N[x, 1, 4] <- mean(exp(out_mod3$sims.list$beta0 + out_mod3$sims.list$beta1 * dat_pred$salt[x] +
                                out_mod3$sims.list$beta2 * mean(dat_pred$turb)))
  pred_N[x, 1, 5] <- mean(exp(out_mod3$sims.list$beta0 + out_mod3$sims.list$beta1 * mean(dat_pred$salt) +
                                out_mod3$sims.list$beta2 * dat_pred$turb[x]))
  
  pred_N[x, 2, 4] <- ci05(exp(out_mod3$sims.list$beta0 + out_mod3$sims.list$beta1 * dat_pred$salt[x] +
                                out_mod3$sims.list$beta2 * mean(dat_pred$turb)))
  pred_N[x, 2, 5] <- ci05(exp(out_mod3$sims.list$beta0 + out_mod3$sims.list$beta1 * mean(dat_pred$salt) +
                                out_mod3$sims.list$beta2 * dat_pred$turb[x]))
  
  pred_N[x, 3, 4] <- ci95(exp(out_mod3$sims.list$beta0 + out_mod3$sims.list$beta1 * dat_pred$salt[x] +
                                out_mod3$sims.list$beta2 * mean(dat_pred$turb)))
  pred_N[x, 3, 5] <- ci95(exp(out_mod3$sims.list$beta0 + out_mod3$sims.list$beta1 * mean(dat_pred$salt) +
                                out_mod3$sims.list$beta2 * dat_pred$turb[x]))
  
  pred_p[x, 1, 3] <- mean(plogis(out_mod3$sims.list$alpha0 + out_mod3$sims.list$alpha1 * dat_pred$effort[x]))
  pred_p[x, 2, 3] <- ci05(plogis(out_mod3$sims.list$alpha0 + out_mod3$sims.list$alpha1 * dat_pred$effort[x]))
  pred_p[x, 3, 3] <- ci95(plogis(out_mod3$sims.list$alpha0 + out_mod3$sims.list$alpha1 * dat_pred$effort[x]))
  
  pred_ret[x, 1, 3] <- mean(plogis(out_mod3$sims.list$gamma0 + out_mod3$sims.list$gamma2 * dat_pred$density[x]))  
  pred_ret[x, 2, 3] <- ci05(plogis(out_mod3$sims.list$gamma0 + out_mod3$sims.list$gamma2 * dat_pred$density[x]))  
  pred_ret[x, 3, 3] <- ci95(plogis(out_mod3$sims.list$gamma0 + out_mod3$sims.list$gamma2 * dat_pred$density[x]))  
  pred_ret[x, 1, 4] <- mean(plogis(out_mod3$sims.list$gamma0 + out_mod3$sims.list$gamma1 +
                                     out_mod3$sims.list$gamma2 * dat_pred$density[x]))  
  pred_ret[x, 2, 4] <- ci05(plogis(out_mod3$sims.list$gamma0 + out_mod3$sims.list$gamma1 +
                                     out_mod3$sims.list$gamma2 * dat_pred$density[x]))  
  pred_ret[x, 3, 4] <- ci95(plogis(out_mod3$sims.list$gamma0 + out_mod3$sims.list$gamma1 +
                                     out_mod3$sims.list$gamma2 * dat_pred$density[x]))  
}

pred_N[, 4, ] <-  pred_p[, 5, ] <- dat_pred$turb
pred_N[, 5, ] <- dat_pred$salt
pred_p[, 4, ] <- dat_pred$effort
pred_ret[, 4, ] <- dat_pred$density


#<> Abundance Plot <>#
pred_N <- as.data.frame.table(pred_N) %>%
  spread(., Var2, Freq) %>%
  separate(Var3, c("Ref", "Model", "Covariate")) %>%
  mutate(Model = recode(Model, "1" = "Model 1", "2" = "Model 2", "3" = "Model 3"))

pred.N.turb <- ggplot(data = filter(pred_N, Covariate == "Turbidity"), aes(x = Turbidity, y = Mean)) + 
  geom_ribbon(aes(ymin = Low, ymax = Up, fill = as.factor(Model)), alpha = 0.4) + 
  geom_line(aes(y = Mean, colour = Model), lwd = 2) +
  guides(colour = FALSE) +
  labs(x = "Standardized turbidity", y = "Abundance") +
  scale_fill_manual(values = c("#CC79A7", "#009E73")) +   #c("blue", "darkturquoise")) +
  scale_colour_manual(values = c("#CC79A7", "#009E73")) + #c("blue", "darkturquoise")) +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(vjust = -2.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.8), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        text = element_text(size = 16)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 

pred.N.salt <- ggplot(data = filter(pred_N, Covariate == "Salinity"), aes(x = Salinity, y = Mean)) + 
  geom_ribbon(aes(ymin = Low, ymax = Up, fill = as.factor(Model)), alpha = 0.4) + 
  geom_line(aes(y = Mean, colour = Model), lwd = 2) +
  guides(colour = FALSE) +
  scale_fill_manual(values = c("#CC79A7", "#0072B2", "#009E73")) + # c("blue", "red", "darkturquoise")) +
  scale_colour_manual(values = c("#CC79A7", "#0072B2", "#009E73")) + #c("blue", "red", "darkturquoise")) +
  labs(x = "Standardized salinity", y = "Abundance") +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(vjust = -2.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.8), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        text = element_text(size = 16)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 


#<> Availability Plot <>#
pred_p <- as.data.frame.table(pred_p) %>%
  spread(., Var2, Freq) %>%
  separate(Var3, c("Ref", "Model", "Covariate")) %>%
  mutate(Model = recode(Model, "2" = "Model 2", "3" = "Model 3"))

p.avail.turb <- ggplot(data = filter(pred_p, Covariate == "Turbidity"), aes(x = Turbidity, y = Mean)) + 
  geom_ribbon(aes(ymin = Low, ymax = Up, fill = as.factor(Model)), alpha = 0.4) + 
  geom_line(aes(y = Mean, colour = Model), lwd = 2) +
  guides(colour = FALSE) +
  scale_fill_manual(values = c("#0072B2")) +
  scale_colour_manual(values = c("#0072B2")) +
  labs(x = "Standardized turbidity", y = expression(q[Availability])) +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(vjust = -2.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.8), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        text = element_text(size = 16)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 

p.avail.vol <- ggplot(data = filter(pred_p, Covariate == "Effort"), aes(x = Effort, y = Mean)) + 
  geom_ribbon(aes(ymin = Low, ymax = Up, fill = as.factor(Model)), alpha = 0.4) + 
  geom_line(aes(y = Mean, colour = Model), lwd = 2) +
  guides(colour = FALSE) +
  scale_fill_manual(values = c("#0072B2", "#009E73")) +
  scale_colour_manual(values = c("#0072B2", "#009E73")) +
  labs(x = "Standardized volume sampled", y = expression(q[Availability])) +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(vjust = -2.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.8), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        text = element_text(size = 16)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 

grid.arrange(pred.N.turb, pred.N.salt, # Figure 4.10
             ncol = 1, nrow = 2)

grid.arrange(p.avail.turb, p.avail.vol, # Figure 4.9
             ncol = 1, nrow = 2)

#<> Retention Plot <>#
pred_ret <- as.data.frame.table(pred_ret) %>%
  spread(., Var2, Freq) %>%
  separate(Var3, c("Ref", "Model", "Covariate")) %>%
  mutate(Model = recode(Model, "2" = "Model 2", "3" = "Model 3"),
         Covariate = recode(Covariate, "sizesmall" = "Small fish", 
                            "sizelarge" = "Large fish")) %>%
  transform(Covariate = factor(Covariate, levels = c("Small fish", "Large fish")))

p.ret.small <- ggplot(data = filter(pred_ret, Covariate == "Small fish"), aes(x = Biomass, y = Mean)) + 
  geom_ribbon(aes(ymin = Low, ymax = Up, fill = as.factor(Model)), alpha = 0.4) + 
  geom_line(aes(y = Mean, colour = Model), lwd = 2) +
  guides(colour = FALSE) +
  ylim(0, 1) +
  scale_fill_manual(values = c("#0072B2", "#009E73")) +
  scale_colour_manual(values = c("#0072B2", "#009E73")) +
  annotate("text", x = 13, y = 0.05, label = "A)", size = 8) +
  labs(x = "Standardized cumulative fish length (mm)", y = expression(q[Retained])) +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(vjust = -2.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.8), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        text = element_text(size = 16)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 

p.ret.large <- ggplot(data = filter(pred_ret, Covariate == "Large fish"), aes(x = Biomass, y = Mean)) + 
  geom_ribbon(aes(ymin = Low, ymax = Up, fill = as.factor(Model)), alpha = 0.4) + 
  geom_line(aes(y = Mean, colour = Model), lwd = 2) +
  guides(colour = FALSE) +
  ylim(0, 1) +
  scale_fill_manual(values = c("#0072B2", "#009E73")) +
  scale_colour_manual(values = c("#0072B2", "#009E73")) +
  annotate("text", x = 13, y = 0.05, label = "B)", size = 8) +
  labs(x = "Standardized cumulative fish length (mm)", y = expression(q[Retained])) +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(vjust = -2.5),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "none", 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        text = element_text(size = 16)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 

grid.arrange(p.ret.small, p.ret.large, # Figure 4.8
             ncol = 1, nrow = 2)


#<><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
#<> Workflow 4.4: Application of the Integrated Abundance Model to the FMWT <>#
# Analysis is based on Walker et al. (2017) and Moriarty et al. (2020). 
# Moriarty, M., Sethi, S.A., Pedreschi, D., Smeltz, T.S., McGonigle, C., Harris, B.P., Wolf, N., and 
#   Greenstreet, S.P.R. 2020. Combining fisheries surveys to inform marine species distribution modelling. 
#   ICES J. Mar. Sci. 77(2): 539-552. doi:10.1093/icesjms/fsz254.
# Walker, N.D., Maxwell, D.L., Le Quesne, W.J.F., and Jennings, S. 2017. Estimating efficiency of survey 
#   and commercial trawl gears from comparisons of catch-ratios. ICES J. Mar. Sci. 
#   74(5): 1448-1457. doi:10.1093/icesjms/fsw250.
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
#<><><><><><><><><><><><><><><><><><>#
#<><><><><> Manipulate Data <><><><><>#
#<><><><><><><><><><><><><><><><><><>#
# 100 = South bay, 200 = central bay, 300's = San Pablo, 400, 600 = Suisun, 500 = Honker, 
# 600 = Suisun Marsh, 700 = Lower Sac, 800 = Lower San Joaq, 900 = Upper San Joaq, 
# Suisun tow: 1.5 x 4.3m open, length = 5.3m, 35mm mesh stretch in body and 6mm stretch at code-end. Towed at 4km/hour
# Fish: SL measured
# Bay midwater: 9 section mesh from 20.3cm stretch at mouth to 1.3cm at codend. 5:1 scope
# Fish: FL besides some species with TL (without forked tails)
# Bay otter: 4.9m headrope, 2.5cm stretch mesh body, 1.3 cm stretch mesh codend, 5:1 scope
# Fish: FL besides some species with TL (without forked tails)
# FMWT: 13.4m2 mouth, 8" stretch mesh to 0.5" stretch in codend
# Fish: FL
#<><> Raw Data Loading <><>#
require(LTMRdata)

ltmr_gam <- LTMRdata::LTMRpilot(quiet = FALSE,
                                convert_length = TRUE,
                                remove_unconverted_lengths = FALSE,
                                size_cutoff = NULL,
                                measured_lengths = FALSE) %>%
  filter(!is.na(Length), # remove any samples where length was not measured
         !is.na(Latitude), # remove sites without coordinates
         !is.na(Longitude)) %>% # remove sites without coordinates
  mutate(Taxa = gsub(" ", "_", Taxa), # replace space with _
         year = year(Date),
         Month = month(Date),
         Season = ifelse(Month > 2 & Month < 6, "Spring",
                         ifelse(Month > 5 & Month < 9, "Summer",
                                ifelse(Month > 8 & Month < 12, "Fall", "Winter"))),
         Year = ifelse(Month > 11, year + 1, year), # california water year runs from Oct1 to Sept 30. However only use December because it is the start of the next season (Sept is still fall, which would force fall to be split between water years)
         Size = ifelse(Length < 50.1, "Small", "Large"),
         TaxaFull = paste(Taxa, Size, sep = "_"), # put taxa and size class into one name
         Gear = as.factor(paste(Source, Method, sep = "_")),
         Tow_volume = ifelse(is.na(Tow_volume), 1.5 * Tow_area, # tow area for Suisun is width * dist * 70% open * minutes towed. I added 1.5 to account for trawl height
                             ifelse(Tow_volume < 1, 1.5 * Tow_area, Tow_volume)), # tow volume for some data was recorded as 0 instead of NA
         Region = as.numeric(Station),
         RegionR = ifelse(Region < 200 & Region > 99, "South Bay",
                          ifelse(Region > 199 & Region < 300, "Central Bay",
                                 ifelse(Region > 299 & Region < 400, "San Pablo",
                                        ifelse(Region < 600 & Region > 499, "Honker",
                                               ifelse(Region < 712 & Region > 699, "Lower Sacramento",
                                                      ifelse(Region < 800 & Region > 711, "Upper Sacramento",
                                                             ifelse(Region < 900 & Region > 799, "Lower San Joaquin",
                                                                    ifelse(Region > 899, "Upper San Joaquin", 
                                                                           ifelse(Region < 100, "Upper Sacramento", "Suisun")))))))))) %>% 
  filter(!is.na(Tow_volume), # remove samples without tow volume information
         Year < 2019) %>% # remove 2019 and after year because full sample has yet to be collected
  transform(RegionR = ifelse(is.na(RegionR), "Suisun", 
                             ifelse(Region == 433, "Honker",
                                    ifelse(Gear == "Bay Study_Midwater trawl" & Region > 700 & Region < 800, "Lower Sacramento",
                                           ifelse(Gear == "Bay Study_Otter trawl" & Region > 700 & Region < 800, "Lower Sacramento", RegionR))))) %>%
  select(Source:Date, Year, Month:Season, Region = RegionR, Gear, Tow_volume, SampleID:TaxaFull, Taxa:Count) %>%
  group_by(Gear, Region, Station, Latitude, Longitude, Date, Year, Season, Month, Tow_volume, TaxaFull) %>%
  summarize(Count = sum(Count)) %>%
  ungroup() %>%
  spread(., TaxaFull, Count) %>% # fill in NA's for missing taxa
  gather(., TaxaFull, Count, Acanthogobius_flavimanus_Large:UnID_Small) %>% # Gather taxa to replace NA with 0
  mutate(Count = replace(Count, which(is.na(Count)), 0)) %>%
  spread(., TaxaFull, Count) %>%
  mutate(Season = as.numeric(as.factor(Season)),
         Yr_Season = as.numeric(as.factor(paste(Year, Season, sep = ("_")))),
         Day = as.numeric(as.factor(julian(Date))),
         JDay = yday(Date),
         FromDay = as.numeric((Date - min(Date))) / (60 * 60 * 24)) %>% # convert seconds to days
  select(Gear:Month, Yr_Season:FromDay, Tow_volume:UnID_Small) %>%
  data.frame(.) %>%
  gather(., Taxa, Counts, Acanthogobius_flavimanus_Large:UnID_Small) %>%
  transform(Gear = factor(Gear, levels = c("FMWT_Midwater trawl",  # make FMWT reference because it has covered cod end experiment
                                           "Bay Study_Midwater trawl", 
                                           "Bay Study_Otter trawl",
                                           "Suisun_Otter trawl")))

#<><><><><><><><><><><><><><><><><><>#
#<> GAM's Analysis <>#
#<><><><><><><><><><><><><><><><><><>#
taxa_ref <- unique(ltmr_gam$Taxa) # to reference all taxa by size breakdowns
coef_fix <- list() # list to store of catch ratio data

for (x in taxa_ref){
  print(x)
  dat <- filter(ltmr_gam, Taxa == x)
  
  pres_gear <- names(which(table(dat[which(dat$Counts > 0), 'Gear']) > 0)) # ID only gears with species present
  gear_ref <- pres_gear[1] # reference gear
  pres <- length(which(table(dat[which(dat$Counts > 0), 'Gear']) > 0)) # how many gear-surveys had a species collected 
  dat <- filter(dat, Gear %in% pres_gear) # only fit to gear in which the fish was captured at least once
  
  prop <- length(dat[which(dat$Counts > 0), 'Counts']) / nrow(dat) # Proportion of samples with fish present at least once
  
  # Only Fit the Model if there are greater than 1 Gear combination with fish present and the species is in at least 1% of samples #
  if (pres < 2 | prop < 0.01){
  } # End of if: Model is not fit
  else {
    mod <- bam(as.integer(round(Counts, digits = 0)) ~ 
                 te(Latitude, Longitude, FromDay, bs = "cs") + 
                 Gear +
                 offset(log(Tow_volume)),
               method = "fREML",
               data = dat,
               family = nb,
               discrete = TRUE)
    
    outfile = paste("C:/Users/bhuntsman/Documents/IEP Survey Review/catchability github/CH_4.4_results/",
                    x, ".Rdata", sep = "") # Name of saved file
    save(mod, file = outfile)	 # Save the results
    
    coef_fix[[x]] <- data.frame(Taxa = x,
                                Gear = pres_gear[2:length(pres_gear)],
                                GearReference = gear_ref,
                                Average = exp(mod$coefficients[2:length(pres_gear)]),
                                Lower = (exp(mod$coefficients[2:length(pres_gear)] - 
                                               (2 * summary(mod)$se[2:length(pres_gear)]))),
                                Upper = (exp(mod$coefficients[2:length(pres_gear)] + 
                                               (2 * summary(mod)$se[2:length(pres_gear)]))))
    
  } # End of Else
} #x

coef_fix <- do.call(rbind, coef_fix)
# save(coef_fix, file = "coef.RData")

#<><><><><><><><><><><><><><><><><><>#
#<> GAM's Plots <>#
#<><><><><><><><><><><><><><><><><><>#
### Figure 4.13: Catch Ratios ###
load(file = "coef.RData")

# For Plotting Simplicity #
cluster_ref <- data.frame(Cluster = c(rep("Open-water Year Around", 6), rep("Open-water Fall", 4), 
                                      rep("Demersal Year Around", 7)),
                          Common = c("Staghorn Scuplin", "Delta Smelt", "Shiner Perch", "Pacific Sardine",
                                     "American Shad", "Splittail", "Topsmelt", "Sacramento Sucker",
                                     "Starry Flounder", "Striped Bass", "California Halibut", "Tule Perch", "Bat Ray",
                                     "Pacific Hering", "Pacific Sanddab", "Mississippi Silverside", "White Sturgeon"),
                          Taxa = c("Leptocottus armatus", "Hypomesus transpacificus", "Cymatogaster aggregata", 'Sardinops sagax',
                                   "Alosa sapidissima", "Pogonichthys macrolepidotus", "Atherinops affinis",
                                   "Catostomus occidentalis", "Platichthys stellatus", "Morone saxatillis", "Paralichthys californicus",
                                   "Hysterocarpus traskii", "Myliobatis californica", "Clupea pallasii", "Citharichthys sordidus",
                                   "Menidia audens", "Acipenser transmontanus"))

sp_code <- read.csv("SpeciesInfo.csv", header = T) %>%
  select(Code = SpeciesCode3, CommonName, Genus, Species, Habitat) %>%
  filter(Code != "CSSR", Code != "CSFR", Code != "CSWR", Code != "CSLFR")# remove the extra chinook salmon species codes

coef_fix <- coef_fix %>%
  separate(Taxa, c("Genus", "Species", "Size"), sep = "_") %>%
  merge(., sp_code, all.x = T) %>%
  mutate(Taxa = paste(Genus, Species, sep = " "),
         Gear = recode(Gear, "Bay Study_Otter trawl" = "Bay Study Otter Trawl",
                       "Suisun_Otter trawl" = "Suisun Otter Trawl",
                       "Bay Study_Midwater trawl" = "Bay Study Midwater Trawl",
                       "FMWT_Midwater trawl" = "FMWT Midwater Trawl")) %>%
  merge(., cluster_ref, all.x = T) %>%
  transform(Gear = factor(Gear, levels = c("FMWT Midwater Trawl",
                                           "Bay Study Midwater Trawl",
                                           "Bay Study Otter Trawl",
                                           "Suisun Otter Trawl")),
            Size = factor(Size, levels = c("Small", "Large")))

coef_values <- filter(coef_fix, !is.na(Average))

nrow(distinct(select(coef_values, Taxa, Size)))
nrow(distinct(select(coef_values, Taxa)))
table(distinct(select(coef_values, Taxa, Size))$Size)

species_list <- distinct(select(coef_values, Taxa, CommonName, Code, Habitat))

## All Species Catch-Ratio FMWT reference ##
figure_4.13 <- ggplot(data = filter(coef_values, GearReference == "FMWT_Midwater trawl"),
       aes(x = Code, y = Average, colour = Habitat)) +
  geom_linerange(aes(xmin = Code, ymin = Lower, xmax = Code, ymax = Upper, colour = Habitat), 
                 lwd = 1.75) +
  geom_point(size = 5) +
  scale_colour_manual(values = c("#CC79A7", "#0072B2", "#009E73")) +
  ylim(0, 10) +
  labs(x = "", y = "Catch-Ratio") +
  geom_hline(yintercept = 1, lty = 2) +
  facet_grid(Gear ~ Size, scale = "free") +
  theme( legend.position = c(0.39, 0.25),
         legend.title = element_blank(),
         legend.key = element_blank(),
         legend.box.background = element_rect(colour = "black", size = 1),
         axis.text = element_text(size = 12),
         axis.title.x = element_text(vjust = -2.5),
         axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.7),
         panel.background = element_blank(),
         panel.border = element_rect(colour = "black", fill = NA, size = 2),
         text = element_text(size = 16),
         plot.title = element_text(vjust = -18, hjust = 0.95))

