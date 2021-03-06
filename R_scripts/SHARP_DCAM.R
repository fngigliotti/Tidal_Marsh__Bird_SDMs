###############################################################
###Dynamic Community Abundance Model - SHARP 10 Year PC Data###
###############################################################

#################
###Description###
#################

#This model analyzes 10 years of point-count data collected by SHARP and examines
#variation in rates of population change over time for all members of the 
#tidal marsh breeding bird community.Explanatory variables (covariates) included 
#in link-functions relate to variability in habitat suitability over space and time,
#and site and survey-specific likelihoods of detection, which influences the 
#observation process. Much of this model code is based off of code presented in 
#Applied Hierarchical Modeling Volumes 1 (2016) and 2 (2020) by Marc Kery and Andy Royle, 
#as well as the supporting R package 'AHMbook'. Future models will expand upon this preliminary
#model to incorporate spatial autocorrelation in data, integrate distance-sampling, CMR, and/or 
#time-to-detction data, and incorporate additional covariates such as vegetative community
#information, climatic variables, or marsh hydrology. CMR data will especially be useful as it
#can help inform the distinction between emigration and death in the model, particularly when
#considering a spatial structure to the data (spatial autocorrelation extension).

#The goal of this code is to create a dynamic community abundance model.
#Breaking down each adjective:
###Dynamic: Explicitly predicting survival and recruitment of each species at sites between years. 
####More importantly, calculating trends in abundance over time of each species.
###Community: Using a community abundance modeling approach. Using distribution information of
####common species to inform distribution information of rare species. Needs to be for species that 
####generally have fairly similar niches (i.e. Tidal Marsh Obligates).
###Abundance: Needs to be abundances or densities of individuals. Important component is trends at
####sites over time, rather than percentage of occupied patches over time.

#######################################
###Load Required Packages and Set WD###
#######################################

library(AHMbook)
library(jagsUI)
library(ggplot2)  
setwd(choose.dir()) #because lazy 
getwd() #where it's at
source("homebrewed_simulation_functions.R") #Load simulation function (requires AHMbook)
#simDCAM() function can be loaded for a 'simple' dynamic CAM with covariates

########################################
###Simulate Dataset for Model Testing###
########################################

#Simple, farily species-poor community for model testing. Observation-level parameters
#specified for model testing. Default parameters employed for covariate relationships. 

sim.data<-simDCAM(nspecies = 20, nsites = 30, nsurveys = 3, nyears = 10, 
                  mean.lambda = 3, mean.phi = 0.7, mean.gamma = 0.4, mean.p = 0.5)

jags.data<-list(nsites = sim.data$nsites, nspecies = sim.data$nspecies, nsurveys = sim.data$nsurveys,
                nyears = sim.data$nyears, cov.lambda = sim.data$habitat, cov.phi = sim.data$habitat, 
                cov.gamma = sim.data$habitat, cov.p = sim.data$wind, y = sim.data$y.obs)

################################
###Specify Model Code in JAGS###
################################

cat(file = "SHARP_DCAM.txt", "
model {
  
  ##########################
  #Priors and Hyperpriors###
  ##########################
  
  # *** Priors and hyperpriors for model on lambda (t = 1 abundance) ***
  # Priors
  for(k in 1:nspecies){ # Loop over # of species 
    alpha.lambda[k] ~ dnorm(mu.alpha.lambda, tau.alpha.lambda) # Intercepts
    beta.lambda[k] ~ dnorm(mu.beta.lambda, tau.beta.lambda) #Coeffs
  }
  
  # Hyperpriors
    #Intercept
    mu.alpha.lambda <- log(mean.lambda)
    mean.lambda ~ dunif(0, 100) #initial abundance (mean across species)
    sd.alpha.lambda ~ dnorm(0, 0.1)I(0,)
    tau.alpha.lambda <- 1/(sd.alpha.lambda*sd.alpha.lambda)
    
    #Coefficient
    mu.beta.lambda ~ dnorm(0,0.1)         
    sd.beta.lambda ~ dnorm(0, 0.1)I(0,)
    tau.beta.lambda <- 1/(sd.beta.lambda*sd.beta.lambda)
  
  # *** Priors and hyperpriors for model on phi (survival) ***
  # Priors
  for(k in 1:nspecies){ # Loop over # of species 
    alpha.phi[k] ~ dnorm(mu.alpha.phi, tau.alpha.phi) #Intercepts
    beta.phi[k] ~ dnorm(mu.beta.phi, tau.beta.phi) #Coeffs
  }
  
  # Hyperpriors
    #Intercept
    mu.alpha.phi <- logit(mean.phi) 
    mean.phi ~ dunif(0, 1) #mean yearly survival across species
    sd.alpha.phi ~ dnorm(0, 0.1)I(0,)
    tau.alpha.phi <- 1/(sd.alpha.phi*sd.alpha.phi)
    
    #Coefficient
    mu.beta.phi ~ dnorm(0, 0.1) #No yearly variation in hyperprior for coeff
    sd.beta.phi ~ dnorm(0, 0.1)I(0,)
    tau.beta.phi <- 1/(sd.beta.phi*sd.beta.phi)
  
  # *** Priors and hyperpriors for model on gamma (recruitment) ***
  # Priors
  for(k in 1:nspecies){ # Loop over # of species 
    alpha.gamma[k] ~ dnorm(mu.alpha.gamma, tau.alpha.gamma) #Intercepts
    beta.gamma[k] ~ dnorm(mu.beta.gamma, tau.beta.gamma) #Coeffs
  }
  
  # Hyperpriors
    #Intercept
    mu.alpha.gamma <- log(mean.gamma) 
    mean.gamma ~ dunif(0, 10) #mean yearly recruitment across species
    sd.alpha.gamma ~ dnorm(0, 0.1)I(0,)
    tau.alpha.gamma <- 1/(sd.alpha.gamma*sd.alpha.gamma)
    
    #Coefficient
    mu.beta.gamma ~ dnorm(0, 0.1) #No yearly variation in hyperprior for coeff
    sd.beta.gamma ~ dnorm(0, 0.1)I(0,)
    tau.beta.gamma <- 1/(sd.beta.gamma*sd.beta.gamma)
  
  # *** Priors and hyperpriors for model on p (detection) ***
  # Priors
  for(k in 1:nspecies){ # Loop over # of species 
    alpha.p[k] ~ dnorm(mu.alpha.p, tau.alpha.p) # Intercepts
    beta.p[k] ~ dnorm(mu.beta.p, tau.beta.p) #Coeffs
  }
  
  # Hyperpriors
    #Intercept
    mu.alpha.p <- logit(mean.p)
    mean.p ~ dunif(0, 1) #Mean detection across species
    sd.alpha.p ~ dnorm(0, 0.1)I(0,)
    tau.alpha.p <- 1/(sd.alpha.p*sd.alpha.p)
    
    #Coefficient
    mu.beta.p ~ dnorm(0, 0.1) #No yearly variation in hyperprior for coeff
    sd.beta.p ~ dnorm(0, 0.1)I(0,)
    tau.beta.p <- 1/(sd.beta.p*sd.beta.p)
    
  ######################
  ###Model Likelihood###
  ######################
  
  ###Ecological submodel: Define state conditional on parameters
  
  ##Abundance
  #Initial Abundance (t=1)
  for(i in 1:nsites){
    for(k in 1:nspecies){
      N[i,k,1] ~ dpois(lambda[i,k]) #initial abundance of species k at site i (year 1)
      log(lambda[i,k]) <- alpha.lambda[k] + beta.lambda[k] * cov.lambda[i]
      #Link function for site/spp.-cov influence on abundance of spp. k at site i
      for(t in 1:(nyears-1)){
        N[i,k,t+1] <- S[i,k,t+1] + R[i,k,t+1] #t+1 abundance of species k at site i
      }#t
    }#k
  }#i
  
  ##Survival (constant across years)
  for(i in 1:nsites){
    for(k in 1:nspecies){
      for(t in 1:(nyears-1)){
        S[i,k,t+1] ~ dbin(phi[i,k], N[i,k,t]) #survival of spp. k at site i (constant across years)
      }#t
      logit(phi[i,k]) <- alpha.phi[k] + beta.phi[k] * cov.phi[i]
      #Link function for site/spp.-cov influence on survival of spp. k at site i
    }#k
  }#i
  
  ##Recruitment (constant across years)
  for(i in 1:nsites){
    for(k in 1:nspecies){
      for(t in 1:(nyears-1)){
        tmp[i,k,t] <- N[i,k,t]*gamma[i,k]
        R[i,k,t+1] ~ dpois(tmp[i,k,t]) #relative recruitment of spp. k at site i (constant across years)
      }#t
      log(gamma[i,k]) <- alpha.gamma[k] + beta.gamma[k] * cov.gamma[i]
      #Link function for site/spp.-cov influence on relative recruitment of spp. k at site i
    }#k
  }#i
 
 
  ###Observation model (incl. covariate effects)
  
  ##Detection
  for (i in 1:nsites){
    for(k in 1:nspecies){
      for (j in 1:nsurveys){
        for (t in 1:nyears){
          # Observed dataset as a function of abundance and detection
          y[i,k,j,t] ~ dbin(p[i,k,j,t], N[i,k,t]) 
          #Detection of sp. k at site i during replicate j and year t 
          logit(p[i,k,j,t]) <- alpha.p[k] + beta.p[k] * cov.p[i,j,t] 
          #Link function for site/spp./rep-cov influence on detection of spp. k at site i during rep. j and year t
        }#t
      }#j
    }#k
  }#i
  
  ########################
  ###Derived Parameters###
  ########################
  
  #???
  
  
}#END OF MODEL CODE
")

######################
###Specify Initials###
######################

Rst <- apply(jags.data$y, c(1,2,4), max) + 10
Rst[,,1] <- NA #No recruitment in year 1
Nst <- apply(jags.data$y, c(1,2,4), max) + 2
Nst[,,2:dim(jags.data$y)[4]] <- NA

inits <- function() list(R=Rst, N=Nst)
                         
#########################
####Parameters to Save###
#########################

#List of saved parameters informs accuracy of base model (not covariate effects)
#Can modify if desired

params<-c("mean.lambda", "mean.phi", "mean.gamma", "mean.p", "mu.alpha.lambda",
          "sd.alpha.lambda", "sd.alpha.phi", "sd.alpha.gamma", "sd.alpha.p", "N")

###################
###Running Model###
###################

#Approximate run time: 220 min
out<-jags(jags.data, inits, params, "SHARP_DCAM.txt",
          n.chains = 3, n.adapt=3000, n.burnin = 20000, 
          n.iter =50000, n.thin =3, parallel=TRUE)

##################################
###Evaluating Model Performance###
##################################

#Percent Error of posteriors

#Community Hyperparameters
((out$mean$mean.lambda-sim.data$mean.lambda)/sim.data$mean.lambda)*100 #Abundance accuracy
((out$mean$sd.alpha.lambda-sim.data$sd.lambda)/sim.data$sd.lambda)*100 #Abundance variance accuracy
((out$mean$mean.phi-sim.data$mean.phi)/sim.data$mean.phi)*100 #Survival accuracy
((out$mean$sd.alpha.phi-sim.data$sd.phi)/sim.data$sd.phi)*100 #Survival variance accuracy
((out$mean$mean.gamma-sim.data$mean.gamma)/sim.data$mean.gamma)*100 #Recruitment accuracy
((out$mean$sd.alpha.gamma-sim.data$sd.gamma)/sim.data$sd.gamma)*100 #Recruitment variance accuracy
((out$mean$mean.p-sim.data$mean.p)/sim.data$mean.p)*100 #Detection accuracy
((out$mean$sd.alpha.p-sim.data$sd.p)/sim.data$sd.p)*100 #Detection variance accuracy


#Species Abundances 
#Abundances of each species across sampled area each year (spp. x year)
Ntot.t<-apply(sim.data$N,c(2,3),sum) #Latent (true) abundances
Ntot.e<-apply(out$mean$N,c(2,3),sum) #Estimated abundances
((Ntot.e-Ntot.t)/Ntot.t)*100 #Abundance accuracy

#Calculating Growth Rates (Continuous)
r.t<-rep(NA, jags.data$nspecies)
r.e<-rep(NA, jags.data$nspecies)
for(k in 1:jags.data$nspecies){
  r.t[k]<-(log(Ntot.t[k,jags.data$nyears]/Ntot.t[k,1])/(jags.data$nyears-1))*100
  r.e[k]<-(log(Ntot.e[k,jags.data$nyears]/Ntot.e[k,1])/(jags.data$nyears-1))*100
}
((r.e-r.t)/r.t)*100 #Growth rate accuracy

###Simple plots to observe results###

#Abundance Accuracy
ggplot(data = data.frame(x = c(0, 8)), aes(x)) + #change limits as needed
  stat_function(fun = dnorm, n = 101, args = list(mean = sim.data$mean.lambda, sd = sim.data$sd.lambda), aes(color = "black")) +
  stat_function(fun = dnorm, n = 101, args = list(mean = out$mean$mean.lambda, sd = out$mean$sd.alpha.lambda), aes(color = "orange")) +
  labs(x = "Abundance Hyperprior Distribution", y = "Density") +
  theme_classic() +
  scale_color_manual(name = NULL, 
                     values =c('black'='black','orange'='orange'), labels = c('True','Estimated'), guide = 'legend')

#Survival Accuracy
ggplot(data = data.frame(x = c(-1, 2)), aes(x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = sim.data$mean.phi, sd = sim.data$sd.phi), aes(color = "black")) +
  stat_function(fun = dnorm, n = 101, args = list(mean = out$mean$mean.phi, sd = out$mean$sd.alpha.phi), aes(color = "orange")) +
  labs(x = "Survival Hyperprior Distribution", y = "Density") +
  theme_classic() +
  scale_color_manual(name = NULL, 
                     values =c('black'='black','orange'='orange'), labels = c('True','Estimated'), guide = 'legend')

#Recruitment Accuracy
ggplot(data = data.frame(x = c(-1, 2)), aes(x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = sim.data$mean.gamma, sd = sim.data$sd.gamma), aes(color = "black")) +
  stat_function(fun = dnorm, n = 101, args = list(mean = out$mean$mean.gamma, sd = out$mean$sd.alpha.gamma), aes(color = "orange")) +
  labs(x = "Recruitment Hyperprior Distribution", y = "Density") +
  theme_classic() +
  scale_color_manual(name = NULL, 
                     values =c('black'='black','orange'='orange'), labels = c('True','Estimated'), guide = 'legend')

#Detection Accuracy
ggplot(data = data.frame(x = c(-1, 2)), aes(x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = sim.data$mean.p, sd = sim.data$sd.p), aes(color = "black")) +
  stat_function(fun = dnorm, n = 101, args = list(mean = out$mean$mean.p, sd = out$mean$sd.alpha.p), aes(color = "orange")) +
  labs(x = "Detection Hyperprior Distribution", y = "Density") +
  theme_classic() +
  scale_color_manual(name = NULL, 
                     values =c('black'='black','orange'='orange'), labels = c('True','Estimated'), guide = 'legend')
#END