#######################################################################
###Spatial Dynamic Community Abundance Model - SHARP 10 Year PC Data###
#######################################################################

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
#model to integrate distance-sampling, CMR, and/or time-to-detction data,and incorporate additional 
#covariates such as vegetative community information, climatic variables, or marsh hydrology.
#CMR data will especially be useful as it can help inform the distinction between 
#emigration and death in the model, particularly when considering a spatial
#structure to the data (spatial autocorrelation extension).

#The goal of this code is to create a spatial dynamic community abundance model.
#Breaking down each adjective:
###Spatial: Explicitly accounting for spatial autocorrelation present in the dataset. This allows
####us to explicitly model immigration and emigration between patches (in a BIDE framework)
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
setwd(choose.dir()) #because lazy 
getwd() #where it's at
source("homebrewed_simulation_functions.R") #Load simulation function (requires AHMbook)
#simSpDCAM() function can be loaded for a 'simple' spatial dynamic CAM with covariates

########################################
###Simulate Dataset for Model Testing###
########################################

#Simple, farily species-poor community for model testing. Observation-level parameters
#specified for model testing. Default parameters employed for covariate relationships. 
#nsites must be an integer^2

sim.data<-simSpDCAM(nspecies = 10, nsites = 25, nsurveys = 3, nyears = 3, mean.lambda = 3, 
                    mean.phi = 0.7, mean.gamma = 0.3, mean.kappa = 0.4, 
                    mean.theta = 1.5, mean.p = 0.5)

jags.data<-list(nsites = sim.data$nsites, nspecies = sim.data$nspecies, nsurveys = sim.data$nsurveys,
                nyears = sim.data$nyears, cov.lambda = sim.data$habitat, cov.phi = sim.data$habitat, 
                cov.gamma = sim.data$habitat, cov.kappa = sim.data$habitat, cov.p = sim.data$wind,
                y = sim.data$y.obs, dist = sim.data$dist)

################################
###Specify Model Code in JAGS###
################################

#Two formulations: First, immigration/emigration can only occur between adjacent cells.
#This greatly decreases computing power but might reduce realism of the model if some
#species in focal community exhibit long distance juvenile or post-breeding dispersal
#or low breeding site fidelity across years. 
#Second formulation allows immigration/emigration to occur between any two cells in the
#spatial matrix. This may increase realism of model for long-distance dispersers, but
#greatly increases computing power requirements because immigration/emigration likelihoods
#need to be calculated for every cell pair during every model iteration. If the 
#focal community, number of cells, and/or number of years is large, the model is 
#unlikely to converge in any 'reasonable' amount of time.

#Nearest Neighbor (Formulation 1)
cat(file = "SHARP_Spatial_NN_DCAM.txt", " 
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
  
  # *** Priors and hyperpriors for model on kappa (emigration) ***
  # Priors
  for(k in 1:nspecies){ # Loop over # of species 
    alpha.kappa[k] ~ dnorm(mu.alpha.kappa, tau.alpha.kappa) #Intercepts
    beta.kappa[k] ~ dnorm(mu.beta.kappa, tau.beta.kappa) #Coeffs
  }
  
  # Hyperpriors
    #Intercept
    mu.alpha.kappa <- logit(mean.kappa) 
    mean.kappa ~ dunif(0, 1) #mean yearly emigration rate across species (conditional on survival)
    sd.alpha.kappa ~ dnorm(0, 0.1)I(0,)
    tau.alpha.kappa <- 1/(sd.alpha.kappa*sd.alpha.kappa)
    
    #Coefficient
    mu.beta.kappa ~ dnorm(0, 0.1) #No yearly variation in hyperprior for coeff
    sd.beta.kappa ~ dnorm(0, 0.1)I(0,)
    tau.beta.kappa <- 1/(sd.beta.kappa*sd.beta.kappa)
  
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
        N[i,k,t+1] <- S[i,k,t+1] - E[i,k,t+1] + R[i,k,t+1] + I[i,k,t+1] #t+1 abundance of species k at site i
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
 
  ##Emigration (constant across years, conditional on survival)
  for(i in 1:nsites){
    for(k in 1:nspecies){
      for(t in 1:(nyears-1)){
        E[i,k,t+1] ~ dbin(kappa[i,k], S[i,k,t+1]) #emigration of spp. k from site i (constant across years)
      }#t
      logit(kappa[i,k]) <- alpha.kappa[k] + beta.kappa[k] * cov.kappa[i]
      #Link function for site/spp.-cov influence on relative emigration of spp. k from site i
    }#k
  }#i
  
  ##Immigration (constant across years, conditional on survival and emigration)
  for(i in 1:nsites){
    for(k in 1:nspecies){
      for(t in 1:(nyears-1)){
        I[i,k,t+1] ~ dpois(Ilambda[i,k,t+1]) #immigration of spp. k into site i (constant across years)
        Ilambda[i,k,t+1] <- sum(E[1:nsites,k,t+1]/nadj[1:nsites]*adj[1:nsites,i])
      }#t
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

#Distance Decay (Formulation 2)
cat(file = "SHARP_Spatial_DD_DCAM.txt", " 
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
  
  # *** Priors and hyperpriors for model on kappa (emigration) ***
  # Priors
  for(k in 1:nspecies){ # Loop over # of species 
    alpha.kappa[k] ~ dnorm(mu.alpha.kappa, tau.alpha.kappa) #Intercepts
    beta.kappa[k] ~ dnorm(mu.beta.kappa, tau.beta.kappa) #Coeffs
  }
  
  # Hyperpriors
    #Intercept
    mu.alpha.kappa <- logit(mean.kappa) 
    mean.kappa ~ dunif(0, 1) #mean yearly emigration rate across species (conditional on survival)
    sd.alpha.kappa ~ dnorm(0, 0.1)I(0,)
    tau.alpha.kappa <- 1/(sd.alpha.kappa*sd.alpha.kappa)
    
    #Coefficient
    mu.beta.kappa ~ dnorm(0, 0.1) #No yearly variation in hyperprior for coeff
    sd.beta.kappa ~ dnorm(0, 0.1)I(0,)
    tau.beta.kappa <- 1/(sd.beta.kappa*sd.beta.kappa)
    
  # *** Priors and hyperpriors for model on theta (distance decay rate for emigrants) ***
  # *** Note the absence of a coefficient, site / survey conditions not predicted to affect ***
  # *** distance decay function (species only expected to impact distance's effect on emigration) ***
  # Priors
  for(k in 1:nspecies){ # Loop over # of species 
    alpha.theta[k] ~ dnorm(mu.alpha.theta, tau.alpha.theta) #Intercepts
  }
  
  # Hyperpriors
    #Intercept
    mu.alpha.theta <- log(mean.theta) 
    mean.theta ~ dunif(0, 10) #mean yearly distance decay rate across species (conditional on survival)
    sd.alpha.theta ~ dnorm(0, 0.1)I(0,)
    tau.alpha.theta <- 1/(sd.alpha.theta*sd.alpha.theta)
    
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
        N[i,k,t+1] <- S[i,k,t] - E[i,k,t] + R[i,k,t] + I[i,k,t] #t+1 abundance of species k at site i
      }#t
    }#k
  }#i
  
  ##Survival (constant across years)
  for(i in 1:nsites){
    for(k in 1:nspecies){
      for(t in 1:(nyears-1)){
        S[i,k,t] ~ dbin(phi[i,k], N[i,k,t]) #survival of spp. k at site i (constant across years)
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
        R[i,k,t] ~ dpois(tmp[i,k,t]) #relative recruitment of spp. k at site i (constant across years)
      }#t
      log(gamma[i,k]) <- alpha.gamma[k] + beta.gamma[k] * cov.gamma[i]
      #Link function for site/spp.-cov influence on relative recruitment of spp. k at site i
    }#k
  }#i
 
  ##Emigration (constant across years, conditional on survival)
  for(i in 1:nsites){
    for(k in 1:nspecies){
      for(t in 1:(nyears-1)){
        E[i,k,t] ~ dbin(kappa[i,k], S[i,k,t]) #emigration of spp. k from site i (constant across years)
      }#t
      logit(kappa[i,k]) <- alpha.kappa[k] + beta.kappa[k] * cov.kappa[i]
      #Link function for site/spp.-cov influence on relative emigration of spp. k from site i
    }#k
  }#i
  
  ##Distance Decay function (constant across years, conditional on emigration)
  for(i in 1:nsites){
    for(j in 1:nsites){
      for(k in 1:nspecies){
      w_ori[i,j,k] <- ifelse(i==j, 1e-6, exp(-1 * alpha.theta[k] * dist[i,j]))
      w[i,j,k] <- w_ori[i,j,k] / sum(w_ori[i,1:nsites,k])
      }#k
    }#j
  }#i

  ##Immigration (constant across years, conditional on survival and emigration)
  for(i in 1:nsites){
    for(k in 1:nspecies){
      for(t in 1:(nyears-1)){
        I[i,k,t] ~ dpois(sum(M[i,1:nsites,k,t])) #immigration of spp. k into site i (constant across years)
        for(j in 1:nsites){
          M[i,j,k,t] <- E[i,k,t] * w[i,j,k]
        }#j
      }#t
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
nyears<-jags.data$nyears
Rst <- apply(jags.data$y, c(1,2,4), max) + 10
Rst[,,nyears] <- NA #No recruitment in final year 
Rst<-Rst[,,-nyears]
Sst <- apply(jags.data$y, c(1,2,4), max) + 10
Sst[,,nyears] <- NA #No survival in final year 
Sst<-Sst[,,-nyears]
Ist <- apply(jags.data$y, c(1,2,4), max) + 5
Ist[,,nyears] <- NA #No immigration in final year 
Ist<-Ist[,,-nyears]
Est <- apply(jags.data$y, c(1,2,4), max) + 5
Est[,,nyears] <- NA #No emmigration in final year 
Est<-Est[,,-nyears]
Nst <- apply(jags.data$y, c(1,2,4), max) + 30
Nst[,,-1] <- NA #Not initializing after year 1

inits <- function() list(R=Rst, N=Nst, I=Ist)
                         
#########################
####Parameters to Save###
#########################

params<-c("mean.lambda", "mean.phi", "mean.gamma", "mean.kappa", "mean.theta", "mean.p", 
          "sd.alpha.lambda", "sd.alpha.phi", "sd.alpha.gamma", "sd.alpha.kappa","sd.alpha.theta",
          "sd.alpha.p", "mu.beta.lambda", "mu.beta.phi", "mu.beta.gamma", "mu.beta.kappa", "mu.beta.p", 
          "N")

###################
###Running Model###
###################

#Test run for syntax examination
#out<-jags(jags.data, inits, params, "SHARP_Spatial_DD_DCAM.txt",
#          n.chains = 3, n.adapt=500, n.burnin = 1000, 
#          n.iter =2000, n.thin =3, parallel=FALSE)

#Full run (ART = 4.5 hrs)
out<-jags(jags.data, inits, params, "SHARP_Spatial_DD_DCAM.txt",
          n.chains = 3, n.adapt=5000, n.burnin = 50000, 
          n.iter =500000, n.thin =3, parallel=TRUE)

##################################
###Evaluating Model Performance###
##################################

#Percent Error of posteriors

#Community Hyperparameters
((out$mean$mean.lambda-sim.data$mean.lambda)/sim.data$mean.lambda)*100 #Abundance accuracy (slightly overestimating)
((out$mean$sd.alpha.lambda-sim.data$sd.lambda)/sim.data$sd.lambda)*100 #Abundance variance accuracy (slightly underestimating)
((out$mean$mean.phi-sim.data$mean.phi)/sim.data$mean.phi)*100 #Survival accuracy (underestimating)
((out$mean$sd.alpha.phi-sim.data$sd.phi)/sim.data$sd.phi)*100 #Survival variance accuracy (greatly overestimating)
((out$mean$mean.gamma-sim.data$mean.gamma)/sim.data$mean.gamma)*100 #Recruitment accuracy (greatly overestimating)
((out$mean$sd.alpha.gamma-sim.data$sd.gamma)/sim.data$sd.gamma)*100 #Recruitment variance accuracy (overestimating)
((out$mean$mean.kappa-sim.data$mean.kappa)/sim.data$mean.kappa)*100 #Emigration accuracy (accurate)
((out$mean$sd.alpha.kappa-sim.data$sd.kappa)/sim.data$sd.kappa)*100 #Emigration variance accuracy (greatly overestimating)
((out$mean$mean.theta-sim.data$mean.theta)/sim.data$mean.theta)*100 #Distance decay accuracy (greatly overestimating)
((out$mean$sd.alpha.theta-sim.data$sd.theta)/sim.data$sd.theta)*100 #Distance decay variance accuracy (greatly overestimating)
((out$mean$mean.p-sim.data$mean.p)/sim.data$mean.p)*100 #Detection accuracy (underestimating)
((out$mean$sd.alpha.p-sim.data$sd.p)/sim.data$sd.p)*100 #Detection variance accuracy (greatly overestimating)

#Species Abundances 
#Abundances of each species across sampled area each year (spp. x year)
Ntot.t<-apply(sim.data$N,c(2,3),sum) #Latent (true) abundances
Ntot.e<-apply(out$mean$N,c(2,3),sum) #Estimated abundances
((Ntot.e-Ntot.t)/Ntot.t)*100 #Abundance accuracy (overestimating all, some greatly)

#Calculating Growth Rates (Continuous)
r.t<-rep(NA, jags.data$nspecies)
r.e<-rep(NA, jags.data$nspecies)
for(k in 1:jags.data$nspecies){
  r.t[k]<-(log(Ntot.t[k,jags.data$nyears]/Ntot.t[k,1])/(jags.data$nyears-1))*100
  r.e[k]<-(log(Ntot.e[k,jags.data$nyears]/Ntot.e[k,1])/(jags.data$nyears-1))*100
}
((r.e-r.t)/r.t)*100 #Growth rate accuracy (pretty accurate, one spp. greatly underestimating)

#END