####################################################################################
###Dynamic Community Abundance Model w/ Temporal Variation- SHARP 10 Year PC Data###
####################################################################################

#################
###Description###
#################

#This model analyzes 10 years of point-count data collected by SHARP and examines
#variation in rates of population change over time for all members of the focal 
#tidal marsh breeding bird community. Focal parameters (abundance, survival,
#recruitment, and detection) are specified as random effects and thus are allowed 
#to vary temporally. Density-dependence in survival and recruitment is also explicit
#in this model. Much of this model code is based off of code presented in 
#Applied Hierarchical Modeling Volumes 1 (2016) and 2 (2020) by Marc Kery and Andy Royle, 
#as well as the supporting R package 'AHMbook', in addition to code presented in the 
#following publications: Bellier et al. 2016, Zhao et al. 2019, and Zhao and Royle 2019.
#Future models will expand upon this preliminary model to incorporate spatial
#autocorrelation in data, integrate distance-sampling, CMR, and/or 
#time-to-detction data, and incorporate additional covariates such as vegetative community
#information, climatic variables, or elevation. CMR data will especially be useful as it
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
setwd(choose.dir()) #because lazy 
getwd() #where it's at
source("homebrewed_simulation_functions.R") #Load simulation functions (requires AHMbook)

########################################
###Simulate Dataset for Model Testing###
########################################

#Simple, farily species-poor community for model exploration. Observation-level parameters
#specified.  

sim.data<-simStochDCAM(nsites = 10, nspecies = 5, nsurveys = 3, nyears = 10, 
                  mean.lambda = 5, mean.phi = 0.4, mean.gamma = 0.6, mean.p = 0.6)

Nm<-apply(apply(sim.data$N,c(1,2),max),2,mean) #Nm for scaling (normally would use y-hat with true data)

jags.data<-list(nsites = sim.data$nsites, nspecies = sim.data$nspecies, nsurveys = sim.data$nsurveys,
                nyears = sim.data$nyears, y = sim.data$y.obs, Nm = Nm)

str(jags.data) #should be sites,species,surveys,years, obs array, and Nm

################################
###Specify Model Code in JAGS###
################################

#No covariates; density-dependence in survival and recruitment; survival, recruitment,
#and detection vary stochastically through time

cat(file = "SHARP_DCAM_No_Covs_Random.txt", " 
model {
  
  ##########################
  #Priors and Hyperpriors###
  ##########################
  
  # *** Priors and hyperpriors for model on lambda (t = 1 abundance) ***
  # Priors
  for(k in 1:nspecies){ # Loop over # of species 
    alpha.lambda[k] ~ dnorm(mu.alpha.lambda, tau.alpha.lambda) # Intercepts
  }
  
  # Hyperpriors
    #Intercept
    mu.alpha.lambda <- log(mean.lambda)
    mean.lambda ~ dunif(0, 50) #initial abundance (mean across species)
    sd.alpha.lambda ~ dnorm(0, 0.1)I(0,)
    tau.alpha.lambda <- 1/(sd.alpha.lambda*sd.alpha.lambda)
    
    
  # *** Priors and hyperpriors for model on phi (survival) ***
  # Priors
  for(k in 1:nspecies){ # Loop over # of species 
    alpha.phi[k] ~ dnorm(mu.alpha.phi, tau.alpha.phi) #Intercepts
    beta.phi[k] ~ dnorm(mu.beta.phi, tau.beta.phi)T(-1,1) #Coeffs
    for(i in 1:nsites){
      for(t in 1:nyears){
        eps.phi[i,k,t] ~ dnorm(0, tau.eps.phi) #Stochasticity
      }
    }
  }
  
  # Hyperpriors
    #Intercept
    mu.alpha.phi <- logit(mean.phi) 
    mean.phi ~ dunif(0, 1) #mean yearly survival across species
    sd.alpha.phi ~ dnorm(0, 0.1)I(0,)
    tau.alpha.phi <- 1/(sd.alpha.phi*sd.alpha.phi)
    
    #Coefficient
    mu.beta.phi ~ dnorm(0, 1) #No yearly variation in hyperprior for coeff w/dens. dep.
    sd.beta.phi ~ dnorm(0, 1)I(0,)
    tau.beta.phi <- 1/(sd.beta.phi*sd.beta.phi)
    
    #Random effect yearly stochasticity
    sd.eps.phi ~ dnorm(0, 0.1)I(0,)
    tau.eps.phi <- 1/(sd.eps.phi*sd.eps.phi)
    
  # *** Priors and hyperpriors for model on gamma (recruitment) ***
  # Priors
  for(k in 1:nspecies){ # Loop over # of species 
    alpha.gamma[k] ~ dnorm(mu.alpha.gamma, tau.alpha.gamma) #Intercepts
    beta.gamma[k] ~ dnorm(mu.beta.gamma, tau.beta.gamma)T(-1,1) #Coeffs
    for(i in 1:nsites){
      for(t in 1:nyears){
        eps.gamma[i,k,t] ~ dnorm(0, tau.eps.gamma) #Stochasticity
      }
    }
  }
  
  # Hyperpriors
    #Intercept
    mu.alpha.gamma <- log(mean.gamma) 
    mean.gamma ~ dunif(0, 5) #mean yearly recruitment across species
    sd.alpha.gamma ~ dnorm(0, 0.1)I(0,)
    tau.alpha.gamma <- 1/(sd.alpha.gamma*sd.alpha.gamma)
    
    #Coefficient
    mu.beta.gamma ~ dnorm(0, 1) #No yearly variation in hyperprior for coeff w/dens. dep.
    sd.beta.gamma ~ dnorm(0, 1)I(0,)
    tau.beta.gamma <- 1/(sd.beta.gamma*sd.beta.gamma)
    
    #Random effect yearly stochasticity
    sd.eps.gamma ~ dnorm(0, 0.1)I(0,)
    tau.eps.gamma <- 1/(sd.eps.gamma*sd.eps.gamma)
    
  # *** Priors and hyperpriors for model on p (detection) ***
  # Priors
  for(k in 1:nspecies){ # Loop over # of species 
    alpha.p[k] ~ dnorm(mu.alpha.p, tau.alpha.p) # Intercepts
    for(i in 1:nsites){
      for(t in 1:nyears){
        eps.p[i,k,t] ~ dnorm(0, tau.eps.p) #Stochasticity
      }
    }
  }
  
  # Hyperpriors
    #Intercept
    mu.alpha.p <- logit(mean.p)
    mean.p ~ dunif(0, 1) #Mean detection across species
    sd.alpha.p ~ dnorm(0, 0.1)I(0,)
    tau.alpha.p <- 1/(sd.alpha.p*sd.alpha.p)
    
    #Random effect yearly stochasticity
    sd.eps.p ~ dnorm(0, 0.1)I(0,)
    tau.eps.p <- 1/(sd.eps.p*sd.eps.p)

  ######################
  ###Model Likelihood###
  ######################
  
  ###Ecological submodel: Define state conditional on parameters
  
  ##Abundance
  #Initial Abundance (t=1)
  for(i in 1:nsites){
    for(k in 1:nspecies){
      N[i,k,1] ~ dpois(lambda[k]) #initial abundance of species k at site i (year 1)
      for(t in 1:(nyears-1)){
        N[i,k,t+1] <- S[i,k,t+1] + R[i,k,t+1] #t+1 abundance of species k at site i
      }#t
    }#k
  }#i
  
  #Log-link abundance (no covariates)
  for(k in 1:nspecies){
    log(lambda[k]) <- alpha.lambda[k] 
  }
  
  ##Survival (density dependence and yearly stochasticity)
  for(i in 1:nsites){
    for(k in 1:nspecies){
      for(t in 1:(nyears-1)){
        S[i,k,t+1] ~ dbin(phi[i,k,t], N[i,k,t]) #survival of spp. k at site i 
        logit(phi[i,k,t]) <- alpha.phi[k] + beta.phi[k]*(N[i,k,t]-Nm[k]) + eps.phi[i,k,t]
      }#t
    }#k
  }#i
  
  ##Recruitment (constant across years)
  for(i in 1:nsites){
    for(k in 1:nspecies){
      for(t in 1:(nyears-1)){
        tmp[i,k,t] <- N[i,k,t]*gamma[i,k,t]
        R[i,k,t+1] ~ dpois(tmp[i,k,t]) #relative recruitment of spp. k at site i (constant across years)
        log(gamma[i,k,t]) <- alpha.gamma[k] + beta.gamma[k]*(N[i,k,t]-Nm[k]) + eps.gamma[i,k,t]
      }#t
    }#k
  }#i
 
 
  ###Observation model (incl. covariate effects)
  
  ##Detection
  for (i in 1:nsites){
    for(k in 1:nspecies){
      for (t in 1:nyears){
        for (j in 1:nsurveys){
          # Observed dataset as a function of abundance and detection
          y[i,k,j,t] ~ dbin(p[i,k,t], N[i,k,t]) 
        }#j
        logit(p[i,k,t]) <- alpha.p[k] + eps.p[i,k,t]
      }#t
    }#k
  }#i
  
  ########################
  ###Derived Parameters###
  ########################
  
  #Total population of each species each year
  for(k in 1:nspecies){
    for(t in 1:nyears){
      N.tot[k,t]<-sum(N[,k,t])
    }#t
  }#k
  
  #Rates of range-wide population change for each species (exponential growth eq.)
  for(k in 1:nspecies){
    r.e[k]<-((N.tot[k,nyears]/N.tot[k,1])^(1/nyears))-1
  }
  
  #SDs of hyperpriors on natural scale
  log(sd.lambda)<-sd.alpha.lambda
  logit(sd.phi)<-sd.alpha.phi
  log(sd.gamma)<-sd.alpha.gamma
  logit(sd.p)<-sd.alpha.p
  
}#END OF MODEL CODE
")

######################
###Specify Initials###
######################

#Finding the correct values is very challenging and requires much trial 
#and error...
Rst <- apply(jags.data$y, c(1,2,4), max, na.rm =TRUE) + 10
Rst[Rst=='-Inf']<-10
Rst[,,1] <- NA #No recruitment in year 1
Nst <- apply(jags.data$y, c(1,2,4), max, na.rm = TRUE) + 20
Nst[Nst=='-Inf']<-20
Nst[,,2:dim(jags.data$y)[4]] <- NA
inits <- function() list(R=Rst, N=Nst)

#########################
####Parameters to Save###
#########################

params<-c("mean.lambda","sd.lambda", "mean.phi", "sd.phi", "mu.beta.phi", "sd.beta.phi",
          "sd.eps.phi", "mean.gamma", "sd.gamma", "mu.beta.gamma", 
          "sd.beta.gamma", "sd.eps.gamma", "mean.p", "sd.p", "sd.eps.p")

###################
###Running Model###
###################

#Test run for syntax examination (compilation may take a while...)
out<-jags(jags.data, inits, params, "SHARP_DCAM_No_Covs_Random.txt",
          n.chains = 3, n.adapt=100, n.burnin = 500, 
          n.iter =1000, n.thin =2, parallel=FALSE)

#Full model run (ART=105 hours)
#out<-jags(jags.data, inits, params, "SHARP_DCAM_No_Covs_Random.txt",
#          n.chains = 3, n.adapt=1000, n.burnin = 20000, 
#          n.iter =60000, n.thin =3, parallel=TRUE)

#Save model output in case of R crash (large-ish file output)
#save(out, file="SHARP_DCAM_No_Covs_Out.R")
#load(file="SHARP_DCAM_No_Covs_Out.R")

####################################
###Evaluating Parameter Estimates###
####################################

#Percentage of estimated latent variables that converged (should be 100% ideally)
sum(out$summary[,8]<1.1)/length(out$summary[,8])*100

#Traceplots to observe mixing
traceplot(out, params) 

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

##END