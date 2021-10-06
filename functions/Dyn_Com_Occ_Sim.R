###Simulating Data to Test Out a SIMPLE Dynamic Occupancy Model
###For later use in SHARP survey analyses
###Code available at: https://github.com/mikemeredith/AHM_code/tree/master/AHM2_ch05

###Load Required Packages###
library(AHMbook)
library(jagsUI)


###Simulate Required Dataset###
##Parameters selected based on thinned SHARP long-term dataset
str(data <- simDCM(#Species, sites, surveys, years
                   nspec = 100, nsites = 50, nsurveys = 2, nyears = 10,
                   #Mean initial occupancy and variance in occupancy among species (intercept)
                   mean.psi1 = 0.4, sig.lpsi1 = 1, 
                   #Species-specific variance in occupancy among sites (site-cov coefficient)
                   mu.beta.lpsi1 = 0, sig.beta.lpsi1 = 1,
                   #Annual variation in persistence
                   range.mean.phi = c(0.6, 0.9), sig.lphi = 1, 
                   #Species-specific variance in persistence among sites (site-cov coefficient)
                   mu.beta.lphi = 0, sig.beta.lphi = 1, 
                   #Annual variation in colonization
                   range.mean.gamma = c(0.1, 0.3), sig.lgamma = 1,
                   #Species-specific variance in colonization among sites (site-cov coefficient)
                   mu.beta.lgamma = 0, sig.beta.lgamma = 1, 
                   #Annual variation in detection
                   range.mean.p = c(0.4, 0.6), sig.lp = 1, 
                   #Species-specific variance in detection among sites 
                   mu.beta.lp = 0, sig.beta.lp = 1, 
                   #Annual variation in detection also?
                   range.beta1.survey = c(0, 0), range.beta2.survey = c(0, 0), 
                   #Temporal trends in detection heterogeneity at the site or survey level
                   trend.sd.site = c(0, 0), trend.sd.survey = c(0, 0), 
                   #Shows plots of data?
                   show.plot = TRUE) )

###Prepare Datasets for Model Input###
#The actual dataset that needs to be input into the model
str(data$y) 
y<-data$y #save as separate object
str(jags.data <- list(y = y, nsite = dim(y)[1], nsurvey = dim(y)[2],
                  nyear = dim(y)[3], nspec = dim(y)[4]))

###Model Specification###
# Specify model in BUGS language
cat(file = "DCM.txt", "
model {
  # Specify priors: Declare species-level effects as random (reduces RMSE for parameters vs. fixed effects)
  for(k in 1:nspec){               # Loop over species
    logit(psi1[k]) <- lpsi1[k]     # Initial occupancy
    lpsi1[k] ~ dnorm(mu.lpsi1, tau.lpsi1)
    logit(phi[k]) <- lphi[k]       # Persistence
    lphi[k] ~ dnorm(mu.lphi, tau.lphi)
    logit(gamma[k]) <- lgamma[k]   # Colonization
    lgamma[k] ~ dnorm(mu.lgamma, tau.lgamma)
    logit(p[k]) <- lp[k]           # Detection
    lp[k] ~ dnorm(mu.lp, tau.lp)
  }
  # Specify hyperpriors: Priors for the hyperparameters
  mu.lpsi1 <- logit(mean.psi1)     # Initial occupancy
  mean.psi1 ~ dunif(0, 1)
  tau.lpsi1 <- pow(sd.lpsi1, -2)
  sd.lpsi1 ~ dunif(0, 10)
  mu.lphi <- logit(mean.phi)       # Persistence
  mean.phi ~ dunif(0, 1)
  tau.lphi <- pow(sd.lphi, -2)
  sd.lphi ~ dunif(0, 10)
  mu.lgamma <- logit(mean.gamma)   # Colonization
  mean.gamma ~ dunif(0, 1)
  tau.lgamma <- pow(sd.lgamma, -2)
  sd.lgamma ~ dunif(0, 10)
  mu.lp <- logit(mean.p)           # Detection
  mean.p ~ dunif(0, 1)
  tau.lp <- pow(sd.lp, -2)
  sd.lp ~ dunif(0, 10)
  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsite){              # Loop over sites
    for(k in 1:nspec){             # Loop over species
      # Initial conditions of system
      z[i,1, k] ~ dbern(psi1[k])   # Presence/absence at start of study
      # State transitions
      for (t in 2:nyear){          # Loop over years
        z[i,t,k] ~ dbern(z[i,t-1,k] * phi[k] + (1-z[i,t-1, k]) * gamma[k])
      }
    }
  }
  # Observation model
  for (i in 1:nsite){              # Loop over sites
    for(k in 1:nspec){             # Loop over species
      for (j in 1:nsurvey){        # Loop over surveys
        for (t in 1:nyear){        # Loop over years
          y[i,j,t,k] ~ dbern(z[i,t,k] * p[k])
        }
      }
    }
  }
  # Derived parameters: Number of occupied sites and population occupancy
  for(k in 1:nspec){               # Loop over species
    n.occ[1, k] <- sum(z[,1,k])    # Number of occupied sites
    psi[1, k] <- psi1[k]           # Population occupancy
    for (t in 2:nyear){            # Loop over years
      n.occ[t, k] <- sum(z[,t,k])
      psi[t, k] <- psi[t-1, k] * phi[k] + (1-psi[t-1, k]) * gamma[k]
    }
  }
}
")

###Initials/Parameters/MCMC Settings###
zst <- apply(y, c(1,3,4), max) # Observed occurrence as inits for z
inits <- function(){ list(z = zst)}

##Parameters monitored (could also add "z")
params <- c("mu.lpsi1", "sd.lpsi1", "mu.lphi", "sd.lphi", "mu.lgamma",
            "sd.lgamma", "mu.lp", "sd.lp", "psi1", "phi", "gamma", "p", "n.occ", "psi")

##MCMC settings
# na <- 1000 ; ni <- 6000 ; nt <- 6 ; nb <- 3000 ; nc <- 3
na <- 1000 ; ni <- 600 ; nt <- 1 ; nb <- 300 ; nc <- 3  # ~~~~ for testing, 7 mins


###Running Model

##Call JAGS (ART 188 min), check convergence and summarize posteriors
out <- jags(jags.data, inits, params, "DCM.txt", n.adapt = na, n.chains = nc,
             n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T)

##Observe results of estimated parameters vs. true values
op <- par(mfrow = c(4,2))
plot(density(plogis(out$sims.list$mu.lpsi1)), lwd = 2, col = 'gray',
     main = 'Community mean of lpsi1', xlab = 'mu.lpsi1', ylab = 'Density',
     xlim = c(0, 0.5), frame = F)
abline(v = data$mean.psi1, col = 'black', lwd = 2)
plot(density(out2$sims.list$sd.lpsi1), lwd = 2, col = 'gray',
     main = 'Community SD of lpsi1', xlab = 'sd.lpsi1', ylab = 'Density',
     xlim = c(0, 4), frame = F)
abline(v = dat$sig.lpsi1, col = 'black', lwd = 2)

plot(density(plogis(out2$sims.list$mu.lphi)), lwd = 2, col = 'gray',
     main = 'Community mean of lpsi1', xlab = 'mu.lphi', ylab = 'Density',
     xlim = c(0.4, 1), frame = F)
abline(v = dat$range.mean.phi[1], col = 'black', lwd = 2)
plot(density(out2$sims.list$sd.lphi), lwd = 2, col = 'gray',
     main = 'Community SD of lphi', xlab = 'sd.lphi', ylab = 'Density',
     xlim = c(0, 2), frame = F)
abline(v = dat$sig.lphi, col = 'black', lwd = 2)

plot(density(plogis(out2$sims.list$mu.lgamma)), lwd = 2, col = 'gray',
     main = 'Community mean of lgamma', xlab = 'mu.lgamma', ylab = 'Density',
     xlim = c(0, 0.3), frame = F)
abline(v = dat$range.mean.gamma[1], col = 'black', lwd = 2)
plot(density(out2$sims.list$sd.lgamma), lwd = 2, col = 'gray',
     main = 'Community SD of lgamma', xlab = 'sd.lgamma', ylab = 'Density',
     xlim = c(0, 2), frame = F)
abline(v = dat$sig.lgamma, col = 'black', lwd = 2)

plot(density(plogis(out2$sims.list$mu.lp)), lwd = 2, col = 'gray',
     main = 'Community mean of lp', xlab = 'mu.lp', ylab = 'Density',
     xlim = c(0, 1), frame = F)
abline(v = dat$range.mean.p[1], col = 'black', lwd = 2)
plot(density(out2$sims.list$sd.lp), lwd = 2, col = 'gray',
     main = 'Community SD of lp', xlab = 'sd.lp', ylab = 'Density',
     xlim = c(0, 4), frame = F)
abline(v = dat$sig.lp, col = 'black', lwd = 2)
par(op)


##Simulating data to test out a much more complicated iteration of a Dynamic Community Occupancy Model
###This model is a modified version of the code presented in Chapter 5.7 of AHM vol. 2

##############################################
###Load Required Packages and Set Directory###
##############################################

library(AHMbook)
library(jagsUI)
setwd("C:\\Users\\12152\\Documents\\UConn_Files\\Project_Files\\Data_Analysis\\Survey_Data_Analyses")

###############################
###Simulate Required Dataset###
###############################

##Parameters selected based on thinned SHARP long-term dataset

str(data <- simDCM(#Species, sites, surveys, years
              nspec = 100, nsites = 20, nsurveys = 2, nyears = 3,
              #Mean initial occupancy and variance in occupancy among species (intercept)
              mean.psi1 = 0.6, sig.lpsi1 = 3, #large variance in occ
              #Species-specific variance in occupancy among sites (site-cov coefficient)
              mu.beta.lpsi1 = 0.3, sig.beta.lpsi1 = 1, #slight positive association
              #Annual variation in persistence
              range.mean.phi = c(0.6, 0.9), sig.lphi = 1, 
              #Species-specific variance in persistence among sites (site-cov coefficient)
              mu.beta.lphi = -0.2, sig.beta.lphi = 1, #slight negative association
              #Annual variation in colonization
              range.mean.gamma = c(0.1, 0.3), sig.lgamma = 1,
              #Species-specific variance in colonization among sites (site-cov coefficient)
              mu.beta.lgamma = 0.3, sig.beta.lgamma = 1, #slight positive association
              #Annual variation in detection
              range.mean.p = c(0.1, 0.2), sig.lp = 3, #large variance in det 
              #Species-specific variance in detection among sites (site-cov coefficient)
              mu.beta.lp = -0.3, sig.beta.lp = 1,  #slight negative association
              #Annual variation in detection also?
              range.beta1.survey = c(0, 0), range.beta2.survey = c(0, 0), 
              #Temporal trends in detection heterogeneity at the site or survey level
              trend.sd.site = c(0, 0), trend.sd.survey = c(0, 0), 
              #Shows plots of data (should usually be FALSE)
              show.plot = FALSE))

###Manipulate data output so ready to be input into model code###

# Species detected in 0 years (missed)
missed <- which(data$nyears.det == 0)      

# Toss out species never detected
y <- data$y                  # copy 4D array
y <- y[,,,-missed]           # Drop data from 19 missed species
str(y)                       # 20 sites x 2 reps x 3 years x # species seen

# Get sample sizes
nsites <- dim(y)[1]
nsurveys <- dim(y)[2]
nyears <- dim(y)[3]
nspec <- dim(y)[4]

# Augment the data set with 100 additional potential species
# Create zero 4D array for M species
nz <- 100                        # Number of 'zero species'
M <- nspec + nz                  # Size of augmented data set
yaug <- array(0, dim = c(nsites, # Prefill with zeroes
            nsurveys, nyears, M)) 
dim(yaug)                        # Check if it went well: and it did !

# Fill in the observed data into this larger array
yaug[,,,1:dim(y)[4]] <- y
str(yaug)
sum(y) ; sum(yaug)               # Quick sum check: should give same sum

# Get covariates (already scaled and standardized from simulation function)
high.marsh <- as.vector(data$Xpsi1) #occupancy covariate representing amount of high marsh in year 1
tide <- data$Xphi #persistence covariate representing increase in tidal innundation
elev <- data$Xgamma #colonization covariate representing increase in elevation
wind <- data$Xp #detection covariate representing wind

# Bundle and summarize data set
str(jags.data <- list(yaug = yaug, nsites = nsites, nsurveys = nsurveys,
                  nyears = nyears, M = M, psi.cov = high.marsh,
                  phi.cov = tide, gamma.cov = elev, p.cov = wind))

###########################
###Specify model in JAGS### 
###########################

cat(file = "DCM_DA_Cov.txt", "
model {

  ##########################
  #Priors and Hyperpriors###
  ##########################
  
  # *** Priors and hyperpriors for model on psi1 (initial occupancy) ***
  # Priors
  for(k in 1:M){ # Loop over # of species in augmented list (nspec + nz)
    alpha.lpsi1[k] ~ dnorm(mu.alpha.lpsi1, tau.alpha.lpsi1) # Intercept
    beta.lpsi1[k] ~ dnorm(mu.beta.lpsi1, tau.beta.lpsi1) #Cov coeffs (spp. specific)
  }
  # Hyperpriors
  mu.alpha.lpsi1 <- logit(mean.alpha.psi1) #Intercept mean
  mean.alpha.psi1 ~ dunif(0, 1) #Intercept mean
  tau.alpha.lpsi1 <- pow(sd.alpha.lpsi1, -2) #Intercept tau
  sd.alpha.lpsi1 ~ dunif(0, 10) #Intercept SD
  
  mu.beta.lpsi1 ~ dnorm(0, 0.1) #Coeff mean
  tau.beta.lpsi1 <- pow(sd.beta.lpsi1, -2)  #Coeff Tau
  sd.beta.lpsi1 ~ dnorm(0, 0.1)I(0,) # Half-Normal prior for Coeff SD

  # *** Priors and hyperpriors for model on phi (colonization) ***
  # Priors
  for(k in 1:M){ # Loop over # of species in augmented list (nspec + nz)
    for(t in 1:(nyears-1)){ # Loop over t-1 intervals (colonization varies by year)
      alpha.lphi[t,k] ~ dnorm(mu.alpha.lphi[t], tau.alpha.lphi) #Intercept
      # phi intercept for 3 intervals, different mean, same variance
    }
    beta.lphi[k] ~ dnorm(mu.beta.lphi, tau.beta.lphi) #Cov coeffs (spp. specific)
  }
  # Hyperpriors
  for(t in 1:(nyears-1)){ # Loop over t-1 intervals (yearly variation in mean)
    mu.alpha.lphi[t] <- logit(mean.alpha.phi[t]) 
    mean.alpha.phi[t] ~ dunif(0, 1)
  }
  tau.alpha.lphi <- pow(sd.alpha.lphi, -2) #No yearly variation in variance 
  sd.alpha.lphi ~ dnorm(0, 0.1)I(0,)
  
  mu.beta.lphi ~ dnorm(0, 0.01) #No yearly variation in hyperprior for coeffs
  tau.beta.lphi <- pow(sd.beta.lphi, -2)
  sd.beta.lphi ~ dnorm(0, 0.1)I(0,)
  
  # *** Priors and hyperpriors for model on gamma (persistence)***
  # Priors
  for(k in 1:M){ # Loop over # of species in augmented list (nspec + nz)
    for(t in 1:(nyears-1)){ # Loop over t-1 intervals (persistence varies by year)
      alpha.lgamma[t,k] ~ dnorm(mu.alpha.lgamma[t], tau.alpha.lgamma) #Intercept
      # gamma intercept for 3 intervals, different mean, same variance
    }
    beta.lgamma[k] ~ dnorm(mu.beta.lgamma, tau.beta.lgamma) #Cov coeffs (spp. specific)
  }
  # Hyperpriors
  for(t in 1:(nyears-1)){ # Loop over t-1 intervals (yearly variation in mean)
    mu.alpha.lgamma[t] <- logit(mean.alpha.gamma[t])
    mean.alpha.gamma[t] ~ dunif(0, 1)
  }
  tau.alpha.lgamma <- pow(sd.alpha.lgamma, -2) #No yearly variation in variance
  sd.alpha.lgamma ~ dnorm(0, 0.1)I(0,)
  
  mu.beta.lgamma ~ dnorm(0, 0.1)
  tau.beta.lgamma <- pow(sd.beta.lgamma, -2)
  sd.beta.lgamma ~ dnorm(0, 0.1)I(0,)
  
  # *** Priors and hyperpriors for model on p (detection) ***
  # Priors
  for(k in 1:M){ # Loop over # of species in augmented list (nspec + nz)
    for(t in 1:nyears){ # Loop over t intervals (detection varies by year)
      alpha.lp[t,k] ~ dnorm(mu.alpha.lp[t], tau.alpha.lp) # Intercept
      # p intercept for t years, different mean, same variance
    }
    beta.lp[k] ~ dnorm(mu.beta.lp, tau.beta.lp) #Cov coeffs (spp. specific)
  }
  # Hyperpriors
  for(t in 1:nyears){ # Loop over t years (yearly variation in mean)
    mu.alpha.lp[t] <- logit(mean.alpha.p[t])
    mean.alpha.p[t] ~ dunif(0, 1)
  }
  tau.alpha.lp <- pow(sd.alpha.lp, -2) #No yearly variation in variance
  sd.alpha.lp ~ dnorm(0, 0.1)I(0,)
  
  mu.beta.lp ~ dnorm(0, 0.1)
  tau.beta.lp <- pow(sd.beta.lp, -2)
  sd.beta.lp ~ dnorm(0, 0.1)I(0,)
  
  ######################
  ###Model Likelihood###
  ######################
  
  ###Data augmentation submodel likelihood
  
  omega ~ dunif(0, 1) # Prior for data augmentation parameter
  # omega ~ dbeta(0.001, 1) # Scale prior (Link, 2013) as an alternative
  for(k in 1:M){ # Loop over # of species in augmented list (nspec + nz)
    w[k] ~ dbern(omega) #Presence/absence of species in metacommunity 
  }
  
  ###Ecological submodel: Define state conditional on parameters
  
  for (i in 1:nsites){ # Loop over # of sites
    for(k in 1:M){ # Loop over # of species in augmented list
      # Initial conditions of system (includes covariate effects)
      z[i,1, k] ~ dbern(psi1[i, k]) #Initial occupancy for year 1
      
      #Logit link for year 1 occupancy probability with cov
      logit(psi1[i,k]) <- alpha.lpsi1[k] + beta.lpsi1[k] * psi.cov[i] 
      
      # State transitions (includes covariate effects)
      for (t in 2:nyears){ # Loop over # of years
        # Transition in occupancy probability between years due 
        # to prior colonization and persistence probabilities 
        z[i,t,k] ~ dbern(z[i,t-1,k]*phi[i,t-1,k] + (1-z[i,t-1, k])*gamma[i,t-1,k])
        
        #Prior colonization probability (intercept varies but relationship with cov is constant)
        logit(phi[i,t-1,k]) <- alpha.lphi[t-1,k] + beta.lphi[k] * phi.cov[i,t-1]  
        
        #Prior persistence probability (intercept varies but relationship with cov is constant)
        logit(gamma[i,t-1,k]) <- alpha.lgamma[t-1,k] + beta.lgamma[k] * gamma.cov[i,t-1]
      }
    }
  }
  
  ###Observation model (incl. covariate effects)
  
  # Multiplication with w is now here (improves mixing if incorporated in det rather than occ)
  for (i in 1:nsites){
    for(k in 1:M){
      for (j in 1:nsurveys){
        for (t in 1:nyears){
          # Observed dataset as a function of occupancy and detection (and multiplied by w here)
          yaug[i,j,t,k] ~ dbern(w[k] * z[i,t,k] * p[i,j,t,k]) 
          
          # Logit-link for detection with covariate
          logit(p[i,j,t,k]) <- alpha.lp[t,k] + beta.lp[k] * p.cov[i,j,t] 
        }
      }
    }
  }
  
  ########################
  ###Derived Parameters###
  ########################
  
  # Number of occupied sites per year
  # (multiplication with w here because z not modified by w earlier)
  for(k in 1:M){
    for (t in 1:nyears){
      n.occ[t, k] <- sum(w[k] * z[,t,k])
    }
  }
  # Species richness: total and per site/year
  Ntotal <- sum(w[]) # Total species richness (community size)
  for(i in 1:nsites){
    for(t in 1:nyears){
      for(k in 1:M){
        tmp[i,t,k] <- w[k] * z[i,t,k]
      }
      Nspec[i,t] <- sum(tmp[i,t,]) # Species richness per site and year (this is the informative value!)
    }
  }
} 
###END OF MODEL CODE
")

#######################################
###Initials/Parameters/MCMC Settings###
#######################################

# Initial values (simply initialize all at 1)
zst <- array(1, dim = c(nsites, nyears, M))
wst <- rep(1, M)
inits <- function(){ list(w = wst, z = zst)}

# Parameters monitored
params <- c("omega", "mu.alpha.lpsi1", "sd.alpha.lpsi1", "mu.beta.lpsi1",
            "sd.beta.lpsi1", "mu.alpha.lphi", "sd.alpha.lphi", "mu.beta.lphi",
            "sd.beta.lphi", "mu.alpha.lgamma", "sd.alpha.lgamma", "mu.beta.lgamma",
            "sd.beta.lgamma", "mu.alpha.lp", "sd.alpha.lp", "mu.beta.lp",
            "sd.beta.lp", "alpha.lpsi1", "beta.lpsi1", "alpha.lphi", "beta.lphi",
            "alpha.lgamma", "beta.lgamma", "alpha.lp", "beta.lp", "Ntotal", "Nspec",
            "n.occ")                             # Could add "z"

# MCMC settings
# na <- 1000 ; ni <- 150000 ; nt <- 50 ; nb <- 100000 ; nc <- 3
na <- 1000 ; ni <- 1500 ; nt <- 1 ; nb <- 1000 ; nc <- 3  # ~~~ for testing, 5 mins

###################
###Running Model###
###################

# Call JAGS (ART 386 min), check convergence and summarize posteriors
out <- jags(jags.data, inits, params, "DCM_DA_Cov.txt", n.adapt = na,
            n.iter = ni, n.thin = nt, n.burnin = nb, n.chains = nc, parallel = TRUE)

op <- par(mfrow = c(3,3)) 
par(op) ; traceplot(out)

out$mean$Ntotal #should be very close to nspecies

################################################################################
cat(file = "DCM_DA_Multiple_Covs.txt", "
model {

  ##########################
  #Priors and Hyperpriors###
  ##########################
  
  # *** Priors and hyperpriors for model on psi1 (initial occupancy) ***
  # Priors
  for(k in 1:M){ # Loop over # of species in augmented list (nspec + nz)
    alpha.lpsi1[k] ~ dnorm(mu.alpha.lpsi1, tau.alpha.lpsi1) # Intercept
      for(g in 1:4){ # Loop over # of coefficients (= number of covariates)
      beta.lpsi1[g, k] ~ dnorm(mu.beta.lpsi1[g], tau.beta.lpsi1[g]) #Cov coeffs (spp. specific)
    }
  }
  # Hyperpriors
  mu.alpha.lpsi1 <- logit(mean.alpha.psi1) #Intercept mean
  mean.alpha.psi1 ~ dunif(0, 1) #Intercept mean
  tau.alpha.lpsi1 <- pow(sd.alpha.lpsi1, -2) #Intercept tau
  sd.alpha.lpsi1 ~ dunif(0, 10) #Intercept SD
  for(g in 1:4){ # Loop over 4 coefficients
    mu.beta.lpsi1[g] ~ dnorm(0, 0.1) #Coeff mean
    tau.beta.lpsi1[g] <- pow(sd.beta.lpsi1[g], -2)  #Coeff Tau
    sd.beta.lpsi1[g] ~ dnorm(0, 0.1)I(0,) # Half-Normal prior for Coeff SD
    # curve(dnorm(x, 0, sqrt(1 / 0.1)), 0, 20) # howsit look like ?
  }
  # *** Priors and hyperpriors for model on phi (colonization) ***
  # Priors
  for(k in 1:M){ # Loop over # of species in augmented list (nspec + nz)
    for(t in 1:(nyears-1)){ # Loop over t-1 intervals (colonization varies by year)
      alpha.lphi[t,k] ~ dnorm(mu.alpha.lphi[t], tau.alpha.lphi) #Intercept
      # phi intercept for 3 intervals, different mean, same variance
    }
    for(g in 1:4){ # Loop over # of coefficients (= number of covariates)
      beta.lphi[g, k] ~ dnorm(mu.beta.lphi[g], tau.beta.lphi[g]) #Cov coeffs (spp. specific)
    }
  }
  # Hyperpriors
  for(t in 1:(nyears-1)){ # Loop over t-1 intervals (yearly variation in mean)
    mu.alpha.lphi[t] <- logit(mean.alpha.phi[t]) 
    mean.alpha.phi[t] ~ dunif(0, 1)
  }
  tau.alpha.lphi <- pow(sd.alpha.lphi, -2) #No yearly variation in variance 
  sd.alpha.lphi ~ dnorm(0, 0.1)I(0,)
  for(g in 1:4){ # Loop over # of coefficients (= number of covariates)
    mu.beta.lphi[g] ~ dnorm(0, 0.01) #No yearly variation in hyperpriors for coeffs
    tau.beta.lphi[g] <- pow(sd.beta.lphi[g], -2)
    sd.beta.lphi[g] ~ dnorm(0, 0.1)I(0,)
  }
  # *** Priors and hyperpriors for model on gamma (persistence)***
  # Priors
  for(k in 1:M){ # Loop over # of species in augmented list (nspec + nz)
    for(t in 1:(nyears-1)){ # Loop over t-1 intervals (persistence varies by year)
      alpha.lgamma[t,k] ~ dnorm(mu.alpha.lgamma[t], tau.alpha.lgamma) #Intercept
      # gamma intercept for 3 intervals, different mean, same variance
    }
    for(g in 1:4){ # Loop over # of coefficients (= number of covariates)
      beta.lgamma[g, k] ~ dnorm(mu.beta.lgamma[g], tau.beta.lgamma[g]) #Cov coeffs (spp. specific)
    }
  }
  # Hyperpriors
  for(t in 1:(nyears-1)){ # Loop over t-1 intervals (yearly variation in mean)
    mu.alpha.lgamma[t] <- logit(mean.alpha.gamma[t])
    mean.alpha.gamma[t] ~ dunif(0, 1)
  }
  tau.alpha.lgamma <- pow(sd.alpha.lgamma, -2) #No yearly variation in variance
  sd.alpha.lgamma ~ dnorm(0, 0.1)I(0,)
  for(g in 1:4){ # Loop over # of coefficients (= number of covariates)
    mu.beta.lgamma[g] ~ dnorm(0, 0.1)
    tau.beta.lgamma[g] <- pow(sd.beta.lgamma[g], -2)
    sd.beta.lgamma[g] ~ dnorm(0, 0.1)I(0,)
  }
  # *** Priors and hyperpriors for model on p (detection) ***
  # Priors
  for(k in 1:M){ # Loop over # of species in augmented list (nspec + nz)
    for(t in 1:nyears){ # Loop over t intervals (detection varies by year)
      alpha.lp[t,k] ~ dnorm(mu.alpha.lp[t], tau.alpha.lp) # Intercept
      # p intercept for t years, different mean, same variance
    }
    for(g in 1:6){ # Loop over # of coefficients (= number of covariates)
      beta.lp[g, k] ~ dnorm(mu.beta.lp[g], tau.beta.lp[g]) #Cov coeffs (spp. specific)
    }
  }
  # Hyperpriors
  for(t in 1:nyears){ # Loop over t years (yearly variation in mean)
    mu.alpha.lp[t] <- logit(mean.alpha.p[t])
    mean.alpha.p[t] ~ dunif(0, 1)
  }
  tau.alpha.lp <- pow(sd.alpha.lp, -2) #No yearly variation in variance
  sd.alpha.lp ~ dnorm(0, 0.1)I(0,)
  for(g in 1:6){ # Loop over # of coefficients (= number of covariates)
    mu.beta.lp[g] ~ dnorm(0, 0.1)
    tau.beta.lp[g] <- pow(sd.beta.lp[g], -2)
    sd.beta.lp[g] ~ dnorm(0, 0.1)I(0,)
  }
  
  ######################
  ###Model Likelihood###
  ######################
  
  ###Data augmentation submodel likelihood
  
  omega ~ dunif(0, 1) # Prior for data augmentation parameter
  # omega ~ dbeta(0.001, 1) # Scale prior (Link, 2013) as an alternative
  for(k in 1:M){ # Loop over # of species in augmented list (nspec + nz)
    w[k] ~ dbern(omega) #Presence/absence of species in metacommunity 
  }
  
  ###Ecological submodel: Define state conditional on parameters
  
  for (i in 1:nsites){ # Loop over # of sites
    for(k in 1:M){ # Loop over # of species in augmented list
      # Initial conditions of system (includes covariate effects)
      z[i,1, k] ~ dbern(psi1[i, k]) #Initial occupancy for year 1
      #Logit link for year 1 occupancy probability
      logit(psi1[i,k]) <- alpha.lpsi1[k] + 
      beta.lpsi1[1,k] * covs[i,1] + beta.lpsi1[2,k] * covs[i,2] +
      beta.lpsi1[3,k] * covs[i,3] + beta.lpsi1[4,k] * covs[i,4]
      # State transitions (includes covariate effects)
      for (t in 2:nyears){ # Loop over # of years
        # Transition in occupancy probability between years due 
        # to prior colonization and persistence probabilities 
        z[i,t,k] ~ dbern(z[i,t-1,k]*phi[i,t-1,k] + (1-z[i,t-1, k])*gamma[i,t-1,k])
        logit(phi[i,t-1,k]) <- alpha.lphi[t-1,k] + #Prior colonization probability 
        beta.lphi[1,k] * covs[i,1] + beta.lphi[2,k] * covs[i,2] +
        beta.lphi[3,k] * covs[i,3] + beta.lphi[4,k] * covs[i,4]
        logit(gamma[i,t-1,k]) <- alpha.lgamma[t-1,k] + #Prior persistence probability 
        beta.lgamma[1,k] * covs[i,1] + beta.lgamma[2,k] * covs[i,2] +
        beta.lgamma[3,k] * covs[i,3] + beta.lgamma[4,k] * covs[i,4]
      }
    }
  }
  
  ###Observation model (incl. covariate effects)
  
  # Multiplication with w is now here (improves mixing if incorporated in det rather than occ)
  # Also note the cosinor function for cyclic timeOfDay (= tod) covariate (not necessary if not using)
  for (i in 1:nsites){
    for(k in 1:M){
      for (j in 1:nsurveys){
        for (t in 1:nyears){
          # Observed dataset as a function of occupancy and detection (and multiplied by w here)
          yaug[i,j,t,k] ~ dbern(w[k] * z[i,t,k] * p[i,j,t,k]) 
          # Logit-link for detection with covariates
          logit(p[i,j,t,k]) <- alpha.lp[t,k] +
          beta.lp[1,k] * covs[i,1] + beta.lp[2,k] * covs[i,2] +
          beta.lp[3,k] * covs[i,3] + beta.lp[4,k] * covs[i,4] +
          beta.lp[5,k] * cos(2*pi*tod[i,j,t]/Tday) +
          beta.lp[6,k] * sin(2*pi*tod[i,j,t]/Tday)
        }
      }
    }
  }
  
  ########################
  ###Derived Parameters###
  ########################
  
  # Number of occupied sites per year
  # (multiplication with w here because z not modified by w earlier)
  for(k in 1:M){
    for (t in 1:nyears){
      n.occ[t, k] <- sum(w[k] * z[,t,k])
    }
  }
  # Species richness: total and per site/year
  Ntotal <- sum(w[]) # Total species richness (community size)
  for(i in 1:nsites){
    for(t in 1:nyears){
      for(k in 1:M){
        tmp[i,t,k] <- w[k] * z[i,t,k]
      }
      Nspec[i,t] <- sum(tmp[i,t,]) # Species richness per site and year (this is the informative value!)
    }
  }
} 
###END OF MODEL CODE
")
