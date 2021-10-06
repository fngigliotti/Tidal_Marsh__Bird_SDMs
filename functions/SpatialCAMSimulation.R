##########################################
###Multispecies Spatial Abundance Model###
##########################################

#Franco Gigliotti 
#03/01/2020

#################
###Description###
#################

#This code executes a simulation of a multispecies abundance dataset 
#and analyses the data using a multispecies spatial abundance model.
#The model is carried out in JAGS and model code is based off of a very
#similar model (single-species spatial abundance model; per AHM Book 9.4.2). 


#######################################
###Load Required Packages and Set WD###
#######################################

library(AHMbook)
library(jagsUI)
library(fields)
setwd(choose.dir()) #because lazy 
getwd() #where it's at
source("homebrewed_simulation_functions.R") #Load simulation function (requires AHMbook)

######################
###Simulate Dataset###
######################

##Use parameters relevant to SHARP data
data<-simSpatialCAM(nspecies = 50, nsites = 50, nsurveys = 3, mean.lambda = 3,
                    sig.loglam = 1, mu.beta.loglam = 1, sig.beta.loglam = 1,
                    mean.p = 0.3, sig.lp = 1, mu.beta.lp = -1, sig.beta.lp = 1,
                    grid.size = 50, variance.RF = 1, theta.RF = 10, show.plots = FALSE,
                    repoducible = TRUE)

############################
###Prepare Data for Model###
############################

# Scale both sets of coordinates
head(coordgrid <- scale(cbind(data$x.coord, data$y.coord))) # All grid cells
head(sitelocs <- coordgrid[data$surveyed.sites,])           # only surveyed cells

#Define number of knots necessary for spline generation
nd <- max(20,min(data$nsites/4,150)) #Ruppert et al. 2003 "best knot number"
#Generate knots (takes a few minutes depending on size of grid)
system.time(knots <- cover.design(R = coordgrid, nd = nd, nruns = 10)) 

#Define the Z matrix for the random effects/knot coefficients (code straight from AHM)
knotlocs <- knots$design
omega <- (e2dist(knotlocs, knotlocs)/10)^3
svd.omega <- svd(omega)
sqrt.omega <- t(svd.omega$v %*% (t(svd.omega$u)*sqrt(svd.omega$d)))
Zk <- (e2dist(coordgrid, knotlocs)/10)^3
Zmat <- t(solve(sqrt.omega, t(Zk)))

#Bundle data
y<-data$y.obs
str(bdata <- list(y = y, nsites = dim(y)[1], nsurveys = dim(y)[2], 
                  nspecies = dim(y)[3], site.cov = data$site.cov, 
                  survey.cov = data$survey.cov, n.knots = nd, Zmat = Zmat))

#########################
###Model Specification###
#########################

# Specify model in JAGS
cat(file = "CAM_GAM.txt", "
model {

##Community priors (with hyperparameters) for species-specific parameters
for(k in 1:nspecies){
  phi[k] ~ dunif(0,1)                              # Zero-inflation
  beta0[k] ~ dnorm(mu.beta0, tau.beta0)            # Abundance intercepts
  alpha0[k] ~ dnorm(mu.alpha0, tau.alpha0)         # Detection intercepts
  beta1[k] ~ dnorm(mu.beta1, tau.beta1)            # Abundance coeff (slope)
  alpha1[k] ~ dnorm(mu.alpha1, tau.alpha1)         # Detection coeff (slope)
}


##Hyperpriors for community hyperparameters
##Abundance model
#Intercept
mu.beta0 <- logit(mu.lbeta0)
mu.lbeta0 ~ dunif(0,1)
tau.beta0 <- pow(sd.beta0, -2)
sd.beta0 ~ dunif(0, 10) #probably should be fairly large?

#Coefficient for site-level (abundance) covariate
mu.beta1 <- logit(mu.lbeta1)
mu.lbeta1 ~ dunif(0,1)
tau.beta1 <- pow(sd.beta1, -2)
sd.beta1 ~ dunif(0, 10) #probably should be larger than sd.beta0?

##Detection model
#Intercept
mu.alpha0 <- logit(mu.lalpha0)
mu.lalpha0 ~ dunif(0,1)
tau.alpha0 <- pow(sd.alpha0, -2)
sd.alpha0 ~ dunif(0, 10) #probably should be fairly small

#Coefficient for survey-level (detection) covariate
mu.alpha1 <- logit(mu.lalpha1)
mu.lalpha1 ~ dunif(0,1)
tau.alpha1 <- pow(sd.alpha1, -2)
sd.alpha1 ~ dunif(0, 10) #probably a small sd value 

##Priors on random effects parameters representing the splines
for (k in 1:n.knots){
  b[k] ~ dnorm(0, tau.b)
}

# Prior on random effects dispersion
tau.b ~ dgamma(0.1, 0.1)
sd.b <- pow(tau.b, -2)

#Smoothing parameter for error term
for (i in 1:nsites){
  smooth2[i] <- inprod(Zmat[i,], b[])
  smooth[i] <- smooth2[i] - mean(smooth2[])
}

##Ecological model for true abundance (process model)
for (i in 1:nsites){
  for(k in 1:nspecies){
    a[i,k] ~ dbern(phi[k])   # zero-inflation
    N[i,k] ~ dpois(a[i,k] * lambda[i,k]) #latent abundance
    log(lambda[i,k]) <- beta0[k] + beta1[k] * site.cov[i] + smooth[i]
  }
}


##Observation model for replicated counts
for (i in 1:nsites){
  for (j in 1:nsurveys){
    for(k in 1:nspecies){
      y[i,j,k] ~ dbin(p[i,j,k], N[i,k])
      logit(p[i,j,k]) <- alpha0[k] + alpha1[k] * survey.cov[i,j] 
    }
  }
}

# Derived parameters: Total population size of each species on grid
for(k in 1:nspecies){
  N.super[k] <- sum(N[,k])
}

# Other derived quantities
for(k in 1:nspecies){
  mean.lambda[k] <- phi[k] * exp(beta0[k]) # Expected abundance on natural scale
  logit(mean.p[k]) <- alpha0[k]            # Mean detection on natural scale
}
}
")

################################
###Inits/Params/MCMC Settings###
################################

##Store variables for init setting
nspecies<-bdata$nspecies
nsites<-bdata$nsites

#Initial values
ast <- matrix(rep(1, nspecies*nsites), nrow = nsites) #set all sites occupied
Nst <- apply(y, c(1,3), max)
Nst[is.na(Nst)] <- 2
Nst[Nst == 0] <- 2
inits <- function(){list(a = ast, N = Nst, b = runif(bdata$n.knots, -0.5, 0.5))}

##Parameters Monitored
params <- c("phi", "mean.lambda", "beta0", "mu.beta0", "sd.beta0", "beta1",
            "mu.beta1", "sd.beta1", "mean.p", "alpha0", "mu.alpha0", 
            "sd.alpha0", "alpha1", "mu.alpha1", "sd.alpha1", "b", "sd.b",
            "smooth", "N.super")

##MCMC settings
 na <- 5000 ; ni <- 100000 ; nt <- 10 ; nb <- 20000 ; nc <- 3 #~10 hours, but not long enough
#na <- 1000 ; ni <- 3000 ; nt <- 2 ; nb <- 1500 ; nc <- 3  # ~~~ for testing (like an hour...)
 
###################
###Running Model###
###################

# Call JAGS (ART maybe like 10 hrs?), gauge convergence and summarize posteriors
out <- jags(bdata, inits, params, "CAM_GAM.txt", n.adapt = na,
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

######################
###Examining Output###
######################

#Results for improving prior determination
out$mean$mu.beta0 #1.4
out$mean$sd.beta0 #1
out$mean$mu.beta1 #1.3
out$mean$sd.beta1 #0.97
out$mean$mu.alpha0 #-0.9
out$mean$sd.alpha0 #1.1
out$mean$mu.alpha1 #-0.99
out$mean$sd.alpha1 #1
out$mean$sd.b #HUGE!

test<-rep(NA,50)
for(k in 1:50){
  test[k]<-sum(data$N.super[,k])
}

#Accuracy of estimated total population sizes of each species
summary(test/out$mean$N.super) #still underestimating some spp. and overestimating others
accuracy<-test/out$mean$N.super
detection<-out$mean$mean.p
abundance<-out$mean$mean.lambda
det.accuracy<-detection/apply(data$p,3,mean)
ab.accuracy<-abundance/apply(data$lambda,2,mean)
library(ggplot2)  
ggplot() +
  geom_density(aes(x=accuracy), fill = "gray") +
  geom_vline(xintercept = 1) +
  labs(x = "Estimated Species Abundance / True Abundance") +
  theme_classic()

ggplot() +
  geom_density(aes(x=det.accuracy), fill = "gray") +
  geom_vline(xintercept = 1) +
  labs(x = "Estimated Species Detection / True Detection") +
  theme_classic()


sum(test)/sum(out$mean$N.super) #underestimating counts overall for the grid

save(out, file = "CAM_GAM_2.R")
#Estimating phi, beta0, beta1, and hyperparams well, other parameters not so much
#Really having trouble with b and smooth estimates, how fix?
load(file = "CAM_GAM_2.R")

