###################################
###Simulation Functions for SDMs###
###################################

#Franco Gigliotti 
#04/08/2021

#################
###Description###
#################

##This file contains functions to simulate datasets to be used in SDMs.
##Code for these functions rely heavily on functions in the 
##R package 'AHMBook' by Kery and Royle. To observe source code for these functions,
##uncomment and run the following two lines of code, or check out the source code
##for the 'AHMBook' package.

#print(getAnywhere('simNmixSpatial')) #outputs source code for S.S.A.M
#print(getAnywhere('simComm')) #outputs source code for C.A.M

#########################################################################
###Function #1: simulate data for a spatial community abundance model ###
#########################################################################

simSpatialCAM<-function (nspecies = 50, nsites = 50, nsurveys = 3, mean.lambda = 2, 
                         sig.loglam = 1, mu.beta.loglam = 1, sig.beta.loglam = 1, 
                         mean.p = 0.25, sig.lp = 1, mu.beta.lp = -1, sig.beta.lp = 1, 
                         grid.size = 50, variance.RF = 1, theta.RF = 10, 
                         show.plots = FALSE, repoducible = FALSE) 
{
  ##Function Arguments:
  #nspecies: True number of species in community  
  #nsites: Number of sampled sites
  #nsurveys: Number of replicated surveys conducted at each site
  #mean.lambda: Mean abundance of each species within each cell in the grid
  #sig.loglam: Variance among species specific abundances within each cell in the grid
  #mu.beta.loglam: Mean of coefficient value for site level covariate across species
  #sig.beta.loglam: Variance in coefficient values for site level covariate across species
  #mean.p: Mean detection probability of each species within each cell in the grid
  #sig.lp: Variance among species specific detection probabilities within each cell in the grid
  #mu.beta.lp: Mean of coefficient value for survey level covariate across species
  #sig.beta.lp: Variance in coefficient values for survey level covariate across species
  #grid.size: length of one dimension of area grid for focal region. If grid size is 50,
  #area available to be sampled is 50x50 = 2500 cells. nsites is then the number of 
  #cells in this grid that were sampled during survey efforts
  #variance.RF: Variance of Gaussian random field for negative exp. correlation between sites 
  #theta.RF: Parameter governing correlation in Gaussian random field; large theta = high correlation
  #show.plots: Logical argument for whether or not to show generated spatial variance plots
  #reproducible: Logical argument for whether or not to make generated dataset reproducible           
 
  total.sites <- grid.size*grid.size #Total area that could be sampled
  seed<-ifelse(repoducible==FALSE, NA, 42) #Reproducible spatial field?
  s <- simExpCorrRF(variance = variance.RF, theta = theta.RF, size = grid.size,
                    show.plots = show.plots, seed = seed) #Simulate spatial field
  
  #Reproducible data generation?
  if(is.na(seed)==FALSE){
    set.seed(seed)
  }else{
    
  }
  
  #site x survey x species arrays to hold observation/detection data (observation process)
  y.super <- y.obs <- y.all<- p <- array(NA, c(total.sites, nsurveys, nspecies)) 
  #Name dimensions
  dimnames(y.all) <- dimnames(y.obs) <- dimnames(p) <- list(paste("site", 
                                                                  1:total.sites, sep = ""), paste("survey", 1:nsurveys, 
                                                                                                  sep = ""), paste("spp", 1:nspecies, sep = "")) 
  
  #site x species arrays to hold abundance data (occupancy process)
  N.super <- lambda <- matrix(NA, total.sites, nspecies) 
  #Name dimensions
  dimnames(N.super) <- dimnames(lambda) <- list(paste("site", 
                                                      1:total.sites, sep = ""), paste("spp", 1:nspecies, sep = ""))
  detected.at.all <- rep(NA, nspecies) #Species every detected during surveys?
  
  #Simulate covariates and parameters
  site.cov <- rnorm(total.sites) #Site specific covariate
  survey.cov <- matrix(rnorm(total.sites * nsurveys), ncol = nsurveys) #Site/survey covariate
  mu.loglam <- log(mean.lambda) #Mean abundance of each species for log-link
  mu.lp <- ifelse(mean.p == "1", 500, qlogis(mean.p)) #Mean detection of each species for logit-link
  beta0 <- rnorm(nspecies, mu.loglam, sig.loglam) #Intercepts for spp. specific abundance
  beta1 <- rnorm(nspecies, mu.beta.loglam, sig.beta.loglam) #Coefficients for spp. specific abundance covariate
  alpha0 <- rnorm(nspecies, mu.lp, sig.lp) #Intercepts for spp. specific detection
  alpha1 <- rnorm(nspecies, mu.beta.lp, sig.beta.lp) #Intercepts for spp. specific detection covariate
  x.coord <- rep(1:grid.size, each = grid.size) #generate x.coords for grid
  y.coord <- rep(1:grid.size, times = grid.size) #generate y coords for grid
    
  #Simulating data from simulated parameters
  for (k in 1:nspecies) {
    lambda[, k] <- exp(beta0[k] + beta1[k] * site.cov + c(s$field)) #species x site abundance
    for (j in 1:nsurveys) {
      p[, j, k] <- plogis(alpha0[k] + alpha1[k] * survey.cov[, j]) #species x site x survey detection
    }
  }
  for (k in 1:nspecies) {
    N.super[, k] <- rpois(total.sites, lambda[, k]) #species x site latent abundance state
  }
  
  #"Observed" data if every grid cell had been sampled (super dataset)
  for (k in 1:nspecies) {
    for (i in 1:total.sites) {
      for (j in 1:nsurveys) {
        y.super[i, j, k] <- rbinom(1, N.super[i, k], p[i, j, k])
      }
    }
  }
  
  #Thinning to just sampled sites
  surveyed.sites <- sort(sample(1:total.sites, size = nsites)) #Only surveyed sites
  y.all<-y.super
  N.all<-N.super
  y.all[-surveyed.sites, , ] <- NA #Now all unsurveyed grid cells are NAs (true dataset)
  N.all[-surveyed.sites, ] <- NA #Now all unsurveyed grid cells are NAs (true dataset)
  tmp <- N.all
  tmp<-ifelse(is.na(tmp),0,tmp)
  tmp <- apply(tmp, 2, sum) #Counts of each species within sampled area
  occurring.in.sample <- as.numeric(tmp > 0) #Was species ever present at sampled sites?
  
  for (k in 1:nspecies) { #Was species ever detected at sampled sites?
    detected.at.all[k] <- if (any(y.all[, , k] > 0)) 
      TRUE
    else FALSE
  }
  
  #Desired datasets
  y.obs <- y.all[, , detected.at.all] #dataframe with only observed species at sampled sites (observation data)
  Ntotal.fs <- sum(occurring.in.sample) #Number of present species at sampled sites
  Ntotal.obs <- sum(detected.at.all) #Number of detected species at sampled sites
  
  #Save desired parameters
  return(list(nspecies = nspecies, nsites = nsites, nsurveys = nsurveys, 
              mean.lambda = mean.lambda, mu.loglam = mu.loglam, 
              sig.loglam = sig.loglam, mu.beta.loglam = mu.beta.loglam, 
              sig.beta.loglam = sig.beta.loglam, mean.p = mean.p, 
              mu.lp = mu.lp, sig.lp = sig.lp, mu.beta.lp = mu.beta.lp, 
              sig.beta.lp = sig.beta.lp, lambda = lambda, p = p, 
              N.super = N.super, y.super = y.super, N.all = N.all,
              y.all = y.all, y.obs = y.obs, Ntotal.fs = Ntotal.fs, 
              Ntotal.obs = Ntotal.obs, grid.size = 50, variance.RF = 1, 
              theta.RF = 10, field = s$field, site.cov = site.cov,
              survey.cov = survey.cov, x.coord = x.coord, y.coord = y.coord,
              surveyed.sites = surveyed.sites))
}

#########
###END###
#########



#######################################################################
###Function #2: Simulate data for a simple community abundance model###
#######################################################################


simCAM<-function (nsites = 30, nreps = 3, nspecies = 100, mean.psi = 0.25, sig.lpsi = 1, 
                  mu.beta.lpsi = 0, sig.beta.lpsi = 0, mean.lambda = 2, sig.loglam = 1, 
                  mu.beta.loglam = 1, sig.beta.loglam = 1, mean.p = 0.25, sig.lp = 1, 
                  mu.beta.lp = 0, sig.beta.lp = 0, show.plot = TRUE){
  
  y.all <- y.obs <- p <- array(NA, c(nsites, nreps, nspecies))
  dimnames(y.all) <- dimnames(y.obs) <- dimnames(p) <- list(paste("site", 
                                                                  1:nsites, sep = ""), paste("rep", 1:nreps, 
                                                                                             sep = ""), paste("sp", 1:nspecies, sep = ""))
  N <- lambda <- matrix(NA, nsites, nspecies)
  dimnames(N) <- dimnames(lambda) <- list(paste("site", 
                                                1:nsites, sep = ""), paste("sp", 1:nspecies, 
                                                                           sep = ""))
  detected.at.all <- rep(NA, nspecies)
  habitat <- sort(rnorm(nsites))
  wind <- matrix(rnorm(nsites * nreps), ncol = nreps)
  mu.loglam <- log(mean.lambda)
  mu.lp <- ifelse(mean.p == "1", 500, qlogis(mean.p))
  beta0 <- rnorm(nspecies, mu.loglam, sig.loglam)
  beta1 <- rnorm(nspecies, mu.beta.loglam, sig.beta.loglam)
  alpha0 <- rnorm(nspecies, mu.lp, sig.lp)
  alpha1 <- rnorm(nspecies, mu.beta.lp, sig.beta.lp)
  for (k in 1:nspecies) {
    lambda[, k] <- exp(beta0[k] + beta1[k] * habitat)
    for (j in 1:nreps) {
      p[, j, k] <- plogis(alpha0[k] + alpha1[k] * wind[ ,j])
    }
  }
  for (k in 1:nspecies) {
    N[, k] <- rpois(nsites, lambda[, k])
  }
  tmp <- apply(N, 2, sum)
  occurring.in.sample <- as.numeric(tmp > 0)
  for (k in 1:nspecies) {
    for (i in 1:nsites) {
      for (j in 1:nreps) {
        y.all[i, j, k] <- rbinom(1, N[i, k], p[i, j, k])
      }
    }
    detected.at.all[k] <- if (any(y.all[, , k] > 0)) 
      TRUE
    else FALSE
  }
  y.obs <- y.all[, , detected.at.all]
  detected.at.site <- apply(y.obs > 0, c(1, 3), any)
  ymax.obs <- apply(y.all, c(1, 3), max)
  Ntotal.fs <- sum(occurring.in.sample)
  Ntotal.obs <- sum(detected.at.all)
  
  return(list(nsites = nsites, nreps = nreps, 
              nspecies = nspecies, mean.lambda = mean.lambda, mu.loglam = mu.loglam, 
              sig.loglam = sig.loglam, mu.beta.loglam = mu.beta.loglam, 
              sig.beta.loglam = sig.beta.loglam, mean.p = mean.p, 
              mu.lp = mu.lp, sig.lp = sig.lp, mu.beta.lp = mu.beta.lp, 
              sig.beta.lp = sig.beta.lp, lambda = lambda, p = p, 
              N = N, y.all = y.all, y.obs = y.obs, ymax.obs = ymax.obs, 
              Ntotal.fs = Ntotal.fs, Ntotal.obs = Ntotal.obs))
}

#########
###END###
#########

###################################################################
###Function #3: Simulate data for a Time-Removal N-Mixture Model###
###################################################################

simNmixTR <- function(nsites = 100, nsurveys = 3, nwin = 5, mean.lambda = 5,
                      beta.llam = 1, mean.theta = 0.7, beta.ltheta = -0.5,
                      mean.p = 0.3, beta.lp = -0.3, reproducible = TRUE){
  ##Function Arguments:
  #nsites: Number of sampled sites
  #nsurveys: Number of replicated surveys conducted at each site
  #nwin: Number of time windows for time-removal sampling
  #mean.lambda: Average abundance of the species at any single site
  #beta.llam: Coefficient value for site level covariate of species abundance
  #mean.theta: Average availability probability of the species at any single site
  #beta.ltheta:Coefficient value for survey level covariate of species availability
  #mean.p: Average detection probability of the species at any single site
  #beta.lp: Coefficient value for survey level covariate of species detection
  #reproducible: Logical argument for whether or not to make generated dataset reproducible           
  
  ###Function Code###
  seed<-ifelse(reproducible==FALSE, NA, 42) #Reproducible?
  
  #Reproducible data generation?
  if(is.na(seed)==FALSE){
    set.seed(seed)
  }else{
    
  }
  
  #Arrays to hold simulated data
  y <- array(NA, c(nsites, nsurveys, nwin)) #observation array
  M <- rep(NA, nsites)                      #Local population size
  N <- array(NA, c(nsites, nsurveys))       #Individuals available for detection
  
  #Name dimensions
  dimnames(y) <- list(paste("site", 1:nsites, sep = ""), paste("survey", 1:nsurveys, 
                                                               sep = ""), paste("time", 1:nwin, sep = "")) 
  dimnames(N) <- dimnames(y)[1:2]
  
  #Simulate covariates and parameters
  site.cov.ab <- rnorm(nsites) #Site specific covariate for abundance
  survey.cov.p <- matrix(rnorm(nsites * nsurveys), ncol = nsurveys) #Site/survey covariate for det
  survey.cov.t <- matrix(rnorm(nsites * nsurveys), ncol = nsurveys) #Site/survey covariate for avail
  lambda <- rep(NA, nsites) #Site level variation in latent abundance due to covariates 
  theta <- matrix(rep(NA,nsites*nsurveys), ncol = nsurveys) #Site and survey level availability
  p <- matrix(rep(NA,nsites*nsurveys), ncol = nsurveys) #Site and survey level detection 
  beta0 <- log(mean.lambda) #Mean abundance for log-link
  alpha0 <- ifelse(mean.p == "1", 500, qlogis(mean.p)) #Mean detection for logit-link
  gamma0 <- ifelse(mean.theta == "1", 500, qlogis(mean.theta))#Mean availability for logit-link
  beta1 <- beta.llam    #Coefficient for abundance covariate
  alpha1 <- beta.lp     #Coefficient for detection covariate
  gamma1 <- beta.ltheta # Coefficient for availability covariate
  
  
  #Simulating abundance and detection data from simulated parameters
  for(i in 1:nsites){
    lambda[i] <- exp(beta0 + beta1 * site.cov.ab[i]) 
    M[i] <- rpois(1, lambda[i]) #latent abundance at each site
    for(j in 1:nsurveys){
      p[i, j] <- plogis(alpha0 + alpha1 * survey.cov.p[i, j]) #site/survey det prob
    }
  }
  
  #Simulating availability data from simulated parameters
  for(i in 1:nsites){
    for(j in 1:nsurveys){
      theta[i, j] <- plogis(gamma0 + gamma1 * survey.cov.t[i, j]) #site/survey avail prob
    }
  }
  mu.theta<-1-(1-theta)^(1/nwin) # mean availability prob for site/survey/time window
  
  #Saving observations as an indexible dataframe
  M.obs<-cbind(seq(1:nsites), M)
  M.obs<-as.data.frame(M.obs)
  colnames(M.obs)<-c("site","count")
  
  #Sampling from closed population using simulated parameters
  M.obs.tmp<-matrix(rep(M.obs$count,nsurveys),ncol = nsurveys)
  t<-array(data = NA, dim = c(nrow(M.obs.tmp),nsurveys,nwin))
  for(i in 1:nrow(M.obs.tmp)){
    for(j in 1:nsurveys){
      for(k in 1:nwin){
        t[i,j,k]<-ifelse(M.obs.tmp[i,j]==0,0,rbinom(1,M.obs.tmp[i,j],mu.theta[i,j])) # capture some in time period 
        M.obs.tmp[i,j]<- M.obs.tmp[i,j] - t[i,j,k] #remove them
      }
    }
  }
  
  #Total counts of available individuals during each site x survey
  N.avail <- apply(t,c(1,2),sum)
  
  #Observation array given available individuals and detection probability
  n.obs.tmp <- array(data = NA, dim = dim(t))
  for(i in 1:nrow(N.avail)){
    for(j in 1:nsurveys){
      for (k in 1:nwin){
        n.obs.tmp[i,j,k] <- rbinom(1, t[i,j,k], p[i,j])
      }
    }
  }
  
  #Total counts of observed individuals during each site x survey
  n.obs <- apply(n.obs.tmp,c(1,2),sum)
  
  #Metadata for site name, survey number, and number of observations
  obs.md<-array(data = NA, dim = c(nsites*nsurveys,3))
  obs.md[,1]<-seq(1:nsites)
  obs.md[,2]<-rep(1:nsurveys, each = nsites)
  obs.md<-as.data.frame(obs.md)
  colnames(obs.md)<-c("site", "survey", "count")
  for(i in 1:nsites){
    for(j in 1:nsurveys){
      obs.md[obs.md$site==i & obs.md$survey==j,3]<-n.obs[i,j]
    }
  }
  
  #Removing non-observations from dataframes
  obs.md<-obs.md[obs.md[,3]!=0,]
  
  #Number of unique detections across all sites and surveys (NOT # of detected individuals!)
  n.dets <- nrow(obs.md)
  
  #Detection x time interval matrix generation for multinomial distribution 
  obs.tmp<-n.obs.tmp[,1,]
  for(i in 2:nsurveys){
    obs.tmp<-rbind(obs.tmp,n.obs.tmp[,i,])
  }
  tint<-obs.tmp[which(apply(obs.tmp,1,sum)!=0),]
  
  #Check to make sure everything is still indexed correctly
  obs.md[,3]==apply(tint,1,sum)
  
  #Saving data to output object
  return(list(nsites = nsites, nsurveys = nsurveys, K = nwin, tint = tint, 
              n.obs = n.obs, n.dets = n.dets, obs.md = obs.md, ab.cov = site.cov.ab, 
              det.cov = survey.cov.p, avail.cov = survey.cov.t, mean.lambda = mean.lambda,
              beta.llam = beta.llam, mean.theta = mean.theta, beta.ltheta = beta.ltheta, 
              mean.p = mean.p, beta.lp = beta.lp, lambda = lambda, theta = theta, p = p,
              beta0 = beta0, alpha0 = alpha0, gamma0 = gamma0, beta1 = beta1, alpha1 = alpha1,
              gamma1 = gamma1, mu.theta = mu.theta, N.avail = N.avail, M.obs=M.obs))
}

#########
###END###
#########

#################################################################################
###Function #4: simulate data for a CAM with Distance sampling and Time of Det###
#################################################################################

simCAMDSTD<-function (nspecies = 50, nsites = 50, nreps = 3, mu.lambda = 2, 
                         sig.loglam = 0.5, mu.beta.loglam = 0.3, sig.beta.loglam = 0.5, 
                         mu.p = 0.4, sig.lp = 0.25, mu.beta.lp = -0.3, sig.beta.lp = 0.5, 
                         mu.avail = 0.8, sig.lavail = 0.8, mu.beta.lavail = -0.3, 
                         sig.beta.lavail = 0.5, mu.pres = 0.7, sig.lpres = 0.8,
                         nSubSamp = 5, B = 100, delta = 10, omega = 0.75,
                         repoducible = FALSE){

seed<-ifelse(repoducible==FALSE, NA, 42) #Reproducible data generation?

if(is.na(seed)==FALSE){
  set.seed(seed)
}else{
  
}

###Function Arguments###

#nspecies: metacommunity species richness 
#nsites: number of sites
#nsurveys: number of repeated surveys

#mu.lambda: Mean site abundance of all species
#sig.loglam: Standard deviation among species specific mean abundances (given on log scale)
#mu.beta.loglam: Mean of coefficient value for abundance covariate across species (given on log scale)
#sig.beta.loglam: Standard deviation in coefficient values for abundance covariate across species (given on log scale)
#mu.p: Mean detection probability of all species
#sig.lp: Standard deviation among species specific detection probabilities (given on log scale)
#mu.beta.lp: Mean of coefficient value for detection covariate across species (given on log scale)
#sig.beta.lp: Standard deviation in coefficient values for detection covariate across species (given on log scale)
#mu.avail: Mean availability probability of all species 
#sig.lavail: Standard deviation among species specific availability probabilities (given on logit scale)
#mu.beta.lavail: Mean of coefficient value for availability covariate across species (given on logit scale)
#sig.beta.lavail: Standard deviation in coefficient values for availability covariate across species (given on logit scale)
#mu.pres: Mean presence probability of all species 
#sig.lpres: Standard deviation among species specific presence probabilities (given on logit scale)

#nSubSamp: time-of-detection sub-intervals per occasion
#B: maximum distance for distance sampling
#delta: width of distance sampling bins
#omega: likelihood of species presence in metacommunity 

detected.at.all <- rep(NA, nspecies)
habitat <- sort(rnorm(nsites))
wind <- matrix(rnorm(nsites * nreps), ncol = nreps)
time <- matrix(rnorm(nsites * nreps), ncol = nreps)
mu.loglam <- log(mu.lambda)

###Calculating mu.lp with given mu.p value
sigmas<-rep(NA,200)
for(i in 1:200){
  sigma<-i
  sigma.1<-integrate(function(d){(2/B^2)*d*exp((-d*d)/(2*sigma^2))},lower = 0,upper = B)
  sigmas[i]<-sigma.1$value
}
mu.lp <- log(which(abs(sigmas-mu.p)==min(abs(sigmas-mu.p))))
###

mu.lavail <- ifelse(mu.avail == "1", 500, qlogis(mu.avail))
mu.lpres <- ifelse(mu.pres == "1", 500, qlogis(mu.pres))
rho0 <- rnorm(nspecies, mu.lpres, sig.lpres)
beta0 <- rnorm(nspecies, mu.loglam, sig.loglam)
beta1 <- rnorm(nspecies, mu.beta.loglam, sig.beta.loglam)
alpha0 <- rnorm(nspecies, mu.lp, sig.lp)
alpha1 <- rnorm(nspecies, mu.beta.lp, sig.beta.lp)
gamma0 <- rnorm(nspecies, mu.lavail, sig.lavail)
gamma1 <- rnorm(nspecies, mu.beta.lavail, sig.beta.lavail)

w<-pPres<-sigma<-rep(NA,nspecies) #create vectors of length nspecies
lambda<-array(NA, dim = c(nspecies,nsites))     #lambda array for species x site 
pAvail<-pAvailPrime<-array(NA, dim = c(nspecies,nsites,nreps)) #availability arrays for species x site x survey
for(k in 1:nspecies){                           #Parameters for each species (species heterogeneity)
  w[k] <- rbinom(1,1,omega)                     #Presence/absence of each species in metacommunity 
  pPres[k] <- plogis(rho0[k])                   #Species-specific probability of presence within sampled area
  #sigma[k] <- runif(1,20,150)                   #scale parameter for half-normal distance sampling (THIS IS WHERE COVARIATES CAN ENTER!)
  for(i in 1:nsites){
    lambda[k,i] <- exp(beta0[k] + beta1[k] * habitat[i])      #heterogeneity in supercommunity abundances across species
    for(j in 1:nreps){
      pAvail[k,i,j] <- plogis(gamma0[k] + gamma1[k] * time[i,j]) #probability of ever being available|present 
      pAvailPrime[k,i,j] <- 1-(1-pAvail[k,i,j])^(1/nSubSamp)     #probability of availability in a subinterval
    }
  }
}

##Distance sampling values
db = seq(0,B, by=delta)                   #dist bin breaks
midpt = seq(.5*delta,B, by=delta)         #midpoints for distance bins
nD = length(db)-1                         #number of distance bins
pix = (2*midpt*delta )/(B*B)              #relative area of each bin
Area = (pi*B*B)/10000                     #area of a single plot in ha
totArea = Area*nsites                     #total area surveyed
sites = 1:nsites                          #site labels

#placeholders
dims<-c(nspecies,nsites,nreps)
Npres<-Npres.true<-sigmaObs<-array(data=NA, dim = dims)
sigma.vec<-obsP<-M.p<-nind<-rep(NA,nspecies)	
M.super<-array(data = NA, dim = c(nspecies,nsites))
PresSite<-surveyPres<-u1Pres<-u2Pres<-dPres<-pInd<-obs<-obsSite<-caphist.tmp<-dObs<-u1Obs<-u2Obs<-dClass<-survey<-counts<-as.list(rep(NA,nspecies))

##Number of Individuals of species k present at site i during survey j

for (k in 1:nspecies){
  for (i in 1:nsites){
    M.super[k,i]<- rpois(1, lambda[k,i]*w[k]) #latent abundance of each species at each site
    for (j in 1:nreps){ #Yes necessary
      Npres[k,i,j]<-rbinom(1, M.super[k,i], pPres[k]) #number of present individuals of each species x site x rep
    }#j
  }#i
}#k


##assign location of individuals of each species present at each site during each survey
angle<-r2<-r<-u1<-u2<-d<-r<-y<-p<-array(NA, c(max(Npres),nspecies,nsites,nreps)) #placeholders
for (k in 1:nspecies){
  for (i in 1:nsites){
    for (j in 1:nreps){
      angle[,k,i,j] <- c(runif(Npres[k,i,j], 0, 2*pi), rep(NA, max(Npres)-Npres[k,i,j]))
      r2[,k,i,j] <- c(runif(Npres[k,i,j], 0, 1), rep(NA, max(Npres)-Npres[k,i,j]))
      r[,k,i,j] <- B * sqrt(r2[,k,i,j])
      u1[,k,i,j] <- r[,k,i,j] * cos(angle[,k,i,j]) + B
      u2[,k,i,j] <- r[,k,i,j] * sin(angle[,k,i,j]) + B
      d[,k,i,j] <- sqrt((u1[,k,i,j] - B)^2 + (u2[,k,i,j] - B)^2)
      Npres.true[k,i,j] <- sum(d[,k,i,j] <= B, na.rm=T)
      sigmaObs[k,i,j] <- exp(alpha0[k] + alpha1[k]*wind[i,j])	
      p[,k,i,j] <- ifelse(d[,k,i,j] < (B), 1, 0) * exp(-d[,k,i,j] * d[,k,i,j]/(2 * (sigmaObs[k,i,j]^2))) 
    }#j 
  }#i
}#k


NpresTot<-apply(Npres.true,c(1,3),sum)      #total number of individuals of each species present during each survey
NpresAll<-apply(NpresTot,1,sum)             #NpresAll[k] should = length(pInd[k])

for (k in 1:nspecies){
  PresSite[[k]]<- as.vector(unlist(apply(Npres.true[k,,],2, function(x) rep(1:nsites, x)))) #site of present individuals
  surveyPres[[k]]<-rep(1:nreps, NpresTot[k,])                            #Survey of present individual
  u1Pres[[k]]<-as.vector(u2[,k,,])[(as.vector(!is.na(u2[,k,,])))]       #x coordinates of present individuals
  u2Pres[[k]]<-as.vector(u2[,k,,])[(as.vector(!is.na(u2[,k,,])))]       #y coordinates of present individuals 
  dPres[[k]]<-as.vector(d[,k,,])[(as.vector(!is.na(d[,k,,])))]          #distance of present individuals
  pInd[[k]]<-as.vector(p[,k,,])[(as.vector(!is.na(p[,k,,])))]           #detection probability of each individual of each species given distance 
}#k


##Availablity if Present (CMR data)
Avail<-newpMat<-array(NA, c(nspecies, max(NpresAll), nSubSamp))  #sub-interval (minute) availability

##create true availability histories
for(k in 1:nspecies){
  if(NpresAll[k]==0){
    for(n in 1:max(NpresAll)){
      Avail[k,n,]<-rep(NA,nSubSamp)
    }#n
  }else{
    for(n in 1:NpresAll[k]){
      Avail[k,n,]<-rbinom(nSubSamp,1,pAvailPrime[k,PresSite[[k]][n],surveyPres[[k]][n]]) #available or not (1/0) during a sub-interval
    }#n
  }
}#k

IndAvail <- apply(Avail,c(1,2),max)            #available (y/n)
IndAvail<-ifelse(is.na(IndAvail),0,IndAvail)   #replacing NAs with 0s
NavailTot<-array(NA, c(nspecies, nreps))

for(k in 1:nspecies){
  if(sum(IndAvail[k,])>0){
    NavailTot[k,]<- tabulate(surveyPres[[k]][which(IndAvail[k,]==1)], nreps) #Total available individuals per survey 
  }else{
    NavailTot[k,]<-0
  } 
}#k

##Detection if available
newp<-y<-matrix(NA, nrow = nrow(IndAvail), ncol = ncol(IndAvail))
for(k in 1:nspecies){
  pInd[[k]]<-append(pInd[[k]], rep(0, ncol(IndAvail) - length(pInd[[k]])))  #Make dims equal
  newp[k,] <- pInd[[k]]*IndAvail[k,]   #new detection probability if individual is available
  y[k,]<-rbinom(IndAvail[k,], 1, newp[k,])                                     #observed y/n
  obs[[k]]<-which(y[k,]>0)                                             #observed individuals
  nind[k]<-length(obs[[k]])                                           #number of individuals observed per species
  caphist.tmp[[k]]<-as.array(Avail[k,obs[[k]],], nrow = length(obs[[k]])) #capture history for observed individuals
  obsSite[[k]]<-PresSite[[k]][obs[[k]]]                                    #observation site (vector of site numbers)
  u1Obs[[k]] <- as.vector(u1[,k,,])[(as.vector(!is.na(u1[,k,,])))][obs[[k]]]   #x coord of observed individuals 
  u2Obs[[k]] <- as.vector(u2[,k,,])[(as.vector(!is.na(u2[,k,,])))][obs[[k]]]   #y coord of observed individuals 
  dObs[[k]] <- as.vector(d[,k,,])[(as.vector(!is.na(d[,k,,])))][obs[[k]]]      #distance of observed individuals
  dClass[[k]] <- c(dObs[[k]]%/%delta + 1)                                      #distance bin	
  survey[[k]] <- surveyPres[[k]][obs[[k]]]                                     #survey (vector of survey numbers)
}#k 

#5-d capture history of species x site x rep x individual x time interval
caphist<-array(NA, c(nspecies, nsites, nreps, max(unlist(lapply(caphist.tmp, function(x) nrow(x)))), nSubSamp)) #making caphist a 5d array
for(k in 1:nspecies){
  for(i in 1:length(obsSite[[k]])){
    if(sum(caphist.tmp[[k]])>0 && length(dim(caphist.tmp[[k]]))>1){
      caphist[k,obsSite[[k]][i],survey[[k]][i],i,1:nSubSamp]<-caphist.tmp[[k]][i,]
    }else if(sum(caphist.tmp[[k]])>0){
      caphist[k,obsSite[[k]][i],survey[[k]][i],i,1:nSubSamp]<-caphist.tmp[[k]]
    }else{
      caphist[k,,,,]<-NA
    }
  }#i
}#k

#4-d array of observations species x site x survey x distance bin 
dClassN<-array(0,c(nspecies, nsites, nreps, nD)) 
for(k in 1:nspecies){ #indexing by site, survey, distance bin
  if(length(obsSite[[k]])>0){
    counts[[k]]<- table(obsSite[[k]], dClass[[k]], survey[[k]])
  }else{
    counts[[k]]<-array(0,dim=c(1,nD,nreps))
  } 
}#k

for(k in 1:nspecies){ #filling 4d array with counts 
  for(i in 1:dim(counts[[k]])[3]){
    dClassN[k,as.numeric(rownames(counts[[k]])),i,as.numeric(colnames(counts[[k]]))]<-counts[[k]][,,i]
  }#j 
}#k

nObs <- apply(dClassN, c(1,2,3), sum) #observations per species, site, and survey
nobsSurvey<-apply(nObs,c(1,3),sum)    #observations per species per survey
numObs<-length(which(nind!=0))        #Number of observed species 
Ntot<-sum(w[])                        #Latent species richness
Mtot<-rep(NA,nspecies)
for(k in 1:nspecies){
  Mtot[k]<-sum(M.super[k,])           #Latent abundances of each species
}

##Store and check data for JAGS
str(return(list(nspec=nspecies, nsites=nsites, nrep=nreps, nD=nD, midpt=midpt, delta=delta, B=B, nObs= nObs,
                  nind=nind, numObs=numObs, nSubSamp=nSubSamp, dClassN=dClassN, caphist=caphist, totArea=totArea,
                  habitat=habitat, wind=wind, time=time, Mtot=Mtot, Ntot=Ntot, mu.lp=mu.lp)))

}#END of Function



########################################################################
###Function #5: Simulate data for a dynamic community abundance model###
########################################################################

simDCAM<-function (nsites = 20, nsurveys = 3, nspecies = 30, nyears = 3,
                   mean.lambda = 4, sd.lambda = 1, mu.beta.lambda = 1, 
                   sd.beta.lambda = 0.5, mean.phi = 0.7, sd.phi = 0.3,
                   mu.beta.phi = 0.3, sd.beta.phi = 0.5, mean.gamma = 0.4, sd.gamma = 0.3,
                   mu.beta.gamma = 0.3, sd.beta.gamma = 0.5, mean.p = 0.5, sd.p = 0.2,
                   mu.beta.p = -0.4, sd.beta.p = 0.5, show.plot = TRUE){
  
  spec <- 1:nspecies
  site <- 1:nsites
  year <- 1:nyears
  survey <- 1:nsurveys
  lambda <- phi <- gamma <- array(NA, dim = c(nsites, nspecies), 
                            dimnames = list(paste("Site", site ,sep = ""),
                            paste("Spec", spec ,sep ="")))
  
  N <- S <- R <- array(NA, dim = c(nsites, nspecies, nyears), 
                 dimnames = list(paste("Site", site, sep = ""), 
                 paste("Spec", spec, sep = ""), paste("Year", year, sep = "")))
  
  y.all <- y.obs <- p <- array(NA, dim = c(nsites, nspecies, nsurveys, nyears),
                               dimnames = list(paste("Site", 
                                          site, sep = ""), paste("Spec", spec, 
                                          sep = ""), paste("Survey", survey, sep = ""),
                                          paste("Year", year, sep = "")))
  
  detected.at.all <- rep(NA, nspecies)
  habitat <- sort(rnorm(nsites))
  wind <- array(data = rnorm(nsites * nsurveys * nyears), dim = c(nsites, nsurveys, nyears))
  
  mu.loglam <- log(mean.lambda)
  alpha.lambda <- rnorm(nspecies, mu.loglam, sd.lambda)
  beta.lambda <- rnorm(nspecies, mu.beta.lambda, sd.beta.lambda)
  
  for (k in 1:nspecies) {
    lambda[, k] <- exp(alpha.lambda[k] + beta.lambda[k] * habitat) #abundance 
  }
  
  mu.lphi <- ifelse(mean.phi == "1", 500, qlogis(mean.phi))
  alpha.phi <- rnorm(nspecies, mu.lphi, sd.phi)
  beta.phi <- rnorm(nspecies, mu.beta.phi, sd.beta.phi)
  
  for (k in 1:nspecies) {
    phi[, k] <- plogis(alpha.phi[k] + beta.phi[k] * habitat) #survival probability
  }
  
  mu.loggam <- log(mean.gamma)
  alpha.gamma <- rnorm(nspecies, mu.loggam, sd.gamma)
  beta.gamma <- rnorm(nspecies, mu.beta.gamma, sd.beta.gamma)
  
  for (k in 1:nspecies) {
    gamma[, k] <- exp(alpha.gamma[k] + beta.gamma[k] * habitat) #recruitment probability
  }
  
  mu.lp <- ifelse(mean.p == "1", 500, qlogis(mean.p))
  alpha.p <- rnorm(nspecies, mu.lp, sd.p)
  beta.p <- rnorm(nspecies, mu.beta.p, sd.beta.p)
  
  for(k in 1:nspecies){
    for(j in 1:nsurveys){
      for(t in 1:nyears){
        p[, k, j, t] <- plogis(alpha.p[k] + beta.p[k] * wind[,j,t]) #detection probability
      }
    }
  }
  
  for (k in 1:nspecies) {
    N[, k, 1] <- rpois(nsites, lambda[, k]) #Abundance year 1
  }
  
  for(i in 1:nsites){
    for(k in 1:nspecies){
      for(t in 1:(nyears-1)){
        S[i,k,t+1] <- rbinom(1, N[i,k,t], phi[i,k]) #survival of spp. k at site i (constant across years)
        R[i,k,t+1] <- rpois(1,N[i,k,t]*gamma[i,k]) #relative recruitment of spp. k at site i (constant across years)
        N[i,k,t+1] <- S[i,k,t+1] + R[i,k,t+1] #latent abundance in year t+1
      }
    }
  }
 
  for (i in 1:nsites) {
    for (k in 1:nspecies) {
      for (j in 1:nsurveys) {
        for (t in 1:nyears){
          y.all[i, k, j, t] <- rbinom(1, N[i, k, t], p[i, k, j, t])
        }
      }
    }
  }
  
  for (k in 1:nspecies) {
  detected.at.all[k] <- if (any(y.all[ , k, , ] > 0)) 
    TRUE
  else FALSE
  }
  
  y.obs <- y.all[ , detected.at.all, , ]
  Ntotal.obs <- sum(detected.at.all)
  
  return(list(nsites = nsites, nsurveys = nsurveys, nspecies = nspecies, nyears=nyears,
              habitat = habitat, wind = wind, mean.lambda = mean.lambda, mu.loglam = mu.loglam, 
              sd.lambda = sd.lambda, mu.beta.lambda = mu.beta.lambda, 
              sd.beta.lambda = sd.beta.lambda, mean.phi = mean.phi,
              sd.phi = sd.phi, mu.beta.phi = mu.beta.phi, sd.beta.phi = sd.beta.phi,
              mean.gamma = mean.gamma, sd.gamma = sd.gamma,
              mu.beta.gamma = mu.beta.gamma, sd.beta.gamma = sd.beta.gamma, 
              mean.p = mean.p, sd.p = sd.p, mu.beta.p = mu.beta.p, sd.beta.p = sd.beta.p, 
              lambda = lambda, phi = phi, gamma = gamma, p = p, 
              N = N, y.all = y.all, y.obs = y.obs, Ntotal.obs = Ntotal.obs))
}

#########
###END###
#########


################################################################################
###Function #6: Simulate data for a spatial dynamic community abundance model###
################################################################################

simSpDCAM<-function (nsites = 25, nsurveys = 3, nspecies = 30, nyears = 3,
                   mean.lambda = 2, sd.lambda = 1, mu.beta.lambda = 0, 
                   sd.beta.lambda = 0.5, mean.phi = 0.7, sd.phi = 0.3,
                   mu.beta.phi = 0, sd.beta.phi = 0.5, mean.gamma = 0.3, sd.gamma = 0.3,
                   mu.beta.gamma = 0, sd.beta.gamma = 0.5, mean.kappa = 0.4, sd.kappa = 0.3, 
                   mu.beta.kappa = 0, sd.beta.kappa = 0.5, mean.theta = 1, sd.theta = 0.3,
                   mean.p = 0.5, sd.p = 0.3, mu.beta.p = 0, sd.beta.p = 0.5){
  stopifnot(sqrt(nsites)%%1==0)
  spec <- 1:nspecies
  site <- 1:nsites
  year <- 1:nyears
  survey <- 1:nsurveys
  
  lat<- rep(1:sqrt(nsites), times = sqrt(nsites))
  long<- rep(1:sqrt(nsites), each = sqrt(nsites))
  dist <- as.matrix(dist(cbind(lat, long), method='euclidean', diag=T, upper=T)) # distance matrix
  w_ori <- array(NA, dim = c(nsites, nsites, nspecies), 
                 dimnames = list(paste("Site", site ,sep = ""),
                                 paste("Site", site ,sep = ""),
                                 paste("Spec", spec ,sep ="")))
    
  lambda <- phi <- gamma <- kappa <- array(NA, dim = c(nsites, nspecies), 
                                  dimnames = list(paste("Site", site ,sep = ""),
                                                  paste("Spec", spec ,sep ="")))
  
  N <- S <- R <- E <- I <- array(NA, dim = c(nsites, nspecies, nyears), 
                            dimnames = list(paste("Site", site, sep = ""), 
                                       paste("Spec", spec, sep = ""), paste("Year", year, sep = "")))
  M <- array(NA, dim = c(nsites, nsites, nspecies, nyears), 
             dimnames = list(paste("Site", site, sep = ""), paste("Site", site, sep = ""),  
                             paste("Spec", spec, sep = ""), paste("Year", year, sep = "")))
  y.all <- y.obs <- p <- array(NA, dim = c(nsites, nspecies, nsurveys, nyears),
                               dimnames = list(paste("Site", 
                                                     site, sep = ""), paste("Spec", spec, 
                                                                            sep = ""), paste("Survey", survey, sep = ""),
                                               paste("Year", year, sep = "")))
  
  detected.at.all <- rep(NA, nspecies)
  habitat <- sort(rnorm(nsites))
  wind <- array(data = rnorm(nsites * nsurveys * nyears), dim = c(nsites, nsurveys, nyears))
  
  mu.loglam <- log(mean.lambda)
  alpha.lambda <- rnorm(nspecies, mu.loglam, sd.lambda)
  beta.lambda <- rnorm(nspecies, mu.beta.lambda, sd.beta.lambda)
  
  for (k in 1:nspecies) {
    lambda[, k] <- exp(alpha.lambda[k] + beta.lambda[k] * habitat) #abundance 
  }
  
  mu.lphi <- ifelse(mean.phi == "1", 500, qlogis(mean.phi))
  alpha.phi <- rnorm(nspecies, mu.lphi, sd.phi)
  beta.phi <- rnorm(nspecies, mu.beta.phi, sd.beta.phi)
  
  for (k in 1:nspecies) {
    phi[, k] <- plogis(alpha.phi[k] + beta.phi[k] * habitat) #survival probability
  }
  
  mu.loggam <- log(mean.gamma)
  alpha.gamma <- rnorm(nspecies, mu.loggam, sd.gamma)
  beta.gamma <- rnorm(nspecies, mu.beta.gamma, sd.beta.gamma)
  
  for (k in 1:nspecies) {
    gamma[, k] <- exp(alpha.gamma[k] + beta.gamma[k] * habitat) #recruitment probability
  }
  
  mu.lkappa <- ifelse(mean.kappa == "1", 500, qlogis(mean.kappa))
  alpha.kappa <- rnorm(nspecies, mu.lkappa, sd.kappa)
  beta.kappa <- rnorm(nspecies, mu.beta.kappa, sd.beta.kappa)
  
  for (k in 1:nspecies) {
    kappa[, k] <- plogis(alpha.kappa[k] + beta.kappa[k] * habitat) #emigration probability
  }    
    
  mu.theta <- log(mean.theta)
  alpha.theta <- rnorm(nspecies, mu.theta, sd.theta)
  theta <- exp(alpha.theta)                                        #distance decay
  for (k in 1:nspecies) {
    w_ori[,,k] <- exp(-1 * theta[k] * dist)
    diag(w_ori[,,k]) <- 0
  }
  
  w <- w_ori / rowSums(w_ori)
  
  mu.lp <- ifelse(mean.p == "1", 500, qlogis(mean.p))
  alpha.p <- rnorm(nspecies, mu.lp, sd.p)
  beta.p <- rnorm(nspecies, mu.beta.p, sd.beta.p)
  
  for(k in 1:nspecies){
    for(j in 1:nsurveys){
      for(t in 1:nyears){
        p[, k, j, t] <- plogis(alpha.p[k] + beta.p[k] * wind[,j,t]) #detection probability
      }
    }
  }
  
  for (k in 1:nspecies) {
    N[, k, 1] <- rpois(nsites, lambda[, k]) #Abundance year 1
  }
  
  for(k in 1:nspecies){
    for(t in 1:(nyears-1)){
      S[,k,t] <- rbinom(nsites, N[,k,t], phi[,k]) #survival of spp. k at site i (constant across years)
      R[,k,t] <- rpois(nsites, N[,k,t]*gamma[,k]) #relative recruitment of spp. k at site i (constant across years)
      E[,k,t] <- rbinom(nsites, S[,k,t], kappa[,k]) #emigration of spp. k from site i (constant across years, conditional on survival)
      for(j in 1:nsites){
        M[,j,k,t] <- E[,k,t] * w[,j,k] #Distance decay function for emigrants
      }
      I[,k,t] <- rpois(nsites,sum(M[,,k,t])) #immigration of spp. k into site i (constant across years)
      N[,k,t+1] <- S[,k,t] - E[,k,t] + R[,k,t] + I[,k,t] #latent abundance in year t+1
    }
  }
  for (i in 1:nsites) {
    for (k in 1:nspecies) {
      for (j in 1:nsurveys) {
        for (t in 1:nyears){
          y.all[i, k, j, t] <- rbinom(1, N[i, k, t], p[i, k, j, t])
        }
      }
    }
  }
  
  for (k in 1:nspecies) {
    detected.at.all[k] <- if (any(y.all[ , k, , ] > 0)) 
      TRUE
    else FALSE
  }
  
  y.obs <- y.all[ , detected.at.all, , ]
  Ntotal.obs <- sum(detected.at.all)
  
  return(list(nsites = nsites, nsurveys = nsurveys, nspecies = nspecies, nyears=nyears,
              habitat = habitat, wind = wind, dist = dist,
              mean.lambda = mean.lambda, mu.loglam = mu.loglam, 
              sd.lambda = sd.lambda, mu.beta.lambda = mu.beta.lambda, 
              sd.beta.lambda = sd.beta.lambda, mean.phi = mean.phi,
              sd.phi = sd.phi, mu.beta.phi = mu.beta.phi, sd.beta.phi = sd.beta.phi,
              mean.gamma = mean.gamma, sd.gamma = sd.gamma,
              mu.beta.gamma = mu.beta.gamma, sd.beta.gamma = sd.beta.gamma, 
              mean.kappa = mean.kappa, sd.kappa = sd.kappa,
              mu.beta.kappa = mu.beta.kappa, sd.beta.kappa = sd.beta.kappa,
              mean.theta = mean.theta, sd.theta = sd.theta,
              mean.p = mean.p, sd.p = sd.p, mu.beta.p = mu.beta.p, sd.beta.p = sd.beta.p, 
              lambda = lambda, phi = phi, gamma = gamma, kappa = kappa, theta = theta, p = p, 
              N = N, y.all = y.all, y.obs = y.obs, Ntotal.obs = Ntotal.obs))
}

#########
###END###
#########


##########################################################################