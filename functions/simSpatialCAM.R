##############################################################
###Simulation Function for Spatial Community Abundance Data###
##############################################################

#Franco Gigliotti 
#03/03/2020

#################
###Description###
#################

##This function file simulates a dataset to be used in a spatial community abundance
##modeling framework. Code for this function relies heavily on two functions in the 
##R package 'AHMBook' by Kery and Royle. To observe source code for these two functions,
##uncomment and run the following two lines of code.

#print(getAnywhere('simNmixSpatial')) #outputs source code for S.S.A.M
#print(getAnywhere('simComm')) #outputs source code for C.A.M

##############
###Function###
##############

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