#################
###Description###
#################

###This code simulates data to test out a Dynamic Spatial Occupancy Model
###For potential later use in SHARP survey analyses (though note single-species!)
###Code is based off of code presented in Ch. 9 of AHM Vol 2
###Code available at: https://github.com/mikemeredith/AHM_code/tree/master/AHM2_ch09

###Load Required Packages###
library(AHMbook)
library(R2OpenBUGS) #Need to use BUGS rather than JAGS for spatial models!
#Unfortunately, neither OpenBUGS nor WinBUGS are provided as modules on the 
#HPC. Maybe can run in another program like Nimble? Or just run locally and
#accept the fact that the spatial models will bog down my laptop/UConn Desktop

#Identify path to OpenBUGS software and set wd
bugs.path<-"C:/Program Files (x86)/OpenBUGS/OpenBUGS323"
setwd(choose.dir())
###Simulate Required Dataset###

##Getting a feel for the spatial autocorrelation (play with theta)
simExpCorrRF(variance = 1, theta = 50, size = 200)
simDynoccSpatial()

##Parameters selected based on thinned SHARP long-term dataset



####Ignore everything below here, this is not the code you are looking for...
# 9.2 Data simulation for spatial N-mixture and occupancy models
# ==============================================================

library(AHMbook)
? simExpCorrRF
str(dat <- simExpCorrRF(variance = 1, theta = 1, size = 50, seed = 1))

str(tmp <- simExpCorrRF(theta = 0.0001, size = 200))
str(tmp <- simExpCorrRF(theta = 1, size = 200))
str(tmp <- simExpCorrRF(theta = 5, size = 200))
str(tmp <- simExpCorrRF(theta = 10, size = 200))
str(tmp <- simExpCorrRF(theta = 100, size = 200))
str(tmp <- simExpCorrRF(theta = 10000, size = 200)) # cool patterns !

data(BerneseOberland)
head(bo <- BerneseOberland)
str(bo)

# ~~~~ extra code for figure 9.1 ~~~~~~~~~~
library(raster)
op <- par(mfrow = c(1, 2), mar = c(3,3,3,5), cex.main = 2, cex.axis = 1.5)
r1 <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = bo$elevation))
plot(r1, col = topo.colors(100), axes = FALSE, box = FALSE,
     main = "Elevation (m a.s.l.)", zlim = c(0, 3000))
r1 <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = bo$forest))
plot(r1, col = topo.colors(100), axes = FALSE, box = FALSE,
     main = "Forest cover (%)", zlim = c(0, 100))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Create Gaussian random field
set.seed(10) # Fig. 9.2
s <- simExpCorrRF(theta = 10, size = 50)

# Choose sample sizes: 2500 sites and 3 surveys
nsites <- 2500 # Number of sites (corresponding to our 50 by 50 grid)
nreps <- 3     # Number of replicate observations

# Scale the real Bernese Oberland covariates
elev <- standardize(bo$elevation)
forest <- standardize(bo$forest)

# Ecological process
beta0 <- 2           # Abundance model: intercept
beta1 <- 2           # Linear effect of elevation: positive
beta2 <- -2          # Quadratic effect of elevation: negative
loglam0 <- beta0 + beta1 * elev + beta2 * elev^2
loglam <- beta0 + beta1 * elev + beta2 * elev^2 + c(s$field)
lam0 <- exp(loglam0) # without spatial autocorrelation
lam <- exp(loglam)   # with spatial autocorrelation

# ~~~~ extra code for figure 9.3 ~~~~
# Plot expected counts (lambda) as a function of covariates only,
#   i.e., excluding spatial field, and including spatial field (Fig. 20-4)
op <- par(mfrow = c(1,2), mar = c(5,8,5,2), cex.lab = 1.5)
plot(bo$elevation, lam0, cex = 1, pch = 16, xlab = "Elevation",
     ylab = "Expected counts (lambda)", main = "Excluding spatial field",
     frame = FALSE, col = rgb(0, 0, 0, 0.3))
plot(bo$elevation, lam, cex = 1, pch = 16, xlab = "Elevation",
     ylab = "Expected counts (lambda)", main = "Including spatial field",
     frame = FALSE, col = rgb(0, 0, 0, 0.3))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Determine actual abundances as Poisson random variables with parameter lam
N <- rpois(n = nsites, lambda = lam)
table(N)            # Distribution of abundances across sites
sum(N > 0) / nsites # Finite-sample occupancy
(totalN <- sum(N))  # Total population size in all 2500 sites

# Create wind speed observational covariate
wind <- matrix(rnorm(nsites*nreps), nrow = nsites, ncol = nreps)

# Observation process
alpha0 <- 0   # logit-linear intercept
alpha1 <- -1  # slope on forest
alpha2 <- -1  # slope on wind speed
p <- array(NA, dim = c(nsites, nreps))
for(j in 1:nreps){
  p[,j] <- plogis(alpha0 + alpha1 * forest + alpha2 * wind[,j])
}

# Count things
y <- array(dim = c(nsites, nreps)) # Array for counts
for (j in 1:nreps){
  y[,j] <- rbinom(n = nsites, size = N, prob = p[,j])
}
str(y)
# int [1:2500, 1:3] 1 0 0 0 0 0 0 3 0 0 ...
summary(N)
summary(c(y))
#  Min. 1st Qu. Median  Mean 3rd Qu.   Max.
# 0.000   0.000  1.000 2.877   3.000 66.000
#  Min. 1st Qu. Median  Mean 3rd Qu.   Max.
# 0.000   0.000  0.000 1.174   1.000 49.000

# Compare true and observed total abundance
(true <- totalN)                  # True
(obs <- sum(apply(y, 1, max)))    # Observed
cat("Underestimation of total abundance:", round(100*(1-obs/true)), "%\n")
# [1] 7192
# [1] 4371
# Underestimation of total abundance: 39 %

# Select a sample of sites for surveys
# set.seed(100)
set.seed(100, sample.kind = "Rounding")
sample.size <- 500
sample.sites <- sort(sample(1:nsites, size = sample.size))

# ~~~~ extra code for figure 9.4 ~~~~
op <- par(mfrow = c(1,3), mar = c(3,3,3,6))
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = N))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
     main = "Abundance (N, truncated at 6)", zlim = c(0, 6))
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = apply(p, 1, mean)))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
     main = "Average detection probability")
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = apply(y, 1, max)))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
     main = "Max count (truncated at 6)", zlim = c(0, 6))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

yobs <- y                    # Make a copy
yobs[-sample.sites,] <- NA   # Turn counts of unsurveyed sites into NAs
head(sample.sites)           # Look at the simulated data set
head(yobs)

simNmixSpatial(nsurveys = 3, mean.lambda = exp(2), beta = c(2, -2),
               mean.p = 0.5, alpha = c(-1, -1), sample.size = 500, variance.RF = 1,
               theta.RF = 10, seeds = c(10, 100), truncN = 6, show.plots = TRUE)

simOccSpatial(nsurveys = 3, mean.psi = 0.6, beta = c(2, -2),
              mean.p = 0.4, alpha = c(-1, -1), sample.size = 500, variance.RF = 1,
              theta.RF = 10, seeds = c(10, 100), show.plots = TRUE)

dat <- simNmixSpatial(nsurveys = 3, mean.lambda = exp(2), beta = c(2, -2),
                      mean.p = 0.5, alpha = c(-1, -1), sample.size = 500, variance.RF = 1,
                      theta.RF = 10,seeds = c(10, 100), truncN = 6, show.plots=FALSE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 9.3 Fitting a nonspatial N-mixture model to the simulated data
# ==============================================================

dat <- simNmixSpatial(nsurveys = 3, mean.lambda = exp(2), beta = c(2, -2),
                      mean.p = 0.5, alpha = c(-1, -1), sample.size = 500, variance.RF = 1,
                      theta.RF = 10,seeds = c(10, 100), truncN = 6, show.plots=FALSE)

# ~~~ extra WinBUGS code for Nmix fitting ~~~~~~~~~~~~~~~~~~~
# Bundle data
bdata <- with(dat, list(y = yobs, nsites = dim(y)[1], nrep = dim(y)[2],
                        elev = elevationS, forest = forestS, wind = wind))
str(bdata)
# List of 6
# $ y     : int [1:2500, 1:3] NA NA 0 NA NA NA NA 3 NA NA ...
# $ nsites: int 2500
# $ nrep  : int 3
# $ elev  : num [1:2500] 1.06 1.836 1.763 1.305 0.268 ...
# $ forest: num [1:2500] 1.146 -0.363 -0.363 0.208 0.493 ...
# $ wind  : num [1:2500, 1:3] 0.534 1.369 -0.426 0.747 -0.414 ...

# Specify model in BUGS language
cat(file = "Nmix.txt", "
model {
  # Priors
  beta0 <- log(mean.lam)
  mean.lam ~ dunif(0, 20)
  alpha0 <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  for(v in 1:2){
    alpha[v] ~ dnorm(0, 0.1)
    beta[v] ~ dnorm(0, 0.1)
  }
  # Model for abundance
  for (i in 1:nsites){
   loglam[i] <- beta0 + beta[1] * elev[i] + beta[2] * pow(elev[i],2)
   loglam.lim[i] <- min(1000, max(-1000, loglam[i]))  # Stabilize log
    lam[i] <- exp(loglam.lim[i])
    N[i] ~ dpois(lam[i])
  }
  # Measurement error model
  for (i in 1:nsites){
    for (j in 1:nrep){
      y[i,j] ~ dbin(p[i,j], N[i])
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(1000, max(-1000, lp[i,j]))  # Stabilize logit
      lp[i,j] <- alpha0 + alpha[1] * forest[i] + alpha[2] * wind[i,j]
    }
  }
  # Derived parameters: total population size in grid
  Ntotal <- sum(N[])
}
")

# Initial values
Nst <- apply(dat$yobs, 1, max)# Max observed abundance as inits for N
Nst[is.na(Nst)] <- 2
Nst[Nst == 0] <- 2
inits <- function(){ list(N = Nst, mean.lam = 1, beta = rep(0, 2),
                          mean.p = 0.5, alpha = rep(0, 2))}

# Parameters monitored
params <- c("mean.lam", "beta0", "beta", "mean.p", "alpha0", "alpha", "Ntotal", "lam")

# MCMC settings
ni <- 2000    ;    nt <- 10    ;    nb <- 1000    ;    nc <- 3  # 13 mins

# Call WinBUGS from R (ART 18 min) and summarize posteriors
library(R2OpenBUGS)
out1 <- bugs(bdata, inits, params, "Nmix.txt", n.chains = nc, n.thin = nt,
             n.iter = ni, n.burnin = nb, debug = TRUE,
             working.directory = getwd())

print(out1, dig = 2)
# save result for comparison with other models
save(out1, file="AHM2_09.3_out1.RData")


# Chapter 9 : SPATIAL MODELS OF DISTRIBUTION AND ABUNDANCE
# ========================================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 20 mins
# Run time with the full number of iterations: 3 hrs

library(AHMbook)
library(R2OpenBUGS)
bugs.dir <- "C:/Program Files (x86)/OpenBUGS/OpenBUGS323" 
# the location of the OpenBUGS.exe file on your machine
setwd(choose.dir()) #don't forget to set wd!
library(fields)
library(raster)

# ~~~~~~~ uses data from 9.2 ~~~~~~~~
data(BerneseOberland)
bo <- BerneseOberland
RNGversion("3.5.3")
dat <- simNmixSpatial(nsurveys = 3, mean.lambda = exp(2), beta = c(2, -2),
                      mean.p = 0.5, alpha = c(-1, -1), sample.size = 500, variance.RF = 1,
                      theta.RF = 10,seeds = c(10, 100), truncN = 6, show.plots=FALSE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 9.4 Descriptive models of spatial autocorrelation
# =================================================

# 9.4.1 Fitting an intrinsic conditional autoregressive model
# -----------------------------------------------------------
?dnearneigh
# Compute the queen's neighborhood
library(spdep)
coordgrid <- cbind(bo$x, bo$y)
neigh <- dnearneigh(coordgrid, d1 = 0, d2 = sqrt(2) * 1000 + 1)
winnb <- nb2WB(neigh)      # Function to get CAR ingredients for BUGS

str(winnb)
# List of 3
# $ adj    : int [1:19404] 2 51 52 1 3 51 ...   # Index of neighbours
# $ weights: num [1:19404] 1 1 1 1 1 1 1 11 ... # Weights
# $ num    : int [1:2500] 3 5 5 5 5 5 5 5 ...   # Size of neighbourhood

# Bundle data
# str(bdata <- list(y = yobs, nsites = dim(y)[1], nrep = dim(y)[2],
# adj = winnb$adj, weights = winnb$weights, num = winnb$num,
# elev = as.numeric(elev), forest = as.numeric(forest), wind = wind))
str(bdata <- with(dat, list(y = yobs, nsites = dim(y)[1], nrep = dim(y)[2],
                            adj = winnb$adj, weights = winnb$weights, num = winnb$num,
                            elev = elevationS, forest = forestS, wind = wind)))
# List of 9
# $ y      : int [1:2500, 1:3] NA NA 0 NA NA NA NA 3 NA NA ...
# $ nsites : int 2500
# $ nrep   : int 3
# $ adj    : int [1:19404] 2 51 52 1 3 51 52 53 2 4 ...
# $ weights: num [1:19404] 1 1 1 1 1 1 1 1 1 1 ...
# $ num    : int [1:2500] 3 5 5 5 5 5 5 5 5 5 ...
# $ elev   : num [1:2500] 1.06 1.836 1.763 1.305 0.268 ...
# $ forest : num [1:2500] 1.146 -0.363 -0.363 0.208 0.493 ...
# $ wind   : num [1:2500, 1:3] 0.534 1.369 -0.426 0.747 -0.414 ...

# Specify model in BUGS language
cat(file = "CAR.Nmix.txt", "
model {
  # Specify priors
  beta0 <- log(mean.lam)
  mean.lam ~ dunif(0, 20)
  alpha0 <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  for(v in 1:2){
    alpha[v] ~ dnorm(0, 0.1)
    beta[v] ~ dnorm(0, 0.1)
  }
  # CAR prior distribution for spatial random effects
  eta[1:nsites] ~ car.normal(adj[], weights[], num[], tau)
  v.eta ~ dnorm(0, 0.01)I(0,)
  tau <- 1/v.eta
  # Model for abundance
  for (i in 1:nsites){
    loglam[i] <- beta0 + beta[1] * elev[i] + beta[2] * pow(elev[i],2) + eta[i]
    loglam.lim[i] <- min(1000, max(-1000, loglam[i])) # 'Stabilize' log
    lam[i] <- exp(loglam.lim[i])
    N[i] ~ dpois(lam[i])
  }
  # Measurement error model
  for (i in 1:nsites){
    for (j in 1:nrep){
      y[i,j] ~ dbin(p[i,j], N[i])
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(1000, max(-1000, lp[i,j])) # 'Stabilize' logit
      lp[i,j] <- alpha0 + alpha[1] * forest[i] + alpha[2] * wind[i,j]
    }
  }
  # Derived parameters: Total population size on grid
  Ntotal <- sum(N[])
}
")

# Initial values
Nst <- apply(dat$yobs, 1, max) # Max observed abundance as inits for N
Nst[is.na(Nst)] <- 2
Nst[Nst == 0] <- 2
inits <- function(){ list(N = Nst, mean.lam = 1, beta = rep(0, 2), mean.p = 0.5,
                          alpha = rep(0, 2), eta = rep(0, nrow(coordgrid)))}

# Parameters monitored
params <- c("mean.lam", "beta0", "beta", "mean.p", "alpha0", "alpha", "v.eta",
            "Ntotal", "eta", "lam")

# MCMC settings
# ni <- 200000 ; nt <- 100 ; nb <- 100000 ; nc <- 3
ni <- 15000 ; nt <- 1 ; nb <- 10000 ; nc <- 3  # ~~~ for testing, 17 mins

# Call WinBUGS from R (ART longish) and summarize posteriors
# out2 <- bugs(bdata, inits, params, "CAR.Nmix.txt", n.chains = nc, n.thin = nt,
# n.iter = ni, n.burnin = nb, bugs.directory = bugs.dir,
# debug = TRUE, working.directory = getwd()) # Close program manually !
out2 <- bugs(bdata, inits, params, "CAR.Nmix.txt", n.chains = nc, n.thin = nt,
             n.iter = ni, n.burnin = nb, working.directory = getwd())  # ~~~~ for testing

print(out2$summary[1:10,], 2)
#               mean      sd    2.5%      25%       50%      75%   97.5% Rhat n.eff
# mean.lam    3.5569   0.386    2.86    3.289    3.5310    3.799    4.34    1  3000
# beta0       1.2630   0.109    1.05    1.191    1.2620    1.335    1.47    1  3000
# beta[1]     1.9796   0.164    1.67    1.865    1.9770    2.091    2.31    1  2100
# beta[2]    -1.9704   0.136   -2.24   -2.061   -1.9690   -1.878   -1.71    1  2100
# mean.p      0.5010   0.026    0.45    0.484    0.5013    0.518    0.55    1  2300
# alpha0      0.0041   0.102   -0.20   -0.065    0.0051    0.071    0.20    1  2300
# alpha[1]   -0.9065   0.088   -1.08   -0.965   -0.9063   -0.848   -0.73    1  2300
# alpha[2]   -0.9531   0.066   -1.08   -0.996   -0.9522   -0.909   -0.82    1  3000
# v.eta       2.3962   0.404    1.69    2.119    2.3670    2.652    3.26    1  2000
# Ntotal   7235.9227 649.527 6113.90 6771.000 7187.5000 7650.250 8659.17    1  3000

# ~~~~ save for comparison with other models ~~~~~~~
save(out2, file="AHM2_09.4.1_out2.RData")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

with(dat, cbind(beta0, beta[1], beta[2], alpha0, alpha[1], alpha[2],
                Ntotal = sum(N), summaxC = sum(apply(y,1,max))))
#      beta0 beta1 beta2 alpha0 alpha1 alpha2 Ntotal summaxC
# [1,]     2     2    -2      0     -1     -1   7192    4371

# ~~~ extra code for figure 9.6 ~~~~~
# Compute average detection probability for each cell
phat <- array(NA, dim = c(2500, 3))
for(j in 1:3){
  phat[,j] <- plogis(out2$mean$alpha0 + out2$mean$alpha[1] * dat$forestS +
                       out2$mean$alpha[2] * dat$wind[,j])
}
pmean <- apply(phat, 1, mean)

op <- par(mfrow = c(3, 2), mar = c(3,3,3,4))
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = dat$lam))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
     main = "Density lambda (true)", zlim = c(0, 10))
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = out2$mean$lam))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
     main = "Density lambda (estimate)", zlim = c(0, 10))
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = c(dat$field)))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
     main = "Spatial effect eta (true)")
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = out2$mean$eta))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
     main = "Spatial effect eta (estimate)")
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = out2$sd$lam))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
     main = "Estimation uncertainty (lambda)", zlim = c(0, 10))
r <- rasterFromXYZ(data.frame(x = bo$x, y = bo$y, z = pmean))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
     main = "Average detection probability)", zlim = c(0, 1))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~ extra code for table of results ~~~~
# needs out1 object from section
load("AHM2_09.3_out1.RData")
truth <- with(dat, c(alpha0 = alpha0, alpha1 = alpha[1], alpha2 = alpha[2],
                     beta0 = beta0, beta1 = beta[1], beta2 = beta[2], Ntotal = sum(N)))
Bayes.est.Nmix0 <- rbind(out1$summary[c('alpha0','alpha[1]', 'alpha[2]',
                                        'beta0', 'beta[1]', 'beta[2]', 'Ntotal'), 1:2])
Bayes.est.CAR <- rbind(out2$summary[c('alpha0', 'alpha[1]', 'alpha[2]',
                                      'beta0', 'beta[1]', 'beta[2]', 'Ntotal'), 1:2])
print(cbind(truth, Bayes.est.Nmix0, Bayes.est.CAR), 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        truth    mean      sd      mean      sd
# alpha0     0    0.31   0.073    0.0041   0.102
# alpha1    -1   -0.90   0.067   -0.9065   0.088
# alpha2    -1   -1.14   0.060   -0.9531   0.066
# beta0      2    1.29   0.052    1.2630   0.109
# beta1      2    2.45   0.102    1.9796   0.164
# beta2     -2   -1.80   0.086   -1.9704   0.136
# Ntotal  7192 5863.50 238.104 7235.9227 649.527


###################################################################################


# Chapter 9 : SPATIAL MODELS OF DISTRIBUTION AND ABUNDANCE
# ========================================================
# Code from proofs dated 2020-08-19

# Approximate execution time for this code: 7 mins
# Run time with the full number of iterations: 30 mins

library(AHMbook)
library(raster)
library(jagsUI) #YAY can use JAGS for GAMs!!

# 9.4 Descriptive models of spatial autocorrelation
# =================================================

# 9.4.2 Fitting a two-dimensional generalized additive model
# ----------------------------------------------------------

# Re-create data set and sample 500 sites (using the same seeds)
library(AHMbook)
RNGversion("3.5.3")
str(dat <- simNmixSpatial(nsurveys = 3, mean.lambda = exp(2),
                          beta = c(2, -2), mean.p = 0.5, alpha = c(-1, -1), sample.size = 500,
                          variance.RF = 1, theta.RF = 10, seeds = c(10, 100), truncN = 6,
                          show.plots = TRUE))

# Scale both sets of coordinates
head(coordgrid <- scale(cbind(dat$xcoord, dat$ycoord))) # All 2500 cells
head(sitelocs <- coordgrid[dat$surveyed.sites,])        # 500 surveyed cells

library(fields)
system.time( knots <- cover.design(R = coordgrid, nd = 125, nruns = 10,
                                   num.nn = 200, max.loop=20) )                        # takes about 4 min

# Define the Z matrix for the random effects/knot coefficients
knotlocs <- knots$design
omega <- (e2dist(knotlocs, knotlocs)/10)^3
svd.omega <- svd(omega)
sqrt.omega <- t(svd.omega$v %*% (t(svd.omega$u)*sqrt(svd.omega$d)))
Zk <- (e2dist(coordgrid, knotlocs)/10)^3
Zmat <- t(solve(sqrt.omega, t(Zk)))

# Visualize basis vectors (as in Fig. 9.7)
head(Zmat) # Look at first couple of values
library(raster)
devAskNewPage(ask = dev.interactive(orNone=TRUE))  # ~~~ better for testing
op <- par(mfrow = c(3, 3), mar = c(2,2,4,11), "ask")
# for(i in 1:125){
for(i in 1:9){  # ~~~ for testing
  r <- rasterFromXYZ(data.frame(x = coordgrid[,1], y = coordgrid[,2],
                                z = Zmat[,i]))
  plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
       main = paste("Basis vector", i), legend = TRUE,
       axis.args = list(cex.axis = 2), legend.width = 2)
  points(knots$design[i,1], knots$design[i,2], col = 'black',
         pch = 16, cex = 1.5)
}
par(op)

# Bundle data
n.knots <- 125 # Fill in the desired spline bases and #knots
y <- dat$yobs # 2500 sites, but not surveyed NA'd out
str(bdata <- list(y = y, nsites = dim(y)[1], sample.size = dat$sample.size,
                  nrep = dim(y)[2], elev = dat$elevationS, forest = dat$forestS,
                  wind = dat$wind, n.knots = n.knots, Zmat = Zmat))
# List of 9
# $ y          : int [1:2500, 1:3] NA NA 0 NA NA NA NA 3 NA NA ...
# $ nsites     : int 2500
# $ sample.size: num 500
# $ nrep       : int 3
# $ elev       : num [1:2500] 1.06 1.836 1.763 1.305 0.268 ...
# $ forest     : num [1:2500] 1.146 -0.363 -0.363 0.208 0.493 ...
# $ wind       : num [1:2500, 1:3] 0.534 1.369 -0.426 0.747 -0.414 ...
# $ n.knots    : num 125
# $ Zmat       : num [1:2500, 1:125] 0.00548 0.00545 0.00544 ...

# Specify model in BUGS language
cat(file = "2dSplines.Nmix.txt", "
model {
  # Specify priors: intercepts and slopes of lambda and p
  beta0 <- log(mean.lam)
  mean.lam ~ dunif(0, 30)
  alpha0 <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  for(v in 1:2){
    alpha[v] ~ dnorm(0, 0.1)
    beta[v] ~ dnorm(0, 0.1)
  }
  # Priors on random effects parameters representing the splines
  for (k in 1:n.knots){
    b[k] ~ dnorm(0, tau.b)
  }
  # Prior on random effects dispersion
  tau.b ~ dgamma(0.1, 0.1)
  sd.b <- pow(tau.b, -2)
  # Model for abundance
  for (i in 1:nsites){
    N[i] ~ dpois(lam[i])
    log(lam[i]) <- beta0 + beta[1] * elev[i] + beta[2] * pow(elev[i],2) + smooth[i]
    smooth[i] <- smooth2[i] - mean(smooth2[])
    smooth2[i] <- inprod(Zmat[i,], b[])
  }
  # Measurement error model
  for (i in 1:nsites){
    for (j in 1:nrep){
      y[i,j] ~ dbin(p[i,j], N[i])
      logit(p[i,j]) <- alpha0 + alpha[1] * forest[i] + alpha[2] * wind[i,j]
    }
  }
  # Derived parameters: Total population size on grid
  Ntotal <- sum(N[])
}
")

# Initial values
Nst <- apply(dat$yobs, 1, max)
Nst[is.na(Nst)] <- 2
Nst[Nst == 0] <- 2
inits <- function(){ list(N = Nst, mean.lam = exp(rnorm(1)), beta = rnorm(2),
                          mean.p = runif(1), alpha = rnorm(2), b = runif(bdata$n.knots, -0.5, 0.5))}

# Parameters monitored
params <- c("mean.lam", "beta0", "beta", "mean.p", "alpha0", "alpha", "Ntotal",
            "sd.b", "b", "lam", "smooth")

# MCMC settings
# na <- 10000 ; ni <- 40000 ; nt <- 20 ; nb <- 20000 ; nc <- 3
na <- 1000 ; ni <- 4000 ; nt <- 2 ; nb <- 2000 ; nc <- 3  # ~~~ for testing, 3 mins

# Call JAGS (ART 40 min), gauge convergence and summarize posteriors
out3 <- jags(bdata, inits, params, "2dSplines.Nmix.txt", n.adapt = na,
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
summary(out3)
which(out3$summary[,8] > 1.1) ; sum(out3$summary[,8] > 1.1)
options(scipen = 10)
print(out3$summary[1:18, -c(4:6, 8:11)], dig = 3) # not shown

# ~~~ extra code to do the table ~~~~~~~
load("AHM2_09.3_out1.RData")
load("AHM2_09.4.1_out2.RData")
truth <- with(dat, c(alpha0 = alpha0, alpha1 = alpha[1], alpha2 = alpha[2], beta0 = beta0,
                     beta1 = beta[1], beta2 = beta[2], Ntotal = sum(N)))
Bayes.est.Nmix0 <- rbind(out1$summary[c('alpha0', 'alpha[1]', 'alpha[2]',
                                        'beta0', 'beta[1]', 'beta[2]', 'Ntotal'), 1:2])
Bayes.est.CAR <- rbind(out2$summary[c('alpha0', 'alpha[1]', 'alpha[2]',
                                      'beta0', 'beta[1]', 'beta[2]', 'Ntotal'), 1:2])
Bayes.est.spline <- rbind(out3$summary[c('alpha0', 'alpha[1]', 'alpha[2]',
                                         'beta0', 'beta[1]', 'beta[2]', 'Ntotal'), 1:2])
print(cbind(truth, Bayes.est.Nmix0, Bayes.est.CAR, Bayes.est.spline), 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Truth, nonspatial model, spatial Nmix with CAR, spatial Nmix with GAM
#        truth    mean      sd      mean      sd    mean      sd
# alpha0     0    0.31   0.073    0.0041   0.102    0.10   0.091
# alpha1    -1   -0.90   0.067   -0.9065   0.088   -0.91   0.072
# alpha2    -1   -1.14   0.060   -0.9531   0.066   -1.01   0.063
# beta0      2    1.29   0.052    1.2630   0.109    1.32   0.092
# beta1      2    2.45   0.102    1.9796   0.164    1.92   0.127
# beta2     -2   -1.80   0.086   -1.9704   0.136   -2.01   0.106
# Ntotal  7192 5863.50 238.104 7235.9227 649.527 6790.76 381.564

# ~~~~ extra code for figures 9.8, 9.9 and 9.10 ~~~~~~~~
# Fig. 9.8
op <- par(mfrow = c(2, 2), mar = c(3,3,3,4))
r <- rasterFromXYZ(data.frame(x = coordgrid[,1], y = coordgrid[,2], z = c(dat$field)))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
     main = "Spatial effect (true)")
r <- rasterFromXYZ(data.frame(x = coordgrid[,1], y = coordgrid[,2], z = out3$mean$smooth))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
     main = "Spatial effect (estimate)")
r <- rasterFromXYZ(data.frame(x = coordgrid[,1], y = coordgrid[,2], z = dat$lam))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
     main = "Density lambda (Truth)", zlim = c(0, 10))
r <- rasterFromXYZ(data.frame(x = coordgrid[,1], y = coordgrid[,2], z = out3$mean$lam))
plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
     main = "Density lambda (Estimate)", zlim = c(0, 10))
par(op)

# Fig. 9.9, a comparison of cell-wise density
xylim <- c(0, 50)
op <- par(mfrow = c(1, 2))
plot(c(dat$lam), out2$mean$lam, xlab = 'True', ylab = 'Estimated',
     main = 'Spatial Nmix (CAR)', frame = FALSE, pch = 16, cex = 1.5,
     col = rgb(0,0,0,0.3), xlim = xylim, ylim = xylim)
abline(0, 1)
abline(lm(out2$mean$lam ~ c(dat$lam)), col = 'blue', lwd = 3)
plot(c(dat$lam), out3$mean$lam, xlab = 'True', ylab = 'Estimated',
     main = 'Spatial Nmix (2d splines)', frame = FALSE, pch = 16, cex = 1.5,
     col = rgb(0,0,0,0.3), xlim = xylim, ylim = xylim)
abline(0, 1)
abline(lm(out3$mean$lam ~ c(dat$lam)), col = 'blue', lwd = 3)
par(op)

# Fig 9.10 Visualize contributions of piece-wise regressions vectors
op <- par(mfrow = c(3, 4), mar = c(1,1,4,1))
for(i in 1:12){
  piecewise.reg <- Zmat[,i] * out3$mean$b[i]
  r <- rasterFromXYZ(data.frame(x = coordgrid[,1], y = coordgrid[,2],
                                z = piecewise.reg))
  plot(r, col = topo.colors(20), axes = FALSE, box = FALSE,
       main = paste("Piecewise regression", i) )
  points(knots$design[i,1], knots$design[i,2], col = 'black', pch = 16, cex = 1)
}
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

