##############################################
###N-mixture Model for SALS with TR (no DS)###
##############################################

##Might be good to check out the model described in AHM Vol. 1 section 10.12.2?
##Uses TTD (TR) data and also looks as sex-ratios!?

#################
###Description###
#################

##This is a single-species static abundance model used to estimate densities
##of saltmarsh sparrow (SALS) among sampled sites in coastal New England. Survey
##data were collected by SHARP. Time-Removal (time-to-detection) data is integrated 
##into the model, but distance sampling data is not used.


#######################################
###Load Required Packages and Set WD###
#######################################
library(jagsUI)
library(AHMbook)
setwd(choose.dir())
source(file="simNmixTR.R") #for data simulation function if desired

#############################
###Read-in Data and Format###
#############################
#Data
data_2011<-read.csv("SHARP_Surveys_2011.csv")
veg_2011<-read.csv("SHARP_Veg_2011.csv")
  
##Quick summary of data to see what's going on here
str(data_2011) ; summary(data_2011)

##For this analysis, only using observations of SALS
data_2011<-data_2011[-which(data_2011$TotalCountN=="#VALUE!"),] #remove error rows
##Only selecting sites that were surveyed all 3 times throughout the season
sites<-unique(data_2011$PointID) #all unique sites

for(i in 1:length(sites)){
  site.index<-which(data_2011$PointID==sites[i])
  if(rje::is.subset(c(1,2,3), data_2011$VisitNum[site.index])==TRUE){
    site.index[i]<-site.index[i]
  }else{
    site.index[i]<-NA
  }
} 
summary(is.na(sites)) #all sites sampled 3 times, so we're all good

data_2011$TotalCountN<-as.numeric(data_2011$TotalCountN) #making numeric
SALS_Obs<-data_2011[grep("SALS",data_2011$AlphaCode),] #SALS observations only
SALS_Obs<-SALS_Obs[-which(rowSums(SALS_Obs[,46:50])==0),] #removing sites w/out T-to-D data
sites<-unique(data_2011$PointID) #names of point count locations (1,652 sites)


#Filling necessary values
nsites<-length(sites)
nsurveys<-3
K<-5 #Number of time intervals

#Making sure number of SALS observed matches number of time removals
counts<-rowSums(SALS_Obs[,46:50])
length(which(SALS_Obs$TotalCountN!=counts)) #too many...
#Modifying counts to make sure that number of unique time removals == number of individuals per site
SALS_Obs$TotalCountN<-counts
n.dets<-nrow(SALS_Obs) #number of unique detection histories

##New approach for modeling: Multinomial probability for number of individuals in each time interval
##No need to make each line one individual now!
##Need two matrices of data for multinomial probabilities (only observations used):

tint<-as.matrix(SALS_Obs[,46:50]) #Observations of SALS per min during each survey
#Metadata for site name, survey number, and number of inds for each row in tint matrix
obs.md<-data.frame(cbind(SALS_Obs$PointID, SALS_Obs$VisitNum, SALS_Obs$TotalCountN))
obs.md$X2<-as.numeric(obs.md$X2)
obs.md$X3<-as.numeric(obs.md$X3)
colnames(obs.md)<-c("site","survey","count")

#Fill observation array (code is ugly but works fine)
n.obs<-array(0,c(nsites,nsurveys))
for(i in 1:nrow(SALS_Obs)){
  n.obs[grep(SALS_Obs$PointID[i],sites), 
        SALS_Obs$VisitNum[i]]<-sum(SALS_Obs[i,46:50])
}

dimnames(n.obs)<-list(sites,c(1:3)) #Naming dimensions

##Next Step: Add Covariates!
#Noise for detection, Temperature for availability, Vegetation data informs abundance, 
#time of survey data once I can figure out how to quickly automate sunrise times
#Survey level covariates (noise)

noise<-array(NA,c(nsites,nsurveys))
dimnames(noise)<-list(sites,c(1:3)) 
for(i in 1:nrow(data_2011)){
  noise[grep(data_2011$PointID[i],sites), data_2011$VisitNum[i]]<-data_2011$Noise[i]
}

#Temperature (quadratic effect)
temps<-array(NA,c(nsites,nsurveys))
dimnames(temps)<-list(sites,c(1:3)) 
for(i in 1:nrow(data_2011)){
  temps[grep(data_2011$PointID[i],sites), data_2011$VisitNum[i]]<-data_2011$TempF[i]
}

#Veg (just using % high marsh right now as a covariate for SALS)
high.marsh<-rep(NA,nsites)
names(high.marsh)<-sites 
for(i in 1:nrow(veg_2011)){
  high.marsh[grep(veg_2011$BirdPtID[i],sites)]<-veg_2011$HighMarshCC[i]
}

#Interpolating veg data for sites missing veg survey w/ mean value (49/1655 = 2.96% of sites)
high.marsh[which(is.na(high.marsh)==TRUE)]<-round(mean(high.marsh, na.rm = TRUE))
high.marsh[which(high.marsh==0.5)]<-1 #replacing 0.5 values with 1 to maintain consistency
#Same for noise and temps
noise[which(is.na(noise)==TRUE)]<-round(mean(noise, na.rm = TRUE))
temps[which(is.na(temps)==TRUE)]<-round(mean(temps, na.rm = TRUE))

##Standardize covariates
noise1<-as.numeric(scale(noise[,1]))
noise2<-as.numeric(scale(noise[,2]))
noise3<-as.numeric(scale(noise[,3]))
noise<-cbind(as.vector(noise1),as.vector(noise2),as.vector(noise3))
temps1<-as.numeric(scale(temps[,1]))
temps2<-as.numeric(scale(temps[,2]))
temps3<-as.numeric(scale(temps[,3]))
temps<-cbind(as.vector(temps1),as.vector(temps2),as.vector(temps3))
rm(noise1,noise2,noise3,temps1,temps2,temps3)
high.marsh<-as.numeric(scale(high.marsh)) #categorical scaled to continuous...

#Finally, rename sites to integers for correct indexing within model
for(i in 1:nrow(obs.md)){
  obs.md$site[i]<-which(sites==obs.md$site[i])
}
obs.md$site<-as.numeric(obs.md$site)
obs.md<-as.matrix(obs.md) #making matrix for JAGS

##Combine data for model
str(jags.data<-list(nsites = nsites, nsurveys = nsurveys, K = K, tint = tint, 
                    n.obs = n.obs, n.dets = n.dets, obs.md = obs.md, noise = noise, 
                    temps = temps, high.marsh = high.marsh))

##Or, if wanting to simulate data for the model instead...obviously can be modified if needed
data<-simNmixTR()
#
#str(jags.data<-list(nsites = data$nsites, nsurveys = data$nsurveys, K = data$K, 
#                    tint = data$tint, n.obs = data$n.obs, n.dets = data$n.dets, 
#                    obs.md = data$obs.md, noise = data$det.cov, 
#                    temps = data$avail.cov, high.marsh = data$ab.cov))

#########################
###Model Specification###
#########################

cat("
model {
#i=site ; j=survey ; k=time-window

###Prior distributions for basic parameters

# Intercepts
beta0 ~ dnorm(0,0.01)      # intercept for lambda
beta.a0 ~ dnorm(0,0.01)    # intercept for availability
alpha0 ~ dnorm(0,0.01)     # intercept for detection
# Coefficients
beta1 ~dnorm(0,0.01)       # slope for lambda covariate
beta.a1 ~ dnorm(0,0.01)    # slope for availability covariate
#beta.a2 ~ dnorm(0,0.01)    # slope for availability covariate (quadratic effect)
alpha1 ~ dnorm(0,0.01)     # slope for detection covariate

###Likelihoods

# Abundance 
for(i in 1:nsites){ 
  M[i] ~ dpois(lambda[i])               # Latent abundance at each site
  # Add site-level covariates to lambda
  log(lambda[i]) <- beta0 + beta1*high.marsh[i] 
}

# Availability
for(i in 1:nsites){
  for(j in 1:nsurveys){ 
    for (k in 1:K){
      pi.pa[i,j,k] <- p.a[i,j] * pow(1-p.a[i,j], (k-1))  
      pi.pa.c[i,j,k] <- pi.pa[i,j,k]/sum(pi.pa[i,j,1:K]) # Conditional probabilities of availability
    }
    pi.pa[i,j,K+1] <- 1 - sum(pi.pa[i,j,1:K]) #residual probability? Not sure if necessary
    phi[i,j] <- sum(pi.pa[i,j,1:K]) # Probability of ever available
    # Add covariates for availability here TIME-REMOVAL (availability)
    logit(p.a[i,j]) <- beta.a0 + beta.a1*temps[i,j]
  }
}

# Conditional observation model for categorical covariates 
# (SALS counts per bin as a function of availability probability)
for (i in 1:n.dets){
  tint[i,1:K] ~ dmulti(pi.pa.c[obs.md[i,1],obs.md[i,2],1:K], obs.md[i,3]) 
}#i
    
# Detection
for(i in 1:nsites){
 for(j in 1:nsurveys){
  logit(p[i,j]) <- alpha0 +  alpha1*noise[i,j] 
 }
}

##Combining Probabilities to Determine Latent States

for(i in 1:nsites){
  for(j in 1:nsurveys){
    # Binomial model for # of captured individuals
    n.obs[i,j] ~ dbin(p[i,j]*phi[i,j], M[i]) # Number of detected individuals at each site
    N.avail[i,j] ~ dbin(phi[i,j], M[i])       # Number of available individuals at each site
  }
}

# Derived quantities
Mtot <- sum(M)           # Total population size
#If wanting to est. density, need to truncate observations to <= 100m
#Dsuper <- Mtot/totArea     # Density of SALS within sampled area (individuals/ha)
p.mean<- mean(p)         # Mean perceptibility across sites and surveys
phi.mean <- mean(phi)    # Mean availability across sites and surveys
}
", fill=TRUE, file="SALS_TR.txt")


################################
###Inits/Params/MCMC Settings###
################################

# Create initial values (including for M and N.avail) and list parameters to save
# n.obs <- data$n.obs #if simulating data
N.avail.st <- n.obs + 1
Mst <- apply(N.avail.st, 1, max)  + 3

inits <- function(){list(M=Mst, N.avail=N.avail.st, alpha0=rnorm(n=1), 
                         beta0=rnorm(n=1), beta.a0=rnorm(n=1), alpha1=rnorm(n=1),
                         beta1=rnorm(n=1), beta.a1=rnorm(n=1))}

params <- c("beta.a0", "beta.a1", "alpha0", "alpha1", "beta0",
            "beta1", "phi.mean", "p.mean", "Mtot")

# MCMC settings
#For testing:
#ni <- 500   ;   nb <- 200   ;   nt <- 2   ;   nc <- 3
#For running:
ni <- 25000   ;   nb <- 5000   ;   nt <- 4   ;   nc <- 3


###################
###Running Model###
###################

# Run JAGS in parallel (ART 30 min), check convergence and summarize posteriors
out2a <- jags(data=jags.data, inits=inits, parameters=params, 
              model.file ="SALS_TR.txt",n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, 
              parallel = TRUE)

# Run JAGS in parallel (ART 52 min) for alternative model specification (same inits/params)
#out2b <- jags(data=jags.data, inits=inits, parameters=params, 
#              model.file ="SALS_TR_2.txt",n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, 
#              parallel = TRUE)

######################
###Examining Output###
######################

traceplot(out2a)   ;   print(out2a, 3)

#Save for later to reanalyze
#save(out2a, file = "SALS_TD_N_Mix_Out.R")
#Load other model specification to compare with new specs.
load(file = "SALS_TD_N_Mix_Out.R")
cbind(out2a$summary[,c(1,2,8)],out2b$summary[,c(1,2,8)])

#Specification B not converging reliably, smaller n.eff. Overall, specification
#A seems like the better approach (also makes more sense intuitively IMO).
############################################
###Trying out another model specification###
############################################

cat("
model {
#i=site ; j=survey ; k=time-window

###Prior distributions for basic parameters

# Intercepts
beta0 ~ dnorm(0,0.01)      # intercept for lambda
beta.a0 ~ dnorm(0,0.01)    # intercept for detection
alpha0 ~ dnorm(0,0.01)     # intercept for availability
# Coefficients
beta1 ~dnorm(0,0.01)       # slope for lambda covariate
beta.a1 ~ dnorm(0,0.01)    # slope for detection covariate
alpha1 ~ dnorm(0,0.01)     # slope for availability covariate

###Likelihoods

# Abundance 
for(i in 1:nsites){ 
  M[i] ~ dpois(lambda[i])               # Latent abundance at each site
  # Add site-level covariates to lambda
  log(lambda[i]) <- beta0 + beta1*high.marsh[i] 
}

# Detection
for(i in 1:nsites){
  for(j in 1:nsurveys){ 
    for (k in 1:K){
      pi.pa[i,j,k] <- p.a[i,j] * pow(1-p.a[i,j], (k-1))  
      pi.pa.c[i,j,k] <- pi.pa[i,j,k]/sum(pi.pa[i,j,1:K]) # Conditional probabilities of detection
    }
    pi.pa[i,j,K+1] <- 1 - sum(pi.pa[i,j,1:K]) #residual probability
    p[i,j] <- sum(pi.pa[i,j,1:K]) # Probability of ever detected
    pmarg[i,j] <- p[i,j]*phi[i,j]
    # Add covariates for detection here TIME-REMOVAL (detection)
    logit(p.a[i,j]) <- beta.a0 + beta.a1*noise[i,j]
  }
}

# Conditional observation model for categorical covariates 
# (SALS counts per bin as a function of detection probability)
for (i in 1:n.dets){
  tint[i,1:K] ~ dmulti(pi.pa.c[obs.md[i,1],obs.md[i,2],1:K], obs.md[i,3]) 
}#i
    
# Availability
for(i in 1:nsites){
 for(j in 1:nsurveys){
  logit(phi[i,j]) <- alpha0 + alpha1*temps[i,j]
 }
}

##Combining Probabilities to Determine Latent States

for(i in 1:nsites){
  for(j in 1:nsurveys){
    # Binomial model for # of captured individuals
    n.obs[i,j] ~ dbin(pmarg[i,j], M[i]) # Number of detected individuals at each site
    N.avail[i,j] ~ dbin(phi[i,j], M[i])     # Number of available individuals at each site
  }
}

# Derived quantities
Mtot <- sum(M)           # Total population size
N.a.tot <- sum(N.avail)  # Total available population size
#If wanting to est. density, need to truncate observations to <= 100m
#Dsuper <- Mtot/totArea     # Density of SALS within sampled area (individuals/ha)
p.mean<- mean(p)         # Mean perceptibility across sites and surveys
phi.mean <- mean(phi)    # Mean availability across sites and surveys
}
", fill=TRUE, file="SALS_TR_2.txt")
####################################

#####EXTRA (PROBABLY UNNECESSARY) CODE
#Duplicating rows for the number of birds observed
repped.counts<-which(SALS_Obs$TotalCountN>1)
sum(SALS_Obs$TotalCountN[repped.counts]) #508 individuals
SALS_Obs_Repped<-SALS_Obs[1,]
for(i in 1:length(repped.counts)){
  temp<-SALS_Obs[repped.counts[i],46:50]
  mat<-matrix(data = NA, nrow = sum(temp), ncol = 5)
  test<-c(rep(1,temp[1]),rep(0,nrow(mat)-temp[1]))
  test<-append(test,c(rep(0,sum(test)),rep(1,temp[2]),rep(0,nrow(mat)-(sum(test) + temp[2]))))
  test<-append(test,c(rep(0,sum(test)),rep(1,temp[3]),rep(0,nrow(mat)-(sum(test) + temp[3]))))
  test<-append(test,c(rep(0,sum(test)),rep(1,temp[4]),rep(0,nrow(mat)-(sum(test) + temp[4]))))
  test<-append(test,c(rep(0,sum(test)),rep(1,temp[5]),rep(0,nrow(mat)-(sum(test) + temp[5]))))
  mat<-matrix(test, nrow = sum(temp), ncol = 5, byrow = F)
  n<-SALS_Obs$TotalCountN[repped.counts[i]]-1
  temp2<-SALS_Obs[repped.counts[i],]
  for(j in 1:n){
    temp2<-rbind(temp2,SALS_Obs[repped.counts[i],])
  }
  temp2[,46:50]<-mat
  SALS_Obs_Repped<-rbind(SALS_Obs_Repped,temp2)
}

#All counts should now be 1 individual per line
SALS_Obs_Repped<-SALS_Obs_Repped[-1,] #removing first (placeholder) row
SALS_Obs_Repped$TotalCountN<-1 #making 1
SALS_Obs<-SALS_Obs[-repped.counts,]
SALS_Obs<-rbind(SALS_Obs,SALS_Obs_Repped) #combining back together

#Double-checking counts again to make sure everything worked
counts<-rowSums(SALS_Obs[,46:50])
length(which(SALS_Obs$TotalCountN!=counts)) #WORKED!!!

#Replacing capture history with time to detection variable (K=5)
SALS_Obs$tint<-NA
for(i in 1:length(SALS_Obs$tint)){
  if(SALS_Obs$Min.1[i]==1){
    SALS_Obs$tint[i]<-1
  }else if(SALS_Obs$Min.2[i]==1){
    SALS_Obs$tint[i]<-2
  }else if(SALS_Obs$Min.3[i]==1){
    SALS_Obs$tint[i]<-3
  }else if(SALS_Obs$Min.4[i]==1){
    SALS_Obs$tint[i]<-4
  }else{
    SALS_Obs$tint[i]<-5
  }
}
SALS_Obs$tint<-as.numeric(SALS_Obs$tint)


# Create the observed encounter frequencies per site (include the zeros! )
n <- rep(0,length(unique(data_2011$PointID))) # The full site vector
names(n) <- unique(data_2011$PointID)
n[names(table(SALS_Obs$PointID))] <- table(SALS_Obs$PointID)  # Put in the counts
names(n) <- seq(1:length(unique(data_2011$PointID)))
site <- match(SALS_Obs$PointID,unique(data_2011$PointID))
nobs <- nrow(SALS_Obs) 

