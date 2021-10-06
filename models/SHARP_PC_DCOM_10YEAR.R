###############################################################
###Dynamic Community Occupancy Model - SHARP 10 Year PC Data###
###############################################################

#################
###Description###
#################

#This model analyzes 10 years of point-count data conducted by SHARP and examines
#spatial variation in rates of population change over time for all members of the 
#tidal marsh breeding bird community. Using an integrated modeling approach, the 
#detection component in the model is parsed into availability and detection components
#to improve the accuracy of occupancy estimation. Explanatory variables (covariates)
#included in link-functions relate to variability in habitat suitability over space
#and time, as well as to site and survey-specific likelihoods of detection or availability.
#Much of this model code is based off of code presented in Applied Hierarchical Modeling 
#Volume 1 (2016) by Marc Kery and Andy Royle, as well as the supporting R package 'AHMbook'. 


#######################################
###Load Required Packages and Set WD###
#######################################

library(AHMbook)
library(jagsUI)
library(abind)
setwd(choose.dir())
#setwd("C:/Users/12152/Documents/UConn_Files/Project_Files/Data_Analysis/Survey_Data_Analyses")

##################
###Read in Data###
##################

survey.data<-read.csv("SHARP_Surveys_All.csv")

veg.data<-read.csv("SHARP_Rapid_Veg_All.csv")

#########################
###Prep Data for Model###
#########################

##Quick summary of data to see what's going on here
str(survey.data) ; summary(veg.data)

##For this analysis, only using observations of tidal marsh breeding birds and 
##only sites that were surveyed more than once per season
#Remove any error rows first
survey.data<-survey.data[-which(survey.data$TotalCountN=="#VALUE!"),]

##Only selecting sites that were surveyed multiple times per year

sites<-unique(survey.data$PointID) #all unique sites
years<-sort(unique(survey.data$Year)) #all sampled years (10)

#data frame to take values
survey.data.thinned<-survey.data[1,] #first row just a placeholder

for(i in 1:length(sites)){ #takes a LONG TIME to run due to large dataset (multiple hours)
  for(j in 1:length(years)){
    site.index<-which(survey.data$PointID==sites[i] & survey.data$Year==years[j])
    if(rje::is.subset(c(1,2), survey.data$VisitNum[site.index])==TRUE){
      survey.data.thinned<-rbind(survey.data.thinned,survey.data[site.index,])
    }else{
      
    }
  } 
}

#getwd()
#write.csv(survey.data.thinned, file = "Survey_Data_Thinned.csv")
#survey.data<-read.csv("Survey_Data_Thinned.csv") #replace full dataset with thinned data
summary(survey.data)

#Columns to save
survey.data<-survey.data[,colnames(survey.data)%in%c("AlphaCode","PointID","VisitNum","Year","TempF","Noise")]
survey.data<-survey.data[-which(survey.data$VisitNum==4),] #remove 4th site visits
species<-unique(survey.data$AlphaCode) #Alpha codes of observed species
v.rare<-which(table(survey.data$AlphaCode)<10) #Ultra-rare species (<10 obs over 10 years)
unk<-species[grep("^UN",species)] #Unknown species
v.rare<-append(names(v.rare),names(unk))
species<-species[!species%in%v.rare] #Ultra-rare and unknown removed
##Modification: Only running models for marsh specialist birds
species<-c("SESP", "NESP", "SALS", "MAWR", "CLRA", "WILL", "SNEG", "ABDU", "GREG", "VIRA", 
       "BTGR", "GLIB", "RWBL", "MALL", "TRES", "GBHE", "SAVS", "SOSP", "YEWA", "BARS", 
       "ALFL", "COYE")
survey.data<-survey.data[survey.data$AlphaCode%in%species,] #ditto for survey data
sites<-unique(survey.data$PointID) #names of point count locations
years<-sort(unique(survey.data$Year)) #names of years

#Filling necessary values
nspec<-length(species)
nsites<-length(sites)
nsurveys<-3 #Some sites surveyed 3 times, others twice (not a problem, JAGS can handle!)
nyears<-length(years)

#Fill observation array (code is ugly but works fine, takes a while to run though!)
y<-array(0,c(nspec,nsites,nsurveys,nyears))
for(i in 1:nrow(survey.data)){
  y[grep(survey.data$AlphaCode[i],species), grep(survey.data$PointID[i],sites), 
    survey.data$VisitNum[i], grep(survey.data$Year[i],years)]<-1
}
#save(y,file = "y_obs_array_SHARP_MSB_DCOM.R")
dimnames(y)<-list(species,sites,c(1:3),years) #Naming dimensions

#Give NAs for third sampling occasion at sites only sampled twice
survey.data.index<-1 #first value just a placeholder
for(i in 1:length(sites)){ #takes a LONG TIME to run due to large dataset (multiple hours)
  for(j in 1:length(years)){
    site.index<-which(survey.data$PointID==sites[i] & survey.data$Year==years[j])
    if(rje::is.subset(c(1,2,3), survey.data$VisitNum[site.index])==TRUE){
      survey.data.index<-append(survey.data.index,site.index)
    }else{
      
    }
  } 
}
survey.data.thinned<-survey.data[-survey.data.index,] #Sites sampled twice

#Fill NAs in observation array (code is ugly but works fine, takes a while to run though!)
for(i in 1:nrow(survey.data.thinned)){
  y[grep(survey.data.thinned$AlphaCode[i],species), grep(survey.data.thinned$PointID[i],sites), 
    survey.data.thinned$VisitNum[i], grep(survey.data.thinned$Year[i],years)]<-NA
}

#save(y,file = "y_obs_array_SHARP_DCOM.R") #resave file

#Adding model covariates
#Noise and Temperature for p, Vegetation data informs psi, phi, and gamma

noise<-array(NA,c(nsites,nsurveys,nyears))
dimnames(noise)<-list(sites,c(1:3),years) 
for(i in 1:nrow(survey.data)){
  noise[grep(survey.data$PointID[i],sites), survey.data$VisitNum[i],
        grep(survey.data$Year[i],years)]<-survey.data$Noise[i]
}

#Temperature 
temps<-array(NA,c(nsites,nsurveys,nyears))
dimnames(temps)<-list(sites,c(1:3),years) 
for(i in 1:nrow(survey.data)){
  temps[grep(survey.data$PointID[i],sites), survey.data$VisitNum[i],
        grep(survey.data$Year[i],years)]<-survey.data$TempF[i]
}

#Veg (just using % high marsh and % low marsh right now as covariates)
#Making a year column in the matrix
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
veg.data$Year<-substrRight(as.character(veg.data$Date),4)
veg.data$Year<-as.numeric(veg.data$Year)

high.marsh<-array(NA,c(nsites,nyears))
dimnames(high.marsh)<-list(sites,years) 

for(i in 1:nrow(veg.data)){
  high.marsh[grep(veg.data$BirdPtID[i],sites), 
             grep(veg.data$Year[i],years)]<-veg.data$HighMarshCC[i]
}

low.marsh<-array(NA,c(nsites,nyears))
dimnames(low.marsh)<-list(sites,years) 

for(i in 1:nrow(veg.data)){
  low.marsh[grep(veg.data$BirdPtID[i],sites), 
             grep(veg.data$Year[i],years)]<-veg.data$LowMarshCC[i]
}

#Interpolating veg data for sites missing veg survey w/ mean value (49/1655 = 2.96% of sites)
high.marsh[which(is.na(high.marsh)==TRUE)]<-round(mean(high.marsh, na.rm = TRUE))
high.marsh[which(high.marsh==0.5)]<-1 #replacing 0.5 values with 1 to maintain consistency
low.marsh[which(is.na(low.marsh)==TRUE)]<-round(mean(low.marsh, na.rm = TRUE))
low.marsh[which(low.marsh==0.5)]<-1 #replacing 0.5 values with 1 to maintain consistency

#Same for noise and temps
noise[which(is.na(noise)==TRUE)]<-round(mean(noise, na.rm = TRUE))
temps[which(is.na(temps)==TRUE)]<-round(mean(temps, na.rm = TRUE))

##Standardize covariates
for(i in 1:nsurveys){
  for(j in 1:nyears){
    noise[,i,j]<-as.numeric(scale(noise[,i,j]))
    temps[,i,j]<-as.numeric(scale(temps[,i,j]))
  }
}

for(t in 1:nyears){
  high.marsh[,t]<-as.numeric(scale(high.marsh[,t]))
  low.marsh[,t]<-as.numeric(scale(low.marsh[,t]))
}

#save first year of veg data as initial psi cov
high.marsh.1<-high.marsh[,1]
low.marsh.1<-low.marsh[,1]

##Save to arrays
p.covs<-abind(noise,temps, along=4) # site x survey x year x cov
psi.covs<-abind(high.marsh.1,low.marsh.1, along = 2) #site x cov (year 1 only)
phi.covs<-abind(high.marsh,low.marsh, along = 3)  #site x year x cov
gamma.covs<-abind(high.marsh,low.marsh, along = 3)  #site x year x cov
 

##Combine data for model
str(jags.data<-list(nspec = nspec, nsites = nsites, nsurveys = nsurveys, 
                    nyears = nyears, p.covs = p.covs, psi.covs = psi.covs,
                    phi.covs = phi.covs, gamma.covs = gamma.covs, y = y))

#save(jags.data, file = "SHARP_PC_DCOM_10YR_Model_Data.R")

################################
###Specify Model Code in JAGS###
################################

cat(file = "DCM_Few_Covs.txt", "
model {

  ##########################
  #Priors and Hyperpriors###
  ##########################
  
  # *** Priors and hyperpriors for model on psi1 (initial occupancy) ***
  # Priors
  for(k in 1:nspec){ # Loop over # of species 
    alpha.lpsi1[k] ~ dnorm(mu.alpha.lpsi1, tau.alpha.lpsi1) # Intercept
      for(g in 1:2){ # Loop over # of coefficients (= number of covariates)
      beta.lpsi1[g, k] ~ dnorm(mu.beta.lpsi1[g], tau.beta.lpsi1[g]) #Cov coeffs (spp. specific)
    }
  }
  
  # Hyperpriors
  mu.alpha.lpsi1 <- logit(mean.alpha.psi1) #Intercept mean
  mean.alpha.psi1 ~ dunif(0, 1) #Intercept mean
  tau.alpha.lpsi1 <- pow(sd.alpha.lpsi1, -2) #Intercept tau
  sd.alpha.lpsi1 ~ dunif(0, 10) #Intercept SD
  for(g in 1:2){ # Loop over 2 coefficients
    mu.beta.lpsi1[g] ~ dnorm(0, 0.1) #Coeff mean
    tau.beta.lpsi1[g] <- pow(sd.beta.lpsi1[g], -2)  #Coeff Tau
    sd.beta.lpsi1[g] ~ dnorm(0, 0.1)I(0,) # Half-Normal prior for Coeff SD
    # curve(dnorm(x, 0, sqrt(1 / 0.1)), 0, 20) # howsit look like ?
  }
  
  # *** Priors and hyperpriors for model on phi (colonization) ***
  # Priors
  for(k in 1:nspec){ # Loop over # of species 
    for(t in 1:(nyears-1)){ # Loop over t-1 intervals (colonization varies by year)
      alpha.lphi[t,k] ~ dnorm(mu.alpha.lphi[t], tau.alpha.lphi) #Intercept
      # phi intercept for t years, different mean, same variance
    }
    for(g in 1:2){ # Loop over 2 coefficients
      beta.lphi[g,k] ~ dnorm(mu.beta.lphi[g], tau.beta.lphi[g]) #Cov coeffs (spp. specific)
    }
  }
  
  # Hyperpriors
  for(t in 1:(nyears-1)){ # Loop over t-1 intervals (yearly variation in mean)
    mu.alpha.lphi[t] <- logit(mean.alpha.phi[t]) 
    mean.alpha.phi[t] ~ dunif(0, 1)
  }
  tau.alpha.lphi <- pow(sd.alpha.lphi, -2) #No yearly variation in variance 
  sd.alpha.lphi ~ dnorm(0, 0.1)I(0,)
  for(g in 1:2){ # Loop over 2 coefficients
    mu.beta.lphi[g] ~ dnorm(0, 0.01) #No yearly variation in hyperpriors for coeffs
    tau.beta.lphi[g] <- pow(sd.beta.lphi[g], -2)
    sd.beta.lphi[g] ~ dnorm(0, 0.1)I(0,)
  }
  
  # *** Priors and hyperpriors for model on gamma (persistence)***
  # Priors
  for(k in 1:nspec){ # Loop over # of species 
    for(t in 1:(nyears-1)){ # Loop over t-1 intervals (persistence varies by year)
      alpha.lgamma[t,k] ~ dnorm(mu.alpha.lgamma[t], tau.alpha.lgamma) #Intercept
      # gamma intercept for t years, different mean, same variance
    }
    for(g in 1:2){ # Loop over 2 coefficients
      beta.lgamma[g,k] ~ dnorm(mu.beta.lgamma[g], tau.beta.lgamma[g]) #Cov coeff (spp. specific)
    }
  }
  
  # Hyperpriors
  for(t in 1:(nyears-1)){ # Loop over t-1 intervals (yearly variation in mean)
    mu.alpha.lgamma[t] <- logit(mean.alpha.gamma[t])
    mean.alpha.gamma[t] ~ dunif(0, 1)
  }
  tau.alpha.lgamma <- pow(sd.alpha.lgamma, -2) #No yearly variation in variance
  sd.alpha.lgamma ~ dnorm(0, 0.1)I(0,)
  for(g in 1:2){ # Loop over 2 coefficients
    mu.beta.lgamma[g] ~ dnorm(0, 0.1)
    tau.beta.lgamma[g] <- pow(sd.beta.lgamma[g], -2)
    sd.beta.lgamma[g] ~ dnorm(0, 0.1)I(0,)
  }
  
  # *** Priors and hyperpriors for model on p (detection) ***
  # Priors
  for(k in 1:nspec){ # Loop over # of species 
    for(t in 1:nyears){ # Loop over t intervals (detection varies by year)
      alpha.lp[t,k] ~ dnorm(mu.alpha.lp[t], tau.alpha.lp) # Intercept
      # p intercept for t years, different mean, same variance
    }
    for(g in 1:2){ # Loop over # of coefficients (= number of covariates)
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
  for(g in 1:2){ # Loop over # of coefficients (= number of covariates)
    mu.beta.lp[g] ~ dnorm(0, 0.1)
    tau.beta.lp[g] <- pow(sd.beta.lp[g], -2)
    sd.beta.lp[g] ~ dnorm(0, 0.1)I(0,)
  }
  
  ######################
  ###Model Likelihood###
  ######################
  
  ###Ecological submodel: Define state conditional on parameters
  
  for (i in 1:nsites){ # Loop over # of sites
    for(k in 1:nspec){ # Loop over # of species in augmented list
      # Initial conditions of system (includes covariate effects)
      z[i,1, k] ~ dbern(psi1[i, k]) #Initial occupancy for year 1
      #Logit link for year 1 occupancy probability
      logit(psi1[i,k]) <- alpha.lpsi1[k] + 
      beta.lpsi1[1,k] * psi.covs[i,1] + beta.lpsi1[2,k] * psi.covs[i,2] 
      
      # State transitions (includes covariate effects)
      for (t in 2:nyears){ # Loop over # of years
        # Transition in occupancy probability between years due 
        # to prior colonization and persistence probabilities 
        z[i,t,k] ~ dbern(z[i,t-1,k]*phi[i,t-1,k] + (1-z[i,t-1, k])*gamma[i,t-1,k])
        logit(phi[i,t-1,k]) <- alpha.lphi[t-1,k] + #Prior colonization probability 
        beta.lphi[1,k] * phi.covs[i,t-1,1] + beta.lphi[2,k] * phi.covs[i,t-1,2] 
        logit(gamma[i,t-1,k]) <- alpha.lgamma[t-1,k] + #Prior persistence probability 
        beta.lgamma[1,k] * gamma.covs[i,t-1,1] + beta.lgamma[2,k] * gamma.covs[i,t-1,2] 
      }
    }
  }
  
  ###Observation model (incl. covariate effects)
  
  for (i in 1:nsites){
    for(k in 1:nspec){
      for (j in 1:nsurveys){
        for (t in 1:nyears){
          # Observed dataset as a function of occupancy and detection
          y[k,i,j,t] ~ dbern(z[i,t,k] * p[k,i,j,t]) 
          # Logit-link for detection with covariates
          logit(p[k,i,j,t]) <- alpha.lp[t,k] +
          beta.lp[1,k] * p.covs[i,j,t,1] + beta.lp[2,k] * p.covs[i,j,t,2] 
        }
      }
    }
  }
  
  ########################
  ###Derived Parameters###
  ########################
  
  # Number of occupied sites per year
  for(k in 1:nspec){
    for (t in 1:nyears){
      n.occ[t, k] <- sum(z[,t,k])
    }
  }
  # Species richness: total and per site/year
  for(i in 1:nsites){
    for(t in 1:nyears){
      Nspec[i,t] <- sum(z[i,t,]) # Species richness per site and year (this is the informative value!)
    }
  }
  
}#END OF MODEL CODE
")


###Initials/Parameters/MCMC Settings###
zst <- apply(y, c(2,4,1), max) # Observed occurrence as inits for z
inits <- function(){ list(z = zst)}

##Parameters monitored (could also add "z")
params <- c("mu.lpsi1", "sd.lpsi1", "mu.lphi", "sd.lphi", "mu.lgamma",
            "sd.lgamma", "mu.lp", "sd.lp", "psi1", "phi", "gamma", "p", "n.occ",
            "Nspec")

##MCMC settings
# na <- 1000 ; ni <- 6000 ; nt <- 6 ; nb <- 3000 ; nc <- 3
na <- 1000 ; ni <- 600 ; nt <- 1 ; nb <- 300 ; nc <- 3  # ~~~~ for testing, 7 mins


###Running Model

##Call JAGS (ART 188 min), check convergence and summarize posteriors
out <- jags(jags.data, inits, params, "DCM_Few_Covs.txt", n.adapt = na, n.chains = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = F)


####Separate Analysis to figure out where Rapid Demo points are, ignore#####
setwd(choose.dir())
data<-read.csv("SHARP_Rapid_Demo_Locs.csv")
points<-unique(data$PointID)
unique.points<-data[which(data$PointID==points[1])[1],c(5,12,13)]
for(i in 2:length(points)){
  unique.points<-rbind(unique.points,data[which(data$PointID==points[i])[1],c(5,12,13)])
}
write.csv(unique.points,file="SHARP_Rapid_Demo_Locs.csv")
