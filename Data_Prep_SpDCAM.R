###SHARP Data Preparation for Spatial DCAM###

##############################
###Set WD and Load Packages###
##############################

setwd("C:/Users/12152/Documents/UConn_Files/Project_Files/Data_Analysis/Survey_Data_Analyses")

##################
###Read in Data###
##################

survey.data<-read.csv("SHARP_Surveys_All_Years.csv")

veg.data<-read.csv("SHARP_Rapid_Veg_All_Years.csv")

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
survey.data<-survey.data[,colnames(survey.data)%in%c("AlphaCode","TotalCount", "PointID", "Point_X", "Point_Y", "VisitNum","Year","Noise")]
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

#Fill count array (code is ugly but works fine, takes a while to run though!)
y<-array(0,c(nspec,nsites,nsurveys,nyears))
for(i in 1:nrow(survey.data)){
  y[grep(survey.data$AlphaCode[i],species), grep(survey.data$PointID[i],sites), 
    survey.data$VisitNum[i], grep(survey.data$Year[i],years)]<-survey.data$TotalCount[i]
}
#save(y,file = "y_count_array_SHARP_DCAM.R")
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
#Noise for p, Vegetation data informs other params

noise<-array(NA,c(nsites,nsurveys,nyears))
dimnames(noise)<-list(sites,c(1:3),years) 
for(i in 1:nrow(survey.data)){
  noise[grep(survey.data$PointID[i],sites), survey.data$VisitNum[i],
        grep(survey.data$Year[i],years)]<-survey.data$Noise[i]
}

#Veg (using % cover of five most ubiquitous marsh species now as covariates)
#START FROM HERE!!
#Making a year column in the matrix
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
veg.data$Year<-substrRight(as.character(veg.data$Date),4)
veg.data$Year<-as.numeric(veg.data$Year)

#Determining 5 most ubiquitous marsh species
top.5<-names(summary(veg.data$DomSp1)[1:5])
#Making veg array to store % cover of each of 5 dominant species
sites<-unique(veg.data$BirdPtID) #all unique sites
years<-sort(unique(veg.data$Year)) #all sampled years (10)
nsites<-length(sites)
nyears<-length(years)
veg<-array(NA,c(nsites,nyears,5))
dimnames(veg)<-list(sites,years,top.5) 
####Code is ugly, can likely streamline with a function but having trouble writing...####
for(i in 1:nrow(veg.data)){
  if(veg.data$DomSp1[i]%in%top.5==TRUE){
    veg[grep(veg.data$BirdPtID[i],sites), 
             grep(veg.data$Year[i],years), 
             grep(veg.data$DomSp1[i],top.5)]<-veg.data$Sp1Percent[i]
  }else {}
}

for(i in 1:nrow(veg.data)){
  if(veg.data$DomSp2[i]%in%top.5==TRUE){
    veg[grep(veg.data$BirdPtID[i],sites), 
        grep(veg.data$Year[i],years), 
        grep(veg.data$DomSp2[i],top.5)]<-veg.data$Sp2Percent[i]
  }else {}
}

for(i in 1:nrow(veg.data)){
  if(veg.data$DomSp3[i]%in%top.5==TRUE){
    veg[grep(veg.data$BirdPtID[i],sites), 
        grep(veg.data$Year[i],years), 
        grep(veg.data$DomSp3[i],top.5)]<-veg.data$Sp3Percent[i]
  }else {}
}

for(i in 1:nrow(veg.data)){
  if(veg.data$DomSp4[i]%in%top.5==TRUE){
    veg[grep(veg.data$BirdPtID[i],sites), 
        grep(veg.data$Year[i],years), 
        grep(veg.data$DomSp4[i],top.5)]<-veg.data$Sp4Percent[i]
  }else {}
}


for(i in 1:nrow(veg.data)){
  if(veg.data$DomSp5[i]%in%top.5==TRUE){
    veg[grep(veg.data$BirdPtID[i],sites), 
        grep(veg.data$Year[i],years), 
        grep(veg.data$DomSp5[i],top.5)]<-veg.data$Sp5Percent[i]
  }else {}
}


for(i in 1:nrow(veg.data)){
  if(veg.data$DomSp6[i]%in%top.5==TRUE){
    veg[grep(veg.data$BirdPtID[i],sites), 
        grep(veg.data$Year[i],years), 
        grep(veg.data$DomSp6[i],top.5)]<-veg.data$Sp6Percent[i]
  }else {}
}


for(i in 1:nrow(veg.data)){
  if(veg.data$DomSp7[i]%in%top.5==TRUE){
    veg[grep(veg.data$BirdPtID[i],sites), 
        grep(veg.data$Year[i],years), 
        grep(veg.data$DomSp7[i],top.5)]<-veg.data$Sp7Percent[i]
  }else {}
}


for(i in 1:nrow(veg.data)){
  if(veg.data$DomSp8[i]%in%top.5==TRUE){
    veg[grep(veg.data$BirdPtID[i],sites), 
        grep(veg.data$Year[i],years), 
        grep(veg.data$DomSp8[i],top.5)]<-veg.data$Sp8Percent[i]
  }else {}
}


for(i in 1:nrow(veg.data)){
  if(veg.data$DomSp9[i]%in%top.5==TRUE){
    veg[grep(veg.data$BirdPtID[i],sites), 
        grep(veg.data$Year[i],years), 
        grep(veg.data$DomSp9[i],top.5)]<-veg.data$Sp9Percent[i]
  }else {}
}


for(i in 1:nrow(veg.data)){
  if(veg.data$DomSp10[i]%in%top.5==TRUE){
    veg[grep(veg.data$BirdPtID[i],sites), 
        grep(veg.data$Year[i],years), 
        grep(veg.data$DomSp10[i],top.5)]<-veg.data$Sp10Percent[i]
  }else {}
}
veg[is.na(veg)]<-0
####End of Veg fill section####

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
