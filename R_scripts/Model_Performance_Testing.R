###############################
###Model Performance Testing###
###############################

#################
###Description###
#################

#This code is used to test the accuracy and precision of competing models used
#to analyze variation in the composition of a tidal marsh breeding bird community 
#through time. Running the code requires multiple scripts located within the 
#"Tidal_Marsh_Bird_SDMs" GitHub repository. Competing models are compared using the
#same simulated data with changes in each focal parameter to observe model sensitivity
#to parameter inputs. Each section in this code examines two competing models. Note
#that no environmental model covariates are incorporated into any model beyond what is 
#required for model specification (i.e. distances of observed birds are used to specify
#a distance sampling model, which affects detection probability, but noise, an environmental
#covariate which can also impact detection, is not incorporated). Model specifications
#and descriptions can be found in the R code for each model, as well as the code used
#to simulate data for each model. Code relies heavily on code presented in 
#Applied Hierarchical Modeling Volumes 1 (2016) and 2 (2020) by Marc Kery and
#Andy Royle, as well as their supporting R package 'AHMbook'. 


############################
###Load Required Packages###
############################
library(AHMbook)
library(jagsUI)
source("~functions/homebrewed_simulation_functions.R") #Load simulation function (requires AHMbook)



