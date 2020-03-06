### Last editted 06/03/2020

#This code file was written by Ruth Kelly in April 2017 updated January and April 2018
# to derive demographic metrics for species with no seedbanks from annual matrices 
# in Compadre (i.e. annual periodicity of measurements = 1). 
# Matrices have already been checked for various criteria ,
# see details in "overview_code_25_07_2018.rmd"

### Thanks to Kevin Healy (https://github.com/healyke) for help with this code 

### The functions 'lifeTimeRepEvents' and 'makeLifeTable_mx0' were adapted from code 
# in the package 'Mage' (https://github.com/jonesor/compadreDB/tree/master/Mage see notes below)
#

# 
# The following metrics are derived for each matrix
# * Generation Time
# * Mean Progressive Growth (SSD weighted)
# * Mean Retrogressive Growth (SSD weighted)
# * Mean Reproduction
# * Net reproductive rate (SSD weighted)
# * Mature life expectancy
# * Mean life expectancy
# * Survival as mean life/life span
# * Gini Index for degree of iteroparity


#

### clear history

rm(list = ls())


########### required files 
### "annual_1yr_Seedbank_16_04_2018.RData"

#### load dataset ####

S1_data <- readRDS("annual_1yr_Seedbank_08_06_2018.RData")

## add libraries
# 
# library("popbio")
# library("popdemo")
# library("MASS")
# library("ineq")
# library("expm")



### a brief check to make sure data looks sensible

names(S1_data)

# [1] "metadata"  "mat" "matrixClass" "version"  


### assign each unit to an object

metadata_S1 <- S1_data$metadata
dim(metadata_S1)
# 102 52
length(unique(metadata_S1$SpeciesAccepted))
# 32
metadata_S1$UID1 <- seq(1, nrow(metadata_S1),1)

mat_S1 <- S1_data$mat
#length(mat_S1)
#mat_S1[[1]]

matrixClass_S1 <- S1_data$matrixClass
length(matrixClass_S1)
matrixClass_S1[[1]]

#######################
S1_data$version

### check for survivalIssue

metadata_S1$SurvivalIssue[which(metadata_S1$SurvivalIssue > 1)]
### 1
### This equals 1 . very odd. 

length(which(metadata_S1$SurvivalIssue == 1))
#30

#### ---- Calculate survival curves ----
# 

# The 'makeLifeTable_mx0' function is adapted from 'makelifetable' from the 
# package 'Mage' to calculate "survival probabilty to time x or % still alive at time x"
# and "mx (reproduction per capita at stage x)."  
# The ammendment to this function means that it now also calculates fecundity in year
# 0 rather than assuming no reproduction in this time step.  
#  Time step's in our matrices = 1 year throughout. 

### 
source("makelifetable_mx0_function_18_04_2018.R")  

#### store results in : 
lxmx_curve_S1 <- list()

### Calculate from first above-ground stage. 
startLife <- 2

### 
for(i in 1:nrow(metadata_S1)) {
  lxmx_curve_S1[[i]] <-  makeLifeTable_mx0(matU =  mat_S1[[i]]$matU, matF =  mat_S1[[i]]$matF, 
                                       startLife = startLife, nSteps = 10000)
}

#### ---- Calculate life span ---- 


### Calculate from this how long things actually live by seeing how many are 
## still alive at time point x.  
## The function 'exceptionalLife', by Kevin Healy does this.
### I've modified the function so the arbitary maximum life is 10000, instead of 100
source("exceptionalLife_function_KHealy_RK_2018.R")


### Here, I calculate the age at which only 0.1% remain alive and treat this as 99.9% lifespan, 
### 99% dead, 95% dead and the age at which 50% are dead and use this as mean lifespan. 

S_dead <- list()
dead_999 <- c()
dead_99 <- c()
dead_95 <- c()
dead_50 <- c()


for(i in 1:nrow(metadata_S1)) {
  S_dead[[i]] <- exceptionalLife(mat_S1[[i]]$matU)
  dead_999[i]<- S_dead[[i]][[4]]
  dead_99[i] <- S_dead[[i]][[3]]
  dead_95[i] <- S_dead[[i]][[2]]
  dead_50[i] <- S_dead[[i]][[1]]
}



# warnings()

which(is.na(dead_50))
# integer(0)

which(is.na(dead_95))
# integer(0)

which(is.na(dead_99))
#[1] 80 82

which(is.na(dead_999))
# [1] 80 82

metadata_S1$SurvivalIssue[which(is.na(dead_99))]
# [1] 1 1

metadata_S1$SpeciesAccepted[which(is.na(dead_999))]
#[1] "Sapium sebiferum" "Sapium sebiferum"



### these issues arise for "Sapium sebiferum" where survival = 1 in the last stage
## i.e. death was not observed in the oldest plants. 
###  


#### ---- Life time reproductive events ----

### Age at maturity - Calculated using 'lifeTimeRepEvents' function adapted from 
#  package "Mage". This version ignores clonality matrices as there are none
# in this dataset. 

# ### Calculating 'Mature Life Expectancy' using the function lifeTimeRepevents. 
# 
# Here we use a slightly older version of the function lifeTimeRepEvents  
## from the Mage package by Owen Jones which calculates - 
#   * pRep: probability of achiving sexual maturity
# * 'La': mean age at maturity (in the same units as the matrix population model).
# * 'meanLifeExpectancy': mean life expectancy conditional on entering the
# life cycle in the first reproductive stage
# * 'remainingMatureLifeExpectancy': Life expectancy from mean maturity. This
# is mean life expectancy - mean age at maturity ('La' above). This value can
# be negative because both mean life expectancy and mean age at maturity are
# means of their respective distributions. 

##### We require - 'La': mean age at maturity and 'meanLifeExpectancy' conditional on reaching maturity


lifeTimeRepEvents <- function(matU, matF, startLife = 1){
  #Function to determine probability of reaching reproduction, age at 
  #maturity and reproductive lifespan (Code adapted from H. Caswell's 
  #matlab code, and Morris & Doak):
  
  if(missing(matU)){stop('matU missing')}
  # if(missing(matF) & missing(matC)){stop('matF or matC missing. You must provide at least one')}
  # if(sum(matF + matC,na.rm=T)==0){stop('matF and matC contains only 0 values')}

    
    matDim <- dim(matU)[1]
    surv <- colSums(matU)
    
    out = NULL
    
    if (sum(matF)>0){
      
      fecLifeStages <- colSums(matF)
      fecLifeStages[which(fecLifeStages>0)] <- 1
      
      #Probability of survival to first sexual reprod event
      Uprime <- matU
      Uprime[ , which(fecLifeStages == 1)] <- 0
      Mprime <- matrix(0, 2, matDim)
      for (p in 1:matDim) {
        if (fecLifeStages[p] == 1){ Mprime[2,p] <- 1 }
        else{Mprime[1, p] = 1 - surv[p]}
      }
      Bprime <- Mprime %*% (MASS::ginv(diag(matDim) - Uprime))
      out$pFec <- Bprime[2, startLife]
      
      #Age at first sexual reproduction (LaFec; Caswell 2001, p 124)
      D <- diag(c(Bprime[2,]))
      Uprimecond <- D %*% Uprime %*% MASS::ginv(D)
      expTimeReprod <- colSums(MASS::ginv(diag(matDim) - Uprimecond))
      out$LaFec <- LaFec <- expTimeReprod[startLife]
      
      #Mean life expectancy conditional on entering the life cycle in the first reproductive stage
      firstFecLifeStage <- min(which(fecLifeStages == 1))
      N <- solve(diag(matDim[1]) - matU)
      out$meanLifeExpectancyFec <- colSums(N)[firstFecLifeStage]
      
    
      
      
      #Life expectancy from mean maturity
      out$remainingMatureLifeExpectancyFec <- colSums(N)[startLife] - LaFec
    }

  return(out)
}

La1 <- c()
RepLifeExp <- c()


for(i in 1:nrow(metadata_S1)) {
  tryCatch({
    La1[i] <- lifeTimeRepEvents(mat_S1[[i]]$matU, mat_S1[[i]]$matF, startLife = startLife)[[2]]
    RepLifeExp[i] <- lifeTimeRepEvents(mat_S1[[i]]$matU, mat_S1[[i]]$matF, startLife = startLife)[[3]]
  }
  , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


summary(La1)


summary(RepLifeExp)


which(is.na(La1))

## same as previous issues
# metadata_S1[which(is.na(La1)),] # Sapium_sebiferum

######

#### ---- Calculate survival index 3 ways!----

# Calculate survival parameter based 
# on mean life expectancy/life span (99%)
# mean life expectancy/life span (99.9%)
# mean life expectancy/(mean age at reproduction + length of reproductive life)
# This parameter varies between close to 0 (but 0 is impossible) meaning average life expectancy is very short compared to life span (i.e. a lot of juvenile mortality), and 1 mean life expectancy = life span, (i.e. all individuals life for the same length of time).

Survival_99 <- dead_50/dead_99
Survival_999 <- dead_50/dead_999
Survival_LaRepLife <- dead_50/(La1 + RepLifeExp)
## summary(Survival_NS)

#### ---- Net reproductive rate ----

### 

net_rep <- c()


for(i in 1:nrow(metadata_S1)) {
  net_rep[i] <-  net.reproductive.rate(mat_S1[[i]]$matA, r = mat_S1[[i]]$matF)
}


summary(net_rep)


### 

which(is.na(net_rep))
# [1] 82
### same issue as previous.. 

#### ---- mean annual per capita reproductive rates ---- 

#### mean annual reproductive rate.. weighting seed by their probability of 
## becoming seedlings

source("pSeedling_function_15_06_2017.r")

mean_rep <- c()

startLife <- 2

for(i in 1:nrow(metadata_S1)) {
  
  ## use function pSeedling to calculate prop survival to seedling stage, from first stage 
  
  p_germ <- pSeedling(mat_S1[[i]]$matU, seed1 = 1, seedling1 = 2)
  
  if(p_germ > 1.0001){stop("Error germination rate > 1")}
  
  ### Adjust first row of the F matrix by the probability of surviving to seedling ##stage from seedbank 
  
  matFtemp <- mat_S1[[i]]$matF
  matFtemp[1,] <- mat_S1[[i]]$matF[1,]*p_germ
  
  ##Then sum up the reproductive rates across the F matrix
  
  repo_sum <- vector()
  rep1 <- nrow(matFtemp)
  
  for(j in 1:rep1){
    repo_sum[j] <- sum(matFtemp[,j])
  }
  
  #### Calculate stable stage distribution based on full matA
  SSD <- eigen.analysis(mat_S1[[i]]$matA)$stable.stage
  
  ### cut1 is an automatic trimmer for SSD based on startLife
  cut1 <- seq(1:(startLife-1))
  
  ## reweight SSD it so it is the stable stage as if the seeds had not been measured.  
  SSD2 <- SSD[startLife:length(SSD)]/(1-sum(SSD[cut1]))
  
  
  ##Then weight the reproduction from non-seed stages, 
  ## (i.e. weight repo_sum[-1]) by the stable stage distribtion of 
  ## germinated stages
  mean_rep[i] <- repo_sum[startLife:length(repo_sum)] %*% SSD2
  
  
}

summary(mean_rep)


######
###### ---- Generation time ---- 

#Generation time

#Basically some type of measure of how long it takes for a cohort to relace itself.
#Turns out there is a bunch of ways to do this. I used the popbio package which 
#calculates the time it takes for a pop to grow by a factor of R0 the net reprodictive 
#rate (how long it would for everyone to replace themselves) and is calulated  
#log(net_repo_rate)/log(population growth rate)


gen_time <- c()

for(i in 1:nrow(metadata_S1)) {
  
  
  gen_time[i] <- generation.time(mat_S1[[i]]$matA, r = mat_S1[[i]]$matF)
  
}

summary(gen_time)



which(is.na(gen_time))


#### ---- Progressive growth (SSD weighted) ----
######

#A measure of how quickly you move through the growth stages of the transition matrix,
# excluding seed to seedling as this is not comparable. 

#length of the matrix

## Here I am using matU directly, this is okay because I have already removed matrices where clonality was measured.  i.e. matC == 0

prog_growth <- c()

for(i in 1:nrow(metadata_S1)) {
  
  #Progression
  ## Put zeros in the upper triangle, leaving only the growth transitions in the matrix U
  ## 
  
  matU <- mat_S1[[i]]$matU
  
  
  matU[upper.tri(matU,diag = TRUE )] <- c(0)
  ###sum up the proggression from each stage
  

  
  prog_sum <- vector()
  
  for(j in 1:nrow(matU)){
    prog_sum[j] <- sum(matU[,j])
  }
  
  ### Here keep only seedling stages and onwards
  
  prog_sum <- prog_sum[startLife:length(prog_sum)]
  
  
  SSD <- eigen.analysis(mat_S1[[i]]$matA)$stable.stage
  ### now cut SSD and reweight it so it is the stable stage if seeds were not 
  ## measured.  
  
  ### cut1 is an automatic trimmer for SSD based on startLife
  cut1 <- seq(1:(startLife-1))
  
  ## reweight it so it is the stable stage if seeds were not measured.  
  SSD2 <- SSD[startLife:length(SSD)]/(1-sum(SSD[cut1]))
  
  ##and weight it against the stable state distribution
  prog_growth[i] <- prog_sum %*% SSD2
  
} 


summary(prog_growth) 


### Zeros here are annual species, i.e. they can never progress upwards in an annual
## time step from when they are above ground already. 


### Calculating retrogression (SSD weighted)

#Retrogression, calculated as per progression, but we blank out the upper triangle 
# instead of the lower.  Again, I use matU directly as clonal matrices are zero in 
# this dataset.


#### ---- Retrogressive growth (SSD weighted) ----

startLife = 2 

retro_growth <- c()


for(i in 1:nrow(metadata_S1)) {
  
  #Progression
  ## Put zeros in the upper triangle, leaving only the growth transitions in the matrix U
  ## 
  
  matU <- mat_S1[[i]]$matU
  
  matU[lower.tri(matU,diag = TRUE )] <- c(0)
  ###sum up the progression from each stage
  
  retro_sum <- vector()
  
  for(j in 1:nrow(matU)){
    retro_sum[j] <- sum(matU[,j])
  }
  
  
  ### Here keep only seedling stages and onwards
  
  retro_sum <- retro_sum[startLife:length(retro_sum)]
  
  
  SSD <- eigen.analysis(mat_S1[[i]]$matA)$stable.stage
  ### now cut SSD and reweight it so it is the stable stage if seeds were not 
  ## measured.  
  
  ### cut1 is an automatic trimmer for SSD based on startLife
  cut1 <- seq(1:(startLife-1))
  
  ## reweight it so it is the stable stage if seeds were not measured.  
  SSD2 <- SSD[startLife:length(SSD)]/(1-sum(SSD[cut1]))
  
  ##and weight it against the stable state distribution
  retro_growth[i] <- retro_sum %*% SSD2
  
} 


summary(retro_growth) 

### Note from looking at the papers it seems that a lot of the 'retrogression' in this database
## is actually caused by grazing. 


#### ---- "Gini index" ----  

##Variation in fecundity accross lifespan.
#### Here with 2 versions one across whole lifespan - Gini_whole_life
#### and one across the reproductive life, i.e. after the start of reproduction

## Here we need mx from the lxmx curve, i.e. the reproduction per capita at stage x 
## We then chop this to the lifespan of 99% and 99.9% lifespan, 
## and calculate the Gini index using the function Gini in the package 'ineq'
## corr = TRUE corrects for small samples (i.e. sort lifespans)


### need to use Trycatch to stop the errors here, because NS_dead_99 = NA from some species, primarily due to survival issues ...  Also, NA's are arising in the Gini index where no reproduction is happening within the 99% lifespan. 


GiniF_life_99 <- c()
for(i in 1:nrow(metadata_S1)) {
  tryCatch({
    repro1 <- lxmx_curve_S1[[i]]$mx[1:(dead_99[i])]
    
    GiniF_life_99[i] <- Gini(repro1, corr = TRUE)
    
  }
  
  , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

summary(GiniF_life_99)
# 


#### calculate with lifespan = 99.9%
GiniF_life_999 <- c()
for(i in 1:nrow(metadata_S1)) {
  tryCatch({
    repro1 <- lxmx_curve_S1[[i]]$mx[1:(dead_999[i])]
    
    GiniF_life_999[i] <- Gini(repro1, corr = TRUE)
    
  }
  
  , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

summary(GiniF_life_999)

GiniF_life_99[which(metadata_S1$OrganismType == "Annual")]
### always 1 - which means totally unequal. 

### Many of these are annual species.  Where there is no survival into the second year 
### and therefore mx is zero throughout.  These will be set to Gini = 0 in the in 
### error checking as they have no variation. 

### GiniF_repro = Gini accross life starting at average age at maturity (La1) (rounded down)
### 

GiniF_repro <- c()
for(i in 1:nrow(metadata_S1)) {
  tryCatch({
    repro1 <- lxmx_curve_S1[[i]]$mx[floor(La1[i]):(dead_999[i])]
    
    GiniF_repro[i] <- Gini(repro1, corr = TRUE)
    
    ### this if statement sets Gini value to NA where the La1 estimate
    ## is higher than the dead999 estimate. Otherwise these are counted in reverse order
    ## in the above steps ie 8:5 = 8,7,6,5
    if(floor(La1[i]) > dead_999[i]){
      GiniF_repro[i] <- NA
    }
    
  }
  
  , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

summary(GiniF_repro)

setdiff(which(is.na(GiniF_life_99)), which(is.na(GiniF_repro)))
## integer(0)


#######################


#### ---- now join everything into one dataset and export ----

demo_data_S1 <- cbind(metadata_S1, mean_rep, net_rep, dead_50, dead_99, dead_999, 
                      Survival_99, Survival_999,Survival_LaRepLife, prog_growth,
                      retro_growth, RepLifeExp, La1, 
                      GiniF_repro,  GiniF_life_99, GiniF_life_999, gen_time )

demo_data_S1$NSeedStages <- 1
 
ncol(demo_data_S1)
# # [1] 52
# ncol(demo_data_S1)
# # [1] 72


names(demo_data_S1)[53:68]  <- c("Mean_repro", "Net_repro", "Age_50_dead", "Age_99_dead", "Age_999_dead",
                                 "survival_index99", "Survival_index999", "Survival_index_La_RepExp",
                                 "Growth", "Retrogression", 
                                 "Repro_life_expectancy", "Age_at_Maturity", "GiniF_repro",
                                 "Gini_life99", "Gini_life999","Gen_time")

write.csv(demo_data_S1, "demography_1yr_seedbank_08_06_2018.csv")

###########################################################

