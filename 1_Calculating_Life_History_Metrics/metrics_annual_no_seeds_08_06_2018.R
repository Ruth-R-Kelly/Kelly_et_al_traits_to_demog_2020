### Last editted 06/03/2020

#This code file was written by Ruth Kelly in April 2017 
# to derive demographic metrics for species with no seedbanks from annual matrices 
# in Compadre (i.e. annual periodicity of measurements = 1). 
# Matrices have already been checked for various criteria ,
# see details in "overview_code_06_06_2018.rmd"

### Thanks to Kevin Healy (https://github.com/healyke) for help with this code. 

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

########### required files 
### "annual_no_Seedbank_16_04_2018.RData"

### delete background files

rm(list = ls())

#### load dataset ####

NS_data <- readRDS("annual_no_Seedbank_08_06_2018.RData")

## add libraries
# 
# library("popbio")
# library("popdemo")
# library("expm")
# library("MASS")
# library("ineq")


### a brief check to make sure data looks sensible

names(NS_data)

# [1] "metadata"    "mat"         "matrixClass" "version"  


### assign each unit to an object

metadata_NS <- NS_data$metadata
dim(metadata_NS)
# 137 52

length(unique(metadata_NS$SpeciesAccepted))
#59

metadata_NS$UID1 <- seq(1, nrow(metadata_NS),1)

mat_NS <- NS_data$mat
#length(mat_NS)
#mat_NS[[1]]

matrixClass_NS <- NS_data$matrixClass

#######################
NS_data$version


#### Make a life-history table for each population matrix #### 
# 

# The 'makeLifeTable_mx0' function is adapted from 'makelifetable' from the 
# package 'Mage' to calculate "survival probabilty to time x or % still alive at time x"
# and "mx (reproduction per capita at stage x)."  
# The ammendment to this function means that it now also calculates fecundity in year
# 0 rather than assuming no reproduction in this time step.  
#  Time step's in our matrices = 1 year throughout. 

### 
source("makelifetable_mx0_function_18_04_2018.R")  

#### ---- Calculate survival curves ----

#### store results in : 
lxmx_curve_NS <- list()

### 
for(i in 1:nrow(metadata_NS)) {
  lxmx_curve_NS[[i]] <-  makeLifeTable_mx0(matU =  mat_NS[[i]]$matU, matF =  mat_NS[[i]]$matF, 
                                       startLife = 1, nSteps = 10000)
}


#### ---- Calculate life span ----

### Calculate from this how long things actually live by seeing how many are 
## still alive at time point x.  


## The function 'exceptionalLife', by Kevin Healy does this.
### I've modified the function so the arbitary maximum life is 10000, instead of 100
source("exceptionalLife_function_KHealy_RK_2018.R")


### Here, I calculate the age at which only 0.1% remain alive and treat this as 99.9% lifespan, 
### 99% dead, 95% dead and the age at which 50% are dead and use this as mean lifespan. 

NS_dead <- list()
NS_dead_999 <- c()
NS_dead_99 <- c()
NS_dead_95 <- c()
NS_dead_50 <- c()


for(i in 1:nrow(metadata_NS)) {
  NS_dead[[i]] <- exceptionalLife(mat_NS[[i]]$matU)
  NS_dead_999[i]<- NS_dead[[i]][[4]]
  NS_dead_99[i] <- NS_dead[[i]][[3]]
  NS_dead_95[i] <- NS_dead[[i]][[2]]
  NS_dead_50[i] <- NS_dead[[i]][[1]]
}


which(is.na(NS_dead_50))

summary(NS_dead_50)

which(is.na(NS_dead_95))

summary(NS_dead_95)


which(is.na(NS_dead_99))


unique(metadata_NS$SpeciesAccepted[which(is.na(NS_dead_99))])
#[1] "Molinia caerulea" "Primula elatior"  "Clidemia hirta"  
metadata_NS$SurvivalIssue[which(is.na(NS_dead_99))]
## [1] 1 1 1 1 1

summary(NS_dead_999)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    5.00   41.75  111.00  347.20  404.00 6981.00       5 

which(is.na(NS_dead_999))
# [1] 13 15 16 25 62


#### These species have a survival issue of 1 in the last stage. 
## I.e. no death was observed in the last stage. 
###  


#### ---- Life time reproductive events ----

### Age at maturity - Calculated using 'lifeTimeRepEvents' function adapted from 
#  package "Mage". This version ignores clonality matrices as there are none
# in this dataset. 


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


### La1 is age at first reproduction
La1 <- c()

### RepLifeExpectancy is lifespan from mean maturity
RepLifeExp <- c()


for(i in 1:nrow(metadata_NS)) {
  tryCatch({
    La1[i] <- lifeTimeRepEvents(mat_NS[[i]]$matU, mat_NS[[i]]$matF, startLife = 1)[[2]]
    RepLifeExp[i] <- lifeTimeRepEvents(mat_NS[[i]]$matU, mat_NS[[i]]$matF, startLife = 1)[[3]]
  }
  , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


summary(La1)


summary(RepLifeExp)

metadata_NS$SurvivalIssue[which(is.na(La1))]
# [1] 1 1 1 1

## mat_NS[which(is.na(La1))]

### same matrices that have issues for lifespan

#### ---- Calculate survival index 3 ways!----

# Calculate survival parameter based 
    # on mean life expectancy/life span (99%)
    # mean life expectancy/life span (99.9%)
    # mean life expectancy/(mean age at reproduction + length of reproductive life)
# This parameter varies between close to 0 (but 0 is impossible) meaning average life expectancy is very short compared to life span (i.e. a lot of juvenile mortality), and 1 mean life expectancy = life span, (i.e. all individuals life for the same length of time).

Survival_NS_99 <- NS_dead_50/NS_dead_99
Survival_NS_999 <- NS_dead_50/NS_dead_999
Survival_NS_LaRepLife <- NS_dead_50/(La1 + RepLifeExp)
## summary(Survival_NS)



#### ---- Net reproductive rate ----

net_rep <- c()


for(i in 1:nrow(metadata_NS)) {
  net_rep[i] <-  net.reproductive.rate(mat_NS[[i]]$matA, r = mat_NS[[i]]$matF)
}


summary(net_rep)

### 

### which(is.na(net_rep))

#### same species as previous.  


#### ---- mean annual per capita reproductive rates ---- 

### Note count is of new seedlings per capita at stable stage distribution

mean_rep <- c()

for(i in 1:nrow(metadata_NS)) {
  
  ##first sum up the reproductive rates across the F matrix
  N <- nrow(mat_NS[[i]]$matA)
  
  repo_sum <- vector()
  
  for(j in 1:(N)){
    repo_sum[j] <- sum(mat_NS[[i]]$matF[,j])
  }
  
  SSD <- eigen.analysis(mat_NS[[i]]$matA)$stable.stage
  ##then weight it against the stable state distribtion
  mean_rep[i] <- repo_sum %*% SSD
  
}

summary(mean_rep)


###### ---- Generation time ---- 

#A measure of how long it takes for a cohort to relace itself.
#Turns out there is a bunch of ways to do this. I used the popbio package which 
#calculates the time it takes for a pop to grow by a factor of R0 the net reprodictive 
#rate (how long it would for everyone to replace themselves) and is calulated  
#log(net_repo_rate)/log(population growth rate)


gen_time <- c()

for(i in 1:nrow(metadata_NS)) {
  
  
  gen_time[i] <- generation.time(mat_NS[[i]]$matA, r = mat_NS[[i]]$matF)
  
}

summary(gen_time)


#### ---- Calculating progressive growth (SSD weighted) ----

#A measure of how quickly you move through the growth stages of the transition matrix. 
## Here I am using matU directly, this is okay because I have already removed matrices
## where clonality was measured.  i.e. matC == 0

prog_growth <- c()

for(i in 1:nrow(metadata_NS)) {
  
  #Progression
  ## Put zeros in the upper triangle, leaving only the growth transitions in the matrix U
  ## 
  
  matU <- mat_NS[[i]]$matU
  
  matU[upper.tri(matU,diag = TRUE )] <- c(0)
  ###sum up the proggression from each stage
  
  prog_sum <- vector()
  
  for(j in 1:nrow(matU)){
    prog_sum[j] <- sum(matU[,j])
  }
  
  SSD <- eigen.analysis(mat_NS[[i]]$matA)$stable.stage
  ##and weight it against the stable state distribution
  prog_growth[i] <- prog_sum %*% SSD
  
} 

summary(prog_growth)

#### ---- Calculating retrogression (SSD weighted) ---- 

#Retrogression, calculated as per progression, but we blank out the upper triangle 
# instead of the lower.  Again, I use matU directly as clonal matrices are zero in 
# this dataset.


retro_growth <- c()


for(i in 1:nrow(metadata_NS)) {
  
  #Progression
  ## Put zeros in the upper triangle, leaving only the growth transitions in the matrix U
  ## 
  
  matU <- mat_NS[[i]]$matU
  
  matU[lower.tri(matU,diag = TRUE )] <- c(0)
  ###sum up the progression from each stage
  
  retro_sum <- vector()
  
  for(j in 1:nrow(matU)){
    retro_sum[j] <- sum(matU[,j])
  }
  
  SSD <- eigen.analysis(mat_NS[[i]]$matA)$stable.stage
  ##and weight it against the stable state distribution
  retro_growth[i] <- retro_sum %*% SSD
  
} 

summary(retro_growth) 


#### ---- "Gini index" ---- 

##Variation in fecundity accross lifespan.
#### Here with 2 versions one across whole lifespan - Gini_whole_life
#### and one across the reproductive life, i.e. after the start of reproduction

## Here we need mx from the lxmx curve, i.e. the reproduction per capita at stage x 
## We then chop this to the lifespan of 99% or 99.9% of the population + 1 year to account for 
## repreduction in the year of death (e.g. annuals etc), 
## and calculate the Gini index using the function Gini in the package 'ineq'
## corr = TRUE corrects for small samples (i.e. sort lifespans)


### need to use Trycatch to stop the errors here, because NS_dead_99 = NA from some species, primarily due to survival issues ...  Also, NA's are arising in the Gini index where no reproduction is happening within the 99% lifespan. 

startlife <- 1

GiniF_life_99 <- c()

for(i in 1:nrow(metadata_NS)) {
  tryCatch({
    repro1 <- lxmx_curve_NS[[i]]$mx[startlife:(NS_dead_99[i])]
    
    GiniF_life_99[i] <- Gini(repro1, corr = TRUE)
    
  }
  
  , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

## for life span = 99.9%
GiniF_life_999 <- c()
for(i in 1:nrow(metadata_NS)) {
  tryCatch({
    repro1 <- lxmx_curve_NS[[i]]$mx[startlife:(NS_dead_999[i])]
    
    GiniF_life_999[i] <- Gini(repro1, corr = TRUE)
    
  }
  
  , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


### GiniF_repro = Gini accross life starting at average age at maturity (La1) ending at
### end of average reproductive lifespan La1 + RepLife
### 

GiniF_repro <- c()
for(i in 1:nrow(metadata_NS)) {
  tryCatch({
 
    repro1 <- lxmx_curve_NS[[i]]$mx[floor(La1[i]):(NS_dead_999[i])]
    
    GiniF_repro[i] <- Gini(repro1, corr = TRUE)

    ### this if statement sets Gini value to NA where the La1 estimate
    ## is higher than the dead999 estimate. Otherwise these are counted in reverse order
    ## in the above steps ie 8:5 = 8,7,6,5
    if(floor(La1[i]) > NS_dead_999[i]){
      GiniF_repro[i] <- NA
    }
    
    
        
  }

  , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

summary(GiniF_repro)


## Note: correlation between the two metrics is very high. 

# cor.test(GiniF_life_999, GiniF_repro)

#### now join everything into one dataset and export

demo_data_NS <- cbind(metadata_NS, mean_rep, net_rep, NS_dead_50, NS_dead_99, NS_dead_999, 
                      Survival_NS_99, Survival_NS_999,Survival_NS_LaRepLife, prog_growth,
                      retro_growth, RepLifeExp, La1, 
                      GiniF_repro, GiniF_life_99, GiniF_life_999,gen_time)
names(demo_data_NS)
demo_data_NS$NSeedStages <- 0

dim(metadata_NS)
dim(demo_data_NS)

names(demo_data_NS)[53:69] <- c("Mean_repro", "Net_repro", "Age_50_dead", "Age_99_dead", "Age_999_dead",
                                 "survival_index99", "Survival_index999", "Survival_index_La_RepExp",
                                 "Growth", "Retrogression", 
                                 "Repro_life_expectancy", "Age_at_Maturity", "GiniF_repro",
                                 "Gini_life99", "Gini_life999",  "Gen_time", "NSeedStages")


write.csv(demo_data_NS, "demography_non_seedbank_08_06_2018.csv")

###########################################################
