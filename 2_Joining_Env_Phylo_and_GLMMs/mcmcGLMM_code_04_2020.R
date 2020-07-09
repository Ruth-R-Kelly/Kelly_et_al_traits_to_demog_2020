### MCMCglmm code rewritten 18/06/2018 
## Datasets include 99.9 lifespan calculations and averaging of populations within 5km to 
## match the scale of the climate data. 

### load libraries

### If you run this code without running the overview file
# "overview_analysis_07_04_2020.Rmd" you will need to load these libraries here. 
# library(MCMCglmm)
# library(MCMCpack)
# library(ape)
# library("caper")
# library("phangorn")
# library("vegan")
# #######

### clear previous files for safety :)
rm(list = ls())

### read in plant trait data in order of phylo 
plant_data <- read.csv("data_with_needed_env_vars_03_2020.csv", row.names = "X")

# General data checking. 
# str(plant_data)
# head(plant_data)
# row.names(plant_data)
# 
# nrow(plant_data)
# summary(plant_data)

####
phy_tree <- read.nexus("reduced_tree_for_traits_03_2020.nex")  

plot(phy_tree)


## This is an ultrametric tree but R isn't recognising it as such
# because of pre-processing therefore we make it ultrametric following: 
#http://blog.phytools.org/2016/08/fixing-ultrametric-tree-whose-edges-are.html
phy_tree2 <- nnls.tree(cophenetic(phy_tree),phy_tree,rooted=TRUE)

#### check trees similarity
# 
 tips<-phy_tree2$tip.label
 cor(as.vector(cophenetic(phy_tree)[tips,tips]),
     as.vector(cophenetic(phy_tree2)[tips,tips]))

### 1

### check species names align with phylogeny
phy_tree2$tip.label == row.names(plant_data)

names(plant_data)

#### assemble the trait data ####

func_trs <- plant_data[,c(2:5)]

# Checking and for reference 
# summary(func_trs)
# 
# #    Height_mean         LA_mean            LMA_mean          SM_mean        
# # Min.   : 0.0242   Min.   :     4.0   Min.   :0.02048   Min.   :    0.00  
# # 1st Qu.: 0.2574   1st Qu.:   117.4   1st Qu.:0.04615   1st Qu.:    0.74  
# # Median : 0.6187   Median :  1404.6   Median :0.06053   Median :    2.94  
# # Mean   : 8.5384   Mean   :  9583.6   Mean   :0.08167   Mean   :  797.31  
# # 3rd Qu.: 8.8448   3rd Qu.:  4648.0   3rd Qu.:0.08927   3rd Qu.:   26.36  
# # Max.   :61.5852   Max.   :485704.8   Max.   :0.31301   Max.   :32600.00  
# 
# plant_data[which(func_trs$LA_mean > 40000),]

# hist(func_trs[,1])
# hist(func_trs[,2])
# hist(func_trs[,3])
# hist(func_trs[,4])

### all very left skewed .. try log transform

log_func_trs <- log(func_trs)

hist(log_func_trs[,1])
hist(log_func_trs[,2])
hist(log_func_trs[,3])
hist(log_func_trs[,4])
##### better...

s_log_func_trs <- scale(log_func_trs)



#### use PCA to get trait axes.. ####
# Note this is primarily to reduce the number of variables in the final model. 
# If more data were available, fitting the 4 separately could 
# arguably provide more nuance. 

pca1 <- rda(s_log_func_trs)

#### 
plot(pca1) 
## check the direction of the Axes and make amendments accordingly
# for prefered interpretation of directionality . 
# In my case, I like high scores on Axes 1 to be tall with large seeds, 
# and high scores on axis 2 to be high LA (matching previous literature).
### scores below (i.e. multiply by -1 as desired. )


Axis1 <- -scores(pca1, choices = 1, display = "sites")
Axis2 <- -scores(pca1, choices = 2, display = "sites")

summary(pca1)


#### Get climate data ####

clim <- plant_data[,37:38] 

s_clim <- scale(clim)

#### Join trait axes to climate ####
# I call these 'fixed' as these would be the 'fixed factors' in a frequentist GLMM

fixed1 <- as.data.frame(cbind(Axis1, Axis2, s_clim))

#### Select demographic responses 

demo_trs <- plant_data[,c(11,12,18,22,23,25)] 
#summary(demo_trs)

###########


#######################################
### data transforms for normality (ish...)

### log transform all except GiniF which is sqrt transformed

log_demo <- demo_trs
log_demo[,c(1:4,6)] <- log(demo_trs[,c(1:4,6)])
log_demo[,c(5)] <- sqrt(demo_trs[,c(5)])

### scale demo_trs to units sd 

s_log_demo <- scale(log_demo)

#### Notes on scaling to facilitate back transform for plotting ####
# You can skip this, it's notes for me for plotting model predictions.


# Numbers needed to back transform are 
#                 :sd of the logged dataset
#                 :"meanSDX" of the dataset after dividing by sd. 

#  Steps are :  Add the meanSDX to the vector
#               Then multiply by the sd
#               Then take the exponential to get back to the non-logged form 

######## Numbers for reference will be

#### standard deviations for each column
# 
# sds <- as.numeric(sapply(log_demo, sd))
# sds
# 
# ### mean of each column after dividing by sd
# meanx <- as.numeric()
# for(i in c(1:6)) {
#    meanx[i] <- mean(log_demo[,i]/sds[i])
# }
# meanx
# 
# ### transfrom used for normality as a chr
# names(log_demo)
# trans1 <- c(rep("log", times = 4), "sqrt", "log")
# trans1

### arrange in a dataframe with variable names 
# 
# full_trans_scaling <- as.data.frame(cbind(trans1, sds, meanx))
# row.names(full_trans_scaling) <- names(log_demo)
# 
# write.csv(full_trans_scaling, "scaling_pre_analysis_14_09_2018.csv", 
#           row.names = T)

#### join response vars and fixed vars into one dataframe ####

dataset1 <- as.data.frame(cbind(s_log_demo, fixed1, species = row.names(demo_trs)))


#### We need to name our nodes in the phylogeny for the next step to work. 
### strangely most of them currently have no names.. 

phy_tree2$node.label <- c(1:length(phy_tree2$node.label))

### make this into a comparative data object using the function comparative.data in the
## package 'caper'



comp_data <- comparative.data(phy = phy_tree2, 
                              data = dataset1, 
                              names.col = species, 
                              vcv=TRUE)

head(comp_data$data)
## nothing has been dropped. 
head(comp_data$dropped)


# Add animal column - mcmcGLMM needs the species name column to be called 'animal'
mcmc_data <- comp_data$data
#As MCMCglmm requires a colume named animal for it to identify it
#as a phylo model we include an extra colume with the species names in it.
mcmc_data <- cbind(animal = rownames(mcmc_data), mcmc_data)
mcmc_tree <- comp_data$phy

## saving file for later use, you can ignore this. 
#write.csv(mcmc_data, "mcmc_data.csv") 


#### Bayesian multivariate model #### 

# Here I fit the null models and final models described in the paper, compare DIC values
# and check model outputs, QC plots etc. 
# The original process actually involved a lot more model checking, but I've simplified it 
# to the essentials here.

# Note: Model results will differ slightly on each run due to the stochastic
# nature of the monte-carlo markov chain process. However, overall results, 
# order of model selection and direction of results should be robust. 



#### Notes on choice of priors and co-variance structure ####


# Here I am fitted the phylogeny separately for each response traits,
# but not calculating off-diagonals of the matrix which would 
# represent co-variance of response traits with respect to phylogeny. 
# This choice was a combination of concerns about interpretability, 
# over parameterisation of the model and sample size concerns.  
# G structure sets priors for this, in the code below. 

#### Covariance of response traits ###
# I am using -  rcov = ~ us(trait):units - this includes the variance within 
# traits and the covariance between all trait pairs (off-diagonals in the response matrix)
# As recommended in Hadfield 2010 - Hadfield, J.D., 2010. MCMC methods for multi-response 
# generalized linear mixed models: the MCMCglmm R package. Journal 
# of Statistical Software, 33(2), pp.1-22.
# This also meets the necessary condition that residual structure 'allows each linear predictor 
# to have a unique residual."


# The fixed factor specification trait:PC1 allows each response trait to 
# vary in it's response to each predictor.  

# -1 suppresses the global intercept making interpretation of fixed factors
# easier. Hadfield 2010

#### Modelling ####

#Reduce number of interactions during practice
nitt <- c(1000000)
#length of burnin
burnin <- c(200000)
#amount of thinning
thin <- c(1000)

##### run and save null model for future reference

priorNull <- list(R = list(V = diag(6)*0.02, nu = 7), 
               G = list(G1 = list(V = diag(6), nu = 0.002)))

mod_Null <- MCMCglmm(fixed =  cbind(Age_at_Maturity,Mean_repro,survival_index999,
                                    Repro_life_expectancy,GiniF_999, Gen_time)
                     ~ 1,
                     random= ~ idh(trait):animal, 
                     rcov = ~ us(trait):units,
                     family = c("gaussian", "gaussian",
                                "gaussian", "gaussian",
                                "gaussian", "gaussian"),
                     pedigree = mcmc_tree,
                     data = mcmc_data,
                     scale = TRUE,
                     nitt = nitt,
                     burnin = burnin,
                     thin = thin, 
                     prior = priorNull)
###
dim(mod_Null$VCV)


#saveRDS(object = mod_Null, file = "mod_Null1_13022019.RData")

sum_null <- summary(mod_Null)
sum_null
PhyloN_ests <- as.data.frame(sum_null$Gcovariances)
PhyloN_ests
RcovN_ests <- as.data.frame(sum_null$Rcovariances)
RcovN_ests

mod_Null$DIC

### priorNull no phylo 
priorNull2 <- list(R = list(V = diag(6)*0.02, nu = 7))

mod_Null_no_phy <- MCMCglmm(fixed =  cbind(Age_at_Maturity,Mean_repro,survival_index999,
                                    Repro_life_expectancy,GiniF_999, Gen_time)
                     ~ 1,
                     rcov = ~ us(trait):units,
                     family = c("gaussian", "gaussian",
                                "gaussian", "gaussian",
                                "gaussian", "gaussian"),
                     
                     data = mcmc_data,
                     scale = TRUE,
                     nitt = nitt,
                     burnin = burnin,
                     thin = thin, 
                     prior = priorNull2)
###

#saveRDS(object = mod_Null_no_phy, file = "mod_Null_no_phylo_13022019.RData")


##### Candidate models 

prior2 <- list(R = list(V = diag(6)*0.02, nu = 7), 
                  G = list(G1 = list(V = diag(6), nu = 0.002)))

# 1. Functional traits predict life history in absence of climate.

mod_H1 <- MCMCglmm(fixed = cbind(Age_at_Maturity,Mean_repro,survival_index999,
                                    Repro_life_expectancy,GiniF_999, Gen_time) 
                      ~ trait:PC1 + trait:PC2 - 1,
                      random= ~ idh(trait):animal, 
                      rcov = ~us(trait):units,
                      family = c("gaussian", "gaussian",
                                 "gaussian", "gaussian",
                                 "gaussian", "gaussian"),
                      pedigree = mcmc_tree,
                      data = mcmc_data,
                      scale = TRUE,
                      nitt = nitt,
                      burnin = burnin,
                      thin = thin, 
                      prior = prior2)



#2. Climate predicts life history in absence of functional traits.
mod_H2 <- MCMCglmm(fixed = cbind(Age_at_Maturity,Mean_repro,survival_index999,
                                    Repro_life_expectancy,GiniF_999, Gen_time) 
                      ~ trait:Temp_PCA1 + trait:logAridityIndex - 1,
                      random= ~ idh(trait):animal, 
                      rcov = ~us(trait):units,
                      family = c("gaussian", "gaussian",
                                 "gaussian", "gaussian",
                                 "gaussian", "gaussian"),
                      pedigree = mcmc_tree,
                      data = mcmc_data,
                      scale = TRUE,
                      nitt = nitt,
                      burnin = burnin,
                      thin = thin, 
                      prior = prior2)

#3. Traits and climate predict life history but do not interact.
mod_H3 <- MCMCglmm(fixed = cbind(Age_at_Maturity,Mean_repro,survival_index999,
                                    Repro_life_expectancy,GiniF_999, Gen_time) 
                      ~ trait:PC1 + trait:PC2 +
                        trait:Temp_PCA1 + trait:logAridityIndex - 1,
                      random= ~ idh(trait):animal, 
                      rcov = ~us(trait):units,
                      family = c("gaussian", "gaussian",
                                 "gaussian", "gaussian",
                                 "gaussian", "gaussian"),
                      pedigree = mcmc_tree,
                      data = mcmc_data,
                      scale = TRUE,
                      nitt = nitt,
                      burnin = burnin,
                      thin = thin, 
                      prior = prior2)

#4. The interaction between traits and climate predict life history
mod_H4 <- MCMCglmm(fixed = cbind(Age_at_Maturity,Mean_repro,survival_index999,
                                    Repro_life_expectancy,GiniF_999, Gen_time) 
                      ~ trait:PC1 + trait:PC2 +
                        trait:Temp_PCA1 + trait:logAridityIndex +
                        trait:PC1:Temp_PCA1 +  trait:PC1:logAridityIndex + 
                     trait:PC2:Temp_PCA1 +  trait:PC2:logAridityIndex - 1,
                      random= ~ idh(trait):animal, 
                      rcov = ~us(trait):units,
                      family = c("gaussian", "gaussian",
                                 "gaussian", "gaussian",
                                 "gaussian", "gaussian"),
                      pedigree = mcmc_tree,
                      data = mcmc_data,
                      scale = TRUE,
                      nitt = nitt,
                      burnin = burnin,
                      thin = thin, 
                      prior = prior2)


###### Model 1 - Full model with interactions ####



mod_fin_1 <- MCMCglmm(fixed = cbind(Age_at_Maturity,Mean_repro,survival_index999,
                                                     Repro_life_expectancy,GiniF_999, Gen_time) 
                      ~ trait:PC1 + trait:PC2 +
                        trait:Temp_PCA1 + trait:logAridityIndex +
                        trait:PC1:Temp_PCA1 +  trait:PC1:logAridityIndex - 1,
                                       random= ~ idh(trait):animal, 
                                       rcov = ~us(trait):units,
                                       family = c("gaussian", "gaussian",
                                                  "gaussian", "gaussian",
                                                  "gaussian", "gaussian"),
                                       pedigree = mcmc_tree,
                                       data = mcmc_data,
                                       scale = TRUE,
                                       nitt = nitt,
                                       burnin = burnin,
                                       thin = thin, 
                                       prior = prior2)



###

mod_fin_2 <- MCMCglmm(fixed = cbind(Age_at_Maturity,Mean_repro,survival_index999,
                                    Repro_life_expectancy,GiniF_999, Gen_time) 
                      ~ trait:PC1 + trait:PC2 +
                        trait:Temp_PCA1 + trait:logAridityIndex +
                        trait:PC1:Temp_PCA1 +  trait:PC1:logAridityIndex - 1,
                      random= ~ idh(trait):animal, 
                      rcov = ~us(trait):units,
                      family = c("gaussian", "gaussian",
                                 "gaussian", "gaussian",
                                 "gaussian", "gaussian"),
                      pedigree = mcmc_tree,
                      data = mcmc_data,
                      scale = TRUE,
                      nitt = nitt,
                      burnin = burnin,
                      thin = thin, 
                      prior = prior2)

####


mod_fin_3 <- MCMCglmm(fixed = cbind(Age_at_Maturity,Mean_repro,survival_index999,
                                    Repro_life_expectancy,GiniF_999, Gen_time) 
                      ~ trait:PC1 + trait:PC2 +
                        trait:Temp_PCA1 + trait:logAridityIndex +
                        trait:PC1:Temp_PCA1 +  trait:PC1:logAridityIndex - 1,
                      random= ~ idh(trait):animal, 
                      rcov = ~us(trait):units,
                      family = c("gaussian", "gaussian",
                                 "gaussian", "gaussian",
                                 "gaussian", "gaussian"),
                      pedigree = mcmc_tree,
                      data = mcmc_data,
                      scale = TRUE,
                      nitt = nitt,
                      burnin = burnin,
                      thin = thin, 
                      prior = prior2)

mod_fin_1$DIC

mod_fin_2$DIC

mod_fin_3$DIC

###### checking model assumptions
#### check chains

#plot the fixed effects
plot(mod_fin_1$Sol)
plot(mod_fin_2$Sol)
plot(mod_fin_3$Sol)


#plot the variance and rcov terms.
# First terms are the phylogenetic contributions to the 
# trait values.  The following (7:42) are the covariance of traits matrix
plot(mod_fin_1$VCV)
plot(mod_fin_2$VCV)
plot(mod_fin_3$VCV)


summary(effectiveSize(mod_fin_1$Sol))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   716.4   800.0   800.0   821.3   800.0  1144.3 

summary(effectiveSize(mod_fin_2$Sol))
#   Min.  1st Qu.  Median    Mean 3rd Qu.    Max. 
#  711.1   800.0   800.0   811.9   800.0  1076.8 

summary(effectiveSize(mod_fin_3$Sol))
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    545.5   800.0   800.0   804.8   800.0  1105.0 

summary(effectiveSize(mod_fin_1$VCV))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   482.9   800.0   800.0   787.3   800.0  1150.0 

summary(effectiveSize(mod_fin_2$VCV))
#      Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     576.6   800.0   800.0   792.0   800.0   934.5 

summary(effectiveSize(mod_fin_3$VCV))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   540.2   800.0   800.0   828.5   800.0  1667.7 


####  Check autocorrelation #### 
x1 <- autocorr.diag(mod_fin_1$Sol)
summary(abs(as.numeric(x1[2:5,])))
#       Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#   0.0002893 0.0148855 0.0238632 0.0292117 0.0410361 0.0886487 

x2 <- autocorr.diag(mod_fin_2$Sol)
summary(abs(as.numeric(x2[2:5,])))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#  0.0002613 0.0118391 0.0267354 0.0284334 0.0410570 0.0833777

x3 <- autocorr.diag(mod_fin_3$Sol)
summary(abs(as.numeric(x3[2:5,])))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#  0.0001341 0.0119662 0.0283125 0.0306972 0.0455025 0.0890535 

xV1 <- autocorr.diag(mod_fin_1$VCV)
summary(abs(as.numeric(xV1[2:5,])))
#      Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#   0.000148 0.015718 0.027798 0.032053 0.043054 0.138649 
### 

xV2 <- autocorr.diag(mod_fin_2$VCV)
summary(abs(as.numeric(xV2[2:5,])))
#       Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0003456 0.0104798 0.0219000 0.0268810 0.0377658 0.1616338

xV3 <- autocorr.diag(mod_fin_3$VCV)
summary(abs(as.numeric(xV3[2:5,])))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#  0.0008992 0.0123570 0.0283815 0.0374153 0.0518767 0.1932159


autocorr.plot(mod_fin_1$Sol)
autocorr.plot(mod_fin_2$Sol)
autocorr.plot(mod_fin_3$Sol)

autocorr.plot(mod_fin_1$VCV)
autocorr.plot(mod_fin_2$VCV)
autocorr.plot(mod_fin_3$VCV)

#####

### test convergence

gelman.diag(mcmc.list(mod_fin_1$Sol, mod_fin_2$Sol, mod_fin_3$Sol))

# Potential scale reduction factors:

# Point est. Upper C.I.
# traitAge_at_Maturity:PC1                            1.000      1.001
# traitMean_repro:PC1                                 1.000      1.001
# traitsurvival_index999:PC1                          1.000      1.002
# traitRepro_life_expectancy:PC1                      1.000      1.000
# traitGiniF_999:PC1                                  1.002      1.005
# traitGen_time:PC1                                   0.999      1.000
# traitAge_at_Maturity:PC2                            0.999      1.000
# traitMean_repro:PC2                                 1.000      1.003
# traitsurvival_index999:PC2                          1.001      1.006
# traitRepro_life_expectancy:PC2                      1.001      1.006
# traitGiniF_999:PC2                                  1.004      1.016
# traitGen_time:PC2                                   1.001      1.005
# traitAge_at_Maturity:Temp_PCA1                      1.000      1.000
# traitMean_repro:Temp_PCA1                           1.002      1.006
# traitsurvival_index999:Temp_PCA1                    0.999      1.000
# traitRepro_life_expectancy:Temp_PCA1                0.999      1.000
# traitGiniF_999:Temp_PCA1                            1.000      1.004
# traitGen_time:Temp_PCA1                             0.999      1.000
# traitAge_at_Maturity:logAridityIndex                1.001      1.006
# traitMean_repro:logAridityIndex                     1.001      1.002
# traitsurvival_index999:logAridityIndex              1.002      1.009
# traitRepro_life_expectancy:logAridityIndex          1.003      1.011
# traitGiniF_999:logAridityIndex                      1.001      1.005
# traitGen_time:logAridityIndex                       1.003      1.013
# traitAge_at_Maturity:PC1:Temp_PCA1                  1.005      1.018
# traitMean_repro:PC1:Temp_PCA1                       1.002      1.009
# traitsurvival_index999:PC1:Temp_PCA1                1.002      1.006
# traitRepro_life_expectancy:PC1:Temp_PCA1            1.001      1.007
# traitGiniF_999:PC1:Temp_PCA1                        0.999      0.999
# traitGen_time:PC1:Temp_PCA1                         1.002      1.009
# traitAge_at_Maturity:PC1:logAridityIndex            1.002      1.008
# traitMean_repro:PC1:logAridityIndex                 1.000      1.000
# traitsurvival_index999:PC1:logAridityIndex          1.000      1.001
# traitRepro_life_expectancy:PC1:logAridityIndex      1.000      1.001
# traitGiniF_999:PC1:logAridityIndex                  1.001      1.003
# traitGen_time:PC1:logAridityIndex                   0.999      1.000
# 
# Multivariate psrf
# 
# 1.02

combchainsSol <- mcmc.list(mod_fin_1$Sol, mod_fin_2$Sol, mod_fin_3$Sol)
plot(combchainsSol)

combchainsVCV <- mcmc.list(mod_fin_1$VCV, mod_fin_2$VCV, mod_fin_3$VCV)
plot(combchainsVCV)

gelman.diag(mcmc.list(mod_fin_1$VCV[,1:6], mod_fin_2$VCV[,1:6], mod_fin_3$VCV[,1:6]))

# Potential scale reduction factors:
#   
#   Point est. Upper C.I.
# traitAge_at_Maturity.animal            1.001       1.00
# traitMean_repro.animal                 1.016       1.03
# traitsurvival_index999.animal          0.999       1.00
# traitRepro_life_expectancy.animal      1.001       1.00
# traitGiniF_999.animal                  1.003       1.01
# traitGen_time.animal                   1.029       1.04
# 
# Multivariate psrf
# 
# 1.01
# 
# 1


###

### for reference
saveRDS(object = mod_fin_1, file = "mod_fin_1_2020.RData")
# saveRDS(object = mod_fin_2, file = "mod_fin_2_24072018.RData")
# saveRDS(object = mod_fin_3, file = "mod_fin_3_24072018.RData")
# ###

### Run final model without the phylogenetic structure for variance
### comparisons later 
