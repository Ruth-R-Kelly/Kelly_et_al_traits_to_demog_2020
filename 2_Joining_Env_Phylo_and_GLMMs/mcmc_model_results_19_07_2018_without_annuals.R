### MCMCglmm code rewritten 18/06/2018 
## Datasets include 99.9 lifespan calculations and averaging of populations within 5km to 
## match the scale of the climate data. 

### load libraries

library(MCMCglmm)
library(MCMCpack)
library(ape)
library(caper)
library(vegan)
library(pheatmap)
library(phangorn)
#######

setwd("C:/R/traits_to_demo_vr2/New_version_21_07_2018/Final_models_27_07_2018")


### clear previous files 
rm(list = ls())
dir()
### read in plant trait data in order of phylo 
plant_data <- read.csv("data_with_all_env_vars_20_07_2018.csv", row.names = "X")

# str(plant_data)
# head(plant_data)
# row.names(plant_data)

# nrow(plant_data)
summary(plant_data)

####
phy_tree <- read.nexus("reduced_tree_for_traits_08_06_2018.nex")  

plot(phy_tree)

#### force tree to be ultrametric following: 
#http://blog.phytools.org/2016/08/fixing-ultrametric-tree-whose-edges-are.html
phy_tree2 <- nnls.tree(cophenetic(phy_tree),phy_tree,rooted=TRUE)

summary(phy_tree2$edge.length)

is.ultrametric(phy_tree2)
plot(phy_tree2)

#### check trees similarity

tips<-phy_tree2$tip.label
cor(as.vector(cophenetic(phy_tree)[tips,tips]),
    as.vector(cophenetic(phy_tree2)[tips,tips]))

### 1

### check species names align with phylogeny
phy_tree2$tip.label == row.names(plant_data)

names(plant_data)

####

func_trs <- plant_data[,c(2:5)]

summary(func_trs)

#    Height_mean         LA_mean            LMA_mean          SM_mean        
# Min.   : 0.0242   Min.   :     4.0   Min.   :0.02048   Min.   :    0.00  
# 1st Qu.: 0.2574   1st Qu.:   117.4   1st Qu.:0.04615   1st Qu.:    0.74  
# Median : 0.6187   Median :  1404.6   Median :0.06053   Median :    2.94  
# Mean   : 8.5384   Mean   :  9583.6   Mean   :0.08167   Mean   :  797.31  
# 3rd Qu.: 8.8448   3rd Qu.:  4648.0   3rd Qu.:0.08927   3rd Qu.:   26.36  
# Max.   :61.5852   Max.   :485704.8   Max.   :0.31301   Max.   :32600.00  

plant_data[which(func_trs$LA_mean > 40000),]

hist(func_trs[,1])
hist(func_trs[,2])
hist(func_trs[,3])
hist(func_trs[,4])
### all very left skewed .. try log transform

log_func_trs <- log(func_trs)

hist(log_func_trs[,1])
hist(log_func_trs[,2])
hist(log_func_trs[,3])
hist(log_func_trs[,4])
##### better though far from perfect!

s_log_func_trs <- scale(log_func_trs)

##############################

### use PCA to get trait axes.. 

pca1 <- rda(s_log_func_trs)

#### 
plot(pca1) ### check the direction of the Axes and make amendments accordindly to the 
### scores below (i.e. multiply by -1 as desired. )
text(pca1, display = "species")

Axis1 <- -scores(pca1, choices = 1, display = "sites")
Axis2 <- -scores(pca1, choices = 2, display = "sites")

summary(pca1)


#### get climate data 
ncol(plant_data)
clim <- plant_data[,44:45] 

### standardise 
s_clim <- scale(clim)
#### join trait axes to climate 

fixed1 <- as.data.frame(cbind(Axis1, Axis2, s_clim))

#### Select demographic responses 
names(plant_data)
demo_trs <- plant_data[,c(11,12,18,22,23,25)] 
summary(demo_trs)

###########


#######################################
######### data transforms for normality (ish)
names(demo_trs)
hist(log(demo_trs$Age_at_Maturity))
hist(log(demo_trs$Mean_repro))
hist(log(demo_trs$survival_index999))
hist(log(demo_trs$Repro_life_expectancy))
hist(sqrt(demo_trs$GiniF_999))
hist(log(demo_trs$Gen_time))

### log transform all except GiniF which is sqrt transformed

log_demo <- demo_trs
log_demo[,c(1:4,6)] <- log(demo_trs[,c(1:4,6)])
log_demo[,c(5)] <- sqrt(demo_trs[,c(5)])

### scale demo_trs to units sd 

s_log_demo <- scale(log_demo)

#### Notes on scaling to facilitate back transform for plotting ####

names(log_demo)

sd(log_demo$Age_at_Maturity)
### 1.274008
mean(log_demo$Age_at_Maturity)
### 1.645934

summary(s_log_demo)

##########
###### numbers needed to back transform are 
######                 :sd of the logged dataset
######                 :"meanSDX" of the dataset after dividing by sd. 

# ###### Steps are :  Add the meanSDX to the vector
#                     Then multiply by the sd
#                     Then take the exponential to get back to the non-logged form 

######## Numbers for reference will be

#### standard deviations for each column

sds <- as.numeric(sapply(log_demo, sd))
sds

### mean of each column after dividing by sd
meanx <- as.numeric()
for(i in c(1:6)) {
   meanx[i] <- mean(log_demo[,i]/sds[i])
}
meanx

### transfrom used for normality as a chr
names(log_demo)
trans1 <- c(rep("log", times = 4), "sqrt", "log")
trans1

### arrange in a dataframe with variable names 

full_trans_scaling <- as.data.frame(cbind(trans1, sds, meanx))
row.names(full_trans_scaling) <- names(log_demo)
full_trans_scaling
write.csv(full_trans_scaling, "scaling_pre_analysis_14_09_2018.csv", 
          row.names = T)


##############################
#pdf(file = "heatmap_of_demography_cors.pdf")
pheatmap(cor(demo_trs), display_numbers = TRUE)
#dev.off()

###
#pdf(file = "heatmap_of_trait_cors.pdf")
pheatmap(cor(fixed1), display_numbers = TRUE)
#dev.off()

cor(plant_data$bio1, plant_data$Temp_PCA1)
# [1] 0.9581546

### Temp_PCA1 - higher = hotter. Very strong correlation.

cor(plant_data$bio4, plant_data$Temp_PCA1)
#[1] -0.9290489

### Temp_PCA1 - higher = less variable. Very strong correlation.

### Note high correlation between PC1 and TempPCA -0.75

cor.test(fixed1$PC1, fixed1$Temp_PCA1)
####
## 0.6930723 

#### join response vars and fixed vars into one dataframe

dataset1 <- as.data.frame(cbind(s_log_demo, fixed1, species = row.names(demo_trs)))

###### subset dataset to remove annuals. ####

names(plant_data)

dataset1 <- dataset1[which(plant_data$OrganismType != "Annual"),]


head(dataset1)

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

###### time to try MCMCglmm!!

# Add animal column 
mcmc_data <- comp_data$data
#As MCMCglmm requires a colume named animal for it to identify it
#as a phylo model we include an extra colume with the species names in it.
mcmc_data <- cbind(animal = rownames(mcmc_data), mcmc_data)
mcmc_tree <- comp_data$phy

#####
write.csv(mcmc_data, "mcmc_data_27_07_2018.csv") 

##### set iterations, burnin time and in thinning to be shortish for testing


####
names(mcmc_data)

##### To see notes on error structure and selection of random vars see: 
## mcmc_test_code_normalised_model_01_06_2017

#Reduce number of interactions during practice
nitt <- c(1000000)
#length of burnin
burnin <- c(200000)
#amount of thinning
thin <- c(1000)


##### Candidate models 



###### Model 1 - Full model with interactions ####

prior2 <- list(R = list(V = diag(6)*0.02, nu = 7), 
               G = list(G1 = list(V = diag(6), nu = 0.002)))

prior2

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

######

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

#dim(mcmc_data)
###### checking model assumptions
#### check chains

#plot the fixed effects
plot(mod_fin_1$Sol)

#plot the variance and rcov terms.
# First terms are the phylogenetic contributions to the 
# trait values.  The following (7:42) are the covariance of traits matrix
plot(mod_fin_1$VCV)

summary(effectiveSize(mod_fin_1$Sol))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   716.4   800.0   800.0   821.3   800.0  1144.3 

####  Check autocorrelation #### 
x1 <- autocorr.diag(mod_fin_1$Sol)
summary(abs(as.numeric(x1[2:5,])))
#       Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#   0.0002893 0.0148855 0.0238632 0.0292117 0.0410361 0.0886487 

xV1 <- autocorr.diag(mod_fin_1$VCV)
summary(abs(as.numeric(xV1[2:5,])))
#      Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#   0.000148 0.015718 0.027798 0.032053 0.043054 0.138649 
### 
autocorr.plot(mod_fin_1$Sol)

autocorr.plot(mod_fin_1$VCV)

#####

### test convergence

gelman.diag(mcmc.list(mod_fin_1$Sol, mod_fin_2$Sol))

# Potential scale reduction factors:
#   
#   Point est. Upper C.I.
# traitAge_at_Maturity:PC1                            0.999       1.00
# traitMean_repro:PC1                                 1.002       1.01
# traitsurvival_index999:PC1                          1.009       1.02
# traitRepro_life_expectancy:PC1                      1.001       1.01
# traitGiniF_999:PC1                                  1.003       1.02
# traitGen_time:PC1                                   1.000       1.00
# traitAge_at_Maturity:PC2                            1.001       1.00
# traitMean_repro:PC2                                 1.000       1.00
# traitsurvival_index999:PC2                          1.001       1.00
# traitRepro_life_expectancy:PC2                      1.001       1.01
# traitGiniF_999:PC2                                  0.999       1.00
# traitGen_time:PC2                                   1.001       1.01
# traitAge_at_Maturity:Temp_PCA1                      1.001       1.01
# traitMean_repro:Temp_PCA1                           1.004       1.01
# traitsurvival_index999:Temp_PCA1                    1.003       1.01
# traitRepro_life_expectancy:Temp_PCA1                1.002       1.01
# traitGiniF_999:Temp_PCA1                            1.000       1.00
# traitGen_time:Temp_PCA1                             1.000       1.00
# traitAge_at_Maturity:logAridityIndex                1.002       1.00
# traitMean_repro:logAridityIndex                     0.999       1.00
# traitsurvival_index999:logAridityIndex              1.004       1.02
# traitRepro_life_expectancy:logAridityIndex          1.001       1.00
# traitGiniF_999:logAridityIndex                      0.999       1.00
# traitGen_time:logAridityIndex                       1.002       1.01
# traitAge_at_Maturity:PC1:Temp_PCA1                  1.000       1.00
# traitMean_repro:PC1:Temp_PCA1                       1.004       1.02
# traitsurvival_index999:PC1:Temp_PCA1                1.005       1.02
# traitRepro_life_expectancy:PC1:Temp_PCA1            1.000       1.00
# traitGiniF_999:PC1:Temp_PCA1                        1.003       1.01
# traitGen_time:PC1:Temp_PCA1                         1.001       1.01
# traitAge_at_Maturity:PC1:logAridityIndex            1.009       1.05
# traitMean_repro:PC1:logAridityIndex                 1.010       1.02
# traitsurvival_index999:PC1:logAridityIndex          1.001       1.01
# traitRepro_life_expectancy:PC1:logAridityIndex      1.011       1.05
# traitGiniF_999:PC1:logAridityIndex                  1.004       1.02
# traitGen_time:PC1:logAridityIndex                   1.003       1.02
# 
# Multivariate psrf
# 
# 1.05

combchainsSol <- mcmc.list(mod_fin_1$Sol, mod_fin_2$Sol)
plot(combchainsSol)

combchainsVCV <- mcmc.list(mod_fin_1$VCV, mod_fin_2$VCV)
plot(combchainsVCV)

gelman.diag(mcmc.list(mod_fin_1$VCV[,1:6], mod_fin_2$VCV[,1:6]))

#                                     Point est. Upper C.I.
# traitAge_at_Maturity.animal             1.00       1.01
# traitMean_repro.animal                  1.10       1.12
# traitsurvival_index999.animal           1.01       1.03
# traitRepro_life_expectancy.animal       1.01       1.04
# traitGiniF_999.animal                   1.03       1.06
# traitGen_time.animal                    1.04       1.05
# 
# Multivariate psrf
# 
# 1.02

###

summary(mod_fin_1)


setwd("C:/R/traits_to_demo_vr2/New_version_21_07_2018/Final_models_27_07_2018")

saveRDS(object = mod_fin_1, file = "models_without_annual_plants.RData")
###


