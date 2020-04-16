#### code to plot MCMCglmm outputs for 06/07/2018 ####

### remove objects from workspace

rm(list = ls())

### set working directory

setwd("C:/R/trait_demo_analysis_2018/New_version_08_06_2018/plotting_06_07_2018")

library("ggplot2")
library("MCMCglmm")

####
dir()

## "data_with_all_env_vars_08_06_2018.csv"
## "MCMCpredictors_06_07_2018.csv"
## "mod_fin_1_06072018.RData"
## "mod_fin_2_06072018.RData"
## "mod_mcmc_null.RData" 
mod_fin_1 <- readRDS("mod_fin_1_24072018.RData")
mod_mcmc_null <- readRDS("mod_Null1_13022019.RData")
mod_no_phylo <- readRDS("mod_fin_no_phylo_18_12_2019.RData")
explanatory_data <- read.csv( "mcmc_data_27_07_2018.csv"  , row.names = "X")


sum1 <- summary(mod_fin_1)

sum1$DIC

Phylo_ests <- as.data.frame(sum1$Gcovariances) 
Rcov_ests <- as.data.frame(sum1$Rcovariances)

#### separating residual covariance into between life-history traits
#### i.e. one explaining another 
#### and fully residual (within traits) 

sum(Rcov_ests$post.mean) # 2.448038
sum(Phylo_ests$post.mean) # 4.387432

mVarF <- var(as.vector(apply(mod_fin_1$Sol,2,mean) %*% t(mod_fin_1$X)))
mVarF
 # 0.2197406

percfixedVar <- mVarF/(mVarF + sum(Rcov_ests$post.mean) + sum(Phylo_ests$post.mean))
percfixedVar ## 0.03114587

percRcovVar <- sum(Rcov_ests$post.mean)/(mVarF + sum(Rcov_ests$post.mean) + sum(Phylo_ests$post.mean))
percPhyloVar <- sum(Phylo_ests$post.mean)/(mVarF + sum(Rcov_ests$post.mean) + sum(Phylo_ests$post.mean))

marginalF <- (mVarF + sum(Phylo_ests$post.mean))/(mVarF + sum(Rcov_ests$post.mean) + sum(Phylo_ests$post.mean))

percfixedVar # 0.03114587
percRcovVar #  0.3469829
marginalF  # 0.6530171

####  model without phylogeny all other variables the same. 

sum_no_phy <-summary(mod_no_phylo)

########
sum_no_phy

PhyloN_ests <- as.data.frame(sum_no_phy$Gcovariances) # there are none
RcovN_ests <- as.data.frame(sum_no_phy$Rcovariances) # 

RcovN_ests

sum(RcovN_ests$post.mean) #  2.946352

mVarF <- var(as.vector(apply(mod_no_phylo$Sol,2,mean) %*% t(mod_fin_1$X)))
mVarF
## [1] 0.1935091

    

percfixedVar <- mVarF/(mVarF + sum(RcovN_ests$post.mean))
percfixedVar ## 6.1%

percRcovVar <- sum(RcovN_ests$post.mean)/(mVarF + sum(RcovN_ests$post.mean))

marginalF <- (mVarF)/(mVarF + sum(RcovN_ests$post.mean))

percfixedVar # 0.06162982
percRcovVar #  0.9383702
marginalF  # 0.06162982

########################


mod_fin_1$DIC
mod_no_phylo$DIC
