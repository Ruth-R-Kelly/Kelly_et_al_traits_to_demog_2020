### Variance calculation 09/06/2020

### remove objects from workspace

rm(list = ls())

### set working directory if needed 


# Load the library if you haven't done so in the overview.rmd file. 
# library("MCMCglmm")


mod_fin_1 <- readRDS("mod_fin1_2020.RData")

### 
#summary(mod_fin_1)

###
des_mat <- as.matrix(mod_fin_1$X)

#### select columns with data for each predictor variable from the design matrix. 
# these are repeated but shifted one each 80 columns so only the first is needed. 
reduced_mat <- des_mat[1:80,c(1,7,13,19,25,31)]
# Repro_mat <- des_mat[81:160,c(2,8,14,20,26,32)]
# Age_mat == Repro_mat

# Get mean values of fixed effects for each trait. 

sum1 <- summary(mod_fin_1)
sols1 <- sum1$solutions

ageMsols <- sols1[c(1,7,13,19,25,31),1]
Reprosols <- sols1[c(2,8,14,20,26,32),1] 
SurvSols <- sols1[c(3,9,15,21,27,33),1] 
ReproLifeSols <- sols1[c(4,10,16,22,28,34),1] 
Gini_Sols <- sols1[c(5,11,17,23,29,35),1] 
Gen_time_Sols <- sols1[c(6,12,18,24,30,36),1] 

Var_ageM <- var(as.vector(ageMsols %*% t(reduced_mat)))
Var_ageM

Var_repro <- var(as.vector(SurvSols %*% t(reduced_mat)))
Var_repro

Var_surv <- var(as.vector(Reprosols %*% t(reduced_mat)))
Var_surv

Var_reprolife <- var(as.vector(ReproLifeSols %*% t(reduced_mat)))
Var_reprolife

Var_gini <- var(as.vector(Gini_Sols %*% t(reduced_mat)))
Var_gini

Var_gentime <- var(as.vector(Gen_time_Sols %*% t(reduced_mat)))
Var_gentime
####

Total_fixed_var <- Var_ageM + Var_repro + Var_surv + Var_reprolife + Var_gini + Var_gentime


### variance of phylo effects.. 


Phylo_ests <- as.data.frame(sum1$Gcovariances) 
Phylo_var <- sum(Phylo_ests$post.mean)

### variance of residuals 

Rcov_ests <- as.data.frame(sum1$Rcovariances)
Rcov_var <- sum(Rcov_ests$post.mean)

####  Variance components.. 

marginalF <- (Phylo_var+ Total_fixed_var)/(Total_fixed_var + Rcov_var + Phylo_var)
marginalF 
# 0.6966757
#
##### 

perc_fixed <- Total_fixed_var/(Total_fixed_var + Rcov_var + Phylo_var)
#[1] 0.1530505
perc_phylo <- Phylo_var/(Total_fixed_var + Rcov_var + Phylo_var)
#[1] 0.5436252


