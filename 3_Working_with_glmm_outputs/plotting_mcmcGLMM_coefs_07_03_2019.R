#### code to plot MCMCglmm coefficients for 03/09/2018 ####

#### remove previous objects

rm(list = ls())

### set working directory

setwd("C:/R/traits_to_demo_vr2/New_version_21_07_2018/Final_models_27_07_2018/plotting_28_08_2018")
library("ggplot2")
library("MCMCglmm")
library("pheatmap")
library("viridis")
####

# [1] "data_with_all_env_vars_20_07_2018.csv"       
# [2] "mcmc_data_27_07_2018.csv"                    
# [3] "mcmc_model_results_19_07_2018.R"             
# [4] "mod_fin_1_24072018.RData"                    
# [5] "mod_fin_2_24072018.RData"                    
# [6] "mod_fin_3_24072018.RData"                    
# [7] "plotting_27_07_2018"                         
# [8] "plotting_mcmcGLMM_coefs06_07_2018.R"         
# [9] "plotting_mcmcGLMM_outputs.R"                 
# [10] "plotting_mcmcGLMM_outputs06_07_2018.R"       
# [11] "plotting_mcmcGLMM_variance_calcs07_07_2018.R"
# [12] "reduced_tree_for_traits_08_06_2018.nex" 

mod_fin_1 <- readRDS("mod_fin_1_24072018.RData")
explanatory_data <- read.csv("mcmc_data_27_07_2018.csv", row.names = "X")

sum1 <- summary(mod_fin_1)

sum1


str(sum1)

#### get fixed effects estimates from the summary object
fixed_ests <- as.data.frame(sum1$solutions)

names(fixed_ests) <- c("post.mean", 
                       "Lower_95_CI", "Upper_95_CI", "eff.samp",  "pMCMC")


### Make an index of the rows in fixed ests for each life 
### history trait

row.names(fixed_ests)
dim(fixed_ests)

##### Make Index column to life history traits
LH_Index <-  rep(c("Age at maturity", "Annual reproduction",
               "Distribution of mortality", "Mature lifespan",
               "Inequality of reproduction", "Generation time"),6)

fixed_ests$Life_history <- LH_Index

##### Make Index column to predictor vars. #### 

Predictor_Index <-  c(rep("Structure and size",6),
                      rep("Leaf economics",6),
                      rep("Temperature PCA", 6),
                      rep("Moisture availability",6),
                      rep("Size * TemperaturePCA", 6),
                      rep("Size * Moisture", 6))


fixed_ests$Predictor <- Predictor_Index

##### Make order for neater plotting ####
## sort(levels(as.factor(LH_Index)), decreasing = TRUE)

fixed_ests$Life_history <- ordered(fixed_ests$Life_history, 
                            levels = c( 
                              "Inequality of reproduction", 
                              "Annual reproduction",
                            "Mature lifespan",  
                            "Generation time",  
                            "Distribution of mortality",
                            "Age at maturity"
                            ))

####

levels(as.factor(fixed_ests$Predictor))

fixed_ests$Predictor <- ordered(fixed_ests$Predictor, 
                                   levels = c( 
                                     "Size * Moisture" ,
                                     "Size * TemperaturePCA", 
                                     "Moisture availability",
                                     "Temperature PCA",
                                     "Leaf economics",
                                     "Structure and size"
                                   ))

### add column for transparancy in plotting based on significance

TransP <- c()


for(i in 1:nrow(fixed_ests)) {
  x <- fixed_ests$pMCMC[i]
if(x > 0.061) {
  TransP[i] <- "Non-Sig"
} else {
  TransP[i] <- "Sig"
}
  
}

fixed_ests$TransP <- TransP

####  old palette option
# 
# palette1 <- (c("#6baed6", "#08519c",
#                "#d95f0e", "#993404",
#                "#74c476", "#006d2c"))
# 

#pdf("afixedcoefs_with_shading_viridis_04_06_2019.pdf", height = 10, width = 8)
#png("afixedcoefs_with_shading_04_06_2019.png", height = 10, width = 8, units = "in", res = 300)
zp1 <- ggplot(fixed_ests)
zp1 <- zp1 + geom_hline(yintercept = 0, colour = gray(2/3), lty = 2)

zp1 <- zp1 + geom_pointrange(aes(x = Life_history, y = post.mean, 
                                 ymin = Lower_95_CI, 
                                 ymax = Upper_95_CI, 
                                 colour = Predictor, 
                                 fill = Predictor),
                             lwd = 1, position = position_dodge(width = 3/4),
                             shape = 21)

zp1 <- zp1 + geom_rect(mapping = aes(xmin =0.5, xmax = 1.5, ymin = -1.5, ymax = 2),  fill = "grey80", alpha = 0.02)
zp1 <- zp1 + geom_rect(mapping = aes(xmin =2.5, xmax = 3.5, ymin = -1.5, ymax = 2),  fill = "grey80", alpha = 0.02)
zp1 <- zp1 + geom_rect(mapping = aes(xmin =4.5, xmax = 5.5, ymin = -1.5, ymax = 2),  fill = "grey80", alpha = 0.02)

zp1 <- zp1 + geom_pointrange(aes(x = Life_history, y = post.mean, 
                                 ymin = Lower_95_CI, 
                                 ymax = Upper_95_CI, 
                                 colour = Predictor, 
                                 fill = Predictor),
                             lwd = 1, position = position_dodge(width = 3/4),
                             shape = 21)

zp1 <- zp1 + coord_flip() + theme_bw() + ylim(-1.5, 2)
zp1 <- zp1 + theme(legend.title = element_blank())
zp1 <- zp1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank())
zp1 <- zp1 + theme(axis.title.y=element_blank(), axis.title.x = element_blank())
zp1 <- zp1 + theme(axis.text.x = element_text(colour = "blue")) 
zp1 <- zp1 + ggtitle("")
zp1 <- zp1 + scale_color_viridis(discrete = TRUE)
zp1 <- zp1 + scale_fill_viridis(discrete = TRUE)
zp1 <- zp1 + theme(axis.text.y = element_text(colour = "black"))
zp1 <- zp1 + guides(fill = guide_legend(reverse = TRUE), colour = guide_legend(reverse = TRUE))
#print(zp1)
#dev.off()


#### phylogeny coefs ####


phylo_ests <- as.data.frame(sum1$Gcovariances)
names(phylo_ests)
summary(phylo_ests)
row.names(phylo_ests)
#library(dplyr)

#### scale phylo estimates by full model variance. 


P_ests <- as.data.frame(sum1$Gcovariances) 
Rcov_ests <- as.data.frame(sum1$Rcovariances)

##
mVarF <- var(as.vector(apply(mod_fin_1$Sol,2,mean) %*% t(mod_fin_1$X)))
mVarF
## [1] 0.2197406

Model_tot_var <- sum(sum(Rcov_ests$post.mean) + sum(P_ests$post.mean) + mVarF)


names(phylo_ests) <- c("post.mean", 
                       "Lower_95_CI","Upper_95_CI",  "eff.samp")

row.names(phylo_ests)

phylo_ests$Trait <- c(      "Age at maturity",
                            "Annual reproduction",
                            "Median/max lifespan",
                            "Mature lifespan",
                            "Inequality of reproduction",
                            "Generation time"
                            )


phylo_ests$Trait <- ordered(phylo_ests$Trait, 
                            levels = c( 
                              "Generation time",
                              "Annual reproduction",
                              "Mature lifespan",
                              "Inequality of reproduction",
                              "Median/max lifespan",
                              "Age at maturity"
                            ))

summary(phylo_ests$Trait)

phylo_ests$post.mean.scaled <- (phylo_ests$post.mean/Model_tot_var)*100
phylo_ests$Lower_95_CI.scaled <- (phylo_ests$Lower_95_CI/Model_tot_var)*100
phylo_ests$Upper_95_CI.scaled <- (phylo_ests$Upper_95_CI/Model_tot_var)*100

phylo_ests

#pdf("aPhylocoefs.pdf", height = 3.5, width = 5)
#png("aPhylocoefs.png", height = 3.5, width = 5, units = "in", res = 300)
zp1 <- ggplot(phylo_ests)
zp1 <- zp1 + geom_hline(yintercept = 0, colour = gray(2/3), lty = 2)
zp1 <- zp1 + geom_pointrange(aes(x = Trait, y = post.mean.scaled, 
                                 ymin = Lower_95_CI.scaled,
                                 ymax = Upper_95_CI.scaled),
                             lwd = 1, position = position_dodge(width = 1/2),
                             shape = 21, fill = "#08519c", colour = "#08306b")
zp1 <- zp1 + coord_flip() + theme_bw()
zp1 <- zp1 + theme(axis.title.y=element_blank(), axis.title.x = element_blank())
zp1 <- zp1 + theme(axis.text.x = element_text(colour = "black"))
zp1 <- zp1 + theme(axis.text.y = element_text(colour = "black"))
zp1 <- zp1 + ggtitle("Phylogenetic signal")
print(zp1)
#dev.off()


##### Output estimates tables ####
write.csv(fixed_ests,"Final_model_fixed_ests_04_09_2018.csv", row.names = TRUE)
write.csv(phylo_ests, "Final_model_phylo_ests_04_09_2018.csv", row.names = TRUE)

RCov_ests <- as.data.frame(sum1$Rcovariances)
write.csv(RCov_ests, "Final_model_RCov_ests_04_09_2018.csv", row.names = TRUE)
#####

names(RCov_ests)
names(RCov_ests) <- c("post.mean", 
                       "Lower_95_CI","Upper_95_CI",  "eff.samp")

RCov_ests$Trait <- row.names(RCov_ests)
dim(RCov_ests)
row.names(RCov_ests)
Rows_self_cov <- c(1,8,15,22,29,36)
diff_lh_pairs <- RCov_ests[-Rows_self_cov,]
diff_lh_pairs

#png("Resid_covars.png", height = 10, width = 12, units = "in", res = 300)
zp1 <- ggplot(diff_lh_pairs)
zp1 <- zp1 + geom_hline(yintercept = 0, colour = gray(2/3), lty = 2)
zp1 <- zp1 + geom_pointrange(aes(x = Trait, y = post.mean, 
                                 ymin = Lower_95_CI,
                                 ymax = Upper_95_CI),
                             lwd = 1, position = position_dodge(width = 1/2),
                             shape = 21, fill = "grey", colour = "grey20")
zp1 <- zp1 + coord_flip() + theme_bw()
zp1 <- zp1 + theme(axis.title.y=element_blank(), axis.title.x = element_blank())
zp1 <- zp1 + theme(axis.text.x = element_text(colour = "black")) 
zp1 <- zp1 + theme(axis.text.y = element_text(colour = "black"))
zp1 <- zp1 + ggtitle("Residual covariance between life-history traits")
print(zp1)
dev.off()

#### Residual correlation ####

RCov_ests
rcovMean <- RCov_ests$post.mean
summary(rcovMean)
cov_mat <- matrix(rcovMean, ncol = 6)
rcovMean
cov_mat
row.names(cov_mat) <- c("Age at maturity", "Annual reproduction", 
                        "Mean/max life span", "Mature lifespan",
                        "Variation in reproduction", "Generation time")
cov_mat

###
cov_mat <- as.data.frame(cov_mat)
cor_mat <- as.data.frame(cor_mat)  

names(cov_mat) <- row.names(cov_mat)
names(cor_mat) <- row.names(cov_mat)
####




#png("Resid_correlations.png", height = 6, width = 6, units = "in", res = 300)
pheatmap(cor_mat, display_numbers = TRUE, fontsize = 10,
         cluster_cols = TRUE, cluster_rows = TRUE,
         number_color = "grey10" )
#dev.off()

pdf(file = "Resid_correlations.pdf", height = 6, width = 6)
pheatmap(cor_mat, display_numbers = TRUE, fontsize = 10,
         cluster_cols = TRUE, cluster_rows = TRUE,
         number_color = "grey10" )
dev.off()


pdf(file = "heatmap_of_demography_correlations_res2.pdf")
pheatmap(cor_mat2, display_numbers = TRUE, cex = 1.5)
dev.off()



pdf(file = "heatmap_of_demography_covariance_res.pdf")
pheatmap(cov_mat, display_numbers = TRUE, cex = 1.5)
dev.off()


