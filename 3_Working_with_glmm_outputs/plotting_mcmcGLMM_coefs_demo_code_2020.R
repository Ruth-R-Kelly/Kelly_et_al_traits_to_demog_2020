#### code to plot MCMCglmm coefficients for 03/09/2018 ####

#### remove previous objects

rm(list = ls())

### set working directory if needed

### load libraries if you haven't done this using the overview code.
# library("ggplot2")
# library("MCMCglmm")
# library("pheatmap")
# library("viridis")
####

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

#row.names(fixed_ests)
#dim(fixed_ests)

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

### Hashed out lines here call plot to a pdf/phg file in the current folder options. 
# otherwise plot will appear in the rmd output, or r plot window if you run the code in the console.
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
print(zp1)
#dev.off()


#### phylogeny coefs ####


phylo_ests <- as.data.frame(sum1$Gcovariances)
names(phylo_ests)
summary(phylo_ests)
row.names(phylo_ests)


#### scale phylo estimates by full model variance, 
# to make their values interprettable outside of the model context 

# value of full model variance was calculated in the code
# 'model_variance_09_06_2020.R' as .. 8.070693  - ie. variance of fixed, phylo and residuals
Model_tot_var <- 8.070693

### get outputted estimates for phylogeny
phylo_ests <- as.data.frame(sum1$Gcovariances) 

names(phylo_ests) <- c("post.mean", 
                       "Lower_95_CI","Upper_95_CI",  "eff.samp")

row.names(phylo_ests)

phylo_ests$Trait <- c(      "Age at maturity",
                            "Annual reproduction",
                            "Distribution of mortality",
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
                              "Distribution of mortality",
                              "Age at maturity"
                            ))

summary(phylo_ests$Trait)

phylo_ests$post.mean.scaled <- (phylo_ests$post.mean/Model_tot_var)*100
phylo_ests$Lower_95_CI.scaled <- (phylo_ests$Lower_95_CI/Model_tot_var)*100
phylo_ests$Upper_95_CI.scaled <- (phylo_ests$Upper_95_CI/Model_tot_var)*100

sum(phylo_ests$post.mean.scaled)


#### Again hashed out lines here make the outputs print to pdf/png in the current folder
#pdf("Phylocoefs_variance_scaled.pdf", height = 3.5, width = 5)
#png("aPhylocoefs.png", height = 3.5, width = 5, units = "in", res = 300)
zp1 <- ggplot(phylo_ests)
zp1 <- zp1 + geom_hline(yintercept = 0, colour = gray(2/3), lty = 2)
zp1 <- zp1 + geom_pointrange(aes(x = Trait, y = post.mean.scaled, 
                                 ymin = Lower_95_CI.scaled,
                                 ymax = Upper_95_CI.scaled),
                             lwd = 0.75, position = position_dodge(width = 1/2),
                             shape = 21, fill = "grey40", colour = "grey40")
zp1 <- zp1 + coord_flip() + theme_bw()
zp1 <- zp1 + theme(axis.title.y=element_blank(), axis.title.x = element_blank())
zp1 <- zp1 + theme(axis.text.x = element_text(colour = "black"))
zp1 <- zp1 + theme(axis.text.y = element_text(colour = "black"))
print(zp1)
#dev.off()

#pdf("Phylocoefs_variance_unscaled.pdf", height = 3.5, width = 5)
#png("aPhylocoefs.png", height = 3.5, width = 5, units = "in", res = 300)
zp1 <- ggplot(phylo_ests)
zp1 <- zp1 + geom_hline(yintercept = 0, colour = gray(2/3), lty = 2)
zp1 <- zp1 + geom_pointrange(aes(x = Trait, y = post.mean, 
                                 ymin = Lower_95_CI,
                                 ymax = Upper_95_CI),
                             lwd = 0.75, position = position_dodge(width = 1/2),
                             shape = 21, fill = "grey40", colour = "grey40")
zp1 <- zp1 + coord_flip() + theme_bw()
zp1 <- zp1 + theme(axis.title.y=element_blank(), axis.title.x = element_blank())
zp1 <- zp1 + theme(axis.text.x = element_text(colour = "black"))
zp1 <- zp1 + theme(axis.text.y = element_text(colour = "black"))
print(zp1)
#dev.off()

