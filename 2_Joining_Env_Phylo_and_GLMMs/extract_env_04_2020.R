### Extracting environmental variables 
###          20/06/2018
###

### This code matches trait and demographic data with environmental information, 
### as per Coutts et al. 2016 "Extrapolating demography with climate, proximity 
### and phylogeny: approach with caution. Ecology letters, 19(12), pp.1429-1438."

### Specifically, this joins information on Aridity Index, Temperature, and 
### standard deviation of precipitation. See details 'overview_analysis_21_04_2018'

# ###  load packages
# library("sp") # for creating spatial points object
# library("raster") # for loading and manipulating climate and aridity rasters
# library("maptools") # for loading world map with country outlines
# library("vegan") # for calculating PCA
# 

#### remove old files from R environment
rm(list = ls())

#### read in trait demog data 

tr_data <- read.csv("phylo_ordered_trait_demog_03_2020.csv",
                    row.names = "X")



## Edit location data for Lathyrus vernus.  Current coordinate does not match description of the site
## in the paper and needs to be moved about 20km west to fall onto the land and not sea. 

tr_data$Lon[tr_data$Labels == "Lathyrus_vernus"] <- 16.7

### convert this to a spatial object using in 'coordinates' function

coordinates(tr_data) <- ~ Lon + Lat

## Tell R that the coordinate system is WGS 1984
proj4string(tr_data) <- '+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs' 

## load environmental datasets bioclim
## these are stored in a separate folder to keep the storage pressure off git hub
setwd(Clim_files)

###
files <- list.files(Clim_files)
files[1:7]
env_temp <- stack(files[1:7]) #makes the raster stack of bioclim temperature files
# and variation in precipitation

# Extracts the average value for these 7 variables within a 2.5k buffer around these points (i.e. 5k diameter). 
bc_5km_mean <- extract(x = env_temp, y = tr_data, method = 'simple', buffer = 2500, fun = mean, sp = TRUE) 

#summary(bc_5km_mean)

names(bc_5km_mean@data)
####  Read Aridity Index - stored in separate folder and is an adf file. 

setwd(Arid_files)
ai_image <- raster('w001001.adf')

### rescale ai_values to correct scale as per metadata 

ai_image2 <- ai_image/10000

### values indicate. 
# Aridity Index Value 		Climate Class
#  < 0.03			           Hyper Arid
# 0.03 – 0.2			           Arid
# 0.2 – 0.5			            Semi-Arid
# 0.5 – 0.65			        Dry sub-humid
# > 0.65			                 Humid

### extract data on aridity
ai_5km_mean <- extract(x = ai_image2, y = bc_5km_mean, method = 'simple', buffer = 2500, fun = mean, sp = TRUE) 

### rename the aridity variable because it's got a terrible name!
names(ai_5km_mean@data)[36] <- "AridityIndex"
#

ai_5km_mean@data$logAridityIndex <- log(ai_5km_mean@data$AridityIndex)

#### return to main working directory
setwd(Main_files)

### save spatial points file as an Robject 
save(ai_5km_mean, file = 'spatialPoints_withclim_03_2020.Rdata')

#dim(ai_5km_mean@data)

####

### save dataframe from points layer back into a sensible points object. 

data_env <- ai_5km_mean@data 
#class(data_env)
 
### give sensible names to the bioclim vars
names(data_env)[29:35] <- paste("bio", c(1:7), sep = "")

# create the PCA of tempature to make the biplot. Each site should only be used once in each 

cor(data_env[,29:35], method = "spearman")

### strong correlations use PCA axis instead as per Coutts et al. 


temp_pca <- rda(data_env[,29:35], scale = T)
plot(temp_pca, display = "sp")

scores(temp_pca)$species

#                PC1         PC2
# bio1  1.71859594  0.53068737
# bio2 -0.09963527  1.69174870
# bio3  1.74582074  0.15934699
# bio4 -1.70282816  0.56932903
# bio5  1.10457137  1.37459191
# bio6  1.82449038 -0.04826026
# bio7 -1.54810444  0.95761342

summary(temp_pca)

# Importance of components:
#                         PC1    PC2     PC3      PC4       PC5       PC6
# Eigenvalue            4.7567 1.8611 0.31310 0.065306 0.0025098 0.0013096
# Proportion Explained  0.6795 0.2659 0.04473 0.009329 0.0003585 0.0001871
# Cumulative Proportion 0.6795 0.9454 0.99012 0.999454 0.9998129 1.0000000

#### PCA1 explains 68% of variance

### Add temperature PCA axis 1 to dataset

data_env$Temp_PCA1 <- scores(temp_pca)$sites[,1]

cor.test(data_env$Temp_PCA1, data_env$logAridityIndex, method = "spearman" )

# Spearman's rank correlation rho
# 
# data:  data_env$Temp_PCA1 and data_env$logAridityIndex
# S = 80786, p-value = 0.6397
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.05314092 

# No correlation 

#### create dataframe with only env_vars of interest and output
names(data_env)


write.csv(data_env, "data_with_needed_env_vars_03_2020.csv", row.names = TRUE)
####

