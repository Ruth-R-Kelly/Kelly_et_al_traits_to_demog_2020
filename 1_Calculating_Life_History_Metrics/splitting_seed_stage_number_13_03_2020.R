### Last edit 13/03/2020

# This code file was first written by Ruth Kelly in May 2017 to split annual matrices 
# in Compadre (computed now as per 'overview_code_03_2020'), into separate objects based on 
# number of seed stages. And remove species they do not have at least two stages other than 
# seed-bank.

#Prior to this code matrices have already been checked for various criteria including 
#Ergodicity, Irreducibility, Primitivity, Divisibility, and further criteria detailed 
#in the overview Rmarkdown file.  .

## The input file for this file is :
##      "all_selected_matrices_merged_per_population_13_03_2020.RData"

#The new R objects are outputted as: 
# No seedbank recorded = "annual_no_Seedbank_03_2020.RData"
# 1 seedbank stage recorded = "annual_1yr_Seedbank_03_2020.RData"
# 2 seedbank stage recorded = "annual_2yr_Seedbank_03_2020.RData"
# 3 seedbank stage recorded = "annual_3yr_Seedbank_03_2020.RData"


### remove any files in the background
rm(list=ls())

#### read in data #### 
an_data <- readRDS("all_selected_matrices_merged_per_population_13_03_2020.RData")

#### A brief check that the compadre object looks sensible! #### 

names(an_data)

# [1] "metadata"    "mat"         "matrixClass" "version"  

#### make separate objects for metadata, mat and matrixClass ####

metadata1 <- an_data$metadata
length(unique(metadata1$SpeciesAccepted))

metadata1$UID1 <- seq(1, nrow(metadata1), 1)
##

mat1 <- an_data$mat
matrixClass1 <- an_data$matrixClass

### create new UID for metadata to store new row order
metadata1$UID1 <- seq(1, nrow(metadata1), 1)

#####
dim(metadata1)
## 254  52

length(unique(metadata1$SpeciesAccepted))
# 94

#### Next split dataset into groups based on the number of seed stages. ####

### make a vector to count number of seed stages

prop1 <-c()

for(i in 1:length(matrixClass1)) {
prop1[i] <- sum(matrixClass1[[i]]$MatrixClassOrganized %in% "prop")
}

summary(as.factor(prop1))
#  0   1    2   3 
# 137 102   8   7 

#### take metadata of those with no seed stages, and check matrix dimensions are a minimum of 2 

metadata_NS <- metadata1[which(prop1 == 0 ),]
summary(metadata_NS$MatrixDimension)
#      Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     2.000   4.000   5.000   6.058   7.000  17.000

mat_NS <- mat1[metadata_NS$UID1]
matrixClass_NS <- matrixClass1[metadata_NS$UID1]

###### Take metadata of those with 1 seed stages

metadata_N1 <- metadata1[which(prop1 == 1 ),]
summary(metadata_N1$MatrixDimension)
# 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  2.000   3.000   6.000   5.784   7.000  18.000 

mat_N1 <- mat1[metadata_N1$UID1]
matrixClass_N1 <- matrixClass1[metadata_N1$UID1]


#### check that the seed stage is listed first in the matrixClass

prop2 <- c()

for(i in 1:length(matrixClass_N1)) {
prop2[i] <- matrixClass_N1[[i]]$MatrixClassOrganized[[1]] == "prop"
}

summary(as.factor(prop2))

# FALSE  TRUE 
#   20    82 

unique(metadata_N1$SpeciesAccepted[which(prop2 == FALSE)])

# [1] "Arenaria serpyllifolia"  "Myosotis ramosissima"    "Saxifraga tridactylites"
# [4] "Veronica arvensis"  

### check matrix class object to see the the 'prop' stage is a seed stage in all cases.
checkmclass <- matrixClass_N1[which(prop2 == FALSE )] 
# checkmclass
### Select matrices where the propagule stage is second

mat_check_N1 <-mat_N1[which(prop2 == FALSE )] 

#### create a function to sort these matrices back around the correct way

flip_mat <- function(mat1) {mat1[nrow(mat1):1, ncol(mat1):1,drop=FALSE] }

#### use a loop to correct these matrices.. 

mat_flip_t1 <- mat_check_N1

for(j in 1:length(mat_check_N1)) {
  
  for(i in 1:length(mat_check_N1[[j]])) {
    mat_flip_t1[[j]][[i]] <- flip_mat(mat_check_N1[[j]][[i]])
  }
}

#### flip the rows in the matrix class object to reflect this..

flip_row <- function(matr) {matr[nrow(matr):1, ,drop=FALSE] }

matrixClass_flip <- checkmclass

for(i in 1:length(checkmclass)) {
  matrixClass_flip[[i]] <- flip_row(checkmclass[[i]])
}

##########  now add these back to the matrix dataset for species with 1 seed stage.. 

x1 <- which(prop2 == F)

for(i in x1) {
  matrixClass_N1[[i]] <- matrixClass_flip[[which(x1 == i)]] 
} 

####### 

for(i in x1) {
  mat_N1[[i]] <- mat_flip_t1[[which(x1 == i)]] 
} 

###########  final version with 1 seed stage, all at the start of matrix = 
# mat_N1, matrixClass_N1, and metadata_N1




#### 
############  Now 2 seed stages

metadata_N2 <- metadata1[which(prop1 == 2 ),]

summary(metadata_N2$MatrixDimension)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  5.00    5.75    6.00    5.75    6.00    6.00 

### make matching matrixclass and matrix datasets

matrixClass_N2 <- matrixClass1[metadata_N2$UID1]
mat_N2 <- mat1[metadata_N2$UID1]

#### check whether seed stages are in the first two rows of the matrix 

prop3 <- c()
prop4 <- c()

for(i in 1:length(matrixClass_N2)) {
prop3[i] <- matrixClass_N2[[i]]$MatrixClassOrganized[[1]] == "prop"
prop4[i] <- matrixClass_N2[[i]]$MatrixClassOrganized[[2]] == "prop"
}

summary(as.factor(prop3))
# TRUE 
# 8 
summary(as.factor(prop4))

# TRUE 
#   8 

##### three seed stages

metadata_N3 <- metadata1[which(prop1 == 3 ),]
summary(metadata_N3$MatrixDimension)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    7       7       7       7       7       7  

metadata_N3 <- metadata_N3[metadata_N3$MatrixDimension >= 3,]
matrixClass_N3 <- matrixClass1[metadata_N3$UID1]

### Check that seed stages are in first three positions 
propa <- c()
propb <- c()
propc <- c()

for(i in 1:length(matrixClass_N3)) {
propa[i] <- matrixClass_N3[[i]]$MatrixClassOrganized[[1]] == "prop"
propb[i] <- matrixClass_N3[[i]]$MatrixClassOrganized[[2]] == "prop"
propc[i] <- matrixClass_N3[[i]]$MatrixClassOrganized[[3]] == "prop"
}

summary(as.factor(propa))
summary(as.factor(propb))
summary(as.factor(propc))

###all true create matching matrix object

mat_N3 <- mat1[metadata_N3$UID1]

#### Save matrices for future use.. ####

#### Outputting data no seed stages ####

NS1 <- list()
NS1[[1]] <- metadata_NS
NS1[[2]] <- mat_NS
NS1[[3]] <- matrixClass_NS
NS1[[4]] <- an_data$version


names(NS1) <- c("metadata", "mat", "matrixClass", "version")
saveRDS(NS1,
file =  "annual_no_Seedbank_03_2020.RData")

############# one seed stage #####

SS1 <- list()
SS1[[1]] <- metadata_N1
SS1[[2]] <- mat_N1
SS1[[3]] <- matrixClass_N1
SS1[[4]] <- an_data$version


names(SS1) <- c("metadata", "mat", "matrixClass", "version")
saveRDS(SS1,
file = "annual_1yr_Seedbank_03_2020.RData")


########### two seed stages ################

SS2 <- list()
SS2[[1]] <- metadata_N2
SS2[[2]] <- mat_N2
SS2[[3]] <- matrixClass_N2
SS2[[4]] <- an_data$version


names(SS2) <- c("metadata", "mat", "matrixClass", "version")

saveRDS(SS2,
file =  "annual_2yr_Seedbank_03_2020.RData")

################## three seed stages #######################

SS3 <- list()
SS3[[1]] <- metadata_N3
SS3[[2]] <- mat_N3
SS3[[3]] <- matrixClass_N3
SS3[[4]] <- an_data$version


names(SS3) <- c("metadata", "mat", "matrixClass", "version")

saveRDS(SS3,
file =  "annual_3yr_Seedbank_03_2020.RData")

#### fini #### 

