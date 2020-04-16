### R last edit script 13_03_2020 - Ruth Kelly

### This script checks (and where necessary converts the seasonal matrices to annual 
### transition matrices) for species in Compadre vr. 4.0.1 where AnnualPeriodicity > 1 and 
### for which trait data is available to us for at least 4 traits, 
### (Height, seed mass, LA and LMA).

### These species are: 

# [1] "Acacia_suaveolens"     "Agropyron_repens"      "Asarum_canadense"     
# [7] "Bromus_tectorum"       "Calluna_vulgaris"      "Daucus_carota"        
#       "Fagus_grandifolia"     "Fragaria_vesca"       
#   "Lactuca_serriola"      "Melampyrum_pratense"  
# [16] "Mimulus_guttatus"      "Salsola_australis"     "Tragopogon_dubius"   

### Only one of these species "Fagus_grandifolia"  turns out to be usable (with corrections)

### this is then appended to the merged annual matrice which were outputted in step 2. 

### and saved as 
## "selected_annual_matrices_merged_per_population_13_03_2020.RData"



####---- Begin code proper ----


###  remove R objects from memory
rm(list = ls())

### load packages
# library("expm")
# library("pracma")
# library("dplyr")
# library("popdemo")

### Compadre

load("COMPADRE_v.4.0.1.RData")


metadata1 <- compadre$metadata
metadata1$UID1 <- seq(1:nrow(metadata1))

specieslist <- c("Acacia_suaveolens", "Agropyron_repens",  "Asarum_canadense" ,   
                 "Bromus_tectorum"  ,
                    "Calluna_vulgaris" ,  "Daucus_carota",     
                      "Fagus_grandifolia", "Fragaria_vesca",       
                         "Lactuca_serriola", "Melampyrum_pratense",  
                            "Mimulus_guttatus", "Salsola_australis","Tragopogon_dubius" )

metadata1 <- metadata1[metadata1$SpeciesAuthor %in% specieslist, ]

###
unique(metadata1$SpeciesAuthor)
unique(metadata1$MatrixTreatment)
# [1] "Unmanipulated" "Shrub removal" "Shelter"     
## select Unmanipulated

metadata1 <- metadata1[metadata1$MatrixTreatment == "Unmanipulated", ]

### select species

mat1 <- compadre$mat[metadata1$UID1]
matrixclass1 <- compadre$matrixClass[metadata1$UID1]

table(metadata1$SpeciesAuthor, metadata1$AnnualPeriodicity)

# #                       2   3   4 4-5; 2   5   6
# Acacia_suaveolens       0   7   0      0   0   0 ## don't use
# Agropyron_repens        0   0   0      0   0   8 ## don't use
# Asarum_canadense       27   0   0      0   0   0 ## don't use
# Bromus_tectorum         0  30   0      0   0   0 ## don't use
# Calluna_vulgaris        0   0   0      0   1   0 ## don't use
# Daucus_carota           0   1   0      0   0   0 ## don't use treatment should be annual mowing
# Fagus_grandifolia       3   0   0      0   0   0 ## use with correction see below
# Fragaria_vesca          0 111   0      0   0   0 ## don't use
# Melampyrum_pratense     0   0  20      0   0   0 ## don't use
# Mimulus_guttatus        2   0   0      0   0   0 ## don't use
# Salsola_australis       1   0   0      0   0   0 ## don't use - annual mowing
# Tragopogon_dubius       0   3   0      0   0   0 ## don't use 

##### Necessary to reset the UID here so matrices and metadata match

metadata1$UID1 <- seq(1:nrow(metadata1))


#### Check species individually #### 

#Acacia_suaveolens - Matrix given as Annual matrix in Compadre matrix file. Periodicity in original paper was 3 -
#    Measured 1997 April/August, 1997 September, 1997 November and 1998 April/August
#    However, I can't use it in this study because our version only has 2 stages seed and adult. 

# Agropyron_repens  - Can't use this in this study as it has clonal matrices.  Also, digitisation of the
#   matrices is unclear to me as start end months don't match up and original book is unavailable.

# Asarum_canadense - Can't use this in this study as it has clonal matrices.  Matrices are digitized as annual matrices. Although, surveys were biannual. 

##


# Bromus_tectorum - ### seasonal matrices with 3 transitions  ### but oddly constructed and 
# I can't figure out how to convert to annual matrices because multiplication of the season matrices gives a 
# matrix of all zeros.. seems like some transition data is missing?  maybe?

# Calluna_vulgaris - This is an annual matrix. Do not use, includes clonal reproduction via layering.. 

# Daucus_carota - Annual matrix I think from the paper. Treatment should be annual mowing.

# Fagus_grandifolia - Calculations are across two different transition time periods.
### 6 years for closed and 8 years for open canopy models.. 
# online info - "American beech can take up to 40 years to begin producing seeds.
### Large crops are produced by 60 years and the tree's total lifespan may be up to 300
## years"

fg <- metadata1[metadata1$SpeciesAuthor == "Fagus_grandifolia" & metadata1$MatrixPopulation %in% c("Closed canopy", "Open canopy"),]
# fg$DOI.ISBN
# fg$Lat

### Longitude corrected from original paper. 
fg$Lon <- -84.33

fg_mat <- mat1[fg$UID1]
fg1 <- fg_mat
periodicity1 <- c(6,8)

#### 
for(i in 1:length(periodicity1)) {
  x1 <- rootm(fg_mat[[i]]$matA,periodicity1[i])$B  
  x1[x1 < 0] <- 0
  fg1[[i]]$matA <- x1
  x2 <- rootm(fg_mat[[i]]$matU, periodicity1[i])$B
  x2[x2 < 0.001] <- 0
  fg1[[i]]$matU <- x2
  x3 <- fg1[[i]]$matA - fg1[[i]]$matU 
  x3[x3 < 0] <- 0
  fg1[[i]]$matF <- x3
  print(i)
}



#### now use matrices fg1  and metadata fg for this species. Fagus_grandifolia


fg$checkY <- "Y_edits"
fg$checknote <- "fixed coorindates, adjusted matrices to annual periodicity"
fg$habitatAuthor <- "forest"


metadata_keep1 <- fg
mat_keep1 <- fg1

#

## Fragaria_vesca  do not use as only clonal reproduction is included. 

## Melampyrum_pratense  ### seasonal matrices - they don't look like they will cover full annual period
## DOI in metadata is incorrect.. 

## Mimulus_guttatus = these are seasonal matrices that can be used to calculate 
## the annual transition. However, data is a mix of experimental treatments and 
## greenhouse data. Exclude

### Salsola_australis  - annual matrix with errors - treatment should be annual mowing. 
### do not use. 


### Tragopogon_dubius - seasonal matrices - exclude, 
### these multiple to a matrix of all zeros ..

##  rename full datasets to match previous code


### correct survival issue to reflect changes to the matrices. 

SI <- c()

for(i in 1:length(mat_keep1 )){
  SI[i] <- max(colSums(mat_keep1 [[i]]$matU)) }

metadata_keep1$SurvivalIssue <- SI

### remove matrices where SI > 1

metadata_keep2 <- metadata_keep1[metadata_keep1$SurvivalIssue <1,]

mat2 <- mat_keep1[which(metadata_keep1$SurvivalIssue <1)]

###########  check ergodicity, primitivity, irreducibility and seedbank. 


### note Treatment == "Unmanipulated" was already checked above.
full_meta2 <- droplevels(subset(metadata_keep2,
                          MatrixSplit == "Divided" 
                         & MatrixFec == "Yes"
                         & SurvivalIssue <= 1))

nrow(full_meta2)
### 


### now Erg, prim and irred

isErgodic(mat2[[1]]$matA)
isPrimitive(mat2[[1]]$matA)
isIrreducible(mat2[[1]]$matA)

## all good. 

### add this to outputted data from merging step -  "selected_annual_matrices_merged_per_population_27_03_2018.RData"

sel_data <- readRDS("selected_annual_matrices_merged_per_population_13_03_2018.RData")

full_meta2$UID_new <- max(sel_data$metadata$UID_new) + 1

sel_data$metadata <- rbind(sel_data$metadata, full_meta2)
sel_data$mat <- c(sel_data$mat,mat2)

matrixclass1_seas <- matrixclass1[full_meta2$UID1]
sel_data$matrixClass <- c(sel_data$matrixClass, matrixclass1_seas)     

#### write to new object ####

saveRDS(sel_data,
     file =  "selected_matrices_with_seasonal_merged_per_population_13_03_2020.RData")

###