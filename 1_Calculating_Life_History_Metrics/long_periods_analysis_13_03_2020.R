### R script last editted 13_03_2020 - Ruth Kelly

### This script checks (and where necessary converts the longer period matrices to annual 
### transition matrices) for species in Compadre vr. 4.0.1 where AnnualPeriodicity < 1 and 
### for which trait data is available to us for at least 4 traits, 
### (Height, seed mass, LA and LMA).

### these species are: 

## [1] "Cornus_florida"        "Abies_concolor"        "Acer_saccharum"       
## [4] "Pinus_ponderosa"       "Prioria_copaifera"     "Shorea_leprosula_2"      
##  

### then adds the metadata and matrices to the annual and seasonal species 
### compiled in "selected_matrices_with_seasonal_merged_per_population_13_04_2018.RData"

### final dataset is outputted as. . 
### "all_selected_matrices_merged_per_population_16_04_2018.RData"

rm(list = ls())

# library("pracma")
# library("dplyr")
# library("popdemo")


load("COMPADRE_v.4.0.1.RData")


metadata1 <- compadre$metadata
metadata1$UID1 <- seq(1:nrow(metadata1))
mat1 <- compadre$mat


specieslist <- c("Cornus_florida" , "Abies_concolor", "Acer_saccharum",       
                 "Pinus_ponderosa", "Prioria_copaifera" , "Shorea_leprosula_2")    
 
metadata1 <- metadata1[metadata1$SpeciesAuthor %in% specieslist,]

###

unique(metadata1$SpeciesAuthor)
unique(metadata1$MatrixTreatment)
# [1] "Unmanipulated" 


mat1 <- compadre$mat[metadata1$UID1]
matrixclass1 <- compadre$matrixClass[metadata1$UID1]


table(metadata1$SpeciesAuthor, metadata1$AnnualPeriodicity)

# #                     0.1 0.2 0.2; 1 0.333; 0.2 0.4 0.5
# Abies_concolor          0  17      0          0   0   0   ## Correct Annual Periodicity to 0.2
# Acer_saccharum          0   0     10          0   0   0   ## Acer saccharum errors see below
# Cornus_florida          0   1      0          0   0   0   ## OMIT - Cornus florida errors
# Pinus_ponderosa         0   1      0          0   0   0   ## Pinus ponderosa correct to 0.2
# Prioria_copaifera       0   0      0          3   0   0   ## Prioria copaifera errors see below
# Shorea_leprosula_2      3   0      0          0   0   0   ## OMIT - I don't have raw data to check


####


### species checked with original documents, see notes above right.  

AS <- metadata1[metadata1$SpeciesAccepted =="Acer saccharum",]
AS$DOI.ISBN
### Here are some likely errors in this species as matrices in the  original paper are 
### from 1951 to 1988 and 1988-2001 (i.e. annual periodicity 1/37 and 1/13).
### use the matrix from 1988-2001 only and create a new mean matrix. 


AS2 <- droplevels(filter(AS, MatrixStartYear == 1988))
AS2$AnnualPeriodicity <- 1/13

AS2$Lat ## Correct
AS2$Lon ## this should be -88.17
AS2$Lon <- -88.17
## Habitat should be forest. 


AS2$checkY <- "Y_edits"
AS2$checknote <- "fixed longitude and annual periodicity"
AS2$habitatAuthor <- "forest"

### similar for Prioria copaifera , although I can't find original matrices, 
## and can't figure out how individual matrices were calculated from the manuscript.  
## Use only individual matrices from 1985 to 1990 and in set periodicity to 0.2

PC <- metadata1[metadata1$SpeciesAccepted == "Prioria copaifera",]
PC$DOI.ISBN


PC2 <- droplevels(filter(PC, MatrixStartYear == 1985))
PC2$AnnualPeriodicity <- 0.2

PC2$Lat # okay 
PC2$Lon # should be -79.85
PC2$Lon <- -79.85


PC2$checkY <- "Y_edits"
PC2$checknote <- "fixed longitude and annual periodicity"
PC2$habitatAuthor <- "forest"

table(metadata1$SpeciesAccepted, metadata1$MatrixComposite)
#                   Individual Mean Pooled Seasonal
# Abies concolor            12    5      0        0
# Acer saccharum             6    4      0        0
# Cornus florida             0    1      0        0
# Pinus ponderosa            1    0      0        0
# Prioria copaifera          2    1      0        0
# Shorea leprosula           0    3      0        0

metadata1[metadata1$SpeciesAccepted == "Pinus ponderosa",]

### only Pinus ponderosa has only 1 individual study.  This is 12 years and two populations so is really a mean. 

PP <- metadata1[metadata1$SpeciesAccepted == "Pinus ponderosa",]

PP$DOI.ISBN ## should be 10.1111/j.1365-2745.2005.01007.x
PP$Lat

PP$Lon ## should be -119.9
PP$Lon <- -119.9

PP$checkY <- "Y_edits"
PP$checknote <- "fixed longitude"
PP$habitatAuthor <- "forest"

### Get matrices for "Abies concolor" . Use mean matrices as these are calculated per location
### except for the first one which should be omitted.  Also add "matrixPopulation" 
### "Hodgedon Meadows" 

ac <- filter(metadata1, SpeciesAccepted == "Abies concolor"
                    & MatrixComposite == "Mean")

ac <- ac[ac$MatrixPopulation != "Hodgedon Meadows; Crystal Road; Suwanee Creek; SEGI Conifer; Log Creek",]

ac_hog <- filter(metadata1, SpeciesAccepted == "Abies concolor"
                 & MatrixPopulation == "Hodgedon Meadows")

ac <- rbind(ac, ac_hog)

#ac$MatrixPopulation
#ac$Lon

#cbind(ac$MatrixPopulation, ac$Lon)
ac$Lon <- c(-118.8, -118.8, -118.7,-118.7, -119.9)

ac$AnnualPeriodicity <- as.numeric(ac$AnnualPeriodicity)

ac$checkY <- "Y_edits"
ac$checknote <- "fixed longitude, took mean populations at each location"
ac$habitatAuthor <- "forest"

###   I've omitted PC2 because the U matrix won't solve in the rootm 

metadata4 <- rbind(ac, AS2, PP)

#### get matching matrices.

mat_long <- compadre$mat[metadata4$UID1]
matrixClass <- compadre$matrixClass[metadata4$UID1]

###


### check to make sure none of these species have clonality matrices
# for(i in 1:length(mat_long)) {
#   print(sum(mat_long[[i]]$matC))
#   }


periodicity1 <- 1/as.numeric(metadata4$AnnualPeriodicity)




###  create new object same length and size as mat_long


mat_long1 <- mat_long

#### 
for(i in 1:length(periodicity1)) {
x1 <- rootm(mat_long[[i]]$matA,periodicity1[i])$B  
x1[x1 < 0] <- 0
mat_long1[[i]]$matA <- x1
x2 <- rootm(mat_long[[i]]$matU, periodicity1[i])$B
x2[x2 < 0] <- 0
mat_long1[[i]]$matU <- x2
x3 <- mat_long1[[i]]$matA - mat_long1[[i]]$matU 
x3[x3 < 0] <- 0
mat_long1[[i]]$matF <- x3
print(i)
}



## check ergodicity, primitivity, irreducibility and seedbank. 

## create new UID for metadata4
metadata4$UID1 <- seq(1,nrow(metadata4),1)

summary(metadata4)

### recalculate SurvivalIssue based on new matrices

SI <- c()

for(i in 1:length(mat_long1)) {
  SI[i] <- max(colSums(mat_long1[[i]]$matU))
print(i)
  }

### add this to metadata4

metadata4$SurvivalIssue <- SI


### note Treatment == "Unmanipulated" was already checked above.
full_meta2 <- droplevels(subset(metadata4,
                                  MatrixSplit == "Divided" 
                                & MatrixFec == "Yes"
                                & SurvivalIssue <= 1))


nrow(full_meta2)
length(unique(full_meta2$Accepted_species))
### all qualify here.. 

### now check Erg, prim and irreducibility


is_ergodic <- vector()
is_primitive <- vector()
is_irreducible <- vector()
all_true <- vector()

for(i in 1:nrow(full_meta2)){
    is_ergodic[i] <- isErgodic(mat_long1[[i]]$matA)
    is_primitive[i] <- isPrimitive(mat_long1[[i]]$matA)
    is_irreducible[i] <- isIrreducible(mat_long1[[i]]$matA)
  all_true[i] <- all(c(is_ergodic[i],is_primitive[i],is_irreducible[i]) == TRUE)
} 

all_true


######################  keep all ###

### Add matrices metadata and matrixclass info to other annual matrices as saved in the 
### previous code for seasonal matrices  "seasonal_matrices_05_04_2018.R" and 
### outputted from there as "selected_matrices_with_seasonal_merged_per_population_13_04_2018.RData"


### add this to outputted data from merging step -  "selected_annual_matrices_merged_per_population_27_03_2018.RData"

sel_data <- readRDS("selected_matrices_with_seasonal_merged_per_population_13_03_2020.RData")
names(sel_data)

metadata4$UID_new <- max(sel_data$metadata$UID_new) + 1

sel_data$metadata <- rbind(sel_data$metadata, metadata4)
sel_data$mat <- c(sel_data$mat,mat_long1)


sel_data$matrixClass <- c(sel_data$matrixClass, matrixClass)     

#### write to new object

saveRDS(sel_data,
     file = "all_selected_matrices_merged_per_population_13_03_2020.RData" )

###