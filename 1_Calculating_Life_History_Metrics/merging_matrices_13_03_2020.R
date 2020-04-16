### last edit 13/03/2020

### merging matrices from different years in the same population.  
### Then checking ergodicity, primitivity and irreducibility and removing matrices that do not meet these criteria. 

### function 'mean' from popdemo is needed for calculated cell by cell averaging 
## accross matrices within a list.. 

rm(list = ls())

# library("popbio")
# library("dplyr")
# library("popdemo")
### load data

data1 <- readRDS("ammended_annual_compadre_data_13_03_2020.RData")


####  create objects for the metadata and the matrices 
metadata1 <- data1$metadata
metadata1$UID1 <- seq(1, nrow(metadata1),1)

mat1 <- data1$mat

length(unique(metadata1$SpeciesAuthor))
### 99

##########  Get metadata for individual matrices, 
## these are the ones we will merge from

metadata_indiv <- metadata1[metadata1$MatrixComposite == "Individual",]
mat_indiv <- mat1[metadata_indiv$UID1]

length(unique(metadata_indiv$SpeciesAuthor))
##55

### give each population within each study a unique ID

PopID  <- paste(metadata_indiv$SpeciesAuthor, metadata_indiv$MatrixPopulation, sep = "_")
## PopID
PopUnique <- unique(PopID)

length(PopUnique)
## 177

####  use loop to get mean matrices for each population within study combination. 

Allmats <- list()
newSurvivalIssue <- c()
n_indiv_mats <- c()

for(j in 1:length(PopUnique)) {
  
set1 <- which(PopID == PopUnique[j])

matU <- list()
matF <- list()
matC <- list()


for(i in 1:length(set1)) {
matU[[i]] <- mat_indiv[[set1[i]]]$matU
matF[[i]] <- mat_indiv[[set1[i]]]$matF
matC[[i]] <- mat_indiv[[set1[i]]]$matC
}

setMat <- list()

setMat$matA <- mean(matU) + mean(matF) + mean(matC) 
setMat$matU <- mean(matU)
setMat$matF <- mean(matF)
setMat$matC <- mean(matC)

newSurvivalIssue[j] <- max(colSums(setMat$matU))
n_indiv_mats[j] <- length(set1)
Allmats[[j]] <- setMat

}

length(Allmats)
#177

#### "Allmats" is our new vector of population matrices, from individual studies populations

### Next match this back to the metadata with one row for each unique species
### population pair (i.e. PopID)

names(metadata_indiv)

# [1] "SpeciesAuthor"          "SpeciesAccepted"        "CommonName"            
# [4] "Genus"                  "Family"                 "Order"                 
# [7] "Class"                  "Phylum"                 "Kingdom"               
# [10] "OrganismType"           "DicotMonoc"             "AngioGymno"            
# [13] "Authors"                "Journal"                "YearPublication"       
# [16] "DOI.ISBN"               "AdditionalSource"       "StudyDuration"         
# [19] "StudyStart"             "StudyEnd"               "AnnualPeriodicity"     
# [22] "NumberPopulations"      "MatrixCriteriaSize"     "MatrixCriteriaOntogeny"
# [25] "MatrixCriteriaAge"      "MatrixPopulation"       "Lat"                   
# [28] "Lon"                    "Altitude"               "Country"               
# [31] "Continent"              "Ecoregion"              "StudiedSex"            
# [34] "MatrixComposite"        "MatrixTreatment"        "MatrixCaptivity"       
# [37] "MatrixStartYear"        "MatrixStartSeason"      "MatrixStartMonth"      
# [40] "MatrixEndYear"          "MatrixEndSeason"        "MatrixEndMonth"        
# [43] "MatrixSplit"            "MatrixFec"              "Observation"           
# [46] "MatrixDimension"        "SurvivalIssue"          "UID1" 
# [49] "checkY"                 "checknote"              "habitatAuthor" 

##### summarise the variables that will differ, e.g. years for reinsertion into the
### dataframe and choose a single UID from each to link back and select one row 
### per study population from the main metadata_indiv

indiv_summary <- as.data.frame(summarise(group_by(metadata_indiv, SpeciesAuthor, 
                                               MatrixPopulation),
                                         UID_match =  min(UID1),
                                         Max_year = max(MatrixEndYear,na.rm = T),
                                         Min_year = min(MatrixStartYear, na.rm = T)
                                         ))

summary(indiv_summary)
which(is.na(indiv_summary$SpeciesAuthor))
which(is.na(indiv_summary$MatrixPopulation))
               
nrow(indiv_summary)
#### arange this in order of UID
indiv_summary <- arrange(indiv_summary, UID_match)

### now select rows matching the UID from the metadata_indiv

metadata_match1 <- metadata_indiv[metadata_indiv$UID1 %in% indiv_summary$UID_match,]
which(is.na(metadata_match1$SpeciesAuthor))

nrow(metadata_match1)

### check that these UID's all match up 
which(metadata_match1$UID1 != indiv_summary$UID_match)
#####

summary(as.factor(metadata_match1$SpeciesAuthor))

## Substitute in the summarised variables for matrix start and matrix end 
## as appropriate

metadata_match1$MatrixStartYear <- indiv_summary$Min_year
metadata_match1$MatrixEndYear <- indiv_summary$Max_year

#####

### spot checking that the metadata still matches the matrices.. 

# Allmats[[177]]$matU
# metadata_match1[177,]
# metadata1[metadata1$UID1 == metadata_match1$UID1[24],]
# 
# PopUnique[177]
# metadata_match1[177,]
# metadata1$UID[metadata1$SpeciesAuthor == "Bothriochloa_ischaemum" &
#         metadata1$MatrixPopulation == "Pedernales Falls State Park; Doeskin Ranch; Eckhardt; Nagel; Simons; Victoria tracts; Freeman Ranch; Shield Ranch"]


### etc.. these seem good. 

### Add new survival issue column to metadata_match1
summary(newSurvivalIssue)
length(metadata_match1$SurvivalIssue)

metadata_match1$SurvivalIssue <- newSurvivalIssue


######  bind metadata_match1 with mean and pooled matrix info in metadata1

meta_mp <- metadata1[metadata1$MatrixComposite != "Individual",]

mat_mp <- mat1[meta_mp$UID1]

new_survival_other <- c()
for(i in 1:length(mat_mp)) {
  
new_survival_other[i] <-  max(colSums(mat_mp[[i]]$matU)) } 
                              
summary(new_survival_other)

meta_mp$SurvivalIssue <- new_survival_other

### check we haven't lost any rows!
nrow(metadata1)
sum(nrow(meta_mp),nrow(metadata_indiv))

metadata_fin <- rbind(metadata_match1, meta_mp)
#  ## make sure there are no lost species!

length(unique(metadata_fin$SpeciesAuthor))
length(unique(metadata1$SpeciesAuthor))


summary(metadata_fin$SurvivalIssue)

### all good 

### now bind the matrices into one object

mats_fin <- append(Allmats, mat_mp)

metadata_fin$UID_new <- seq(1, nrow(metadata_fin), 1)

 
### quick check that matrices all still match the metadata rows by matching the 
### SurvivalIssue calculated from matrices with the metadata_fin$SurvivalIssue.

check_mats <- c()

for(i in 1:nrow(metadata_fin)) {
check_mats[i] <- max(colSums(mats_fin[[i]]$matU))  == metadata_fin$SurvivalIssue[i]

}

which(check_mats == "FALSE")

####### check ergodicity, primitivity and irreducibility.. 


# #### create vector of NA's in which to store info about ergodicity.. 
# 


erg1 <- rep(NA,length(mats_fin))


#### here try catch forces function to continue even if errors arise for a given 'i', 
### this is why the NA's are needed in the definition of the vector 'erg1'
for(i in 1:length(mats_fin)) {
 tryCatch({
              res1 <- isErgodic(mats_fin[[i]]$matA)
              erg1[i] <- res1
   }
 , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
 }

 length(erg1)
## 261

ergT <- which(erg1 == TRUE)
length(ergT)
## [1] 250


meta_erg <- metadata_fin[ergT,] 
 
length(unique(metadata_fin$SpeciesAccepted))
# 96
# 
length(unique(meta_erg$SpeciesAccepted))
## 94
# 
# #### we lose 2 species to non-ergodicity

# #### check irreducible.. #### 
# 
mat1_ERG1 <- mats_fin[meta_erg$UID_new]
# 
# 

red1 <- rep(NA,length(mat1_ERG1))
length(red1)
# 250


# 
 for(i in 1:length(mat1_ERG1)) {
   res2 <- isIrreducible(mat1_ERG1[[i]]$matA)
   red1[i] <- res2
        }
# 
redT <- which(red1 == TRUE)
length(redT)
# # 246
# 
redF <- which(red1 == FALSE)
redF
### [1] 109 110 199 214

# ### select metadata only for irreducible matrices 
meta_erg_red <- meta_erg[redT,]

# 
# 
# 
length(unique(meta_erg$SpeciesAccepted))
# 94
# 
length(unique(meta_erg_red$SpeciesAccepted))
# 92
# 
# #### we lose 2 species here

# ###################
# 
# #### check primitivity.. #### 
# 
# 
mat1_ERG_RED1 <- mats_fin[meta_erg_red$UID_new]
length(mat1_ERG_RED1)

# 
# 
prim1 <- rep(NA,length(mat1_ERG_RED1))
length(prim1)
# 246
# 
for(i in 1:length(mat1_ERG_RED1)) {
   res3 <- isPrimitive(mat1_ERG_RED1[[i]]$matA)
   prim1[i] <- res3
 }
# 
 primT <- which(prim1 == TRUE)
 length(primT)
# 244
# ### 
 
 primF <- which(prim1 == FALSE)
primF

###[1]  41 231
 
# #####
Metadata_fin2 <-droplevels(meta_erg_red[primT,])
nrow(Metadata_fin2)
# 244
# 
length(unique(Metadata_fin2$SpeciesAccepted))
# 90 species meet all the above criteria
# #####
# 
length(unique(Metadata_fin2$SpeciesAccepted[Metadata_fin2$MatrixComposite == "Individual"]))
#[1] 51
length(unique(Metadata_fin2$SpeciesAccepted[Metadata_fin2$MatrixComposite == "Mean"]))
#  36
length(unique(Metadata_fin2$SpeciesAccepted[Metadata_fin2$MatrixComposite == "Pooled"]))
# 3
# 
# 
# ################ save ergodic, reducible, primitive R matrices for future use.. 
annual_metadata <- Metadata_fin2
annual_mats <- mats_fin[Metadata_fin2$UID_new]
length(annual_mats)

matrixClass <- data1$matrixClass
annual_matrixClass  <- matrixClass[Metadata_fin2$UID1]
# length(annual_matrixClass)

dim(annual_metadata)
# ## [1] 244  52
# 
# ################
length(unique(annual_metadata$SpeciesAccepted))
# ######### spot checking to make sure everything still matches..

# annual_mats[[235]]
# annual_matrixClass[[235]]
# annual_metadata[235,]


#### write matrices, matrixClass and metadata to new object ####

#######

obs1 <- list()
obs1[[1]] <- annual_metadata
obs1[[2]] <- annual_mats
obs1[[3]] <- annual_matrixClass
obs1[[4]] <- data1$version

# 
names(obs1) <- c("metadata", "mat", "matrixClass", "version")
# 
saveRDS(obs1,
         file =  "selected_annual_matrices_merged_per_population_13_03_2018.RData")

#############################

