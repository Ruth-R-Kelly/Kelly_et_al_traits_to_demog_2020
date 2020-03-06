### Code used to extract suitable annual matrices from Compadre, and check data with original sources.
###  Matrices and metadata are ammended to match orignal sources where differences exist between 
### Compadre (VR 4.0.1) and original sources. Matrices are selected for inclusion only where they have
###  annual periodicity, 2 or more dimensions, and conducted in unmanipulated or control environments, 
### are divisible ,have a fecundity matrix and do not include clonality.

### Furthermore, species are only selected for inclusion where information is available for 4 or more traits.

### Here we extract all the individual matrices and average within populations later. 
### Checking for ergodicity, primitivity and irreducibility has been transferred to the next step after merging 
### of individual matrices into temporal mean matrices. 

### load Compadre version COMPADRE_v.4.0.1.Rdata

load("COMPADRE_v.4.0.1.RData")

## add libraries 
# library("popdemo")
# library("vegan")

#### examining compadre object ####

#names(compadre)
#compadre$version
# [1] "metadata"    "matrixClass" "mat"         "version"   

metadata1 <- compadre$metadata
matrixClass1 <- compadre$matrixClass
mat1 <- compadre$mat

metadata1$UID1 <- seq(1,nrow(metadata1),1)


#### number of matrices 
nrow(metadata1)
## 7024 -matrices are in this version

### number of species
length(unique(metadata1$SpeciesAccepted))
## 695 species..


### Below I've made some corrections to errors in the AnnualPeriodicity recorded in Compadre.
### the following species have annual matrices in the original papers, but are listed in Compadre
## as having AnnualPeriodicity between 0.2 and 2. 
### AnnualPeriodicity in Compadre should == "1"


sp1 <- c("Psidium_guajava", "Swietenia_macrophylla", "Shorea_leprosula", "Astrocaryum_mexicanum",
         "Atriplex_vesicaria", "Banksia_ericifolia", "Echium_vulgare",
         "Fumana_procumbens", "Araucaria_cunninghamii")


metadata1[which(metadata1$SpeciesAuthor %in% sp1),21] <- 1


###### ammendments to pinus strobus, all MatrixTreatments should be "Unmanipulated" based
### on original source

metadata1$MatrixTreatment[metadata1$SpeciesAuthor == "Pinus_strobus"] <- "Unmanipulated"

####  Minor ammendments to "Salsola_australis" based on original source
SA <- metadata1[metadata1$SpeciesAuthor =="Salsola_australis",]
SA
mat1[[SA$UID1]]$matA[2,1] <- 0
mat1[[SA$UID1]]$matA[2,7] <- 101
mat1[[SA$UID1]]$matA[4,1] <- 0.146
mat1[[SA$UID1]]$matA[4,2] <- 0.026
mat1[[SA$UID1]]$matF[2,1] <- 0
mat1[[SA$UID1]]$matF[2,7] <- 101
mat1[[SA$UID1]]$matU <- mat1[[SA$UID1]]$matA - mat1[[SA$UID1]]$matF


length(unique(metadata1$SpeciesAccepted))

##  Remove current mixed-sex matrix for "Brosimum alicastrum" and replace with female only 
### matrix, as this mixed-sex matrix is not 'irreducible'
length(unique(metadata1$SpeciesAccepted))

b1 <- metadata1[metadata1$SpeciesAccepted == "Brosimum alicastrum",]


mat1[[b1$UID1]]$matA <- as.matrix(read.csv("Brosimum_alicastrum_27_10_2017_female_matA.csv"))
mat1[[b1$UID1]]$matF <-as.matrix(read.csv("Brosimum_alicastrum_27_10_2017_female_matF.csv"))
mat1[[b1$UID1]]$matU <- mat1[[b1$UID1]]$matA - mat1[[b1$UID1]]$matF
mat1[[b1$UID1]]$matC <- mat1[[b1$UID1]]$matC[1:8,1:8]  

#####################  

matrixClass1[[b1$UID1]]  <- matrixClass1[[b1$UID1]] [c(1:6,8,10),] 
mat1[[b1$UID1]]

#####

matrixClass1[[b1$UID1]]$MatrixClassNumber <- seq(1,8,1)

length(unique(metadata1$SpeciesAccepted))


### ammendments to "Swallenia_alexandrae" matrices say divided 
### but U matrix supplied is empty.  U can be fixed by subtracting F from A.


SA2 <- metadata1[metadata1$SpeciesAuthor == "Swallenia_alexandrae",]

mat1[[SA2$UID1]]$matF[1,4] <- 0.162
mat1[[SA2$UID1]]$matU <- mat1[[SA2$UID1]]$matA -mat1[[SA2$UID1]]$matF
metadata1$SurvivalIssue[SA2$UID1] <- max(colSums(mat1[[SA2$UID1]]$matU))

##### Subsetting use individual matrices where available. 

### Note I have checked that these are okay to include ...

treats <- c("Unmanipulated", "Abandoned",  "Control", 
            "Ambient Condition")

types <-    c("Annual", "Epiphyte", "Fern", "Herbaceous perennial",
              "Liana", "Palm", "Shrub", "Succulent", "Tree") 

meta_sub <- droplevels(subset(metadata1,
                         MatrixSplit == "Divided" 
                         & MatrixFec == "Yes" 
                         & MatrixTreatment %in% treats
                         & MatrixComposite == "Individual"
                         & AnnualPeriodicity == "1"
                         # & SurvivalIssue <= 1.0
                         & OrganismType %in% types))

length(unique(meta_sub$SpeciesAccepted))
## 270

######

### Prunus serotina --- Use MatrixPopulations 'Light conditions'
### and shaded conditions

meta_PS <- droplevels(subset(metadata1,
                               MatrixSplit == "Divided" 
                               & MatrixFec == "Yes" 
                              # & MatrixComposite == "Individual"
                               & AnnualPeriodicity == "1"
                               # & SurvivalIssue <= 1
                              & SpeciesAccepted == "Prunus serotina"
                              & MatrixPopulation %in% c("Light Conditions", "Shaded Conditions")))


meta_sub <- rbind(meta_sub, meta_PS)

####
#### 
### Epipactis atrorubens ---  Use Matrix treatment 'moist site'.


meta_EA <- droplevels(subset(metadata1,
                             MatrixSplit == "Divided" 
                             & MatrixFec == "Yes" 
                             # & MatrixComposite == "Individual"
                             & AnnualPeriodicity == "1"
                             # & SurvivalIssue <= 1.0
                             & SpeciesAccepted == "Epipactis atrorubens"
                             & MatrixTreatment == "Moist site"))

### add to meta_sub
meta_sub <- rbind(meta_sub, meta_EA)

########################################################################################

nrow(meta_sub)
##2434

length(unique(meta_sub$SpeciesAccepted))
#272

########################################################################################

### Use mean matrices for species where no individual matrices exist.. 

x1 <- setdiff(metadata1$SpeciesAccepted, meta_sub$SpeciesAccepted)


length(x1)
# 423


#### make metadata only including species which don't have individual matrices
metadata1_means <- droplevels(metadata1[metadata1$SpeciesAccepted %in% x1,])
### check number of species matches expected.
length(unique(metadata1_means$SpeciesAccepted))

summary(metadata1_means$Lat)

lat_given <- metadata1_means[-which(is.na(metadata1_means$Lat)),]
no_lat_given <- metadata1_means[which(is.na(metadata1_means$Lat)),]

### give preference to datasets with lat and long, i.e. means of single populations 
### with georeferencing over means across populations.. 

meta_sub_means1 <- droplevels(subset(lat_given,
                              MatrixSplit == "Divided" 
                              & MatrixFec == "Yes" 
                              & MatrixTreatment %in% treats
                              & MatrixComposite == "Mean"
                              & AnnualPeriodicity == 1
                            #  & SurvivalIssue <= 1.0
                              & OrganismType %in% types))

length(unique(meta_sub_means1$SpeciesAccepted))
### [1] 144

#### check which species still have no matrices included. 

x2 <- setdiff(no_lat_given$SpeciesAccepted, meta_sub_means1$SpeciesAccepted)


length(x2)
## 115 

#### take mean matrices of suitable species from this set if they are not already 
## included. 


meta_sub_means2 <- droplevels(subset(no_lat_given,
                                     MatrixSplit == "Divided" 
                                     & MatrixFec == "Yes" 
                                     & MatrixTreatment %in% treats
                                     & MatrixComposite == "Mean"
                                     & AnnualPeriodicity == 1
                                     #& SurvivalIssue <= 1.0
                                     & OrganismType %in% types
                                     & SpeciesAccepted %in% x2))


length(unique(meta_sub_means2$SpeciesAccepted))
##35


### combine these two datasets
meta_sub_IM <- rbind(meta_sub, meta_sub_means1, meta_sub_means2)

length(unique(meta_sub_IM$SpeciesAccepted))
### 451

### Here the species "Dipsacus_sylvestris_2", and this is replicated once
### with an age matrix and once with a stage matrix 
## but for the same study and population.
## The stage matrix is the one we will use here. The age matrix has MatrixDimension == 8
## Remove that one. 

meta_sub_IM <- meta_sub_IM[-which(meta_sub_IM$SpeciesAuthor == "Dipsacus_sylvestris_2" & meta_sub_IM$MatrixDimension == 8),]
#####



### add any species where only pooled matrices are available. 

x1 <- setdiff(metadata1$SpeciesAccepted, meta_sub_IM$SpeciesAccepted)

length(x1)
# 244

meta_sub_pooled <- droplevels(subset(metadata1,
                                    MatrixSplit == "Divided" 
                                    & MatrixFec == "Yes" 
                                    & MatrixTreatment %in% treats
                                    & MatrixComposite == "Pooled"
                                    & AnnualPeriodicity == "1"
                                   # & SurvivalIssue <= 1
                                    & OrganismType %in% types
                                    & SpeciesAccepted %in% x1))

length(unique(meta_sub_pooled$SpeciesAccepted))
### [1] 20

### check no intersect in species names here
intersect(meta_sub_IM$SpeciesAccepted, meta_sub_pooled$SpeciesAccepted)

### combine these dataframes

meta_sub_IMP <- rbind(meta_sub_IM, meta_sub_pooled)

###### Remove "Persoonia_glaucescens" based on checking of original source, 
#### this data is from a seed burial experiment on burnt sites.

meta_sub_IMP <- meta_sub_IMP[-which(meta_sub_IMP$SpeciesAuthor == "Persoonia_glaucescens"),]

###################


#### checking species with matched trait data for errors

sp_to_check <- meta_sub_IMP
nrow(sp_to_check)
#3045
length(unique(sp_to_check$SpeciesAuthor))
## 510

### 

name_matches <- read.csv("all_trait_means10_10_2017_added_datasets_vr4.csv")

summary(name_matches)

name_matches2 <- name_matches

name_matches2[is.na(name_matches2)] <- 0

pa<- decostand(name_matches2[,-1], method = "pa")

rownames(pa) <- name_matches2[,1]

####

more4 <- name_matches2[,1][which(rowSums(pa) >3)]


head(more4)

###
x1 <- which(sp_to_check$SpeciesAccepted %in% more4)

### 1011

sp_to_check_in_try <- sp_to_check[x1,]

length(unique(sp_to_check_in_try$SpeciesAuthor))
## 180

study_list <- unique(sp_to_check_in_try$SpeciesAuthor)

study_list <- sort(study_list)


#### "Achillea_millefolium"

study_list[1]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[1],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[1]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[1]] <- "Experimental plot, treatment should = 'weeded'"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor == study_list[1]] <- "Experimental weeded field"
####### Other species in this study are A. odoratum and T. pratense 
 ### Trifolium_pratense

SA <- sp_to_check_in_try$SpeciesAuthor[which(sp_to_check_in_try$DOI.ISBN == CA$DOI.ISBN)]
which(study_list %in% SA)
# [1]   1  12 168

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor %in% SA] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor %in% SA] <- "Experimental plot, treatment should = 'weeded'"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor %in% SA] <- "Experimental weeded field"

###### "Achnatherum_calamagrostis"

study_list[2]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[2],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[2]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[2]] <- "The is a population spread model, life stages are not life stages!"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor == study_list[2]] <- "Badlands of Pyrennees, steep slopes with little vegetation"

#### "Actaea_spicata"

study_list[3]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[3],] 

### fix lat long
siteALon <- rep(17.07, 6)
siteBLon <- rep(17.05, 6)
Lon_new <- c(siteALon, siteBLon)

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[3]]  <- 58.1
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[3]]  <- Lon_new
###

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[3]]  <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[3]]  <- "lat long corrected"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor == study_list[3]] <- as.character(sp_to_check_in_try$Observation[sp_to_check_in_try$SpeciesAuthor == study_list[3]])

####  "Agrimonia_eupatoria"
study_list[4]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[4],] 
CA_full <- metadata1[metadata1$SpeciesAuthor == study_list[4],] 


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[4]]  <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[4]]  <- "Use as is"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor == study_list[4]] <- "road verges, adjacent to gravel road and arable fields"

######## "Agropyron_cristatum"
study_list[5]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[5],] 

CA_full <- metadata1[metadata1$SpeciesAuthor == study_list[5],]

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[5]] <- -107.88

### typos fixed as follows
CA$SurvivalIssue
CA$Authors
CA$MatrixPopulation
CA$Observation

mat1[CA$UID1][[1]]$matU[7,5]
mat1[CA$UID1][[1]]$matU[7,5] <- 0.2
mat1[CA$UID1][[1]]$matU[8,6] <- 0.1
mat1[CA$UID1][[1]]$matA <- mat1[CA$UID1][[1]]$matU + mat1[CA$UID1][[1]]$matF

mat1[CA$UID1][[2]]$matU
mat1[CA$UID1][[2]]$matU[7,5] <- 0.25
mat1[CA$UID1][[2]]$matU[8,6] <- 0.17
mat1[CA$UID1][[2]]$matA <- mat1[CA$UID1][[2]]$matU + mat1[CA$UID1][[2]]$matF

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[5]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[5]] <- "when matrix typos are fixed a surival issue arises, ?splitting of tussocks?, longitude corrected"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor == study_list[5]] <- "mixed-grass grasslands"

#### "Ailanthus_altissima"

study_list[6]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[6],] 

### fix longitude
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[6]] <- -90.5

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[6]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[6]] <- "Fixed longitude, Note only one population measured for one year"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor == study_list[6]] <- "Forest"

#########

####### Other species in this study  

# [1] "Alliaria_petiolata_3" "Allium_vineale"       "Narcissus_poeticus"  
# [4] "Lonicera_maackii"     "Rosa_canina"          "Rosa_multiflora"     
# [7] "Ailanthus_altissima"  "Cerastium_fontanum"   "Cerastium_pumilum"   
# [10] "Iris_germanica" 

SA <- sp_to_check_in_try$SpeciesAuthor[which(sp_to_check_in_try$DOI.ISBN == CA$DOI.ISBN)]
which(study_list %in% SA)
#  [1]   6   7   9  42  43  82  92  99 140 141
### check these now too

#### Alliaria_petiolata_3

study_list[7]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[7],] 

### fix longitude
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[7]] <- -90.5

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[7]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[7]] <- "Fixed longitude. Only one population measured for one year"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor == study_list[7]] <- "Forest"

######  Allium_vineale 
study_list[9]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[9],] 


### fix longitude
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[9]] <- -90.5

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[9]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[9]] <- "Fixed longitude. Only one population measured for one year"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor == study_list[9]] <- "Grassland"

####

######  Cerastium_fontanum 
study_list[42]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[42],] 

### fix longitude
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[42]] <- -90.5

#### does not match original, replace with ..

matA <- mat1[CA$UID1][[1]]$matA[1:2,1:2]
matA[1,1] <- 0.3
matA[1,2] <- 87.64
matA[2,1] <- 0.001
matA[2,2] <- 4.44

matF <- matA
matF[,1] <- 0 

matU <- matA - matF 

mat1[CA$UID1][[1]]$matA <- matA
mat1[CA$UID1][[1]]$matU <- matU
mat1[CA$UID1][[1]]$matF <- matF 
mat1[CA$UID1][[1]]$matC <- matrix(nrow = 2, c(0,0,0,0))
mat1[CA$UID1]

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[42]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[42]] <- "Fixed longitude. Changes made to match original, only one population measured for one year"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor == study_list[42]] <- "Grassland"

#### Cerastium_pumilum
study_list[43]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[43],] 

### fix longitude
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[43]] <- -90.5

#### does not match original, replace with ..

matA <- mat1[CA$UID1][[1]]$matA[1:2,1:2]
matA[1,1:2] <- c(0.66, 71.19)
matA[2,1:2] <- c(0.02, 0.80)

matU <- matA
matU[1:2,2] <- c(0, 0)

matF <- matA - matU
#####

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[43]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[43]] <- "Fixed longitude. Changes made to match original, only one population measured for one year"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor == study_list[43]] <- "Grassland"

###### Iris_germanica 
study_list[82]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor
                         == study_list[82],] 

### fix longitude
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[82]] <- -90.5

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[82]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[82]] <- "Fixed longitude. Incomplete matrix in original paper, omit"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor == study_list[82]] <- "Grassland"

######  Lonicera_maackii
study_list[92]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor
                         == study_list[92],] 

### fix longitude
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[92]] <- -90.5

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[92]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[92]] <- "Fixed longitude. Matches original study" 
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor == study_list[92]] <- "Grassland"

##### Narcissus_poeticus
study_list[99]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor
                         == study_list[99],] 

### fix longitude
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[99]] <- -90.5

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor
                          == study_list[99]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor
                             == study_list[99]] <- "Fixed longitude. Incomplete matrix in original paper, omit"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor == study_list[99]] <- "Grassland"

### Rosa_canina
study_list[140]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor
                         == study_list[140],] 

### fix longitude
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[140]] <- -90.5

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor
                          == study_list[140]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor
              == study_list[140]] <- "Fixed longitude. Incomplete matrix in original paper, omit"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[140]] <- "Grassland"

##### Rosa_multiflora
study_list[141]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor
                         == study_list[141],] 

### fix longitude
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[141]] <- -90.5


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor
                          == study_list[141]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor
                             == study_list[141]] <- "Fixed longitude. Matches original study" 
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[141]] <- "Forest"

###  Continue at 8 skipping  9  42  43  82  92  99 140 141

### "Alliaria_petiolata_4"
study_list[8]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor
                         == study_list[8],] 

### fix lat long to match original paper.. 

Bar <- rep(40.72205,3)
Ed <- rep(41.96451,3)
Farm <- rep(40.67668,3)
Healy <- rep(42.10862,3)
Holl <- rep(42.77765, 3)
Hom <- rep(40.06325,3)
Illi <- rep(40.07925,3)
Ives <- rep(41.98147,3)
John <- rep(42.92559,3)
Rose <- rep(42.81234,3)
Russ <- rep(42.01162,3)
Shia <- rep(42.88546, 2)

new_lat <- c(Bar, Ed, Farm, Healy, Holl, 
             Hom, Illi, 
             Ives, John, Rose, Russ, Shia)


sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor
                          == study_list[8]] <- new_lat

Bar <- rep(-89.5055,3)
Ed <- rep(-85.9962,3)
Farm <- rep(-89.4878,3)
Healy <- rep(-88.2148,3)
Holl <- rep(-86.2025, 3)
Hom <- rep(-87.9787,3)
Illi <- rep(-88.2107,3)
Ives <- rep(-83.932,3)
John <- rep(-85.7699,3)
Rose <- rep(-84.4042,3)
Russ <- rep(-85.9703,3)
Shia <- rep(-84.0491, 2)

new_lon <- c(Bar, Ed, Farm, Healy, Holl, Hom, Illi, 
             Ives, John, Rose, Russ, Shia)

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor
                       == study_list[8]] <- new_lon


#### omit MatrixPopulation "Russ Forest" &  "Rose Lake" lambda values of 
# A matrices 2005 quite different from those given in paper. 

om1 <- which(sp_to_check_in_try$MatrixPopulation %in% c("Russ Forest","Rose Lake"))

sp_to_check_in_try <- sp_to_check_in_try[-om1,]

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor
                          == study_list[8]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor
                             == study_list[8]] <- "coordinates fixed. Omitted Matrix population Russ Forest and Rose Lake. Can't get original matrices, these were sent by email to Compadre. Life-cycle matches the paper and lambda values match except in Russ Forest and Rose Lake 2005 " 

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[8]] <- "Forest"

####### "Anthericum_liliago"
study_list[10]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor
                         == study_list[10],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor
                          == study_list[10]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor
                             == study_list[10]] <- "Has clonality" 

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[10]] <- "open, sunny stony slopes"


###### "Anthericum_ramosum"
study_list[11]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor
                         == study_list[11],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor
                          == study_list[11]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor
                             == study_list[11]] <- "Has clonality" 

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[11]] <- "open, sunny stony slopes and forest"



###### "Anthoxanthum_odoratum"  ### study_list[12] already checked above. 
study_list[12]

#### "Araucaria_cunninghamii"

study_list[13]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor
                         == study_list[13],] 

### Ammended this matrix based on updated version supplied in: 
## Enright and Watson, 1991, Australian Journal of Ecology, 16, 507-520. 
matrixClass1[CA$UID1][[1]] <- matrixClass1[CA$UID1][[1]][2:11,1:3]
 
matA <- matrix(c(0.390,0,0,0,0,198,494,347,314,254,
                 0.012,0.880,0,0,0,0,0,0,0,0,
                 0,0.035,0.952,0,0,0,0,0,0,0,
                 0,0,0.008,0.952,0,0,0,0,0,0,
                 0,0,0,0.010,0.980,0,0,0,0,0,
                 0,0,0,0,0.017,0.935,0,0,0,0,
                 0,0,0,0,0,0.062,0.961,0,0,0,
                 0,0,0,0,0,0,0.036, 0.961,0,0,
                 0,0,0,0,0,0,0,0.029,0.979,0,
                 0,0,0,0,0,0,0,0,0,0.964), 
                 nrow = 10, ncol = 10, byrow = TRUE)
matF <- matrix(c(0,0,0,0,0,198,494,347,314,254,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0), 
               nrow = 10, ncol = 10, byrow = TRUE)

matU <- matA-matF

mat1[CA$UID1][[1]]$matA <- matA
mat1[CA$UID1][[1]]$matU <- matU
mat1[CA$UID1][[1]]$matF <- matF

### 
#######

### added lat long from google maps based on the centroid of the everglades as gps not given in paper.
sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor
                          == study_list[13]] <- -7.25
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor
                          == study_list[13]] <- 146.68

sp_to_check_in_try$MatrixStartYear[sp_to_check_in_try$SpeciesAuthor
                       == study_list[13]] <- 1977

sp_to_check_in_try$MatrixEndYear[sp_to_check_in_try$SpeciesAuthor
                                   == study_list[13]] <- 1982

sp_to_check_in_try$StudyDuration[sp_to_check_in_try$SpeciesAuthor
                                 == study_list[13]] <- 6

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor
                          == study_list[13]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor
          == study_list[13]] <- "edits to lat long and matrices based on 
                          ammended matrices in  Enright and Watson, 1991, Australian Journal of Ecology, 16, 507-520."
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[13]] <- "Forest"


#### "Ardisia_elliptica"

study_list[14]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor
                         == study_list[14],] 


### error in A and you matrices. Correct as follows.

mat1[CA$UID1][[8]]$matU[6,5] <- 0.224
mat1[CA$UID1][[8]]$matU[6,6] <- 0.95
mat1[CA$UID1][[8]]$matA <- mat1[CA$UID1][[8]]$matU +
                               mat1[CA$UID1][[8]]$matF

### added lat long from google maps based on the centroid of the everglades as gps not given in paper.
sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor
                          == study_list[14]] <- 25.5 ### error should have been 25.8
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor
                          == study_list[14]] <- -88.8  ### error should have been -80.822669

### coordinate errors are corrected in a later part of the analysis.

#######

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor
                          == study_list[14]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor
                             == study_list[14]] <- "edits to Transitional Forest 2000 based on original. Also do not matrixPopulation Schinus thicket 2000 as this is a replicate of 1999 and not in the original paper"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[14]] <- "Transitional forest"

### "Arenaria_serpyllifolia"
study_list[15]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[15],] 

### DOI should be : 10.1111/j.1654-1103.2007.tb02519.x
sp_to_check_in_try$DOI.ISBN[sp_to_check_in_try$SpeciesAuthor == study_list[15]] <- "10.1111/j.1654-1103.2007.tb02519.x"
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[15]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[15]] <- "DOI added, excel supplied by author so I can't see it, but transition types match paper"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[15]] <- "grassland and scree" 

#### which other species come from this study.. 

match_sp <- sp_to_check_in_try[which(sp_to_check_in_try$DOI.ISBN == "10.1658/1100-9233(2007)18[91:PDOAIP]2.0.CO;2"),]

which(study_list ==  "Veronica_arvensis" )
#177 

#### "Veronica_arvensis" 
study_list[177]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[177],] 

### DOI should be : 10.1111/j.1654-1103.2007.tb02519.x

sp_to_check_in_try$DOI.ISBN[sp_to_check_in_try$SpeciesAuthor == study_list[177]] <- "10.1111/j.1654-1103.2007.tb02519.x"
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[177]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[177]] <- "DOI added, excel supplied by author so I can't see it, but transition types match paper"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[177]] <- "grassland and scree" 

#### "Armeria_maritima" 
study_list[16]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[16],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[16]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[16]] <- "excel supplied by author so I can't see it, but transition types match paper"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[16]] <- "Mine spoil, recolonisation" 

#### "Artemisia_genipi"

study_list[17]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[17],] 

CA$SurvivalIssue
colSums(mat1[CA$UID1][[3]]$matU)

mat1[CA$UID1][[3]]$matU

### SurvivalIssue in matrix 3 is probably due to rounding, as there is no clonality/splitting in the methods
### and matrix matches original.  Ammend as follows.. 

mat1[CA$UID1][[3]]$matU[4,5] <- 0.3215
mat1[CA$UID1][[3]]$matU[5,5] <- 0.6785

mat1[CA$UID1][[3]]$matA[4,5] <- 0.3215
mat1[CA$UID1][[3]]$matA[5,5] <- 0.6785

### "Use, matches original"
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[17]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[17]] <- "SurvivalIssue in matrix 3 is probably due to rounding, ammended matrix column 5 to fix"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[17]] <- "alpine glacial moraines"

#####  Astrocaryum_mexicanum

study_list[18]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[18],] 

### fix longitude

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[18]] <- -95.12

### error (typo) in cell [1, 10] plot CC - should read 26.86 not 16.86 

mat1[CA$UID1][[7]]$matF[1,10] <- 26.86
mat1[CA$UID1][[7]]$matU
mat1[CA$UID1][[7]]$matA[1,10] <- 26.86

#### Do not use grand mean matrix.. "Plot A; Plot AA; Plot B; Plot BB; Plot C; Plot CC"
om1 <- which(sp_to_check_in_try$SpeciesAuthor == study_list[18] &
               sp_to_check_in_try$MatrixPopulation == "Plot A; Plot AA; Plot B; Plot BB; Plot C; Plot CC")


sp_to_check_in_try <- sp_to_check_in_try[-om1,]

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[18]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[18]] <- "typo corrected in plot CC, average of plot locations removed"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[18]] <- "forest"


# Atriplex_canescens

study_list[19]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[19],] 

### Lat long these don't match the paper 

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[19]] <- 26.65
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[19]] <- -103.93

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[19]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[19]] <- "Coordinates fixed. Matches original"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[19]] <- "grassland"

####

#  "Atriplex_vesicaria"

study_list[20]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[20],] 
CA_full <- metadata1[metadata1$SpeciesAuthor == study_list[20],] 

CA$Authors
CA$MatrixPopulation

### fix lat long based on original paper
sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[20]] <- -32.73
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[20]] <- 137.13

### typo in matrixPopulation A (Unmanipulated), [3,3] in matU and matA 
# should 0.42

mat1[CA$UID][[1]]$matA[3,3] <- 0.42
mat1[CA$UID][[1]]$matU[3,3] <- 0.42

### typo in matrixPopulation G (Unmanipulated), [1,1] in matU and matA 
# should 0.20

mat1[CA$UID][[7]]$matA[1,1] <- 0.2
mat1[CA$UID][[7]]$matU[1,1] <- 0.2

### matrix incorrectly split. #### 

# These matrices are difficult to split because the life cycle description included in the paper
# is for 6 months but the transition matrices are for one year.  

##  the useful paragraph is: 
# Thus, for example, the contributions of reproduction from the 
# 6-monthly matrices are spread out across the yearly matrix rather 
# than being confined to the first two rows in the last column. 
# For this reason, subadults, which were deemed to be non-reproductive,
# are shown with a small reproductive contribution over a 1-year 
# transition interval (the first and second rows of
# the third column of the matrix). Similarly, the entries
# on the first row of columns one and two of the yearly
# matrices include a component of reproduction, but this
# was a small contribution.

### furthermore the paper states that retrogression from adults to 
## juvenile stages is impossible.  i.e. a plant can't go for stages
## 3 or 4 to 1 or 2 except by reproduction. 

## Hence, cells [1:2,3:4] are the fecundity part of A 
## with the exception of some very small component which can not be
## separated in cells [1,1:2].  Because this reproduction is in this 
## column I have separated it into reproduction where it leads to 
## a growth transition of greater than one occuring from this 
## life-stage as this would be mathematical impossible. 


### I have resplit the matrices as follows.. 

matrixClass1[CA$UID1]
####

for(i in 1:nrow(CA)) {

mat1[CA$UID][[i]]$matF <- matrix(rep(0,16), nrow = 4)

mat1[CA$UID][[i]]$matF[1:2,3:4] <- mat1[CA$UID][[i]]$matA[1:2,3:4] 

mat1[CA$UID][[i]]$matU <- mat1[CA$UID][[i]]$matA - mat1[CA$UID][[i]]$matF

}

#####

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[20]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[20]] <- "fixed coorindates, matrix typo and corrected matrix splitting"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[20]] <- "arid open shrublands"

##### Avicennia_germinans

study_list[21]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[21],] 

### fix longitude from original ms
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[21]] <- -71.73

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[21]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[21]] <-
                                      "Longitude fixed"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[21]] <- "tropical rainforest (mangrove)"

###############################

##### Balsamorhiza_sagittata

study_list[22]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[22],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[22]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[22]] <-  "No fecundity given"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[22]] <- "Not using"


#### "Banksia_ericifolia"
study_list[23]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[23],] 

#### can't see original matrices but transition types look reasonable. 
### Matrices from NCEAS and Exeter to compadre.
#### All data are post fire, but fires are natural in the area (Australia)
## and this is a fire adapted tree species. 
## some data on seed and seedling survival < 3 years in from experimental fire rather 
## than wildfire. 


## Latitude longitude coordinates are about 100km from the location given in the paper. 
## Change to -33.48, 151.25 based on google maps.  Accuracy of this is within 10km. 

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[23]] <- -33.48
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[23]] <- 151.25


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[23]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[23]] <-  
   "Can't check original matrices but transitions types look reasonable.  Some seedling data is from experimental fires, but fire is natural in the ecosystem"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[23]] <- "Shrubland/forest"

######################################

#### "Bertholletia_excelsa"

study_list[24]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[24],] 

#### locations are misassigned to the matrices - swop 'El Tigre' and 'El Sena'
## for 'normal year'
cbind(CA$MatrixPopulation, as.character(CA$Observation))
#           [,1]       [,2]         
# [1,] "El Tigre" "Normal year"
# [2,] "El Sena"  "Normal year"
# [3,] "El Tigre" "Dry year"   
# [4,] "El Sena"  "Dry year" 

## 1 should be = "El Sena" "Normal year"
## 2 should be = "El Tigre" "Normal year"

ES <- which(sp_to_check_in_try$SpeciesAuthor == study_list[24] 
            & sp_to_check_in_try$MatrixPopulation == "El Tigre" 
            & sp_to_check_in_try$Observation == "Normal year")

ET <- which(sp_to_check_in_try$SpeciesAuthor == study_list[24] 
            & sp_to_check_in_try$MatrixPopulation == "El Sena" 
            & sp_to_check_in_try$Observation == "Normal year")


sp_to_check_in_try$MatrixPopulation[ES] <- "El Sena"
sp_to_check_in_try$MatrixPopulation[ET] <- "El Tigre"


##### lat long is wrong by by a couple of hundred km!  Fix to 

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[24] 
        & sp_to_check_in_try$MatrixPopulation == "El Tigre"] <- -10.98
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[24] 
      & sp_to_check_in_try$MatrixPopulation == "El Tigre"] <- -65.71

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[24] 
      & sp_to_check_in_try$MatrixPopulation == "El Sena"] <- -11.50
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[24] 
    & sp_to_check_in_try$MatrixPopulation == "El Sena"] <- -67.25

####

### typo in matrix 'El Tigre' 'Normal year'
# [13,13] <- 0.968 , [14,13] <- 0.027

mat1[CA$UID1][[2]]$matA[13,13] <- 0.968
mat1[CA$UID1][[2]]$matA[14,13] <- 0.027

mat1[CA$UID1][[2]]$matU[13,13] <- 0.968
mat1[CA$UID1][[2]]$matU[14,13] <- 0.027
mat1[CA$UID1][[2]]$matF

######


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[24]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[24]] <-  
  "matrices checked, typo and lat long corrected"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[24]] <- "Forest"

############### "Betula_nana"

study_list[25]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
          study_list[25],] 

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[25]] <- 68.62
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[25]] <- -149.28

### 

## Paper describes growth of 'individual branches' not individual plants.
## Can not be used to infer the life cycles of these plants, 
##  only of the branches.


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[25]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[25]] <-  
       "This matrix is for individual branches, not individual plants.  Therefore inappropriate for lifecycle calculations in our case"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[25]] <- "Tundra"

#### Bothriochloa_insculpta

study_list[26]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[26],] 

### Do not use, fragmentation of tussocks is included in matrices, meaning 
## that survival is greater than 1 where individuals multiply by splitting 
## as a result of damage. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[26]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[26]] <-  
  "Do not use, survivalIssue is present in original matrices"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[26]] <- "Savanna grasslands"


### "Bothriochloa_ischaemum"

study_list[27]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[27],] 

sp_to_check_in_try$MatrixComposite[sp_to_check_in_try$SpeciesAuthor == study_list[27]] <- 
             c("Mean", "Individual", "Individual")

### matrices are incorrectly split.  Re-enter U matrices from original source

##### Note here: where the values in the G matrix were greater than those given 
## in the A matrix (within 0.05 in all cases I've changed them to match the values 
## given in the A matrix).  For unclipped 1998-1999 this cells, [2:3,1] and [3:4, 4]

G1 <- matrix(c(0.3641618, 0.3764706, 0.1702128, 0, 0.05,
  0.223, 0.3058824, 0.3617021, 0.0833333, 0,
  0.032, 0.1529412, 0.2340426, 0.239, 0.05,
  0, 0.0235294, 0.1276596, 0.539, 0.35,
  0, 0,  0, 0.0833333, 0.5), nrow = 5, byrow = TRUE)

G1 <- round(G1, 3)

mat1[CA$UID1][[2]]$matU <- G1
mat1[CA$UID1][[2]]$matF <- mat1[CA$UID1][[2]]$matA - G1

mat1[CA$UID1][[2]]

G2 <- matrix(c(0.2444444, 0.2962963, 0.4, 0, 0.125,
              0.0555556, 0.1296296, 0.2, 0, 0.125,
               0, 0.0925926, 0.25, 0.2307692, 0.125,
                 0, 0, 0.05, 0.1538462, 0.125,
                  0, 0, 0, 0.1538462, 0.5), nrow = 5, byrow = TRUE)


mat1[CA$UID1][[3]]$matU <- G2
mat1[CA$UID1][[3]]$matF <- mat1[CA$UID1][[3]]$matA - G2

### Remove mean matrix from list

om1 <- which(sp_to_check_in_try$SpeciesAuthor == study_list[27] &
               sp_to_check_in_try$MatrixComposite == "Mean")

sp_to_check_in_try <- sp_to_check_in_try[-om1,]

###  ## Location added based on centroid of study site, 
## Pedernales Falls State Park, Texas, US

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[27]] <- 30.30
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[27]] <- -98.24

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[27]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[27]] <-  
  "Lat/Lon added, Original matrices were split incorrectly this has been recalculated"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[27]] <- "Grasslands"

###  "Brosimum_alicastrum"

study_list[28]

### this matrix has already been amended above, to be a 'female only' matrix.
### see above.  Checking further details here.

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[28],] 

# Longitude does not match original source. Amend to : -97.233333

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[28]] <- -97.23

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[28]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[28]] <-  
  "matrix editted so only female lifecycle is shown, longitude corrected"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[28]] <- "Tropical forest"


###### "Calamagrostis_canescens"

study_list[29]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[29],] 

### has clonality matrix - omit from study. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[29]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[29]] <-  
  "has clonality"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[29]] <- "NA"

#####"Callitris_intratropica"

study_list[30]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[30],] 

#### Do not use the paper uses seed and seedling data from multiple sources, 
### and fire effects are based on simulations. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[30]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[30]] <-  
  "Frankenstein matrix with simulated treatments"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[30]] <- "NA"

##### "Camellia_japonica"

study_list[31]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[31],] 

#### Fecundity values are not reported.  Reported A matrix is actually the U matrix. 
### I think annualPeriodicity should be 0.2


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[31]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[31]] <-  
  "No fecundity reported.  A matrix is actually the U matrix"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[31]] <- "NA"



#### "Camellia_japonica_2"

study_list[32]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[32],] 

mat1[CA$UID1]
matrixClass1[CA$UID1]
CA$DOI.ISBN
CA$AdditionalSource
CA$Authors

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[32]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[32]] <-  
  "No associated information about location, publication or authors"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[32]] <- "NA"

###  "Carduus_nutans"

study_list[33]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[33],] 

#### these appear correct

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[33]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[33]] <-  
  "Correct, information on matrix splitting was obtained by compadre from author"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[33]] <- "grasslands"

####   "Carduus_nutans_3"

study_list[34]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[34],] 

### fix latitude
sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[34]] <- c(-36.36, -40.25)

### all appear correct
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[34]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[34]] <-  
  "Latitude fixed"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[34]] <- "grasslands"


#### "Carex_bigelowii"

study_list[35]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[35],] 


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[35]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[35]] <-  
  "Do not use, U matrix includes clonal transitions I think.. "
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[35]] <- "NA"

###  "Carex_humilis"

study_list[36]

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[36]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[36]] <-  
  "Do not use, has clonality "
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[36]] <- "NA"

#####"Carlina_vulgaris_2"

study_list[37]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[37],] 


####  Can't find original matrices and it looks like some fecundity has been misassigned? 
### do not use.. 
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[37]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[37]] <-  
  "Can't find original matrices and it looks like some fecundity has been misassigned?"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[37]] <- "NA"

#####  "Castanea_dentata"
study_list[38]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[38],] 

### latitude and longitude do not match the original paper.  

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                     study_list[38] 
                   &sp_to_check_in_try$MatrixPopulation == "Leelanau"] <- -85.759
                                                                       
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[38]
            &sp_to_check_in_try$MatrixPopulation == "Missaukee healthy"] <- -85.15

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[38]
      & sp_to_check_in_try$MatrixPopulation == "Missaukee Diseased"] <- -85.43

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[38]
                       & sp_to_check_in_try$MatrixPopulation == "Missaukee Diseased"] <- 44.44


sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[38]
                & sp_to_check_in_try$MatrixPopulation == "Stivers"] <- 44.48


sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[38]
             & sp_to_check_in_try$MatrixPopulation == "Frankfort"] <- -86.23


sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[38]
                       & sp_to_check_in_try$MatrixPopulation == "County Line"] <- -86.10


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[38]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[38]] <-  
                  "Coordinates fixed. Can't access original matrices, but transition types look correct"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[38]] <- "Forest"

###### "Cecropia_obtusifolia"
study_list[39]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[39],] 

### correct the longitude value.. 

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[39]] <- -95.12


### error in matrixPopulation - "Los Tuxtlas" - 1989
mat1[CA$UID1][[3]]$matA[6,4] <- 0.0769

##########

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[39]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[39]] <-  
                       "error ammended, longitude corrected"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[39]] <- "Forest"

######### "Centaurea_jacea"
study_list[40]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[40],] 

##########

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[40]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[40]] <-  
                              "Has clonality"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[40]] <- "grasslands"

######  "Centaurea_stoebe_2"

study_list[41]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[41],] 

CA_full <- metadata1[metadata1$SpeciesAuthor == study_list[41],]

### Frankenstein matrix.  Do not use. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[41]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[41]] <-  
  "Frankenstein matrix"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[41]] <- "NA"


############ 42 and 43 are already done.. 
## "Cerastium_fontanum" and "Cerastium_pumilum"
### move onto 44.. 

#######  "Chaerophyllum_aureum"

study_list[44]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[44],] 

CA_full <- metadata1[metadata1$SpeciesAuthor == study_list[44],]


### model included density dependent seedling survival, I don't know 
## how this was incorporated into cell [1,2].  Compadre team must have had 
# extra information from somewhere


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[44]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[44]] <-  
  "Survival issue in last stage, also model has density dependent seedling survival, I don't know 
     how this was incorporated into cell [1,2]. Compadre team must have had 
              extra information from somewhere"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[44]] <- "grasslands"

######  Cirsium_acaule

study_list[45]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[45],] 

CA_full <- metadata1[metadata1$SpeciesAuthor == study_list[45],]

### increased spatial accuracy based on paper appendix. 
## Demography was conducted at 3 sites, this is the centroid of that area. 
# 50.52  14.22

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[45]] <- 50.52
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[45]] <- 14.22

##### 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[45]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[45]] <-  
              "edits to improve precision on lat long"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[45]] <- "grasslands"

##### [1] "Cirsium_acaule_2"

study_list[46]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[46],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[46]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[46]] <-  
  "This data is sourced from the same study as Cirsium_acaule.  I don't know why the fecundity differs. Omit."
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[46]] <- "grasslands"


##### Cirsium_arvense

study_list[47]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[47],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[47]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[47]] <-  
  "Frankenstein matrix"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[47]] <- "NA"


#### "Cirsium_dissectum"

study_list[48]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[48],] 

### has clonality

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[48]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[48]] <-  
  "Has clonality"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[48]] <- "NA"

####"Cirsium_palustre"

study_list[49]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[49],] 

### matrix class author for classes 4 and 5 are swopped. 
### change to 

matrixClass1[CA$UID1][[1]]$MatrixClassAuthor[4:5] <- c("≥ 2 years old vegetative","Early flowering")
matrixClass1[CA$UID1][[2]]$MatrixClassAuthor[4:5] <- c("≥ 2 years old vegetative","Early flowering")
matrixClass1[CA$UID1][[3]]$MatrixClassAuthor[4:5] <- c("≥ 2 years old vegetative","Early flowering")


### typo in 3rd matrix.. amend.. 

mat1[CA$UID1][[3]]$matA[6,4] <- 0.070 
mat1[CA$UID1][[3]]$matU[6,4] <- 0.070 
### 59.480727, 18.592338

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[49]] <- 59.480727
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[49]] <-  18.592338
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[49]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[49]] <-  
  "Typo fixed and coordinates refined"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[49]] <- "grasslands"

#####


####  "Cirsium_pannonicum"

study_list[50]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[50],] 

### Lat and long recalculated based on demography centroid in paper appendix.

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[50]] <- 50.54
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[50]] <-  14.16
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[50]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[50]] <-  
  "coordinates refined"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[50]] <- "grasslands"

####  Cirsium_undulatum

study_list[51]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[51],] 

CA_full <- metadata1[metadata1$SpeciesAuthor == 
                       study_list[51],] 

##################

#### fix lat long coordinates

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[51]] <- 38.13
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[51]] <- -99.05

#### error matrix 7 MatrixStartYear 1948  ammend as follows:

mat1[CA$UID1][[7]]$matA[1:2,1:2] <- c(0.3,0.17,0.3,0.17)
mat1[CA$UID1][[7]]$matA[1:2,1:2] 
mat1[CA$UID1][[7]]$matF[1:2,1:2] <- c(0.3,0,0.3,0)
mat1[CA$UID1][[7]]$matU[1:2,1:2] <- c(0,0.17,0,0.17)

###

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[51]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[51]] <-  
  "coordinates refined, typo in matrix corrected"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[51]] <- "grasslands"

####  

x1 <- sp_to_check_in_try[sp_to_check_in_try$DOI.ISBN == unique(CA$DOI.ISBN),]
sp_to_check <- unique(x1$SpeciesAuthor)

### do 134 and 157 whilst data is at hand.. 

#### Ratibida_columnifera

study_list[134]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[134],] 
CA$Lon

###### 
sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[134]] <- 38.13
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[134]] <- -99.05

###
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[134]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[134]] <-  
  "coordinates refined"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[134]] <- "grasslands"

#### "Sphaeralcea_coccinea"

study_list[157]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[157],] 
##################

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[157]] <- 38.13
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[157]] <- -99.05


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[157]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[157]] <-  
  "coordinates refined"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[157]] <- "grasslands"

####  Cirsium_vulgare_3

study_list[52]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[52],] 

CA_full <- metadata1[metadata1$SpeciesAuthor == 
                       study_list[52],] 

### lat incorrect ammend to -35.20

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[52]] <- -35.20

#### use only matrices 2 and 3 here, 1 and 4 are incomplete i.e. MatrixEndYear 1982 and 1983

om1 <- which(sp_to_check_in_try$SpeciesAuthor == study_list[52] &
        sp_to_check_in_try$MatrixEndYear %in% c("1981", "1984"))

sp_to_check_in_try <- sp_to_check_in_try[-om1,]

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[52]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
              == study_list[52]] <- "coordinates fixed, use only matrices MatrixEndYear %=% c(1982, 1983)"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[52]] <- "grasslands"

### "Clidemia_hirta"
study_list[53]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[53],] 

### only mean matrices are included here, individual annual matrices are given in original ms. 

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[53]] <- c(19.64, 19.99)
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[53]] <- c(-155.22, -155.25)

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[53]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[53]] <- "coordinates fixed"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[53]] <- "forest"

######## Colchicum_autumnale_2

study_list[54]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[54],] 

CA$Authors

CA_full <- metadata1[metadata1$SpeciesAuthor == 
                       study_list[54],] 


###  centroid for "Austria" = lat : 48.17, long : 15.71 

### coordinates for 'Germany' don't make sense. e.g. longitude 8°83'2''  ! 
### I've hard to just assume they are right in terms of degrees lat long and average those.. 
### this means precision is about 100km, but hopefully better than that as the locations make
### sense on google maps
### centroid for "Germany" = lat : 50.2, long : 8.5 ,  

#### error cell [6,3] in A and U matrices "matrixPopulation" = Austria, 
### MatrixStartYear = 2010.  [6,3] should be 0.02 not 0.03

mat1[CA$UID1][[3]]$matA[6,3] <- 0.02
mat1[CA$UID1][[3]]$matU[6,3] <- 0.02

### survival issues in matrices caused for mixing clonity in growth matrices. Do not use. 


sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[54] 
          & sp_to_check_in_try$MatrixPopulation == "Austria"] <- 48.17
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[54] 
                       & sp_to_check_in_try$MatrixPopulation == "Austria"] <- 15.71

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[54] 
                       & sp_to_check_in_try$MatrixPopulation == "Germany"] <- 50.2
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[54] 
                       & sp_to_check_in_try$MatrixPopulation == "Germany"] <- 8.5


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[54]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[54]] <- "survival issues in matrices caused for mixing clonity in growth matrices."
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[54]] <- "grasslands"


#### Cucurbita_pepo
study_list[55]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[55],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == 
                            study_list[55]] <- "N"

sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor ==
            study_list[55]] <- "Do not use, errors in assignment of F, also it's a common garden experiment"

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[55]] <- "NA"

#### "Cynoglossum_officinale"
study_list[56]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[56],] 

### I don't think this matrix is correctly transcribed as an annual matrix. 
### for example the seed to seedling transition at 25% happens over 6 months 
### not 12. First transitions need to be compressed such that there
### are only two stages - vegetative and adult  (seed germinates in the same year as deposition) 
### Elsewhere in the paper it states that adult to adult transition 
### also occurs rarely (i.e. ages 2-3 and 3-4)

### change to 
matA <- matrix(c(0,22.5,0.044,0.16), ncol = 2, byrow = TRUE)

matU <- matrix(c(0,0, 0.044,0.16), ncol = 2, byrow = TRUE)

matF <- matrix(c(0,22.5,
                 0,0), ncol = 2, byrow = TRUE)


mat1[CA$UID1][[1]]$matA <- matA
mat1[CA$UID1][[1]]$matU <- matU
mat1[CA$UID1][[1]]$matF <- matF

####  Ammend lat long

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[56]] <- 52.97
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[56]] <- 0.77
      
                
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[56]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
         == study_list[56]] <- "coordinates added, matrices fixed"

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[56]] <- "dune grasslands"

#####

#### "Cytisus_scoparius"
study_list[57]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[57],] 

### 

  #                        Lat      Lon
# 13th Division Prairie   47.10     -122.53
# Discovery Park          47.66     -122.41 
# Johnson Prairie         47.10     -122.53
# Magnuson Park           47.68     -122.25
# Montlake Fill           47.66     -122.29
# Weir Prairie            47.10     -122.53


sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[57]
    & sp_to_check_in_try$MatrixPopulation == "Discovery Park", 27 ] <- 47.66

sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[57]
                   & sp_to_check_in_try$MatrixPopulation == "Discovery Park", 28 ] <- -122.41 

Fort_lewis_base <- c("13th Division Prairie","Johnson Prairie", "Weir Prairie")


sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[57]
                   & sp_to_check_in_try$MatrixPopulation %in% Fort_lewis_base, 27 ] <- 47.10   

sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[57]
                   & sp_to_check_in_try$MatrixPopulation %in% Fort_lewis_base, 28 ] <- -122.53 

sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[57]
                   & sp_to_check_in_try$MatrixPopulation == 
                     "Magnuson Park", 27 ] <- 47.68

sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[57]
                   & sp_to_check_in_try$MatrixPopulation == 
                     "Magnuson Park", 28 ] <- -122.25 

sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[57]
                   & sp_to_check_in_try$MatrixPopulation == 
                     "Montlake Fill", 27 ] <- 47.66

sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == study_list[57]
                   & sp_to_check_in_try$MatrixPopulation == 
                     "Montlake Fill", 28 ] <- -122.29 

#################################

### these appear to be correct, originals are from Parker 2000
### http://www.jstor.org/stable/2641041


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[57]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[57]] <- "coordinates added"

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[57]] <- "urban grasslands and prairie"

####"Dicorynia_guianensis"  

study_list[58]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[58],] 

#### fix coordinates .. 

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[58]] <- 5.27

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                     study_list[58]] <- -52.92


### It's unclear where the matrices were transcribed from
## but they appear reasonable based on described
## lifecycle

####

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[58]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[58]] <- "coordinates added"

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[58]] <- "forest"

####  "Digitalis_purpurea_2"

study_list[59]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[59],] 

#### lat long
sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[59]] <- 60.57

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[59]] <- 6.92

####

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[59]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[59]] <- "coordinates edited to two decimal places"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[59]] <- "grasslands"

#### "Dipsacus_sylvestris_2"

study_list[60]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[60],] 

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[60]] <- -85.40

### error in Matrix Population 'B' cell [6,6] should be 3.67 in A and U matrices
mat1[CA$UID1][[2]]$matA[6,6] <- 0.367
mat1[CA$UID1][[2]]$matU[6,6] <- 0.367


####

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[60]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[60]] <- "Planted population, but unmanipulated thereafter. Longitude fixed"

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[60]] <- "grasslands"

###### Echium_vulgare
study_list[61]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[61],] 

#### lat long

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[61]] <- -43.00

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[61]] <- -76.08


## annual seed production 247
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[61]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[61]] <- "coordinates added, I'm not confident that this
                            matrix is correct, nor that I can correctly make it from the info in the paper"

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[61]] <- "abandoned limestone quarry. vegetation density < 50%"

####### Eperua_falcata

study_list[62]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[62],] 

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[62]] <- 5.30

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[62]] <- -52.92


## typo in matrix A and U, transition [1,1] should = 0.9814
mat1[CA$UID1][[1]]$matA[1,1] <- 0.9814
mat1[CA$UID1][[1]]$matU[1,1] <- 0.9814


## This model assumes equal fecundity accross all age classes, but this is 
## unrealistic and will underestimate age at maturity considerably.
## I have therefore added information from ...

# Effect of reduced-impact logging on seedling recruitment
# in a neotropical forest
# Skye L. Rivett, Jake E. Bicknell ⇑, Zoe G. Davies

##### And adjusted fecundity estimates stating that : 
### species begins reproduction at dbh 0.6-1m, i.e. stage class 6. 
## i.e. assuming all reproduction (6*0.0072) comes from this stage

mat1[CA$UID1][[1]]$matF[1,2:5] <- 0
mat1[CA$UID1][[1]]$matF[1,6] <- 6*0.0072
mat1[CA$UID1][[1]]$matA[1,2:5] <- 0
mat1[CA$UID1][[1]]$matA[1,6] <- 6*0.0072


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[62]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[62]] <- "coordinates fixed, changes made to fecundity estimates.
                              Note: This model assumes equal fecundity accross all age classes, but this is 
                               unrealistic and will underestimate age at maturity considerably.
                               I have therefore added information from Effect of reduced-impact logging on seedling recruitment
                               in a neotropical forest
                              Skye L. Rivett, Jake E. Bicknell ⇑, Zoe G. Davies
                               And adjusted fecundity estimates stating that : 
                               species begins reproduction at dbh 0.6-1m, i.e. stage class 6. 
                               i.e. assuming all reproduction (6*0.0072) comes from this stage"

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[62]] <- "forest"

### Check these two now as well "Oxandra_asbeckii"     "Vouacapoua_americana"
x1 <- sp_to_check_in_try$SpeciesAuthor[sp_to_check_in_try$DOI.ISBN == CA$DOI.ISBN] 

which(study_list %in% x1)

####  62 104 180

# [1] "Oxandra_asbeckii"
study_list[104]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[104],] 



sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[104]] <- 5.30
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[104]] <- -52.92


### DBH at maturity checked in 
##Local disturbance, forest structure and dispersal effects
# on sapling distribution of light-demanding
# and shade-tolerant species in a French Guianian forest
# Olivier Flores *, Sylvie Gourlet-Fleury, Nicolas Picard

### seed production from 10cm has been observed for this species
## at the same site as the transition matrix model. 



sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[104]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[104]] <- "coordinates fixed"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[104]] <- "forest"

##### "Vouacapoua_americana"

study_list[180]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[180],] 

##### lat long ammended
sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[180]] <- 5.30
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[180]] <- -52.92
#### 


####
## This model assumes equal fecundity accross all age classes, but this is 
## unrealistic and will underestimate age at maturity considerably.
## I have therefore added information from ...

# Spatial patterns of two rodent-dispersed rain
# forest trees Carapa procera (Meliaceae) and
# Vouacapoua americana (Caesalpiniaceae) at
# Paracou, French Guiana
# PIERRE-MICHEL FORGET*, FRANCËOIS MERCIER² and FREÂDEÂRIQUE
# COLLINET³

##### And adjusted fecundity estimates based on this assuming that
### species begins reproduction at dbh 0.26m, i.e. stage class 2. 
## i.e. reproduction is resplit as (0.0099*5)/4 in each stage 2:5 

matrixClass1[CA$UID1]
fec1 <- (0.0099*5)/4
fec1
mat1[CA$UID1][[1]]$matF[1,1] <- 0
mat1[CA$UID1][[1]]$matF[1,2:5] <- fec1

mat1[CA$UID1][[1]]$matA <- mat1[CA$UID1][[1]]$matF + mat1[CA$UID1][[1]]$matU

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[180]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[180]] <- "coordinates fixed, changes made to fecundity estimates.
Note: This model assumes equal fecundity accross all age classes, but this is 
unrealistic and will underestimate age at maturity considerably.
I have therefore added information from 
Spatial patterns of two rodent-dispersed rain
forest trees Carapa procera (Meliaceae) and
 Vouacapoua americana (Caesalpiniaceae) at
 Paracou, French Guiana
 PIERRE-MICHEL FORGET*, FRANCËOIS MERCIER² and FREÂDEÂRIQUE
 COLLINET³
And adjusted fecundity estimates given that species begins reproduction 
at dbh 0.26m, i.e. stage class 2"

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[180]] <- "forest"


##### Epipactis_atrorubens

study_list[63]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[63],] 

# lat long are correct!


####  I can't find oiginal matrices . Life stages look correct from diagram. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[63]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[63]] <- "I can't find original matrices . Life stages look correct from diagram"

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[63]] <- "forest"

### Eryngium_maritimum

study_list[64]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[64],] 

### 
### doi: 10.1111/j.1756-1051.2004.tb01647.x
### author name mispelled should be Curie. 


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[64]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
        == study_list[64]] <- "authors report high levels of clonality 
                      and matrices do not take into account seed bank or clonality.
                      Original matrix combines transitions and fecundity not
                      easily split"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[64]] <- NA

### Eupatorium_perfoliatum

study_list[65]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[65],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[65]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[65]] <- "has clonality in matrix"

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[65]] <- NA


##### Euterpe_edulis
study_list[66]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[66],] 

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[66]] <- -22.829

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[66]] <- -47.109


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[66]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[66]] <- "corrections to lat long"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[66]] <- "forest"

###### "Euterpe_precatoria"

study_list[67]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[67],] 

#### lat long fixed
sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[67]] <- -10.98
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == 
                         study_list[67]] <- -65.71




sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[67]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[67]] <- "corrections to lat long"

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[67]] <- "forest"

#### "Fritillaria_meleagris"

study_list[68]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[68],] 


### This is a transplant experiment, do not use. 
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[68]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[68]] <- "This is a transplant experiment, do not use."
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[68]] <- "NA"

######  "Fumana_procumbens"
study_list[69]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[69],] 
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[69]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[69]] <- "Use as is"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[69]] <- "grasslands - low ridges and hillocks in open habitats"

###### "Gentiana_pneumonanthe"

study_list[70]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[70],] 

add1 <- metadata1[metadata1$SpeciesAuthor == study_list[70] 
                  & metadata1$MatrixPopulation == "GEPN_DW3"
                  & metadata1$MatrixStartYear %in% c(1988,1989,1990),]

add1$MatrixComposite <- "Individual"
add1$checkY <- ""
add1$checknote <- ""
add1$habitatAuthor <- ""




sp_to_check_in_try <- rbind(sp_to_check_in_try, add1)

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[70]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[70]] <- "Use all individual year matrices of GEPN_DW3,
                                          these should have been marked as individual 
                                          rather than mean"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[70]] <- "Heath"

###### "Gentianella_campestris"

study_list[71]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[71],] 

####
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[71]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[71]] <- "Unmanipulated condition is a 
                                      managed grazing regime"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[71]] <- "Managed grasslands"
#### "Geranium_sylvaticum"

study_list[72]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[72],] 
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[72]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[72]] <- "Use as is"

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[72]] <- "grasslands"

#### Geum_reptans

study_list[73]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[73],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[73]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[73]] <- "Has clonality, exclude"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[73]] <- "NA"

##### Geum_rivale

study_list[74]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[74],] 
####
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[74]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[74]] <- "Survival issues are in the 
                              original matrix, could be because of clonality not
                              being adequately separated."

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[74]] <- "NA"


#### "Heracleum_mantegazzianum"

study_list[75]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[75],] 

CA$Authors
CA$DOI.ISBN
CA$MatrixPopulation
### Matrixpopulation given as NA.  Give this a name so that later code works correctly! 

sp_to_check_in_try$MatrixPopulation[sp_to_check_in_try$SpeciesAuthor == study_list[75]] <- "West Bohemia, Czech"


######  
sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[75]] <- 50.05
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[75]] <- 12.61


### fecundity is misassigned. 

mat1[CA$UID1][[1]]$matF[2,4] <- mat1[CA$UID1][[1]]$matA[2,4]
mat1[CA$UID1][[1]]$matU <-  mat1[CA$UID1][[1]]$matA - mat1[CA$UID1][[1]]$matF

#####
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[75]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                              == study_list[75]] <- "Fecundity missasigned as retrogression
                                          in cell [2,4], lat long added from previous work DOI:10.1111/j.1365-2664.2005.01092.x"
 sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                  == study_list[75]] <- "grassland"


##### "Heteropogon_contortus"

study_list[76]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[76],] 


### Do not use, fragmentation of tussocks is included in matrices, meaning 
## that survival is greater than 1 where individuals multiply by splitting 
## as a result of damage. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[76]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[76]] <- "Original matrices has SurvivalIssues"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[76]] <- "NA"


#### Hieracium_floribundum
study_list[77]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[77],] 

## correct lon to -80.233
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[77]] <- -80.233

### I'm going to use this is because it's already noted as corrected in the observation column


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[77]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[77]] <- "Longitude fixed"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[77]] <- "grassland"

### Himantoglossum_hircinum
study_list[78]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[78],] 

### lat long match the paper

#### Cell [2,5] in original paper is 0.8227 not 0.813 BUT this causes
### survival > 1 in this column.  I am presuming this was
### previously changed by Compadre staff with reason. No change made here. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[78]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[78]] <- "Use as is"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[78]] <- "grassland"

### ---- Hypochaeris_radicata ---- 
study_list[79]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[79],] 

### Do not use MatrixTreatment should read Annual Mowing 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[79]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[79]] <- "Do not use, MatrixTreatment should read Annual Mowing"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[79]] <- "grassland"

###  "Centaurea_jacea"      "Hypochaeris_radicata" "Succisa_pratensis"  

 x1 <- unique(sp_to_check_in_try$SpeciesAuthor[sp_to_check_in_try$DOI.ISBN == 
                    unique(CA$DOI.ISBN)])
 x1
 ## "Centaurea_jacea"      "Hypochaeris_radicata" "Succisa_pratensis" 
 
 which(study_list == "Succisa_pratensis")
###   [1] 160
 ### Don't use this species for the same reason.. 
 
 sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[160]] <- "N"
 sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                              == study_list[160]] <- "Do not use MatrixTreatment should read Annual Mowing"
 sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                  == study_list[160]] <- "grassland" 
 
 ### "Ipomopsis_aggregata"
 
 study_list[80]
 
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                            study_list[80],] 

# added lat long based on research laboratory location

 sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[80]] <- 38.95
 sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[80]] <- -106.98
 
 ####
 sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[80]] <- "Y_edits"
 sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                              == study_list[80]] <- "lat long added"
 sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                  == study_list[80]] <- "grassland" 
 
 #####  "Iriartea_deltoidea"
 
 study_list[81]
 
 CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                            study_list[81],] 
 
 ### Survival transitions appear to be estimated based on the size class distribution in a single 4 month field season.  
 ## Not annual measurements, do not use. 
 
 ### quote from paper
 # Since all of my data were collected during a four month field season( August-November1989),
 # my life history statistic are time-specific.
 # To estimates survival probabilities, I assumed the distribution of individuals across size classes
 # represents the behavior of a single cohort followed over its lifetime. For each size class,I used the ratio
 # of number of individuals in size class, over the number of individuals in size class, to represent the
 # survival probability for palms in that stage.
 
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[81]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                              == study_list[81]] <- "Survival transitions appear to be estimated based 
                                on the size class distribution in a single 4 month field season.  
                                Not annual transition measurements, do not use."
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                  == study_list[81]] <- "NA" 
 
 #### Iris_germanica
 study_list[82]
 #### this species is already done above..
 
 ### Isatis_tinctoria
 
 study_list[83]
 
 CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                            study_list[83],] 
 
 ### do not use, matrix periodicity is incorrect. Also, this is an experimentally
 ## established population. 
 
 sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[83]] <- "N"
 sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                              == study_list[83]] <- "matrix periodicity is incorrect. Also, this is an experimentally established population."
 sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                  == study_list[83]] <- "NA" 
 
 #####"Lactuca_virosa"

 study_list[84]
 
 CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                            study_list[84],] 
 
 ### I don't think this matrix is correctly transcribed as an annual matrix. 
 ### for example the seed to seedling transition at 25% happens over 6 months 
 ### not 12. First transitions need to be compressed such that there
 ### are only two stages - vegetative and adult  (seed germinates in the same year as deposition) 
 ### Elsewhere in the paper it states that adult to adult transition 
 ### also occurs rarely (i.e. ages 2-3 and 3-4)
 
 ### change to 
 matA <- matrix(c(0,8,0.15,0), ncol = 2, byrow = TRUE)
 
 matU <- matrix(c(0,0, 0.15,0), ncol = 2, byrow = TRUE)
 
 matF <- matrix(c(0,8,
                  0,0), ncol = 2, byrow = TRUE)
 
 
 mat1[CA$UID1][[1]]$matA <- matA
 mat1[CA$UID1][[1]]$matU <- matU
 mat1[CA$UID1][[1]]$matF <- matF
 
 ######### coordinates added
 
 sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[84]] <- 52.97
 sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[84]] <- 0.77
 
 ###
 
 sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[84]] <- "Y_edits"
 sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                              == study_list[84]] <- "coordinates added, matrices fixed (original entry was not an annual transition)"
 sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                  == study_list[84]] <- "dune grasslands"
 
#### Lantana_camara
 
study_list[85]
 
 CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                            study_list[85],] 
 
### Don't use, These aren't really 'unmanipulated' in the sense of being 
# natural populations. Farm is cattle grazed and forest is a plantation forest. 


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[85]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[85]] <- "Don't use, These aren't really 'unmanipulated' in the sense of being 
              natural populations. Farm is cattle grazed and forest is a plantation forest. "
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[85]] <- "NA"


####### "Lantana_camara_2"

study_list[86]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[86],] 

### Don't use, These aren't really 'unmanipulated' in the sense of being 
# natural populations. Farm is cattle grazed and forest is a plantation forest. 


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[86]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[86]] <- "Don't use, These aren't really 'unmanipulated' in the sense of being 
natural populations. Farm is cattle grazed and forest is a plantation forest. "

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[86]] <- "NA"

####### Lathyrus_vernus
study_list[87]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[87],] 

### set habitats for each study site in the paper. 

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[87] &
          sp_to_check_in_try$MatrixPopulation %in% c("T1", "T2", "T3", "T4", "T5", "T6", "N")] <- "grassland"

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[87] &
                                   sp_to_check_in_try$MatrixPopulation %in% c("G", "H", "K", "L")] <- "forest"

#####
### these are slightly different from those in the paper, but I think it's just because of rounding. 
## These ones have more decimal places than those in the paper. 


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[87]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[87]] <- "Habitats added"

#### "Lindera_umbellata"

study_list[88]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[88],] 

### has clonal growth matrix omit.. 
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[88]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[88]] <- "has clonal growth matrix omit"

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[88]] <- "NA"

#####  "Linnaea_borealis"

study_list[89]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[89],] 

### matrix describes dynamics of individual shoots, not whole plants. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[89]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
            == study_list[89]] <- "matrix describes dynamics of individual shoots, not whole plant"

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[89]] <- "NA"


##### "Linum_catharticum"
study_list[90]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[90],] 


### I can't find original matrices but these appear to match the 
# lifecycle described. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[90]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[90]] <- "Treatment should be annual mowing"

sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[90]] <- "grassland"

### Linum_tenuifolium

study_list[91]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[91],] 

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[91]] <- 50.48

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor 
                       == study_list[91]] <- 14.13

####

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[91]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[91]] <- "lat long added"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[91]] <- "grassland"

#####Lonicera_maackii

study_list[92]

### already done above.  

### Lupinus_arboreus
study_list[93]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[93],] 


# Lat long calculated based on study location California Bodega Marine Reserve
# Lat 38.32, -123.06


sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor 
                       == study_list[93]] <- 38.32

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor 
                       == study_list[93]] <- -123.06

###
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[93]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[93]] <- "lat long added.  Note: plants were planted
                            experimentally, and rodent predation manipulated.  Matrices
                            represent control conditions"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[93]] <- "grassland"

## Miconia_albicans
study_list[94]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[94],] 

### Lat long ammended 
sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor 
                       == study_list[94]] <- -15.93

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor 
                       == study_list[94]] <- -47.88

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[94]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[94]] <- "lat long ammended 
                                      Note: this is the unburnt control condition in a region that 
                                     experiences fire regimes"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[94]] <- "shrubland and forest"

### do other species in this paper now too 
sp_to_check_in_try$SpeciesAuthor[sp_to_check_in_try$DOI.ISBN == unique(CA$DOI.ISBN)]
## "Roupala_montana"

which(study_list == "Roupala_montana")
## [1] 142

## "Roupala_montana"
study_list[142]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[142],] 

### clonal do not use

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[142]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[142]] <- "clonal, do not use"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[142]] <- "shrubland and forest"

##  "Miconia_prasina"
study_list[95]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[95],] 


## Matrices are not correctly divided.  The author reports clonality which is not included here.
## Also, although no Fecundity is recorded within quadrats this plant is reproductive at this site, 
## so the all zero fecundity matrix whilst true is potential unrepresentative of the species.  
## Note this in observation column? 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[95]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[95]] <- "clonal, and matrix not correctly divided"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[95]] <- NA

#####################################

##  "Minuartia_obtusiloba"
study_list[96]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[96],] 

# 
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor 
                       == study_list[96]] <- -105.6

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[96]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[96]] <- "Longitude fixed"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[96]] <- "alpine fellfields"

#### Molinia_caerulea

study_list[97]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[97],] 

#### matrices appear to have been supplied by author to compadre

### Use only matrices from Observation "unburned", observation says that this was a 'natural' fire
#  but it is more likely to have been anthropogenic given the habitat and location. 

### remove 'sp_to_check_in_try$Observation == "Natural fire, thus control conditions"'

x1 <- which(sp_to_check_in_try$SpeciesAuthor == study_list[97]
                   & sp_to_check_in_try$Observation == "Natural fire, thus control conditions")


sp_to_check_in_try <- sp_to_check_in_try[-x1,]

#
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[97]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[97]] <- "Use only observation = Unburned"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[97]] <- "heathland"


#### Myosotis_ramosissima

study_list[98]


CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[98],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[98]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[98]] <- "Correct as is"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[98]] <- "grassland and scree"

#####

sp_in_paper <- unique(sp_to_check_in_try$SpeciesAuthor[sp_to_check_in_try$DOI.ISBN == unique(CA$DOI.ISBN)])

which(study_list %in% sp_in_paper)
## [1]  15  98 150 177

### "Saxifraga_tridactylites"
study_list[150]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[150],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[150]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[150]] <- "Correct as is"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[150]] <- "grassland and scree"



##### "Veronica_arvensis"
study_list[177]

### already done above.. 

###  "Narcissus_poeticus"
study_list[99]

### already done above

###  "Narcissus_pseudonarcissus"
study_list[100]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[100],] 

### Do not use, clonal reproduction common in popuations and 
## unclear how it's incorporated in the matrices (or why its omitted)

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[100]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[100]] <- "Do not use, clonal reproduction common in popuations and 
                  unclear how it's incorporated in the matrices (or why its omitted if it is)"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[100]] <- "NA"

### Nothofagus_fusca

study_list[101]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[101],] 

#### additional source does not mention this species I think the 
### correct additional source should be June and Ogden 1975

### MatrixPopulation: This should be 'Mount Colenso'  which is the study site, 
### Antioch Dunes NWR is in California. 


sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor 
                       == study_list[101]] <- -39.75

sp_to_check_in_try$MatrixPopulation[sp_to_check_in_try$SpeciesAuthor 
                       == study_list[101]] <- "Mount Colenso"


sp_to_check_in_try$AdditionalSource[sp_to_check_in_try$SpeciesAuthor 
                  == study_list[101]] <- "June and Ogden, 1975"

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[101]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[101]] <- "Edits to longitude, population name and 
                                   additional source"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[101]] <- "Forest"

#### Orchis_purpurea

study_list[102]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[102],] 

#### coordinates seem fine

### matrices are correctly inputted. There are likely to be problems with lifespan 
# calculated due to low observed mortality last three life stages over time span of study. 


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[102]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[102]] <- "matrices seem to be correctly inputted. There are likely to be problems with lifespan 
  calculated due to low observed mortality last three life stages over time span of study. "
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[102]] <- as.character(sp_to_check_in_try$Observation[sp_to_check_in_try$SpeciesAuthor 
                                                                                         == study_list[102]] )
####  "Oxalis_acetosella"

study_list[103]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[103],] 

### do not use, has clonality matrix
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[103]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[103]] <- "Do not use, clonal"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[103]] <- NA

#### Oxandra_asbeckii

study_list[104]

## already done above

## Paeonia_officinalis
study_list[105]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[105],] 

### do not use CA$MatrixPopulation "managed habitat"

om1 <- which(sp_to_check_in_try$SpeciesAuthor == study_list[105]
                 & sp_to_check_in_try$MatrixPopulation == "Managed habitat")
om1
sp_to_check_in_try <- sp_to_check_in_try[-om1,]


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[105]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[105]] <- "Exclude managed population"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[105]]  <- CA$MatrixPopulation[1:10]


### Pentaclethra_macroloba

study_list[106]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[106],] 
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[106]] <- -83.983


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[106]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[106]] <- "Longitude fixed, matches original matrices in phd thesis"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[106]] <- "Forest"

#### Petrophile_pulchella

study_list[107]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[107],] 

#### can't see original matrices but transition types look reasonable. 
### Matrices from NCEAS and Exeter to compadre.
#### All data are post fire, but fires are natural in the area (Australia)
## and this is a fire adapted tree species. 
## some data on seed and seedling survival < 3 years in from experimental fire rather 
## than wildfire. 


## Latitude longitude coordinates are about 100km from the location given in the paper. 
## Change to -33.48, 151.25 based on google maps.  Accuracy of this is within 10km. 

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[107]] <- -33.48
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[107]] <- 151.25


mat1[CA$UID1]

### matrices digitized from 'old matrices' so I can't check them, however, they seem
## reasonable based on lifecycle

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[107]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[107]] <-  
  "Latitude ammended. Can't check original matrices but transitions types look reasonable.  Some seedling data is from experimental fires, but fire is natural in the ecosystem"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[107]] <- "Shrubland/forest"

#### Picris_hieracioides

study_list[108]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[108],] 

## 
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[108]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[108]] <- "Matrix type listed in Compadre as 'density'. Not enough detail in paper to check matrix"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[108]] <- "abandoned limestone quarry. vegetation density < 50%"

#### Pimpinella_saxifraga

study_list[109]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[109],] 
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[109]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[109]] <- "Use as is"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[109]] <- "unmanaged roadverge"

#### Pinguicula_alpina

study_list[110]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[110],] 
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[110]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[110]] <- "Use as is. Can't access matrices, but lifecycle looks right"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[110]] <- "heath"

#### "Pinguicula_vulgaris"

study_list[111]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[111],] 
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[111]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[111]] <- "Use as is. Can't access matrices, but lifecycle looks right"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[111]] <- "heath"

#### "Pinus_nigra"

study_list[112]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[112],] 
## correct lat
sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[112]] <- -44.266

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[112]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[112]] <- "Latitude fixed, matrices from previous database"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[112]] <- CA$MatrixPopulation

####### "Pinus_palustris"

study_list[113]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[113],] 
## correct long
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[113]] <- -84.21

###
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[113]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[113]] <- "Longitude corrected. But do not use site is managed by annual burning"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[113]] <- "Forest"

#### "Pinus_strobus"
study_list[114]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[114],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[114]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[114]] <- "edits made to MatrixTreatment further up in the code
                             to include all sites as unmanipulated as per paper"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[114]] <- "Forest"


###   "Plantago_coronopus"

study_list[115]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[115],] 

## has clonality do not use. 
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[115]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[115]] <- "has clonality"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[115]] <- NA

#### "Plantago_media"

study_list[116]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[116],] 

### lat long based on locations given in Eriksson and Eriksson 1997
# 10.1111/j.1756-1051.1997.tb00344.x

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor 
                       == study_list[116]] <- 58.96

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor 
                       == study_list[116]] <- 17.58

## matrices correctly divided

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor
                          == study_list[116]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[116]] <- "Lat/long added, 
                                   note sites are grazed by cattle"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[116]] <- "grasslands"

### "Poa_alpina"

study_list[117]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[117],] 

### includes clonality, do not use 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[117]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[117]] <- "has clonality"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[117]] <- NA

#### "Podophyllum_peltatum"

study_list[118]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[118],] 
## This transition matrix relates to transitions between vegetative and sexual states 
## and not fecundity per se. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[118]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[118]] <- "This transition matrix relates to transitions between vegetative and sexual states 
## and not fecundity per se. Not appropriate to include"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[118]] <- NA

### "Potentilla_anserina"

study_list[119]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[119],] 

### includes clonality, do not use 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[119]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[119]] <- "has clonality"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[119]] <- NA

####  "Primula_elatior"

study_list[120]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[120],] 
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[120]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[120]] <- "use as is"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[120]] <- "forest"

#### "Primula_farinosa"
study_list[121]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[121],] 

## locations inferred from population location names in table1 and google maps and 
# https://mapcarta.com/

# Rossholm - 60.58, 17.88
# Karbomosse - 58.03, 13.49
# Flottskär - 59.93, 18.87

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[121]
          & sp_to_check_in_try$MatrixPopulation == "Karbomosse"] <- 58.03

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[121]
       & sp_to_check_in_try$MatrixPopulation == "Karbomosse"] <- 13.49

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[121]
        & sp_to_check_in_try$MatrixPopulation == "Rossholm"] <- 60.58

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[121]
        & sp_to_check_in_try$MatrixPopulation == "Rossholm"] <- 17.88


sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[121]
         & sp_to_check_in_try$MatrixPopulation == "Flottskär"] <- 59.93

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[121]
         & sp_to_check_in_try$MatrixPopulation == "Flottskär"] <- 18.87

### I can't access original matrices but transition types look correct
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[121]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[121]] <- "lat long added"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[121]] <- "Grassland"

#####

#### "Primula_veris_2"
study_list[122]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[122],] 


### Meadow populations were mown for conservation purposes exclude. 
### Pasture populations were grazed, since we have ungrazed populations
## exclude these also. 


### survival issues in matrices are probably caused by the fact that clonal reproduction was
## incorporated in the U matrix and can not be split post-hoc. 

### From text : "Rare events of fission (vegetative reproduction) were handled by enabling 
# two simultaneous transitions from stage i to stage j by the same plant."


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[122]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[122]] <- "clonality included in U matrix, causing survival issues"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[122]] <- NA

##### "Primula_veris_4"

study_list[123]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[123],] 

# At ‘Altenbroek’ will be further referred to as ‘cleared canopy site’ Picea abies trees
# planted in 1990 on formerly intact calcareous grassland
#  were felled and removed from the site in early
# spring of 1999 in an attempt to restore the original
# vegetation. While clear-felling and removing the
# trees, only little soil and vegetation disturbance was
# caused. On the ‘Vrouwenbos’ site ‘closed canopy
# site’, Fraxinus excelsior was planted in 1980 again
# on former grassland and canopy has been closing
# during the last few years.

### whilst both of these sites are categorised as 'unmanipulated'
### one is a plantation forest and the other has been recently clearfelled. 
## Do not use. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[123]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[123]] <- "whilst both of these sites are categorised as 'unmanipulated'
                               one is a plantation forest and the other has been recently clearfelled"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[123]] <- NA

### "Primula_vulgaris"
study_list[124]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[124],] 

### DOI missing should be: 10.1046/j.1365-2745.1998.00280.x
sp_to_check_in_try$DOI.ISBN[sp_to_check_in_try$SpeciesAuthor == study_list[124]] <- "10.1046/j.1365-2745.1998.00280.x"

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[124]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[124]] <- "DOI added"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[124]] <- "forest"


###  "Primula_vulgaris_3"

study_list[125]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[125],] 
## fix longitude 
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[125]] <- -5.5

### Matrices recieved by compadre from author, I can't see the originals but 
## life-cycle looks accurate. 


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[125]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[125]] <- "Longitude corrected"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[125]] <- "fragmented wood-pasture habitat"


#### "Prosopis_glandulosa"

study_list[126]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[126],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[126]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[126]] <- "Survival issues and can not find original matrices"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[126]] <- NA

#### "Prunus_serotina"

study_list[127]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[127],] 


## Matrices constructed using space-for-time substitution.
## I.e. no individuals were repeatedly measured over time and instead 
## transitions were constructed based on ratios on individuals in each stage. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[127]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[127]] <- "matrices constructed using space-for-time substitution.
                             I.e. no individuals were repeatedly measured over time"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[127]] <- NA

#### "Psidium_guajava"

study_list[128]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[128],] 


### This was conducted on a 'farm' where management included - "Pastures are
# extensively managed, no fertilizer is applied, and weeds (including guava
# seedlings) are controlled by manual cutting and with a yearly application of
# 2,4-D at the end of the dry season. ... Mean annual stocking rate is 1 AU/
#  ha/yr (AU = animal unit of 350 kg live weight). Cattle are freely allowed to
# graze under the guavas during periods of fruit production."


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[128]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[128]] <- "This is not unmanipulated, 
                            study is carried out on actively managed farm"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[128]] <- NA

### "Ranunculus_acris"


study_list[129]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[129],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[129]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[129]] <- "has clonality"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[129]] <- NA

### "Ranunculus_bulbosus"
study_list[130]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[130],] 

### not really unmanipulated.. 
# "Sheep grazing takes place chiefly from late autumn to the
# spring, the field is normally grazed by cattle between the sheep-grazing periods and a few
# horses are regularly kept in the field during the late summer. Thistle cutting and manuring
# were done regularly between August and September every year, except in 1970"


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[130]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[130]] <- "This is not unmanipulated, 
study is carried out on actively managed farm"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[130]] <- NA

### "Ranunculus_bulbosus_2"
study_list[131]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[131],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[131]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[131]] <- "Note: this is based on the 
                            same study as Ranunculus_bulbosus. This is not unmanipulated, 
                                study is carried out on actively managed farm"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[131]] <- NA

### "Ranunculus_peltatus"
study_list[132]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[132],] 

## lat long corrected based on details of study location given in text. 
##58.802014, 17.690830

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[132]] <- 58.80
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[132]] <- 17.69

####
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[132]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[132]] <- "lat long corrected based on details of study location given in text"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[132]] <- "Marine"

###  "Ranunculus_repens"

study_list[133]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[133],] 

### not really unmanipulated.. 
# "Sheep grazing takes place chiefly from late autumn to the
# spring, the field is normally grazed by cattle between the sheep-grazing periods and a few
# horses are regularly kept in the field during the late summer. Thistle cutting and manuring
# were done regularly between August and September every year, except in 1970"


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[133]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[133]] <- "This is not unmanipulated, 
study is carried out on actively managed farm"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[133]] <- NA

####### "Ratibida_columnifera"

study_list[134]

### already done above 


#### "Rhizophora_mangle"
study_list[135]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[135],] 
### Ammended Longitude paper on original paper, but this is strangely
## documented in paper.  i.e. it reads 71 degrees 73 minutes W, even 
## though minutes should only go up to 60! 
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[135]] <- -71.74

###
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[135]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[135]] <- "Longitude fixed"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[135]] <- "tropical forest (mangrove)"

#### Rhododendron_maximum

study_list[136]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[136],] 

### This is model based on 'shoots' rather than whole plants.
## i.e. age, death, reproduction etc, refers to individual branches. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[136]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[136]] <- "This is model based on 'shoots' rather than whole plants."
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[136]] <- NA

##### "Rhododendron_maximum_2"

study_list[137]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[137],] 

### This is model based on 'shoots' rather than whole plants.
## i.e. age, death, reproduction etc, refers to individual branches. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[137]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[137]] <- "This is model based on 'shoots' rather than whole plants."
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[137]] <- NA


#### "Rhododendron_maximum_3"

study_list[138]


CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[138],] 

### This is model based on 'shoots' rather than whole plants.
## i.e. age, death, reproduction etc, refers to individual branches. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[138]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[138]] <- "This is model based on 'shoots' rather than whole plants."
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[138]] <- NA



#### "Rhododendron_ponticum"

study_list[139]


CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[139],] 

### survival issues, can not access original files, these are from an MSC thesis
## not on Compadre share drive. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[139]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[139]] <- "Survival issues, original not available"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[139]] <- NA

#### "Rosa_canina"

study_list[140]

### already done above .. move on.. 

#### "Rosa_multiflora"

study_list[141]

### already done above .. move on.. 

#### "Roupala_montana"

study_list[142]

### already done above .. move on.. 

### "Rubus_ursinus"
study_list[143]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[143],] 
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[143]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[143]] <- "Clonality do not use"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[143]] <- NA

### "Salix_arctica"

study_list[144]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[144],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[144]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[144]] <- "Clonality do not use"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[144]] <- NA

##### "Salix_arctica_2"

study_list[145]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[145],] 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[145]] <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[145]] <- "Clonality do not use"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[145]] <- NA

#### "Sanicula_europaea"

study_list[146]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[146],] 

sp_to_check_in_try$Observation <- as.character(sp_to_check_in_try$Observation)


sp_to_check_in_try$Observation[sp_to_check_in_try$SpeciesAuthor == 
                    study_list[146]] <- "Some cattle grazing in study region, 
                     study is conducted in decidous forest I can't tell from the
                     descrpition if the cows are in the forest!"

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[146]] <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[146]] <- "Use as is"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[146]] <- "forest"

##### ---- Sapium_sebiferum ----

study_list[147]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[147],] 

### do not use "Clearcut mixed pine-hardwood forest"
## this is not Unmanipulated

om1 <- which(sp_to_check_in_try$SpeciesAuthor == study_list[147] &
        sp_to_check_in_try$MatrixPopulation == "Clearcut mixed pine-hardwood forest")

sp_to_check_in_try <- sp_to_check_in_try[-om1,]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[147],] 
CA$SurvivalIssue

### Small Survival Issues 0.001 probably induced by rounding. 
## present in original matrices

colSums(mat1[CA$UID1][[1]]$matU)
mat1[CA$UID1][[1]]$matU

### error in matrix A and U cell [4,3]  should be 0.778
## matrixstartyear = 1997, MatrixPopulation = Maritime evergreen forest
mat1[CA$UID1][[1]]$matU[4,3] <- 0.778
mat1[CA$UID1][[1]]$matA[4,3] <- 0.778

### fix to survivalIssue probable rounding error. 
mat1[CA$UID1][[4]]$matU[7:8,7] <- 0.44825
mat1[CA$UID1][[4]]$matA[7:8,7] <- 0.44825
### 

CA$SurvivalIssue
# 1.0010 1.0000 1.0000 1.0001 1.0000
# 1.0001 1.0000 1.0001 1.0001 1.0000 
# 1.0000 1.0001


### fix to survivalIssue probable rounding error. 

mat1[CA$UID1][[6]]$matU


mat1[CA$UID1][[6]]$matU[c(5,8),7] <- 0.05555
mat1[CA$UID1][[6]]$matA[c(5,8),7] <- 0.05555

#####

colSums(mat1[CA$UID1][[8]]$matU)


mat1[CA$UID1][[8]]$matU[9:10,9] <- 0.42855
mat1[CA$UID1][[8]]$matA[9:10,9] <- 0.42855

####

colSums(mat1[CA$UID1][[9]]$matU)
mat1[CA$UID1][[9]]$matU

mat1[CA$UID1][[9]]$matU[c(6,8),7] <- 0.14285
mat1[CA$UID1][[9]]$matA[c(6,8),7] <- 0.14285

### 


colSums(mat1[CA$UID1][[12]]$matU)
mat1[CA$UID1][[12]]$matU

mat1[CA$UID1][[12]]$matU[3:6,5] <- mat1[CA$UID1][[12]]$matU[3:6,5] - 0.000025
mat1[CA$UID1][[12]]$matA[3:6,5] <- mat1[CA$UID1][[12]]$matA[3:6,5] - 0.000025


####
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[147]] <- -79.25

###

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[147]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[147]] <- "minor edits made to matrices to correct survivalIssue of 0.0001, longitude corrected"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[147]] <- "forest"

##### "Sarracenia_purpurea_2" 

study_list[148]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[148],] 

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[148]] <- -73.86

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[148]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[148]] <- "longitude corrected"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[148]] <- "bog/fen wetland"

#### Saxifraga_aizoides

study_list[149]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[149],] 
### survival issue presumably rounding error fix by 

mat1[CA$UID1][[2]]$matU[3:5,4] <- mat1[CA$UID1][[2]]$matU[3:5,4] - 0.00034
mat1[CA$UID1][[2]]$matA[3:5,4] <- mat1[CA$UID1][[2]]$matA[3:5,4] - 0.00034

CA$SurvivalIssue[3]
colSums(mat1[CA$UID1][[3]]$matU)

colSums(mat1[CA$UID1][[3]]$matU)

### survival issue presumably rounding error fix by 

mat1[CA$UID1][[3]]$matU[5:6,6] <- mat1[CA$UID1][[3]]$matU[5:6,6] - 0.0005
mat1[CA$UID1][[3]]$matA[5:6,6] <- mat1[CA$UID1][[3]]$matA[5:6,6] - 0.0005

## 
CA$SurvivalIssue[4]
colSums(mat1[CA$UID1][[4]]$matU)

### previous presumably a rounding error.  

mat1[CA$UID1][[4]]$matU[4:6,6] <- mat1[CA$UID1][[4]]$matU[4:6,6] - 0.00034
mat1[CA$UID1][[4]]$matA[4:6,6] <- mat1[CA$UID1][[4]]$matA[4:6,6] - 0.00034

### following these edits lambda values for A matrix still match those given in paper. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[149]] <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor 
                             == study_list[149]] <- "tiny correction (=< 0.0005) to U and A matrices to correct survival issues"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[149]] <- "alpine moraine"

####  "Saxifraga_tridactylites"
study_list[150]

### done above

### "Scorzonera_hispanica"
study_list[151]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[151],] 

### increase accuracy of lat long from paper. 

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[151]] <- 14.19
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[151]] <- 50.52


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[151]]   <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[151]] <- "accuracy of lat long improved from paper"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[151]] <- "grassland"

#### Scorzonera_humilis

study_list[152]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[152],] 

### not unmanipulated - MatrixTreatment should be 'annual mowing' 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[152]]   <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[152]] <- "not unmanipulated - MatrixTreatment should be annual mowing"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[152]] <- NA

### "Senecio_jacobaea_2"

study_list[153]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[153],] 

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[153]]   <- c(-114.95, -114.97)
sp_to_check_in_try$Observation[sp_to_check_in_try$SpeciesAuthor == study_list[153]] [2] <- "Nutrient poor, fire and salvage logging 11 years pre survey"


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[153]]   <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[153]] <- "Longitude corrected, spreading invasive populations. Can not see how matrices were constructed but life-cycle appears correct"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[153]] <- "Grasslands with some successional pine trees"

##### "Sequoia_sempervirens_2"

study_list[154]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[154],] 

### time steps are 50years AnnualPeriodicity should be 1/50 i.e 0.02

###
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[154]]   <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[154]] <- "annual peroidicity should be 1/50. matrices based on assumptions and inferences about parameters rather than sufficient data"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[154]] <- NA

##### "Shorea_leprosula"

study_list[155]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[155],] 

## add values for start and end years

sp_to_check_in_try$MatrixStartYear[sp_to_check_in_try$SpeciesAuthor == study_list[155]]   <- c(1986, 1990,1995)
sp_to_check_in_try$MatrixEndYear[sp_to_check_in_try$SpeciesAuthor == study_list[155]] <-  c(1990, 1995, 2000)

mat1[CA$UID1][[1]]$matA

### stasis values in matrices do not match originals for 2nd and 3rd matrices

mat1[CA$UID1][[2]]$matA[2,2] <- 0.8340
mat1[CA$UID1][[2]]$matA[3,3] <- 0.8316
mat1[CA$UID1][[2]]$matA[4,4] <- 0.9187
mat1[CA$UID1][[2]]$matA[5,5] <- 0.9336
mat1[CA$UID1][[2]]$matA[6,6] <- 0.9710
mat1[CA$UID1][[2]]$matA[7,7] <- 0.9869

mat1[CA$UID1][[2]]$matU <- mat1[CA$UID1][[2]]$matA - mat1[CA$UID1][[2]]$matF

###
mat1[CA$UID1][[3]]$matA

mat1[CA$UID1][[3]]$matA[2,2] <- 0.8016
mat1[CA$UID1][[3]]$matA[3,3] <- 0.8319
mat1[CA$UID1][[3]]$matA[4,4] <- 0.9138
mat1[CA$UID1][[3]]$matA[5,5] <- 0.9289
mat1[CA$UID1][[3]]$matA[6,6] <- 0.9646
mat1[CA$UID1][[3]]$matA[7,7] <- 0.9790

mat1[CA$UID1][[3]]$matU <- mat1[CA$UID1][[3]]$matA - mat1[CA$UID1][[3]]$matF


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[155]]   <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[155]] <- "Start and end year of matrices added for 2 and 3, annual periodicity corrected, matrices fixed based on original"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[155]] <- "forest"


#### "Silene_acaulis"

study_list[156]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[156],] 

### Survival Issues are most likely caused by round errors, also, original matrices show
## fecundity in earlier stages but this not shown in compadre because it would 
## have rounded to 0 at 5 decimal places. Matrices to been reinputted from Ellis et al. 2012

## Matrix population "Gulch" to be excluded because this is not included in the Ellis et al. 2012
## dataset


om1 <- which(sp_to_check_in_try$SpeciesAuthor == study_list[156] &
               sp_to_check_in_try$MatrixPopulation == "Gulch")

sp_to_check_in_try <- sp_to_check_in_try[-om1,]

#####

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[156],] 


matr1 <- read.csv( "silence_matrices_ellis_2012.csv"  , sep = ",")

### convert matA and matU to character strings
matr1$matA <- as.character(matr1$matA)
matr1$matU <- as.character(matr1$matU)


### remove ; and [ ]  from matA and matU
matr1$matA<- gsub(";", "" , matr1$matA)
matr1$matA<- gsub('^.|.$', '', matr1$matA)

matr1$matU<- gsub(";", "" , matr1$matU)
matr1$matU<- gsub('^.|.$', '', matr1$matU)

### check order of rows matchs in CA and matr1

nrow(matr1)
matr1$site_name
matr1$year
CA$MatrixPopulation
CA$MatrixStartYear

### replace matrices in mat1 with matrices from ellis et al. 2012 
nrow(matr1)
for(i in 1:nrow(CA)) {
x <- as.numeric(unlist(strsplit(matr1$matA[i], c(" "))))
x2 <- matrix(x, nrow = 12, ncol = 12, byrow = T)
mat1[CA$UID1][[i]]$matA <- x2
print(i)
}

#same for matU and get mat F by subtracting U from A 
for(i in 1:nrow(CA)) {
  x <- as.numeric(unlist(strsplit(matr1$matU[i], c(" "))))
  U2 <- matrix(x, nrow = 12, ncol = 12, byrow = T)
  mat1[CA$UID1][[i]]$matU <- U2
  mat1[CA$UID1][[i]]$matF <- mat1[CA$UID1][[i]]$matA - mat1[CA$UID1][[i]]$matU
  
  print(i)
}


##
#for(i in 1:20) {
#print(max(colSums(mat1[CA$UID1][[i]]$matU)))
#}
#### survival issues are solved too : ) 

### lat and long adjusted to match ellis et al. 2012

CA$Lon

CC <- rep(-142.84068,5)
PA <- rep(-142.8215,5)
RG <- rep(-142.84708, 5)
RI <- rep(-142.81453, 5)

Lon1 <- c(CC,PA,RG,RI)


str(sp_to_check_in_try)

sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[156]] <- 61.49
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[156]] <- Lon1

### add habitats

CC <- rep("scree",5)
PA <- rep("tundra meadow",5)
RG <- rep("boulder field",5)
RI <- rep("alpine fell field", 5)

Hab1 <- c(CC,PA,RG,RI)


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[156]]   <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[156]] <- "Lat long improved, matrices reloaded from originals to escape rounding issues, particularly with start of fecundity"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[156]] <- Hab1

#### "Sphaeralcea_coccinea"

study_list[157]

## Already done above.  move on. 

#### "Sporobolus_heterolepis"

study_list[158]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[158],] 


### has clonality omit.. 
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[158]]   <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[158]] <- "Has clonality omit"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[158]] <- NA

###### "Stryphnodendron_excelsum"

study_list[159]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[159],] 

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[159]]   <- -83.983

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[159]]   <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[159]] <- "Longitude corrected"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[159]] <- "forest"


#### "Succisa_pratensis"

study_list[160]

### already done above.. 
#### "Succisa_pratensis_3" 
study_list[161]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[161],] 

CA_full <- metadata1[metadata1$SpeciesAuthor == 
                       study_list[161],]

### I don't have the original matrices but the other details of the study are correct. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[161]]   <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[161]] <- "Use as is"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[161]] <- "grasslands"


#### "Swietenia_macrophylla"

study_list[162]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[162],] 

### lat long corrected from original paper. 
sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[162]]   <- -15.78
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[162]]   <- -62.92

### matrix correct
sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[162]]   <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[162]] <- "Coordinates fixed"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[162]] <- "forest"

##### "Syzygium_jambos"

study_list[163]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[163],] 

### longitude corrected from paper. 
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[163]]   <- -65.83

mat1[CA$UID1][[1]]$matA[1,4:6] <- c(140, 1166, 1519)
mat1[CA$UID1][[1]]$matF[1,4:6] <- c(140, 1166, 1519)

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[163]]   <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[163]] <- "Coordinates ammended, fecundity fixed based on original paper"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[163]] <- "forest"

#####  "Themeda_triandra"

study_list[164]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[164],] 

### Do not use, fragmentation of tussocks is included in matrices, meaning 
## that survival is greater than 1 where individuals multiply by splitting 
## as a result of damage.

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[164]]   <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[164]] <- "Survival issues caused by inclusion of tussock splitting as new individuals in transition matrix"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[164]] <- "grassland (light grazing)"

###  "Themeda_triandra_2"

study_list[165]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[165],] 

### this is the same study as "Themeda_triandra" do not use for the same
## reasons..

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[165]]   <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[165]] <- "Survival issues caused by inclusion of tussock splitting as new individuals in transition matrix"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[165]] <- "grassland (light grazing)"

#### "Thymus_webbianus"

study_list[166]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[166],] 

## I will need to use the in MatrixPopulation Calpe as no death was observed in the final stage in Benidorm. 

om1 <- which(sp_to_check_in_try$SpeciesAuthor == study_list[166]
             & sp_to_check_in_try$MatrixPopulation == "Benidorm")

sp_to_check_in_try <- sp_to_check_in_try[-om1,]

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[166]]   <- "Y"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[166]] <- "Only use population Calpe, as in Benidorm no death was observed"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[166]] <- "rocky escarpment"

##### Tragopogon_pratensis

study_list[167]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[167],] 

### matrices appear to have elements missing?  And I don't have
### the original versions, spreadsheet says matrices are from Ramula.
### A matrices are not ergodic.  


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[167]]   <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[167]] <- "matrices appear to have elements missing?"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[167]] <- NA

##### "Trifolium_pratense"
study_list[168]

### already done above .. ignore.. 

### "Trillium_grandiflorum"
study_list[169]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[169],] 


### fix longitude 

DC <- rep(-80.43, 4)
DH <- rep(-80.03, 3)
DR <- rep(-80.23, 3)
ER <- rep(-80.08, 3)
FX <- rep(-80.2, 3)
GR <- rep(-80, 3)
LR <- rep(-80.05, 3)
RH <- rep(-80.4, 4)
RM <- rep(-79.97, 3)
TW <- rep(-80.35, 4)
WC <- rep(-80.07, 3)
WH <- rep(-79.95, 3)

Lon2 <- c(DC, DH, DR, ER, FX, GR, LR, RH, RM, TW, WC, WH)

sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[169]]   <- Lon2

###
CA$SurvivalIssue

matrixClass1[CA$UID1][[1]]
colSums(mat1[CA$UID1][[1]]$matU)

### matrices appear to be correctly inputted


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[169]]   <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[169]] <- "Longitude fixed, high survival in late stages my cause problems in lifespan calculations"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[169]] <- "forest"

#### "Trillium_grandiflorum_3"

study_list[170]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[170],] 

#### add these from original
sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[170]]   <- 45.58
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[170]]   <- 74.0


### note there's an NA in one matrix. 


#### I don't have original matrices but life-cycle looks to be correctly
## represented in matrices. 

### check matrix 5 which has NA's - MatrixPopulation "Forest C", MatrixStartYear "2003"

CA$MatrixPopulation[5]
CA$MatrixStartYear[5]

### replace NA's first row with 0's.  This transition will always be 0 
## according to the life - cycle diagram

mat1[CA$UID1][[5]]$matU[1,] <- 0
mat1[CA$UID1][[5]]$matU 

CA$MatrixPopulation
CA$SurvivalIssue
CA$Authors

matrixClass1[CA$UID1][[1]]
mat1[CA$UID1][[2]]

colSums(mat1[CA$UID1][[2]]$matU)

hab2 <- c(rep("forest", 6), rep("hedgerow", 5))


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[170]]   <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[170]] <- "One matrix has survivalIssue arising from stage 5 which includes dormancy.  I don't know how this was calculated. Coordinates added from original and matrix error fixed."
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[170]] <- hab2

### "Trillium_ovatum"

study_list[171]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[171],] 

### Longitude incorrect..
### Multiple errors in matrix digitization.  Do not use. 

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[171]]   <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[171]] <- "Multiple errors in matrix digitization, and survival issue with original data in thesis.  Do not use."
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[171]] <- NA

####  "Trollius_europaeus"

study_list[172]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[172],] 

#### original authors added very small values
## to cells with transitions know to exist in the species
## but unobserved in the field.  I have added these
## values here in the appropriate cells according to 
## the supplementary material provided by authors. 
## The is necessary to maintain reducibility of these matrices. 
## value added = 0.00000000000000000001
mat1[CA$UID1][[1]]$matU[4,3] <- 0.00000000000000000001
mat1[CA$UID1][[1]]$matU[6,5] <- 0.00000000000000000001
mat1[CA$UID1][[1]]$matF[1,6] <- 0.00000000000000000001

mat1[CA$UID1][[1]]$matA <- mat1[CA$UID1][[1]]$matF + mat1[CA$UID1][[1]]$matU

#####

mat1[CA$UID1][[2]]$matU[6,5] <- 0.00000000000000000001
mat1[CA$UID1][[2]]$matF[1,6] <- 0.00000000000000000001

#### F matrix, does not match original. cell 1,5 should <- 0

mat1[CA$UID1][[2]]$matF[1,5] <- 0

mat1[CA$UID1][[2]]$matA <- mat1[CA$UID1][[2]]$matF + mat1[CA$UID1][[2]]$matU

####

mat1[CA$UID1][[3]]$matF[1,6] <- 0.00000000000000000001
mat1[CA$UID1][[3]]$matA <- mat1[CA$UID1][[3]]$matF + mat1[CA$UID1][[3]]$matU

#### original matrix  for JAGab 2007-2008 not found in supplement

mat1[CA$UID1][[5]]$matU[2,1] <- 0.00000000000000000001
mat1[CA$UID1][[5]]$matA <- mat1[CA$UID1][[5]]$matF + mat1[CA$UID1][[5]]$matU

#####

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[172]]   <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[172]] <- "small constants (1e-20) added to cells where no transition was observed but must occur, as per original authors.  One fecundity value fixed to match original"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[172]] <- "abandoned wet grassland"

#### "Trollius_laxus"

study_list[173]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[173],] 

### fix longitude
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[173]]   <- -75.78

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[173]]   <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[173]] <- "longitude value fixed"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[173]] <- "forested fen"

##### "Trollius_laxus_2"

study_list[174]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[174],] 

### This is the same study as "Trollius_laxus", but the matrices in 
### "Trollius_laxus" are averaged accross populations and here the populations 
### are maintained separately.  However, these individual populations are quite small (eg. some <50)
### and all in the same site location. Use merged matrices above
## for "Trollius_laxus", and treat site as single location. 


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[174]]   <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[174]] <- "This is the same study as Trollius_laxus, but populations are split into smaller units, leading to small sample sizes for individual transitions"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[174]] <- "forested fen"


##### "Tsuga_canadensis"

study_list[175]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[175],] 


### fix longitude 
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[175]]   <- -78.37 


####

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[175]]   <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[175]] <- "Longitude fixed"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[175]] <- "forest"

##### "Ulex_gallii"

study_list[176]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[176],] 


#####  use matrix from core of range in Devon instead so that
###  lat long can be calculated more accurately

matrixClass1[CA$UID1][[1]] <- matrixClass1[CA$UID1][[1]][1:14,]
matA <- matrix(c(0.333,  0.4, 0.0,  0.0,  10.6,  10.6,  10.6,  10.6,  10.6 , 10.6,  10.6 , 10.6, 10.6, 10.6,
                 0.112, 0.485, 0,    0,    0,      0,     0,     0,     0,     0,     0,     0,    0,   0,   
                    0,   0.5, 0.556, 0,    0,      0,     0,     0,     0,     0,     0,     0,    0,   0,  
                    0,   0,   0.286, 0.511,0,      0,     0,     0,     0,     0,     0,     0,    0,   0, 
                    0,   0,   0.143, 0.368,0.556,  0,     0,     0,     0,     0,     0,     0,    0,   0, 
                    0,   0,     0,   0.053,0.429, 0.485,  0,     0,     0,     0,     0,     0,    0,   0, 
                    0,   0,     0,   0.053,   0 , 0.409,0.485,   0,     0,     0,     0,     0,    0,   0,
                    0,   0,     0,      0,    0,  0.091,0.364,0.595,    0,     0,     0,     0,    0,   0,
                    0,   0,     0,      0,    0,   0,   0.136,0.268,0.344,     0,     0,     0,    0,   0,
                    0,   0,     0,      0,    0,   0,     0,  0.073,0.357, 0.604,     0,     0,    0,   0,
                    0,   0,     0,      0,    0,   0,     0,  0.049,0.214, 0.143, 0.293,     0,    0,   0,
                    0,   0,     0,      0,    0,   0,     0,    0,  0.071, 0.238, 0.692, 0.727,    0,   0,
                    0,   0,     0,      0,    0,   0,     0,    0,     0,    0,     0,   0.258, 0.610,  0,
                    0,   0,     0,      0,    0,   0,     0,    0,     0,    0,     0,     0,   0.375, 0.984),
                  nrow = 14, byrow = T)

matF <- matA
matF[2:14,] <- 0
matU <- matA - matF


### insert these into mat1 object
mat1[CA$UID1][[1]]$matA <- matA
mat1[CA$UID1][[1]]$matF <- matF
mat1[CA$UID1][[1]]$matU <- matU

####
sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[176]]   <- 50.71
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[176]]   <- -3.32

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[176]]   <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[176]] <- "changed to individual matrix for the core range, instead of spatial mean. coordinates added based on paper."
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[176]] <- "heath"

#### "Veronica_arvensis"

study_list[177]

### already done above.. 

###  "Viola_stagnina"

study_list[178]
CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[178],] 


### not 'unmanipulated' by my definitions - author states that these are 'extensively 
# used' and 'annually mown meadows'

sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[178]]   <- "N"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[178]] <- "Not unmaipulated, author states these are extensively 
                                                       used and annually mown meadows"
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[178]] <- "annually mown meadows"

###### "Vochysia_ferruginea"

study_list[179]

CA <- sp_to_check_in_try[sp_to_check_in_try$SpeciesAuthor == 
                           study_list[179],] 

### mismatch between matrix A and U and original (typo).
### correct cell [4,4] from 0.310 to 0.731

mat1[CA$UID1][[1]]$matU[4,4] <- 0.731
mat1[CA$UID1][[1]]$matA[4,4] <- 0.731


#### add latitude and longitude, based on paper.  Note, 
### matrix is a spatial average based on 3 sites of similar elevation
## (10-60m)  within 50km of each other.  I have used the 
## coordinates from the third site 'Fonseca' which is between the other
## two.  
## 


sp_to_check_in_try$Lat[sp_to_check_in_try$SpeciesAuthor == study_list[179]]   <- 12.27
sp_to_check_in_try$Lon[sp_to_check_in_try$SpeciesAuthor == study_list[179]]   <- -83.97

sp_to_check_in_try$Observation[sp_to_check_in_try$SpeciesAuthor == study_list[179]]   <- "spatially averaged matrix, study conducted during post-hurricane recovery"


sp_to_check_in_try$checkY[sp_to_check_in_try$SpeciesAuthor == study_list[179]]   <- "Y_edits"
sp_to_check_in_try$checknote[sp_to_check_in_try$SpeciesAuthor == study_list[179]] <- "Matrix typo corrected, lat long added. Note matrix constructed during post-hurricane recovery of forest" 
sp_to_check_in_try$habitatAuthor[sp_to_check_in_try$SpeciesAuthor 
                                 == study_list[179]] <- "Tropical forest"

#### ---- data checking complete. Export .csv and compadre style object. ---- 

filename1 <- paste("metadata_check_180_species", Sys.Date(), ".csv", sep = "")

write.csv(sp_to_check_in_try, filename1, row.names = F)

##### --- select matrices not marked "N" in column checkY ---- 

summary(as.factor(sp_to_check_in_try$checkY))

sel1 <- which(sp_to_check_in_try$checkY == "N")

sp_selected <- sp_to_check_in_try[-sel1,]

length(unique(sp_selected$SpeciesAccepted))
### 96

summary(as.factor(sp_to_check_in_try$checkY))
##### get matching matrices and matrixClass info

mat_sel <- mat1[sp_selected$UID1]
matrixClass_sel <- matrixClass1[sp_selected$UID1]

####  double check that no species with clonality have slipped in. 

clonal1 <-c()
# 
# ###
# 
for(i in 1:length(mat_sel)) {
  clonal1[i] <- sum(mat_sel[[i]]$matC)
}
# 
# #
length(which(clonal1 == 0))
## [1] 663
length(which(clonal1 > 0))
## [1] 0

#### output as 'compadre style' object, for further use. 

obs1 <- list()
obs1[[1]] <- sp_selected
obs1[[2]] <- mat_sel
obs1[[3]] <- matrixClass_sel
obs1[[4]] <- compadre$version

names(obs1) <- c("metadata", "mat", "matrixClass", "version")

saveRDS(obs1,
        file =  "ammended_annual_compadre_data_08_06_2018.RData")

#####

