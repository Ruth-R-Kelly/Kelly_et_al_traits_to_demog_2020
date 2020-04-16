### Code to assess sensibleness of computed demographic metrics and join these with trait data
## Ruth Kelly 
## last edits 03/04/2020

### Editted to export full dataset i.e. values for each species at 
### each gps location in addition to the a single value per species. 


# Note: Trait data now includes control conditions from experiments, (i.e. nutrient addition and CO2 treatments are excluded throughout). 
# 
# library("dplyr")


### remove old files

rm(list = ls())


data1 <- read.csv("demography_non_seedbank_03_2020.csv")
data2 <- read.csv("demography_1yr_seedbank_03_2020.csv")
data3 <- read.csv("demography_2yr_seedbank_03_2020.csv")
data4 <- read.csv("demography_3yr_seedbank_03_2020.csv")



demogs <- rbind(data1,data2, data3, data4)

# general investigation
length(unique(demogs$SpeciesAccepted))
summary(demogs$OrganismType)
nrow(demogs)

#### go through demographic metrics one by one for plausibility.. 

names(demogs[54:69])

### Columns of main interest to check are:
# [1] "Mean_repro"               "Net_repro"                "Age_50_dead"             
# [4] "Age_99_dead"              "Age_999_dead"             "survival_index99"        
# [7] "Survival_index999"        "Survival_index_La_RepExp" "Growth"                  
# [10] "Retrogression"            "Repro_life_expectancy"    "Age_at_Maturity"         
# [13] "GiniF_repro"              "Gini_life99"              "Gini_life999"            
# [16] "Gen_time"      

########  

summary(demogs[54:69])

### exclude rows where Age_99_dead or Age_999_dead did not calculate correctly.
### These were mainly species with survival Issues (e.g. where no death was observed in the 
### study)

#demogs$SurvivalIssue[which(is.na(demogs$Age_99_dead))]
# # [1] 1 1 1 1 1 1 1
# 
unique(demogs$SpeciesAuthor[which(is.na(demogs$Age_99_dead))])
## # [1] Molinia_caerulea Primula_elatior  Clidemia_hirta   Sapium_sebiferum

#### Remove these.
d_sub <- demogs[-which(is.na(demogs$Age_99_dead)),]

nrow(d_sub)
##247

summary(d_sub[54:69])

#### recheck for errors in other columns... 


### Still some issues with Age_at_Maturity, life_span outliers.. 

##Look at Age_at_Maturity

hist(d_sub$Age_at_Maturity)

summary(d_sub$Age_at_Maturity)

quantile(d_sub$Age_at_Maturity, 0.95)

# view rows where age at maturity > than 95th quantile
sub_late_maturity <- d_sub[which(d_sub$Age_at_Maturity > quantile(d_sub$Age_at_Maturity, 0.95)),]
sub_late_maturity[,c(3,11,48,65)]
# Most of the highest are  Bertholletia excelsa which is known in literature to hav mature ages of > 100
# silene acaulis discussed with tundra biologists and has a very slow life-cycle

### make a cut off at 200 years

d_sub <- d_sub[-which(d_sub$Age_at_Maturity > 200),]
### this removes 1 record for Eperua falcata - 

######  life span

hist(d_sub$Age_99_dead)
sub_long_life <- d_sub[which(d_sub$Age_99_dead > quantile(d_sub$Age_99_dead, 0.95)),]
sub_long_life[,c(3,11,48,57,58)]

d_sub[d_sub$SpeciesAccepted == "Lathyrus vernus",]
### remove "Lathyrus vernus" Age_99_dead more than 100.  This is very unrealistic based on the
## species biology and more than twices as long as the next longest calculation in our dataset. 
# Also, survivalIssue  = 1 for this record
d_sub <- d_sub[-which(d_sub$SpeciesAccepted == "Lathyrus vernus" & d_sub$Age_99_dead > 100),] 

d_sub[d_sub$SpeciesAccepted == "Trillium grandiflorum",] 

### remove this "Trillium grandiflorum" Age_99_dead ==  391.  This is very unrealistic based on the
## species biology (max age in literature = 70, Knight 2006) and 3 times as long as the next longest calculation in our dataset. Also, survivalIssue 
## = 1

### others although very long seem plausible 

d_sub <- d_sub[-which(d_sub$SpeciesAccepted == "Trillium grandiflorum" & d_sub$Age_99_dead == 391),] 

summary(d_sub[,48:69])

#####

## for which species is Gini_life99 = NA

check_gini <- d_sub[which(is.na(d_sub$Gini_life99)),] 
check_gini$OrganismType

## remove species where Gini_life999 = NA

d_sub <- d_sub[-which(is.na(d_sub$Gini_life999)),] 

#### check annual species to make sure estimates are sensible

check_ann <- d_sub[which(d_sub$OrganismType == "Annual"),]
check_ann
##################################################################

#### remove tree species which have life spans calculated 
# < 10 but which are known to be long lived

d_sub_tree <- d_sub[d_sub$OrganismType == "Tree",]
tree_wrong_age <- levels(droplevels(d_sub_tree$SpeciesAuthor[which(d_sub_tree$Age_999_dead < 10)]))
tree_wrong_age

TW <- which(d_sub$SpeciesAuthor %in% tree_wrong_age)
TW

####

d_sub <- d_sub[-TW,]
##[1] 236  76

length(unique(demogs$SpeciesAccepted))
# 94
length(unique(d_sub$SpeciesAccepted))
##[1] 90

####### We've lost 4 species through this process.

# setdiff(demogs$SpeciesAccepted,d_sub$SpeciesAccepted)

######################################################################

### Name matching... 

species_mat <- read.csv("compadre_synonyms_iplant_web_checked27_3_2017.csv")

length(levels(species_mat$Name_submitted))
# 572
length(levels(species_mat$Name_matched))
#569

spe_comp <- unique(d_sub$SpeciesAccepted)
length(spe_comp)
x1 <- intersect(spe_comp, unique(species_mat$Name_matched))

length(x1)
#90

#### All match : )

d_sub$Name_matched <- as.character(d_sub$SpeciesAccepted)

################ 

d_sub$Name_matched <- factor(d_sub$Name_matched)


### some edits to GPS locations are needed based on rechecking of 
# original papers.  

names(d_sub)

# for "Actaea_spicata" location fixed by locating forest site described in paper
d_sub$Lat[d_sub$Name_matched == "Actaea spicata"] <- 58.90
d_sub$Lon[d_sub$Name_matched == "Actaea spicata"] <- 17.92

# Cirsium_palustre location fixed by relocating meadow site described in paper

d_sub$Lat[d_sub$Name_matched =="Cirsium palustre"] <- 58.95
d_sub$Lon[d_sub$Name_matched =="Cirsium palustre"] <- 17.60

# "Scorzonera hispanica" location fixed - lat and long inputted in incorrect order

d_sub$Lat[d_sub$Name_matched == "Scorzonera hispanica"] <- 50.52
d_sub$Lon[d_sub$Name_matched =="Scorzonera hispanica"] <- 14.19

d_sub$Lat[d_sub$Name_matched == "Primula farinosa" &
        d_sub$MatrixPopulation == "Flottskär"] <- 59.93
d_sub$Lon[d_sub$Name_matched == "Primula farinosa" &
            d_sub$MatrixPopulation == "Flottskär"] <-  18.87


# 
summary(d_sub$Lat)
summary(d_sub$Lon)
### in order to merge these we need to fix the absence of Study duration info 
## for  Epipactis_atrorubens, according to the original paper this should be 8
summary(d_sub$StudyDuration)
d_sub[is.na(d_sub$StudyDuration),]
d_sub$StudyDuration[which(is.na(d_sub$StudyDuration))] <- 8

d_sub$gps_loc <- paste(d_sub$Lat, d_sub$Lon)

#### use function summarise 'dplyr' to get mean values of demographic traits 
### per gps location

d_sub$Ecoregion <- as.character(d_sub$Ecoregion)

all_metric_means <- as.data.frame(summarise(group_by(d_sub, Name_matched, 
                          OrganismType, Class, Ecoregion, Lat, Lon, gps_loc, habitatAuthor),
                          Age_at_Maturity = mean(Age_at_Maturity, na.rm = T),
                          Mean_repro = mean(Mean_repro, na.rm = T),
                          Net_repro = mean(Net_repro, na.rm = T),
                          Age_50_dead = mean(Age_50_dead, na.rm = T),
                          Age_99_dead = mean(Age_99_dead, na.rm = T),
                          Age_999_dead = mean(Age_999_dead, na.rm = T),
                          survival_index99 = mean(survival_index99, na.rm = T),
                          survival_index999 = mean(Survival_index999, na.rm = T),
                          survival_index_LaPlusRepExp = mean(Survival_index_La_RepExp, na.rm = T),
                          prog_growth = mean(Growth, na.rm = T),
                          retro_growth = mean(Retrogression, na.rm = T),
                          Repro_life_expectancy = mean(Repro_life_expectancy, na.rm = T),
                          GiniF_999 = mean(Gini_life999, na.rm = T),
                          GiniF_99 = mean(Gini_life99, na.rm = T),
                          Gen_time = mean(Gen_time, na.rm = T),
                          NSeedStages = mean(NSeedStages, na.rm = T),
                          Studylength=  max(StudyDuration, na.rm = T),
                          MatrixDim= max(MatrixDimension, na.rm = T)))

summary(all_metric_means)


write.csv(all_metric_means, "metric_means_for_all_species_locations_by_habitats_13_03_2020.csv")

length(unique(all_metric_means$Name_matched))
#90

#### Select species with multiple populations #### 
# so the populations in the same habitat and within 5km can be merged. 
# Merge population metrics for populations within 5km 


sum1 <- summary(all_metric_means$Name_matched)
sum2 <- sum1[sum1 >1]

sum2
#  Abies concolor        Actaea spicata    Alliaria petiolata 
#       3                     2                    11 
# Bertholletia excelsa        Carduus nutans      Castanea dentata 
#       2                     3                     6 
# Cytisus scoparius       Lathyrus vernus       Orchis purpurea 
#        4                     3                     2 
# Paeonia officinalis       Pinus nigra      Primula farinosa 
#          2                     2                     3 
# Primula vulgaris        Silene acaulis Trillium grandiflorum 
#        4                     4                    11 

### subset d_sub to find these. . 
names(sum2)
# 15 species

multi_pop <- droplevels(all_metric_means[all_metric_means$Name_matched %in% names(sum2),])
multi_pop$Name_hab <- factor(paste(multi_pop$Name_matched, multi_pop$habitatAuthor))

unique(multi_pop$Name_matched)
# 15 species
unique(multi_pop$Name_hab)

## remove these species from all_metric_means dataframe - will be later remerged

keep_names <- setdiff(all_metric_means$Name_matched, names(sum2))

all_metric_means <- all_metric_means[all_metric_means$Name_matched %in% keep_names,]
unique(all_metric_means$Name_matched)
### 
# 75 species

sum3 <- summary(multi_pop$Name_hab)
sum4 <- sum3[sum3>1]
sum4


# a dataframe with species which contains records where multiple matrices 
# exist in the same habitat. We can now search this to see if populations 
# occur within 5km 


multi_pop_hab <- droplevels(multi_pop[multi_pop$Name_hab %in% names(sum4),])



### split into species 
spl1 <- split(multi_pop_hab, multi_pop_hab$Name_hab)

### species with more than one matrix is a given habitat. 
names(spl1)


### calcuate distances between populations for each species 
# in list of species with multiple matrices within a habitat. 
# 
#library(geosphere)

### check lat and long columns make sense. 
summary(multi_pop$Lat)
summary(multi_pop$Lon)


#### ---- Calculate geodesic distances (based on ellipsoid distance) ----

## distGeo - is calculating a 'geodesic' distance based 
## by default on a standard WGS1984 ellipsoid.  
## Resulting measurement is in metres and but I have converted to km
## by dividing by 1000. 

# NOTE COORDINATES MUST BE NUMERIC < THERE IS NO WARNING IF YOU ACCIDENTALLY INPUT OTHERWISE :)

distx <- c()
maxdists <- c()
mindists <- c()
full_dists_comb <- list()
length(spl1)

for(j in 1:length(spl1)) {

  distx <- c()
  combins <- c()
  combins <- t(combn(1:nrow(spl1[[j]]),2)) 

  for(i in 1:nrow(combins)){
 
 
  x1 <- combins[i,1]
  x2 <- combins[i,2]
  
  distx[i] <- distGeo(p1=c(spl1[[j]]$Lon[x1],spl1[[j]]$Lat[x1]),
                      p2 = c(spl1[[j]]$Lon[x2],spl1[[j]]$Lat[x2]))/1000
  }
 
  
maxdists[j] <- max(distx)
mindists[j] <- min(distx)

res1 <- cbind(combins,distx) 
colnames(res1) <- c("pop1", "pop2", "distance")
full_dists_comb[[j]] <- res1

print(j)
}

summary(res1)
### find species where there are populations within 5km

### do any species have all their populations within 5km? 
which(maxdists<5) # integer(0)

# Which have any within 5km
which(mindists<5)
#[1]  6 10

# find species habitat combinations where metrics can be averaged 
names(spl1)[which(mindists<5)]
# [1] "Cytisus scoparius urban grasslands and prairie"
# [2] "Trillium grandiflorum forest"  


### look at these in turn


## [1] "Trillium grandiflorum forest"

which(full_dists_comb[[10]][,3]<5 )

full_dists_comb[[10]]
## [1] 20 46

# 3    4  4.469915
# 7    8  2.489181

### average populations 7 and 8 and remove the others - these are closest
### together

### check all Trillium grandiflorum are in fact forest 

multi_pop[multi_pop$Name_matched == "Trillium grandiflorum",]



### remove this species from the multipop dataset, 
# because we are about to recalculate it and return it 
remove1 <- which(multi_pop$Name_hab == "Trillium grandiflorum forest")
length(remove1)

remove2 <- remove1[c(1:6,9:11)]

multi_pop <- multi_pop[-remove2,] 


#### average the coordinates 
mean_lat <- mean(multi_pop$Lat[which(multi_pop$Name_hab == "Trillium grandiflorum forest")])
mean_long <- mean(multi_pop$Lon[which(multi_pop$Name_hab == "Trillium grandiflorum forest")])

#### stick this centroid in for the location for all populations of Silene
multi_pop$Lat[which(multi_pop$Name_hab == "Trillium grandiflorum forest")] <- mean_lat
multi_pop$Lon[which(multi_pop$Name_hab == "Trillium grandiflorum forest")] <- mean_long
multi_pop$gps_loc[which(multi_pop$Name_hab == "Trillium grandiflorum forest")] <- paste(mean_lat, mean_long)

multi_pop[which(multi_pop$Name_hab == "Trillium grandiflorum forest"),] 

###### 

which(mindists<5)

names(spl1)[6] ##  "Cytisus scoparius urban grasslands and prairie"

full_dists_comb[[6]]

### only two populations are within 5km average these and drop the other populations
### pops 3 and 4


### check all Cytisus scoparius are urban grasslands and prairie

length(which(multi_pop$Name_matched == "Cytisus scoparius"))

remove1 <- which(multi_pop$Name_hab == "Cytisus scoparius urban grasslands and prairie")
length(remove1)

remove2 <- remove1[c(1:2)]
remove2
multi_pop <- multi_pop[-remove2,] 


#### average the coordinates 
mean_lat <- mean(multi_pop$Lat[which(multi_pop$Name_hab == "Cytisus scoparius urban grasslands and prairie")])
mean_long <- mean(multi_pop$Lon[which(multi_pop$Name_hab == "Cytisus scoparius urban grasslands and prairie")])

#### stick this centroid in for the location for all populations of Silene
multi_pop$Lat[which(multi_pop$Name_hab == "Cytisus scoparius urban grasslands and prairie")] <- mean_lat
multi_pop$Lon[which(multi_pop$Name_hab == "Cytisus scoparius urban grasslands and prairie")] <- mean_long
multi_pop$gps_loc[which(multi_pop$Name_hab == "Cytisus scoparius urban grasslands and prairie")] <- paste(mean_lat, mean_long)

multi_pop[which(multi_pop$Name_hab == "Cytisus scoparius urban grasslands and prairie"),] 

####### use summarise again to recalculate averages for this group 


multipop_metric_means <- as.data.frame(summarise(group_by(multi_pop, Name_matched, 
                                                     OrganismType, Class, Ecoregion, Lat, Lon, gps_loc, habitatAuthor),
                                            Age_at_Maturity = mean(Age_at_Maturity, na.rm = T),
                                            Mean_repro = mean(Mean_repro, na.rm = T),
                                            Net_repro = mean(Net_repro, na.rm = T),
                                            Age_50_dead = mean(Age_50_dead, na.rm = T),
                                            Age_99_dead = mean(Age_99_dead, na.rm = T),
                                            Age_999_dead = mean(Age_999_dead, na.rm = T),
                                            survival_index99 = mean(survival_index99, na.rm = T),
                                            survival_index999 = mean(survival_index999, na.rm = T),
                                            survival_index_LaPlusRepExp = mean(survival_index_LaPlusRepExp, na.rm = T),
                                            prog_growth = mean(prog_growth, na.rm = T),
                                            retro_growth = mean(retro_growth, na.rm = T),
                                            Repro_life_expectancy = mean(Repro_life_expectancy, na.rm = T),
                                            GiniF_999 = mean(GiniF_999, na.rm = T),
                                            GiniF_99 = mean(GiniF_99, na.rm = T),
                                            Gen_time = mean(Gen_time, na.rm = T),
                                            NSeedStages = mean(NSeedStages, na.rm = T),
                                            Studylength=  max(Studylength, na.rm = T),
                                            MatrixDim= max(MatrixDim, na.rm = T)))

summary(multipop_metric_means)

unique(multipop_metric_means$Name_matched)
# still 15 



#### I have disabled the next section of code, which selected the random species 
# where duplicates existed. And inserted code which will select those species records
# that were used in the paper. The randomisation code here, is how these were initially 
# selected 
# all_metric_means <- rbind(all_metric_means, multipop_metric_means)

### randomly select one row per species for further analyses

# all_metric_means$UID <- seq(1:nrow(all_metric_means))
# dim(all_metric_means)
# spl1 <- split(all_metric_means,all_metric_means$Name_matched)
# names(all_metric_means)
# unique(all_metric_means$Name_matched)
# ##90

# sample1 <- c()
# 
# for(i in 1:length(spl1)) {
# 
# sample1[i] <-spl1[[i]][sample(nrow(spl1[[i]]), 1), 27]
# 
# }
# 
# 
# sample1


# metric_means <- droplevels(all_metric_means[all_metric_means$UID %in% sample1,])

#### 

multipop_metric_sel <- multipop_metric_means[c(3,4,12,18,20,22,28,31,32,35,36,39,43,46,49),]
dim(multipop_metric_sel)
unique(multipop_metric_sel$Name_matched)

metric_means <- rbind(all_metric_means, multipop_metric_sel)

metric_means <- droplevels(metric_means)

#### 
length(unique(metric_means$Name_matched))
## 90

######
summary(metric_means)

## for reference
summary(as.factor(metric_means$Ecoregion))
# BOR      DES      FGS      MED      TBM TBM; BOR      TCF      TGS      TGV 
# 7        5        3        5       44        1        4        7        2 
#      TMB      TUN     NA's 
#       16        2        1 


summary(as.factor(metric_means$OrganismType))
# Herbaceous perennial     Palm                Shrub                 Tree 
# 53                        3                   13                   23 
# Annual 
# 5 

nrow(metric_means)

write.csv(metric_means,"demo_means_selected_03_04_2020.csv")

###########

