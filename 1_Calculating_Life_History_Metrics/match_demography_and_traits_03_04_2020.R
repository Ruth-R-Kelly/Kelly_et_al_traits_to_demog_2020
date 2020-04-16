### code to merge final demographic metrics dataset with trait data.

## Last edited 08/06/2018

### remove objects from memory 
rm(list = ls())


#### read in trait data
tms <- read.csv("all_trait_means20_04_2018_added_datasets.csv")
names(tms)
summary(tms)


### read in demographic metrics
metric_means <- read.csv("demo_means_selected_03_04_2020.csv")

summary(metric_means)

### merge with trait dataset
try_demog <- merge(tms, metric_means, by.x = "Name_matched", by.y = "Name_matched", all = FALSE)

#######

nrow(try_demog)
## 90 

names(try_demog[1:7])
# [1] "Name_matched" "Height_mean"  "LA_mean"      "LN_mean"      "LMA_mean"     "SSD_mean"    
# [7] "SM_mean"      

nrow(na.omit(try_demog[c(2)]))
### 90 have height
nrow(na.omit(try_demog[c(3)]))
## [1] 81 have leaf area
nrow(na.omit(try_demog[c(4)]))
### 74 have LN
nrow(na.omit(try_demog[c(5)]))
### 90 have LMA
nrow(na.omit(try_demog[c(6)]))
### 64 have stem density
nrow(na.omit(try_demog[c(7)]))
### 88 have seed mass

######## optimal number of species with 3 traits, always with height and seedmass

# [1] "Name_matched" "Height_mean"  "LA_mean"      "LN_mean"      "LMA_mean"     "SSD_mean"    
# [7] "SM_mean"  

nrow(na.omit(try_demog[c(2,3,7)])) # Height, LA, SM
## 80
nrow(na.omit(try_demog[c(2,4,7)])) # Height, LN, SM
##72
nrow(na.omit(try_demog[c(2,5,7)])) # Height, LMA, SM
## 88
nrow(na.omit(try_demog[c(2,6,7)])) # Height, SSD, SM
## 62
nrow(na.omit(try_demog[c(2,3,5)])) # Height, LA, LMA
#81

##### 6 traits
nrow(na.omit(try_demog[,2:7]))
## 52

cor(na.omit(try_demog[,2:7], method = "spearman"))

#              Height_mean      LA_mean     LN_mean    LMA_mean   SSD_mean    SM_mean
# Height_mean  1.000000000 -0.006107886 -0.07682982  0.41349438  0.4151407 0.20432344
# LA_mean     -0.006107886  1.000000000  0.29451235 -0.21276076 -0.1367645 0.02607899
# LN_mean     -0.076829822  0.294512345  1.00000000 -0.54720110 -0.1385651 0.10519988
# LMA_mean     0.413494375 -0.212760760 -0.54720110  1.00000000  0.3345822 0.00565589
# SSD_mean     0.415140749 -0.136764473 -0.13856512  0.33458219  1.0000000 0.38300393
# SM_mean      0.204323443  0.026078990  0.10519988  0.00565589  0.3830039 1.00000000

######## optimal number of species with 4 traits, always with height and seedmass
names(try_demog[,2:7])
## [1] "Height_mean" "LA_mean"     "LN_mean"     "LMA_mean"    "SSD_mean"    "SM_mean"  
nrow(na.omit(try_demog[c(2,3,4,7)])) # Height, LA, LN, Seed size
## 64
nrow(na.omit(try_demog[c(2,3,5,7)])) ## height, LA, LMA, Seed size
##80
nrow(na.omit(try_demog[c(2,3,6,7)])) ## height, LA, SSD, Seed size
## 59 
nrow(na.omit(try_demog[c(2,4,5,7)]))## height, LN, LMA, Seed size
## 72
nrow(na.omit(try_demog[c(2,4,6,7)])) ## height, LN, SSD, Seed size
## 55
nrow(na.omit(try_demog[c(2,5,6,7)]))## height, LMA, SSD, Seed size
## 62

#### ## height, LA, LMA, Seed size, appears to be the best compromise between 
#### sample size and trait combination # 83 at 4 traits.

###### optimal 5 traits.. 


nrow(na.omit(try_demog[c(2,3,4,5,7)]))##  height, LA, LN, LMA, Seed size
## 64
nrow(na.omit(try_demog[c(2,3,4,6,7)]))## height, LA, LN, SSD, Seed size
## 52
nrow(na.omit(try_demog[c(2,4,5,6,7)]))## height, LN, LMA, SSD, Seed size
## 55  
#######
nrow(na.omit(try_demog[c(2,3,5,6,7)]))## height, LA, LMA, SSD, Seed size
## 59
###########################

### write dataset of matching traits and demography for all species.. 
write.csv(try_demog, "matching_trait_demo_data_08_06_2018.csv")

#### ---- write dataset of matching demography with height, LA, LMA, Seed size ---- 

# This is the dataset that feeds into the next step given the folder 2. 

try_4traits <- na.omit(try_demog[c(1,2,3,5,7)]) ## height, LA, LMA, Seed size

try_4traits2 <- merge(try_4traits, metric_means, by.x = "Name_matched", by.y = "Name_matched", all = FALSE)

write.csv(try_4traits2, "matching_4traits_demo_data_03_04_2020.csv")

#### ---- write dataset of 4 traits matching height, LN, LMA, Seed size 

# Note: this was done for some comparisons that may not have made it into the
# final paper. Basically, Leaf Area works better than Leaf Nitrogen, and 
# has some preferable physiological properties too (in my view).

try_4_LN_traits <- na.omit(try_demog[c(1,2,4,5,7)]) ## height, LN, LMA, Seed size

try_4_LN_traits2 <- merge(try_4_LN_traits, metric_means, by.x = "Name_matched", by.y = "Name_matched", all = FALSE)

write.csv(try_4_LN_traits2, "matching_4_H_LN_LMA_SM_traits_demo_data_03_04_2020.csv")

###### 