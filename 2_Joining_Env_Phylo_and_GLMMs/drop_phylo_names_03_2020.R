### code to match phylogenetic tree to trait/demo dataset. 
## 07/04/2020

# library("ape")
# library("dplyr")
#####

### remove uneeded files 
rm(list = ls())


#### ---- read in tree ----
treefile <- read.tree(file ="Vascular_Plants_rooted.dated.tre")

tip_labels <- as.character(treefile$tip.label)
#str(tip_labels)

#### ---- Read in demo and traits ----
### read in .csv trait demography file for 4 traits to figure out which tips are needed.

demog <- read.csv("matching_4traits_demo_data_03_04_2020.csv", row.names = "X.1")

summary(demog)

#### Check which species are and are not already in the phylogeny ####
plants <- as.character(demog$Name_matched)
length(plants)

plants2 <- gsub(' ', '_', plants)
head(plants2)
head(tip_labels)
length(intersect(plants2, tip_labels))
setdiff(plants2, tip_labels)
# [1] "Cirsium_acaule"                "Cirsium_pannonicum"           
# [3] "Epipactis_atrorubens"          "Himantoglossum_hircinum"      
# [5] "Minuartia_obtusiloba"          "Orchis_purpurea"              
# [7] "Petrophile_pulchella"          "Sanicula_europaea"            
# [9] "Sapium_sebiferum"              "Saxifraga_tridactylites"      
# [11] "Sphaeralcea_coccinea"          "Stryphnodendron_microstachyum" 


#### Add species not present according to closest in literature ####


#### ---- Add Cirsiums ---- 

### put these plants into the phylogeny in the correct genera.
### this is a quick fix and should be refined later

## Use plant position of "Cirsium_setosum" for "Cirsium_acuale"
## Use plant position of "Cirsium_helenioides" for "Cirsium_pannonicum"
##  they are  unneeded species in the correct genus. 

plants2[grep("Cirsium",plants2)]
# [1] "Cirsium_acaule"     "Cirsium_palustre"   "Cirsium_pannonicum"  

plants2[grep("Cirsium",plants2)]

tip_labels[grep("Cirsium", tip_labels)]


######ADD IN SPECIES NOT IN TANK TREE based on positions of congeneric species #####

#### I haven't got an up to date phylogeny for cirsium so I'm just 
## placing these in the genus.

#Code for grafting species onto the tank tree
wherePlaced <- c("Cirsium_setosum")
placement <- which(treefile$tip.label == wherePlaced)
treefile$tip.label[placement] <- "Cirsium_acaule" 
########################################

#Code for grafting species onto the tank tree
wherePlaced <- c("Cirsium_helenioides")
placement <- which(treefile$tip.label == wherePlaced)
treefile$tip.label[placement] <- "Cirsium_pannonicum" 
########################################

tip_labels <- as.character(treefile$tip.label)
tip_labels[grep("Cirsium", tip_labels)]

#### ---- "Epipactis_atrorubens" ---- 

plants2[grep("Epipactis",plants2)]
# [1] "Epipactis_atrorubens"

tip_labels[grep("Epipactis", tip_labels)]

# [1] "Epipactis_palustris"   "Epipactis_microphylla" "Epipactis_helleborine"
# [4] "Epipactis_muelleri"    "Epipactis_leptochila" 

### check for known synonyms based on plant list

tip_labels[grep("Amesia", tip_labels)]
### not in phylogeny. 

### match to closest species in genus based in phylotree in 
### Jin et al. 2014 DOI 1471-2229/14/63
### This is "Epipactis_helleborine"

#Code for grafting species onto the tank tree
wherePlaced <- c("Epipactis_helleborine")
placement <- which(treefile$tip.label == wherePlaced)
treefile$tip.label[placement] <- "Epipactis_atrorubens"


## next

setdiff(plants2,treefile$tip.label)
# [1] "Himantoglossum_hircinum"      
# [2] "Minuartia_obtusiloba"         
# [3] "Orchis_purpurea"              
# [4] "Petrophile_pulchella"         
# [5] "Sanicula_europaea"            
# [6] "Sapium_sebiferum"             
# [7] "Saxifraga_tridactylites"      
# [8] "Sphaeralcea_coccinea"         
# [9] "Stryphnodendron_microstachyum"

#### ---- "Himantoglossum_hircinum" ----

plants2[grep("Himantoglossum",plants2)]
# [1]  "Himantoglossum_hircinum"

tip_labels[grep("Himantoglossum", tip_labels)]

# [1] "Himantoglossum_robertianum"

### match to the only other option in genus.

#Code for grafting species onto the tank tree
wherePlaced <- c("Himantoglossum_robertianum")
placement <- which(treefile$tip.label == wherePlaced)
treefile$tip.label[placement] <- "Himantoglossum_hircinum"

### next

#### ---- "Minuartia_obtusiloba"   ----

plants2[grep("Minuartia",plants2)]
# [1]  "Minuartia_obtusiloba"

tip_labels[grep("Minuartia", tip_labels)]

# [1] "Minuartia_picta"        "Minuartia_stricta"     
# [3] "Minuartia_verna"        "Minuartia_rubella"     
# [5] "Minuartia_rupestris"    "Minuartia_biflora"     
# [7] "Minuartia_laricifolia"  "Minuartia_groenlandica"

### check for synonym Lidia obtusiloba

tip_labels[grep("Lidia", tip_labels)]

### not present in tree

### No detailed phylogeny found for genus -
### match to another species in genus  "Minuartia_picta" 
 
#Code for grafting species onto the tank tree
wherePlaced <- c( "Minuartia_picta")
placement <- which(treefile$tip.label == wherePlaced)
treefile$tip.label[placement] <- "Minuartia_obtusiloba"

#### ---- next "Orchis_purpurea"   ----

plants2[grep("Orchis",plants2)]
# [1]  "Orchis_purpurea"

tip_labels[grep("Orchis", tip_labels)]

# [1] "Orchis_quadripunctata" "Orchis_militaris"      
# [3] "Orchis_anthropophora"    

### no non-orchis synonyms

### closest species is "Orchis_militaris"
### according to phylogeny in Jacquemyn et al. 
### New Phytologist (2011) 192: 518–528
### doi: 10.1111/j.1469-8137.2011.03796.x


#Code for grafting species onto the tank tree
wherePlaced <- c("Orchis_militaris")
placement <- which(treefile$tip.label == wherePlaced)
treefile$tip.label[placement] <- "Orchis_purpurea"

### next 

#### ---- "Petrophile_pulchella" ----


plants2[grep("Petrophile",plants2)]
# [1] "Petrophile_pulchella"

tip_labels[grep("Petrophile", tip_labels)]

# [1] "Petrophile_biloba"    "Petrophile_canescens"
# [3] "Petrophile_circinata"   

### no Petrophile pulchella is not recognised on The Plant List
## but is Tropicos

### no detailed phylogeny found for genus,
### use "Petrophile_biloba"  
### also recognised within genus in Tropicos.

#Code for grafting species onto the tank tree
wherePlaced <- c("Petrophile_biloba")
placement <- which(treefile$tip.label == wherePlaced)
treefile$tip.label[placement] <- "Petrophile_pulchella"

#### ---- "Sanicula_europaea" ---- 

### "Sanicula_europaea"

plants2[grep("Sanicula",plants2)]
# [1] "Sanicula_europaea"

### check synonyms

# tip_labels[grep("Astrantia", tip_labels)]
# tip_labels[grep("Caucalis", tip_labels)]

sort(tip_labels[grep("Sanicula", tip_labels)])

# "Sanicula_arctopoides"   "Sanicula_arguta"       
# [3] "Sanicula_bipinnatifida" "Sanicula_canadensis"   
# [5] "Sanicula_chinensis"     "Sanicula_crassicaulis" 
# [7] "Sanicula_elata"         "Sanicula_graveolens"   
# [9] "Sanicula_gregaria"      "Sanicula_odorata"      
# [11] "Sanicula_smallii"      

### use Smallii as closest based on tree in Calvino et al 2007
## doi:10.1016/j.ympev.2007.01.002
### * not all of these species are in this tree

#Code for grafting species onto the tank tree
wherePlaced <- c("Sanicula_smallii")
placement <- which(treefile$tip.label == wherePlaced)
treefile$tip.label[placement] <- "Sanicula_europaea"

setdiff(plants2, treefile$tip.label)

#### ---- "Sapium_sebiferum"  ----


plants2[grep("Sapium",plants2)]
# [1] "Sapium_sebiferum"

### no non-sapium synonyms

sort(tip_labels[grep("Sapium", tip_labels)])
   

### no genus level phylogeny found. Use "Sapium_haematospermum". 
#Code for grafting species onto the tank tree
wherePlaced <- c("Sapium_haematospermum")
placement <- which(treefile$tip.label == wherePlaced)
treefile$tip.label[placement] <- "Sapium_sebiferum"

setdiff(plants2, treefile$tip.label)

#### ---- "Saxifraga_tridactylites" ----

plants2[grep("Saxifraga",plants2)]
# [1] "Saxifraga_aizoides"      "Saxifraga_tridactylites"

### Tridactylites synonyms not present in tree

sort(tip_labels[grep("Tridactylites", tip_labels)])

####
sort(tip_labels[grep("Saxifraga", tip_labels)])

### Can't find a phylogeny with this species. 
### Used phylogeny from Conti et al. 1999 
# Molecular Phylogenetics and Evolution
# Vol. 13, No. 3, December, pp. 536–555, 1999
# Article ID mpev.1999.0673,

### which includes species - S. osloensis Knaben 
## from subsection Tridactylites.  Closest species to 
## this which is in both their tree and this one is 
## "Saxifraga_cespitosa"  

wherePlaced <- c("Saxifraga_cespitosa")
placement <- which(treefile$tip.label == wherePlaced)
treefile$tip.label[placement] <- "Saxifraga_tridactylites" 


#### ---- "Sphaeralcea_coccinea"  ----

plants2[grep("Sphaeralcea",plants2)]
# [1] "Sphaeralcea_coccinea"

### NO non-Sphaeralcea synoymns 

sort(tip_labels[grep("Sphaeralcea", tip_labels)])
# [1] "Sphaeralcea_ambigua"      "Sphaeralcea_angustifolia"
# [3] "Sphaeralcea_palmeri"

## No genus level phylogeny found. Place this species in genus 
## at position of "Sphaeralcea_ambigua" 



####

wherePlaced <- c("Sphaeralcea_ambigua")
placement <- which(treefile$tip.label == wherePlaced)
treefile$tip.label[placement] <- "Sphaeralcea_coccinea"

### ---- "Stryphnodendron_microstachyum" ---- 


plants2[grep("Stryphnodendron",plants2)]
# [1] "Stryphnodendron_microstachyum"

### NO non-Sphaeralcea synoymns 

sort(tip_labels[grep("Stryphnodendron", tip_labels)])
# [1] "Stryphnodendron_adstringens"  "Stryphnodendron_coriaceum"   
# [3] "Stryphnodendron_guianense"    "Stryphnodendron_polystachyum"
# [5] "Stryphnodendron_porcatum"     "Stryphnodendron_pulcherrimum"
# [7] "Stryphnodendron_racemiferum" 

## Closest matches are - 
## "Stryphnodendron_porcatum"     "Stryphnodendron_pulcherrimum"
## based on Phylogeny in Simon et al. 2016
## DOI: 10.1086/684077
###

wherePlaced <- c("Stryphnodendron_pulcherrimum")
placement <- which(treefile$tip.label == wherePlaced)
treefile$tip.label[placement] <- "Stryphnodendron_microstachyum"


#### ---- remove unused species from phylogeny ----

omit_spe <- as.character(setdiff(treefile$tip.label, plants2))

reduced_plant_tree <- drop.tip(treefile, omit_spe)

str(reduced_plant_tree)
class(reduced_plant_tree)
plot(reduced_plant_tree, cex = 0.4)

#### ---- export reduced phylogeny as nexus file ----

write.nexus(reduced_plant_tree, file = "reduced_tree_for_traits_03_2020.nex")
#####

#### --- sort demography dataset to match phylo order ---- ####

# Note this resorting is crucial for later MCMCglmm analyses

head(reduced_plant_tree$tip.label)
test1 <- demog
demog$Name_matched <- plants2

Labels <- reduced_plant_tree$tip.label

merge1 <- merge(as.data.frame(Labels), demog, by.x = "Labels", by.y = "Name_matched",
                sort = FALSE)

merge1$Labels == reduced_plant_tree$tip.label

### set these as rownames

row.names(merge1) <- merge1$Labels

head(merge1)
#### remove spurious X column

names(merge1)

# I have a column called X here which I want to rid of. 
merge2 <- merge1[,-6]

write.csv(merge2, "phylo_ordered_trait_demog_03_2020.csv")

##### 
