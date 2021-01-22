### Short script to add growth forms as per supplementary material and 
## Fig2 of the main paper, these growth forms are based on floras and 
## scientific literature for each species, and stored in "data_growth_form_lookup.csv"

plant_data <- read.csv("data_with_needed_env_vars_03_2020.csv", row.names = "X")
#head(plant_data)
#row.names(plant_data)

### add growth forms 

g_forms <- read.csv("data_growth_form_lookup.csv")
head(g_forms)

### add growth_forms to dataset 
plant_data$GrowthForm <- g_forms$GrowthForm[match(plant_data$Labels, g_forms$Labels)] 

head(plant_data)

### change name of column to reflect that these are the compadre organism types
names(plant_data)[6] <- "OrganismType_Compadre_Database"

write.csv(plant_data, "data_with_needed_env_vars_for_GLMM_analysis.csv", row.names = TRUE)
 
# fin...
