library('vegan')
library('tidyverse')

setwd("C:/Users/15404/Documents/GitHub/research_project/Final_Processing")

getwd()

leaf_data <- read.csv("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/Dataset_Creation/Datasets/leafGenus_plantDNA_leafNutrients.csv", header = TRUE, sep = ",")

###

taxa <- (colnames(leaf_data))

taxa <-taxa[20:length(taxa)]

taxa[0:5]

leaf_data[taxa]

###

nutrients <- (colnames(leaf_data))

nutrients <-nutrients[7:19]

nutrients

###

dim(leaf_data)

leaf_data <- leaf_data %>% drop_na(nutrients)
leaf_data <- leaf_data %>% replace(is.na(.), 0)

dim(leaf_data)

###



leaf_cs <- capscale(leaf_data[, c(taxa)] ~ P_percent_dm + K_percent_dm + 
                    Ca_percent_dm + Mg_percent_dm + Na_percent_dm + S_percent_dm + 
                    Zn_ppm_dm + Mn_ppm_dm + Fe_ppm_dm + Cu_ppm_dm + B_ppm_dm + Al_ppm_dm,
                    data = leaf_data,
                    distance = "bray")


lengths(leaf_cs)

plot(leaf_cs)

anova(leaf_cs)

print(leaf_cs)

leaf_cs <- capscale(leaf_data[, c(taxa)] ~ N_percent_dm + P_percent_dm + K_percent_dm + 
                      Ca_percent_dm + Mg_percent_dm + Na_percent_dm + S_percent_dm + 
                      Zn_ppm_dm + Mn_ppm_dm + Fe_ppm_dm + Cu_ppm_dm + B_ppm_dm + Al_ppm_dm +
                      V2 + site_id + BGR,
                    data = leaf_data,
                    distance = "jaccard")

lengths(leaf_cs)

plot(leaf_cs)

anova(leaf_cs)

print(leaf_cs)
