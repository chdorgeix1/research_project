library('vegan')
library('tidyverse')

setwd("C:/Users/15404/Documents/GitHub/research_project/Final_Processing")

getwd()

leaf_data <- read.csv("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/Dataset_Creation/Datasets/leafGenus_plantDNA_leafNutrients.csv", header = TRUE, sep = ",")
soil_data <- read.csv("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/Dataset_Creation/Datasets/soilGenus_plantDNA_soilNutrients.csv", header = TRUE, sep = ",")

#####

taxa_l <- (colnames(leaf_data))

taxa_l <-taxa_l[20:length(taxa_l)]

taxa_l[0:10]


nutrients_l <- (colnames(leaf_data))

nutrients_l <-nutrients_l[7:19]

nutrients_l

#####

taxa_s <- (colnames(soil_data))

taxa_s <-taxa_s[18:length(taxa_s)]

nutrients_s <- (colnames(soil_data))

nutrients_s <-nutrients_s[7:17]

nutrients_s

#####

df_list <- list(leaf_data, soil_data)

df <- df_list %>% reduce(full_join, by='sample_id')

dim(leaf_data)

dim(soil_data)

dim(df[taxa_s])
dim(df[nutrients_s])
dim(df[nutrients_l])

predictor_cols <- c(names(df[taxa_s]), names(df[nutrients_s]), names(df[nutrients_l]))


formula_str <- paste(". ~ ", paste(predictor_cols, collapse = " + "))

formula_str

taxa_l[(length(taxa_l)-10): length(taxa_l)]

leaf_cs <- capscale(df[, c(taxa_l)] ~ formula_str,
                    data = df,
                    distance = "jaccard")


leaf_cs <- capscale(df[, c(taxa_l)] ~ Phosphorus + Potassium + Calcium + 
                      Magnesium + Sulfur + Sodium + Zinc + Manganese + Iron + 
                      Copper + Boron,
                    data = soil_data,
                    distance = "bray")
