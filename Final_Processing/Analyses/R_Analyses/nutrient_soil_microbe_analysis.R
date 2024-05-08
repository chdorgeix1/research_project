library('vegan')
library('tidyverse')

setwd("C:/Users/15404/Documents/GitHub/research_project/Final_Processing")

getwd()

#soil_data = read.csv("micro_nutrients.csv", header = TRUE)
soil_data <- read.csv("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/Dataset_Creation/Datasets/soilGenus_plantDNA_soilNutrients.csv", header = TRUE, sep = ",")

### Subset Taxa Names

taxa <- (colnames(soil_data))

taxa <-taxa[17:length(taxa)]

### Subset Soil Nutrient Names

nutrients <- (colnames(soil_data))

nutrients <-nutrients[6:16]

###

soil_cs <- capscale(soil_data[, c(taxa)] ~ Phosphorus + Potassium + Calcium + 
                      Magnesium + Sulfur + Sodium + Zinc + Manganese + Iron + 
                      Copper + Boron,
                    data = soil_data,
                    distance = "bray")

plot(soil_cs)

anova(soil_cs)

print(soil_cs)

###

soil_cs <- capscale(soil_data[, c(taxa)] ~ Phosphorus + Potassium + Calcium + 
                      Magnesium + Sulfur + Sodium + Zinc + Manganese + Iron + 
                      Copper + Boron + V2,
                    data = soil_data,
                    distance = "bray")

plot(soil_cs)
anova(soil_cs)

print(soil_cs)

anova_result <- anova(soil_cs)

# Print ANOVA table
print(anova_result)

soil_cs <- capscale(soil_data[, c(taxa)] ~ Phosphorus + Potassium + Calcium + 
                      Magnesium + Sulfur + Sodium + Zinc + Manganese + Iron + 
                      Copper + Boron,
                    data = soil_data,
                    distance = "jaccard")

plot(soil_cs)

anova(soil_cs)

print(soil_cs)

###

soil_cs <- capscale(soil_data[, c(taxa)] ~ Phosphorus + Potassium + Calcium + 
                      Magnesium + Sulfur + Sodium + Zinc + Manganese + Iron + 
                      Copper + Boron + V2,
                    data = soil_data,
                    distance = "jaccard")
soil_cs

anova(soil_cs, by='term')

print(soil_cs)



# Total constrained inertia
total_constrained_inertia <- 2.732931

loadings <- scores(soil_cs)$species

# Calculate proportion of variance explained by each variable
proportions <- (loadings^2) / total_constrained_inertia

# Summarize the results
variables <- c("Phosphorus", "Potassium", "Calcium", "Magnesium", "Sulfur", 
               "Sodium", "Zinc", "Manganese", "Iron", "Copper", "Boron", "V2")

for (i in 1:length(variables)) {
  cat(paste(variables[i], ": ", proportions[i] * 100, "%\n"))
}



soil_cs$colsum

plot(soil_cs, col = as.factor(soil_data$site_id), pch=19)

ordiellipse(soil_cs, groups = soil_data$BGR)

print(soil_cs)

