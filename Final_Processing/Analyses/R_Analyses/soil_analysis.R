library('vegan')
library('tidyverse')

setwd("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/Untitled Folder")

getwd()

soil_data = read.csv("full_soil_microbiome.csv", header = TRUE)

x <- (colnames(soil_data))

x <-x[6:length(x)]
x <-x[0:(length(x)-2)]

x[0:5]
tail(x)

taxa_data = soil_data[, c(x)]

d = vegdist(taxa_data, method="bray", na.rm = TRUE)

soil_PCoA <- wcmdscale(d, eig = TRUE)

ggplot(data = data.frame(soil_PCoA$points),
       aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = soil_data$species)) +
  theme_bw()




d = vegdist(taxa_data, method="jaccard", na.rm = TRUE)

soil_PCoA <- wcmdscale(d, eig = TRUE)

ggplot(data = data.frame(soil_PCoA$points),
       aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = soil_data$species)) +
  theme_bw()

