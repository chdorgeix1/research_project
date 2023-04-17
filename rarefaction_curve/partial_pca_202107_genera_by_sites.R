# Load required libraries
library(vegan)
library(ggplot2)

setwd("C:/Users/15404/Documents/GitHub/research_project/rarefaction_curve/")

# Read data from a file or create a data frame with relative abundance data
data <- read.csv("202107_phylum_with_site.csv", header = TRUE, sep = ",")
# Load the iris dataset

# Perform a PCA on the dataset
pca <- prcomp(data[][5:ncol(data)], center = TRUE, scale. = TRUE)

summary(pca)

# Extract the first two principal components
pc1 <- pca$x[,1]
pc2 <- pca$x[,2]

#pc1 <- pc1[pc1<40]
#pc2 <- pc2[pc2<30]
unique(data$site)

# Create a color vector based on the species of iris
site <- data$site
coldf = data[c('site')]

coldf[coldf == "SLG"] <- 'darkred'
coldf[coldf == "PNR"] <- 'forestgreen'
coldf[coldf == "RF"] <- 'yellow'
coldf[coldf == "FRW"] <- 'blue'
coldf[coldf == "RRL"] <- 'red'

coldf

#coldf <- coldf[-c(50), ]

coldf

pc1

# Plot the PCA with colored points

plot(pc1, pc2, col = coldf$site, pch = 19, 
     xlab = "PC1(13.9)", ylab = "PC2(10.3)", main = "PCA plot of Sites by Phylum")
#legend("topright", legend = levels(unique(data$site)), col = coldf$site, 
#       pch = 19, bty = "n")
legend(-6, 6.5, legend=c(unique(data$site)), cex = 1, 
        fill = c('darkred','forestgreen','yellow','blue','red'))

              
       