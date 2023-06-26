# Load required libraries
library(vegan)
library(ggplot2)

setwd("C:/Users/15404/Documents/GitHub/research_project/rarefaction_curve/")

# Read data from a file or create a data frame with relative abundance data
data <- read.csv("soil_microbe_genus.csv", header = TRUE, sep = ",")
# Load the iris dataset
#data <- subset(data, site == 'FRW' | site == 'RRL')
data <- data[, colSums(data != 0) > 0]


# Perform a PCA on the dataset
pca <- prcomp(data[][4:ncol(data)], center = TRUE, scale. = TRUE)

summary(pca)

# Extract the first two principal components
pc1 <- pca$x[,1]
pc2 <- pca$x[,2]
pc3 <- pca$x[,3]

unique(data$site)

# Create a color vector based on the species of iris
site <- data$site
coldf = data[c('site')]

coldf[coldf == "RRL"] <- 'blue'
coldf[coldf == "FRW"] <- 'blue'
coldf[coldf == "MMP"] <- 'green'
coldf[coldf == "PNR"] <- 'blue'
coldf[coldf == "CMB"] <- 'green'
coldf[coldf == "RF"] <- 'blue'
coldf[coldf == "SLG"] <- 'blue'
coldf[coldf == "LFS"] <- 'green'
coldf[coldf == "SGC"] <- 'green'
coldf[coldf == "RGT"] <- 'green'
coldf[coldf == "HR"] <- 'blue'
coldf[coldf == "PTW"] <- 'green'
coldf[coldf == "GH"] <- 'green'
coldf[coldf == "LM"] <- 'green'
coldf[coldf == "MKP"] <- 'green'
coldf[coldf == "TCP"] <- 'green'
 
coldf[coldf == "WART"] <- 'red'
coldf[coldf == "SBR"] <- 'red'
coldf[coldf == "BTU"] <- 'red'
coldf[coldf == "SKY"] <- 'red'
coldf[coldf == "BLDTTT"] <- 'red'
            
coldf[coldf == "KCKT"] <- 'red'





# Plot the PCA with colored points
#, col = coldf$site
#pch = 19,

#unique(data$site)

plot(pc1, pc2, cex =1,pch = 19, xlim = c(-45, 40), ylim = c(-50,40), col = coldf$site,
     xlab = "PC1(15.2%)", ylab = "PC2(4%)", main = "PCA plot of Soil Microbe Sites by Genus")
#legend("topright", legend = levels(unique(data$site)), col = coldf$site, 
#       pch = 19, bty = "n")
#, 
fill = c('darkred','forestgreen','yellow','blue','red','orange')
legend(8, 10, legend=c(unique(data$site)), cex = 1)

              
       