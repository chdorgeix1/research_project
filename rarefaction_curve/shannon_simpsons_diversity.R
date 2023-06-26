library(vegan)
library(ggplot2)
library(spatialEco)

setwd("C:/Users/15404/Documents/GitHub/research_project/rarefaction_curve/")

df <- read.csv("plant_microbe_phylum.csv", header = TRUE, sep = ",")

df1 = df[][4:ncol(df)]

head(df1)

head(df)

# calculate diversity indices by site
shannon_diversity <- diversity(df1, index = "shannon")
simpson_diversity <- diversity(df1, index = "simpson")

#plot(c(1,1), shannon_diversity)
#plot(c(2,2), simpson_diversity)


df$shannon <- shannon_diversity
df$simpson <- simpson_diversity

unique(df$Site)

boxplot(df$shannon ~ df$Site, df,                                # Change main title and axis labels
        main = "Shannon Diversity of Plant Microbiome Phylum by Site",
        xlab = "Site",
        ylab = "Shannon Diversity",
        col = c('blue', 'red', 'darkgreen', 'yellow', 'orange', 'purple'),
        ylim = c(0.3, 1.8),
        names=unique(df$Site),
        show.names = TRUE)

boxplot(df$simpson ~ df$site, df,                                # Change main title and axis labels
        main = "Simpson Diversity of Phylum by Site",
        xlab = "Site",
        ylab = "Simpson Diversity",
        col = factor(df$site))

