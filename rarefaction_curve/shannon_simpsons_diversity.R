library(vegan)
library(ggplot2)
library(spatialEco)

setwd("C:/Users/15404/Documents/GitHub/research_project/rarefaction_curve/")

df <- read.csv("202107_phylum_with_site.csv", header = TRUE, sep = ",")

df1 = df[][5:ncol(df)]

head(df1)

# calculate diversity indices by site
shannon_diversity <- diversity(df1, index = "shannon")
simpson_diversity <- diversity(df1, index = "simpson")

#plot(c(1,1), shannon_diversity)
#plot(c(2,2), simpson_diversity)


df$shannon <- shannon_diversity
df$simpson <- simpson_diversity

boxplot(df$shannon ~ df$site, df,                                # Change main title and axis labels
        main = "Shannon Diversity of Phylum by Site",
        xlab = "Site",
        ylab = "Shannon Diversity",
        col = c('blue', 'red', 'darkgreen', 'yellow', 'orange', 'purple'))

boxplot(df$simpson ~ df$site, df,                                # Change main title and axis labels
        main = "Simpson Diversity of Phylum by Site",
        xlab = "Site",
        ylab = "Simpson Diversity",
        col = factor(df$site))

