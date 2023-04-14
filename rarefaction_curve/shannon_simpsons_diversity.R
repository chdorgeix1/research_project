library(vegan)
library(ggplot2)

setwd("C:/Users/15404/Documents/GitHub/research_project/rarefaction_curve/")

df <- read.csv("full_202109_genera_with_site.csv", header = TRUE, sep = ",")

head(df[][5:ncol(df)])
df1 = df[][5:ncol(df)]

# calculate diversity indices by site
shannon_diversity <- diversity(df1, index = "shannon")
simpson_diversity <- diversity(df1, index = "simpson")

df$shannon <- shannon_diversity
df$simpson <- simpson_diversity

boxplot(df$shannon ~ df$site, df,                                # Change main title and axis labels
        main = "Shannon Diversity by Site",
        xlab = "Site",
        ylab = "Shannon Diversity")

boxplot(df$simpson ~ df$site, df,                                # Change main title and axis labels
        main = "Simpson Diversity by Site",
        xlab = "Site",
        ylab = "Simpson Diversity")
