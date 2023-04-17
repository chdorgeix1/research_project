library(vegan)
library(ggplot2)
library(spatialEco)

setwd("C:/Users/15404/Documents/GitHub/research_project/rarefaction_curve/")

df <- read.csv("202109_phylum_with_site.csv", header = TRUE, sep = ",")

df2 <- subset(df, site == 'RRL')

df2

sum(df[11:20,][5:ncol(df2)])


sum(df2[0:10,][5:ncol(df2)])

df1 = df[][5:ncol(df)]
head(df1)

df1[df1==0] <- NA

df[df==0] <- NA


shannon <- shannons(df[][5:ncol(df)], counts = FALSE, margin = "row")
shannon
df$shannon <- shannon$H
df$evenness <- shannon$evenness


boxplot(df$shannon ~ df$site, df,                                # Change main title and axis labels
        main = "Shannon Diversity of Phylum by Site",
        xlab = "Site",
        ylab = "Shannon Diversity")

boxplot(df$evenness ~ df$site, df,                                # Change main title and axis labels
        main = "Evenness of Phylum by Site",
        xlab = "Site",
        ylab = "Shannon Diversity")

bact1 <- c(0.1, 0.3, 0.5)
bact2 <- c(0.2, 0.1, 0.1)
bact3 <- c(0.1, 0.1, 0)

df2 <- data.frame(bact1, bact2, bact3)
df2

df2[df2==0] <- NA

shannons(df2, count = FALSE, ens = TRUE, margin = 'row')
