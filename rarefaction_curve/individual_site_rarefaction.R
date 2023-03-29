require(vegan)

setwd("C:/Users/15404/Documents/GitHub/research_project/rarefaction_curve/")

df <- read.csv("202109_phylum_with_site.csv", header = TRUE, sep = ",")

df1 <- subset(df, site == 'KCKT')

#df1 <- df

df2 <- df1
df2[df2 == 0] <- 25
min_df =  min(df2[][5:ncol(df2)])
min_df

df1[][5:ncol(df1)] = df1[][5:ncol(df1)] / min_df
df1[][5:ncol(df1)] = trunc(df1[][5:ncol(df1)])

print(df1[][5:7])

raremax <- min(rowSums(df1[][5:ncol(df1)]))

S <- specnumber(df1)
Srare <- rarefy(df1[][5:ncol(df1)], raremax)

# Plot the rarefaction curve with colored lines
rarecurve(df1[][5:ncol(df1)], step = 20, sample = raremax, col = 'blue', lwd = 2, cex = 0.5)

          