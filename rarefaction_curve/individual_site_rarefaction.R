require(vegan)

setwd("C:/Users/15404/Documents/GitHub/research_project/rarefaction_curve/")

df <- read.csv("full_202109_genera_with_site.csv", header = TRUE, sep = ",")

unique(df$site)

#df1 <- df
df1 <- subset(df, site == 'KCKT')

numsites = (dim(unique(df1[c('site')]))[1])
print(numsites)

coldf = df[c('site')]

coldf[coldf == "WART"] <- 'darkred'
#coldf[coldf == "SBR"] <- 'forestgreen'
#coldf[coldf == "BTU"] <- 'yellow'
#coldf[coldf == "SKY"] <- 'blue'
#coldf[coldf == "BLDTTT"] <- 'red'
#coldf[coldf == "KCKT"] <- 'orange'
  
# Convert the color list to a vector
collist <- as.vector(unlist(coldf))

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
rarecurve(df1[][5:ncol(df1)], step = 20, sample = raremax, col = 'orange', lwd = 2, cex = 0.5)

          