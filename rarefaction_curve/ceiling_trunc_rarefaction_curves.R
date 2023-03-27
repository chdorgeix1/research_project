require(vegan)


setwd("C:/Users/15404/Documents/GitHub/research_project/rarefaction_curve")

df <- read.csv("rarefaction_data_202109_phylum.csv", header = TRUE, sep = ",")

df1 <- df
dim(df1)
sum(df1 != 0) 

(dim(df1))
ncol(df1)
df2 <- df1
df2[df2 == 0] <- 25
#head(df)
#head(df1)
min_df =  min(df2[][4:ncol(df1)])
min_df

df1[][4:ncol(df1)]
df1[][4:ncol(df1)] = df1[][4:ncol(df1)] / min_df
#df1[][4:ncol(df1)] = df1[][4:ncol(df1)] * 7303252
#df1[][4:ncol(df1)] = (df1[][4:ncol(df1)]) / 876
#df[][4:ncol(df)] = df[][4:ncol(df)]
print(df1[][4:ncol(df1)] )

df1[][4:ncol(df1)] = ceiling(df1[][4:ncol(df1)])
print(df1)


raremax <- min(rowSums(df1[][4:ncol(df1)]))
#print(raremax)



S <- specnumber(df1)
print(S)
Srare <- rarefy(df1[][4:ncol(df1)], raremax)


#plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")

rarecurve(df1[][4:ncol(df1)], step = 20, sample = raremax, col = "blue", cex = 0.6)
