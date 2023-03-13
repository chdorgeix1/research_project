require(Rarefy)
require(ade4)
require(adiv)
require(ape)
require(vegan)
require(phyloregion)
require(raster)

setwd("C:/Users/15404/Documents/GitHub/research_project/rarefaction_curve")

df <- read.csv("rarefaction_data.csv", header = TRUE, sep = ",")

ncol(df)
df1 <- df
df1[df1 == 0] <- 25
head(df)
head(df1)
min_df =  min(df1[][4:ncol(df1)])
min_df
  
df[][4:ncol(df)]
df[][4:ncol(df)] = df[][4:ncol(df)] /(min_df)
#df[][4:ncol(df)] = df[][4:ncol(df)]
print(df)

df[][4:ncol(df)] = trunc(df[][4:ncol(df)])
print(df)
raremax <- min(rowSums(df[][4:ncol(df)]))
#print(raremax)



S <- specnumber(df)
print(S)
Srare <- rarefy(df[][4:ncol(df)], raremax)


#plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")

rarecurve(df[][4:ncol(df)], step = 20, sample = raremax, col = "blue", cex = 0.6)

#YAYYYY ^^^





data("duneFVG")

data("duneFVG.xy") #plots geographic coordinates


dist_sp<-dist(duneFVG.xy$tot.xy)

ser_rarefaction<-directionalSAC(duneFVG$total,dist_sp)

plot(1:128,ser_rarefaction$N_Exact,xlab="M",ylab="Species richness",ylim=c(0,71),pch=1)
