install.packages('Rarefy')
install.packages('ade4')
install.packages('adiv')
install.packages('ape')
install.packages('vegan')
install.packages('phyloregion')
install.packages('raster')

library('Rarefy')
library('Rarefy')
library('Rarefy')
library('Rarefy')
library('Rarefy')
library('Rarefy')
library('Rarefy')
library('Rarefy')

require(Rarefy)
require(ade4)
require(adiv)
require(ape)
require(vegan)
require(phyloregion)
require(raster)

Samples <- c("CMB1")
Strain_1 <- c(0.1)
Strain_2 <- c(2)
Strain_3 <- c(0.5)
Strain_4 <- c(12)
Strain_5 <- c(0.2)

df <- data.frame(Samples, Strain_1, Strain_2, Strain_3, Strain_4, Strain_5)

print (df)
print (df[][2:5])

print(df[][2:5] /(min(df[][2:5])))
df[][2:6] = df[][2:6] /(min(df[][2:6]))
print(df)

raremax <- min(rowSums(df[][2:6]))
print(raremax)

S <- 5
Srare <- rarefy(df[][2:6], raremax)


#Srare <- rarefy(BCI, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")

rarecurve(df[][2:6], step = 20, sample = raremax, col = "blue", cex = 0.6)

#YAYYYY ^^^





data("duneFVG")

data("duneFVG.xy") #plots geographic coordinates


dist_sp<-dist(duneFVG.xy$tot.xy)

ser_rarefaction<-directionalSAC(duneFVG$total,dist_sp)

plot(1:128,ser_rarefaction$N_Exact,xlab="M",ylab="Species richness",ylim=c(0,71),pch=1)

