library('vegan')
library('tidyverse')
library('geosphere')

setwd("C:/Users/15404/Documents/GitHub/research_project/Final_Processing")

getwd()

leaf_data <- read.csv("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/Dataset_Creation/Datasets/leafGenus_plantDNA_leafNutrients.csv", header = TRUE, sep = ",")
soil_data <- read.csv("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/Dataset_Creation/Datasets/soilGenus_plantDNA_soilNutrients.csv", header = TRUE, sep = ",")

df_list <- list(leaf_data, soil_data)

df <- df_list %>% reduce(right_join, by='sample_id')

dim(df)

df <- subset(df, select = -c(BGR.x, site_id.x, V2.x, X.x, species.x, X.y)) #Na_percent_dm


dim(df)

df <- df[-c(68, 69), ]  # MMP 16/17 do not have leaf nutrient analysis
df <- df %>% replace(is.na(.), 0)
df <- df[, colSums(df != 0) > 0]
dim(df)
df <- df[vapply(df, function(x) length(unique(x)) > 1, logical(1L))]
dim(df)


#####

leaf_nutrients <- (colnames(df))[2:13]

soil_nutrients <- (colnames(df))[1206:1216]

leaf_taxa <- (colnames(df))[14:1201]

soil_taxa <- (colnames(df))[1218:length(df)-1]

#####

soil_nutrients #named elements

leaf_nutrients #in percents

head(leaf_taxa)
tail(leaf_taxa)

head(soil_taxa)
tail(soil_taxa)

df$latitude <- df$site_id.y
df$longitude <- df$site_id.y

df['latitude'][df['latitude'] == 'CMB'] <- 37.755588
df['latitude'][df['latitude'] == 'LFS'] <- 37.772967
df['latitude'][df['latitude'] == 'LM'] <- 38.162222
df['latitude'][df['latitude'] == 'MKP'] <- 37.591111
df['latitude'][df['latitude'] == 'MMP'] <- 37.766821

df['longitude'][df['longitude'] == 'CMB'] <- -79.19706
df['longitude'][df['longitude'] == 'LFS'] <- -79.184033
df['longitude'][df['longitude'] == 'LM'] <- -79.076388
df['longitude'][df['longitude'] == 'MKP'] <- -80.193055
df['longitude'][df['longitude'] == 'MMP'] <- -79.187052


df$latitude <- as.numeric(df$latitude)
df$longitude <- as.numeric(df$longitude)

geo_df <- df[ , c('longitude', 'latitude')]

geo_matrix <- data.matrix(geo_df)

res <- distm(geo_matrix, fun=distGeo )

df$dist_matrix <- res

##### Example 1 Phyllosphere Microbiome Taxa Capscale 

predictor_cols <- c(leaf_nutrients, soil_nutrients, 'dist_matrix', 'V2.y') #Site, Broad Geographic Location, Soil Nutrients, and Plant Genotype

fmla <- as.formula(paste("df[, c(leaf_taxa)] ~ ", paste(predictor_cols, collapse = "+")))

leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "bray")

leaf_cs

anova(leaf_cs, by='terms')


scores(leaf_cs)['sites']

score_val <- scores(leaf_cs)['sites']

for (cap in scores(leaf_cs)['sites']){
  print(cap[1])
}


class(score_val)




predictor_cols <- c(soil_nutrients, 'latitude', 'longitude', 'V2.y') #Site, Broad Geographic Location, Soil Nutrients, and Plant Genotype

fmla <- as.formula(paste("df[, c(soil_taxa)] ~ ", paste(predictor_cols, collapse = "+")))

soil_cs <- capscale(fmla,
                    data = df,
                    distance = "bray",
                    add = TRUE)

soil_cs

anova(soil_cs, by='terms')


plot(soil_cs)

soil_cs$CA

soil_dist <- vegdist(df[, c(soil_taxa)], method = 'bray')



vp <- varpart(soil_dist, df[, c(leaf_nutrients)], df$dist_matrix)

vp

plot(vp)

test_info <- dbrda(fmla,
      data=df,
      distance = 'bray')
test_info

plot(test_info)

##### Example 2 Leaf Microbiome Taxa Capscale 


predictor_cols <- c(soil_nutrients, leaf_nutrients, 'dist_matrix', 'V2.y') #Site, Broad Geographic Location, Soil Nutrients, and Plant Genotype

fmla <- as.formula(paste("df[, c(leaf_taxa)] ~ ", paste(predictor_cols, collapse = "+")))

leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "bray")

leaf_cs

anova(leaf_cs, by='terms')

plot(leaf_cs)

##### Example 3 Leaf Microbiome Taxa x Soil Microbiome Capscale 

predictor_cols <- c(soil_taxa) #Site, Broad Geographic Location, Soil Nutrients, and Plant Genotype

fmla <- as.formula(paste("df[, c(leaf_taxa)] ~ ", paste(predictor_cols, collapse = "+")))

leaf_soil_cs <- capscale(fmla,
                    data = df,
                    distance = "bray")

leaf_soil_cs






soil_dist <- vegdist(df[, c(soil_taxa)], method = 'bray')

vp <- varpart(soil_dist, df[, c(soil_nutrients)], df$V2.y)

vp

plot(vp)








predictor_cols <- c(leaf_nutrients) #Site, Broad Geographic Location, Soil Nutrients, and Plant Genotype

fmla <- as.formula(paste("df[, c(leaf_taxa)] ~ ", paste(predictor_cols, collapse = "+")))
#fmla <- as.formula(paste("soil_dist ~ ", paste(predictor_cols, collapse = "+")))


soil_cs <- capscale(fmla,
                    data = df,
                    distance = "bray")


soil_cs

leaf_dist <- vegdist(df[, c(leaf_taxa)], method = 'bray')


vp <- varpart(leaf_dist, df[, c(soil_nutrients)], df[, c(soil_taxa)])

vp


plot(vp)

anova(soil_cs, by='terms')











plot(soil_cs)

plot(soil_cs, col = df$site_id.y)

plot(soil_cs)

anova(soil_cs)

scores_soil=scores(soil_cs)

site_scores = scores_soil$sites

sp_scores = scores_soil$species

plot(site_scores)

cbind(site_scores)

spec <- soil_cs$CCA$v
sites <- mod$CCA$u
eig <- mod$CCA$eig
rsum <- soil_cs$rowsum
csum <- mod$colsum
rsum

anova(soil_cs, by = 'onedf')

pdf(file="./Figures/soil_microbiome_PDBRA.pdf", width = 12, height = 12)
plot(soil_cs, main = 'Soil Microbiome PDBRA')
dev.off()
#####

predictor_cols <- c(leaf_nutrients, soil_nutrients, 'site_id.y', 'BGR.y', 'V2.y') #Site, Broad Geographic Location, Soil Nutrients, and Plant Genotype

fmla <- as.formula(paste("df[, c(soil_taxa)] ~ ", paste(predictor_cols, collapse = "+")))


leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "bray")

leaf_cs

plot(leaf_cs)


predictor_cols <- c(soil_nutrients)

fmla <- as.formula(paste("df[, c(leaf_taxa)] ~ ", paste(predictor_cols, collapse = "+")))


leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "bray")

leaf_cs


pdf(file="./Figures/leaf_microbiome_PDBRA.pdf", width = 12, height = 12)
plot(leaf_cs, main = 'Leaf Microbiome PDBRA')
dev.off()

anova(leaf_cs, permutations = 999)

predictor_cols <- c(leaf_nutrients)

fmla <- as.formula(paste("df[, c(leaf_taxa)] ~ ", paste(predictor_cols, collapse = "+")))


leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "bray")

leaf_cs

predictor_cols <- c(soil_nutrients)

fmla <- as.formula(paste("df[, c(leaf_taxa)] ~ ", paste(predictor_cols, collapse = "+")))


leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "bray")

leaf_cs

predictor_cols <- c('V2.y')

fmla <- as.formula(paste("df[, c(leaf_taxa)] ~ ", paste(predictor_cols, collapse = "+")))


leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "bray")

leaf_cs


############## Using Soil Microbes to Predict Leaf Microbes

predictor_cols <- soil_taxa #Site, Broad Geographic Location, Soil Nutrients, and Plant Genotype

fmla <- as.formula(paste("df[, c(leaf_taxa)] ~ ", paste(predictor_cols, collapse = "+")))

length(fmla)

leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "bray")


leaf_cs <- capscale(df[, c(leaf_taxa)] ~ soil_taxa,
                    data = df,
                    distance = "bray")

leaf_cs


leaf_cs

plot(leaf_cs)

















leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "jaccard")

leaf_cs

df[leaf_taxa] <- ceiling(df[leaf_taxa])

predictor_cols <- c(leaf_nutrients, soil_nutrients, 'V2.y')

fmla <- as.formula(paste("df[, c(leaf_taxa)] ~ ", paste(predictor_cols, collapse = "+")))

leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "bray")

leaf_cs

############

predictor_cols <- c('V2.y')


fmla <- as.formula(paste("df[, c(leaf_taxa)] ~ ", paste(predictor_cols, collapse = "+")))

fmla

leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "bray")

leaf_cs


leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "jaccard")

leaf_cs

##############

predictor_cols <- c(leaf_nutrients)


fmla <- as.formula(paste("df[, c(leaf_taxa)] ~ ", paste(predictor_cols, collapse = "+")))

fmla

leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "bray")

leaf_cs


leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "jaccard")

leaf_cs

plot(leaf_cs)

####################


predictor_cols <- c(soil_nutrients)


fmla <- as.formula(paste("df[, c(leaf_taxa)] ~ ", paste(predictor_cols, collapse = "+")))

fmla

leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "bray")

leaf_cs


leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "jaccard")

leaf_cs

plot(leaf_cs)


#################


predictor_cols <- c('V2.y')

df$V2.y

fmla <- as.formula(paste("df[, c(leaf_taxa)] ~ ", paste(predictor_cols, collapse = "+")))

fmla

leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "bray")

leaf_cs


leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "jaccard")

leaf_cs

plot(leaf_cs)















results <- prcomp(df[soil_taxa], scale = FALSE)

results$x

soil_pca <- results$x[,1]

df$soil_pca <- soil_pca


predictor_cols <- c('soil_pca')

fmla <- as.formula(paste("df[, c(leaf_taxa)] ~ ", paste(predictor_cols, collapse = "+")))

fmla

leaf_cs <- capscale(fmla,
                    data = df,
                    distance = "bray")

leaf_cs



taxa_data = df[, c(soil_taxa)]

d = vegdist(taxa_data, method="bray", na.rm = TRUE)

soil_PCoA <- wcmdscale(d, eig = TRUE)

ggplot(data = data.frame(soil_PCoA$points),
       aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = df$BGR)) +
  theme_bw()

ggplot(data = data.frame(soil_PCoA$points),
       aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = df$species.y)) +
  theme_bw()

ggplot(data = data.frame(soil_PCoA$points),
       aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = df$site_id.y)) +
  theme_bw()




