library('tidyverse')
library('vegan')
library('tabula')


setwd("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/")

getwd()

soil_data <- read.csv("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/Dataset_Creation/Datasets/soilGenus_plantDNA.csv", header = TRUE, sep = ",")


taxa <- colnames(soil_data)

taxa[0:20]

taxa <- taxa[7:(length(taxa)-2)]

taxa[(length(taxa)-5): length(taxa)]

soil_data[taxa] <- mutate_if(soil_data[taxa], is.numeric, ~ . * 10000)

soil_data[taxa] <- mutate_if(soil_data[taxa], is.numeric, round)

#df <- data.matrix(soil_data[][taxa])
df <- soil_data[][taxa]

df <- df %>% replace(is.na(.), 0)

head(df[0:5])

#############


H <- diversity(df)

simp <- diversity(df, index = "simpson")

shannon <- diversity(df, index = 'shannon')

invsimp <- diversity(df, "inv")

unbias.simp <- simpson.unb(df)

alpha <- fisher.alpha(df)


#pairs(cbind(H, simp, invsimp, unbias.simp, alpha), col="blue")


###################

S <- specnumber(df) ## rowSums(BCI > 0) does the same...
J <- H/log(S)

S
J


###################

soil_data$alpha <- alpha

soil_data$simp <- simp

soil_data$shannon <- shannon


##################

v2_list <- soil_data$V2

v2_list <- round(v2_list)

soil_data$genetic_spec <- v2_list

soil_data['genetic_spec'][soil_data['genetic_spec'] == 0] <- 'Exaltata'
soil_data['genetic_spec'][soil_data['genetic_spec'] == 1] <- 'Syriaca'

##################

dev.off()
plot(soil_data$V2, soil_data$alpha)


jpeg(file="./Figures/soil_alpha_div_BGR.jpeg")
boxplot(alpha ~ BGR, data = soil_data, main = 'Soil Microbiome Alpha Diversity by BGR',
        xlab = 'BGR', ylab = 'Alpha Diversity', col = rainbow(10))
dev.off()

jpeg(file="./Figures/soil_shannon_div_BGR.jpeg")
boxplot(shannon ~ BGR, data = soil_data, main = 'Soil Microbiome Shannon Index by BGR',
        xlab = 'BGR', ylab = 'Shannon Index', col = rainbow(10))
dev.off()


jpeg(file="./Figures/soil_simp_div_BGR.jpeg")
boxplot(simp ~ BGR, data = soil_data, main = 'Soil Microbiome Simpson Index by BGR',
        xlab = 'BGR', ylab = 'Simpson Index', col = rainbow(10))
dev.off()



jpeg(file="./Figures/soil_alpha_div_site.jpeg")
boxplot(alpha ~ site_id, data = soil_data, main = 'Soil Microbiome Alpha Diversity by Site',
        xlab = 'Site', ylab = 'Alpha Diversity', col = rainbow(10))
dev.off()

jpeg(file="./Figures/soil_shannon_div_site.jpeg")
boxplot(shannon ~ site_id, data = soil_data, main = 'Soil Microbiome Shannon Index by Site',
        xlab = 'Site', ylab = 'Shannon Index', col = rainbow(10))
dev.off()

jpeg(file="./Figures/soil_simp_div_site.jpeg")
boxplot(simp ~ site_id, data = soil_data, xlab = 'Site', ylab = 'Simpson Index', col = rainbow(10))
dev.off()



jpeg(file="./Figures/soil_alpha_div_species.jpeg")
boxplot(alpha ~ genetic_spec, data = soil_data, xlab = 'Plant Species', main = 'Soil Microbiome Alpha Diversity by Plant Host Species',
        ylab = 'Alpha Diversity', col = c('darkgreen', 'yellow'))
dev.off()

jpeg(file="./Figures/soil_shannon_div_species.jpeg")
boxplot(shannon ~ genetic_spec, data = soil_data, main = 'Soil Microbiome Shannon Index by Plant Host Species',
        xlab = 'Plant Species', ylab = 'Shannon Index', col = c('darkgreen', 'yellow'))
dev.off()

jpeg(file="./Figures/soil_simpson_div_species.jpeg")
boxplot(simp ~ genetic_spec, data = soil_data, xlab = 'Plant Species', ylab = 'Simpson Index', col = c('darkgreen', 'yellow'))
dev.off()


t.test(alpha ~ genetic_spec, data = soil_data)

t.test(simp ~ genetic_spec, data = soil_data)

t.test(shannon ~ genetic_spec, data = soil_data)


one.way <- aov(alpha ~ site_id, data = soil_data)

summary(one.way)

one.way <- aov(simp ~ site_id, data = soil_data)

summary(one.way)

one.way <- aov(shannon ~ site_id, data = soil_data)

summary(one.way)



one.way <- aov(alpha ~ BGR, data = soil_data)

summary(one.way)

one.way <- aov(simp ~ BGR, data = soil_data)

summary(one.way)

one.way <- aov(shannon ~ BGR, data = soil_data)

summary(one.way)

# 3 figures arranged in 3 rows and 1 column
dev.off()

pdf(file="./Figures/soil_microbiome_alpha_div.pdf", width = 7, height = 5)

par(mfrow=c(1,3), mar = c(6, 3, 7, 0.5), cex.axis = 0.5, cex.lab = 1.5)
boxplot(alpha ~ site_id, data = soil_data, #main = 'Soil Microbiome Alpha Diversity by Site',
        xlab = 'Site', ylab = '', col = rainbow(10))
boxplot(alpha ~ BGR, data = soil_data, #main = 'Soil Microbiome Alpha Diversity by BGR',
        xlab = 'BGR', ylab = '', col = rainbow(10))
boxplot(alpha ~ genetic_spec, data = soil_data, xlab = 'Plant Species', #main = 'Soil Microbiome Alpha Diversity by Plant Host Species',
        ylab = '', col = c('darkgreen', 'yellow'))
mtext("Soil Microbiome Alpha Diversity is Influenced by Site and BGR",   # Add main title
      side = 3,
      line = -5,
      outer = TRUE)

mtext("Alpha Diversity",   # Add main title
      side = 2,
      line = -1.25,
      outer = TRUE)
dev.off()



pdf(file="./Figures/soil_microbiome_shannon_index.pdf", width = 7, height = 4)

par(mfrow=c(1,3), mar = c(6, 3.5, 7, 0.1), cex.axis = 0.8, cex.lab = 1.5)
boxplot(shannon ~ site_id, data = soil_data, #main = 'Soil Microbiome Alpha Diversity by Site',
        xlab = 'Site', ylab = '', col = rainbow(10), ylim=c(4,5.5), las=2)
boxplot(shannon ~ BGR, data = soil_data, #main = 'Soil Microbiome Alpha Diversity by BGR',
        xlab = 'BGR', ylab = '', col = rainbow(10), ylim=c(4,5.5))
boxplot(shannon ~ genetic_spec, data = soil_data, xlab = 'Plant Species', #main = 'Soil Microbiome Alpha Diversity by Plant Host Species',
        ylab = '', col = c('darkgreen', 'yellow'), ylim=c(4,5.5))
mtext("Soil Microbiome Shannon Index is Influenced by Site and BGR",   # Add main title
      side = 3,
      line = -5,
      outer = TRUE)

mtext("Shannon Index",   # Add main title
      side = 2,
      line = -1.25,
      outer = TRUE)
dev.off()

