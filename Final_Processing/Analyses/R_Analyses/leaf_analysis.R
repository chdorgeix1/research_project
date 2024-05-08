library('vegan')
library('tidyverse')

setwd("C:/Users/15404/Documents/GitHub/research_project/Final_Processing")

getwd()

leaf_data <- read.csv("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/Dataset_Creation/Datasets/leafGenus_plantDNA.csv", header = TRUE)

x <- (colnames(leaf_data))
x <-x[7:length(x)]

x[0:10]

colnames(leaf_data[0:10])

taxa_data = leaf_data[, c(x)]

v2_list <- leaf_data$V2

v2_list <- round(v2_list)

leaf_data$genetic_spec <- v2_list

d = vegdist(taxa_data, method="bray", na.rm = TRUE)

leaf_PCoA <- wcmdscale(d, eig = TRUE)





pdf(file="./Figures/leaf_microbiome_shannon_index.pdf", width = 7, height = 4)

install.packages('gridExtra')
library(gridExtra)

plot1 <- ggplot(data = data.frame(leaf_PCoA$points),
                aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = leaf_data$genetic_spec)) +
  theme_bw() + guide_legend(title='nope')
plot1
plot1 + guides(fill=guide_legend(title='fruit'))

plot1
dev.off()

plot2 <- ggplot(data = data.frame(leaf_PCoA$points),
                aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = leaf_data$site_id)) +
  theme_bw()
plot3 <- ggplot(data = data.frame(leaf_PCoA$points),
       aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = leaf_data$BGR)) +
  theme_bw()
grid.arrange(plot1, plot2, plot3, ncol=3)


par(mfrow=c(1,3), mar = c(6, 3.5, 7, 0), cex.axis = 0.5, cex.lab = 1.5)


dev.off()






d = vegdist(taxa_data, method="jaccard", na.rm = TRUE)

leaf_PCoA <- wcmdscale(d, eig = TRUE)

ggplot(data = data.frame(leaf_PCoA$points),
       aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = leaf_data$genetic_spec)) +
  theme_bw()

ggplot(data = data.frame(leaf_PCoA$points),
       aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = leaf_data$site_id)) +
  theme_bw()

ggplot(data = data.frame(leaf_PCoA$points),
       aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = leaf_data$BGR)) +
  theme_bw()
