library('tidyverse')
library('vegan')
library('tabula')
library('gridExtra')
library('devtools')
library('pairwiseAdonis')

setwd("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/")

getwd()

soil_data <- read.csv("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/Dataset_Creation/Datasets/soilGenus_plantDNA.csv", header = TRUE, sep = ",")

colnames(soil_data[0:25])

soil_taxa <- soil_data[7:length(soil_data)]


v2_list <- soil_data$V2

v2_list <- round(v2_list)

soil_data$genetic_spec <- v2_list

soil_data['genetic_spec'][soil_data['genetic_spec'] == 0] <- 'Exaltata'
soil_data['genetic_spec'][soil_data['genetic_spec'] == 1] <- 'Syriaca'

soil_data$binned_sp <- soil_data$V2

soil_data['binned_sp'][soil_data['binned_sp'] < 0.33] <- 0
soil_data['binned_sp'][soil_data['binned_sp'] > 0.66] <- 1
soil_data['binned_sp'][soil_data['binned_sp'] >= 0.33 & soil_data['binned_sp'] <= 0.66] <- 2

soil_data['binned_sp'][soil_data['binned_sp'] == 0] <- 'Exaltata'
soil_data['binned_sp'][soil_data['binned_sp'] == 1] <- 'Syriaca'
soil_data['binned_sp'][soil_data['binned_sp'] == 2] <- 'Hybrid'

soil_data['BGR'][soil_data['BGR'] == 'HR'] <- 'Cole Mountain'



d_1 = vegdist(soil_taxa, method="bray", na.rm = TRUE)

soil_PCoA <- wcmdscale(d_1, eig = TRUE)

pdf(file="./Figures/soil_micro_rel_abundance_binned_sp_PCoA.pdf", width = 7, height = 4)

p <- ggplot(data = data.frame(soil_PCoA$points),
            aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = soil_data$binned_sp)) + 
  labs(color='Plant Host Species') +
  labs(title="Relative Abundance of soil Microbes by Plant Species PCoA") +
  scale_color_manual(values = c('darkgreen', 'green', 'yellow'))
p

dev.off()

pdf(file="./Figures/soil_micro_rel_abundance_plant_spec_PCoA.pdf", width = 7, height = 4)

p <- ggplot(data = data.frame(soil_PCoA$points),
            aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = soil_data$genetic_spec)) + 
  labs(color='Plant Genotype') +
  labs(title="Relative Abundance of Phyllosphere Microbes by Plant Genotype PCoA") +
  scale_color_manual(values = c('darkgreen', 'yellow'))
p


dev.off()

pdf(file="./Figures/soil_micro_rel_abundance_site_PCoA.pdf", width = 7, height = 4)

p <- ggplot(data = data.frame(soil_PCoA$points),
            aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = soil_data$site_id)) + 
  labs(color='Plant Host Species') +
  labs(title="Relative Abundance of soil Microbes by Site PCoA")
p


dev.off()

pdf(file="./Figures/soil_micro_rel_abundance_bgr_PCoA.pdf", width = 7, height = 4)

p <- ggplot(data = data.frame(soil_PCoA$points),
            aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = soil_data$BGR)) + 
  labs(color='Broad Geographic Location') +
  labs(title="Relative Abundance of soil Microbes by Broad Geographic Region PCoA")
p


dev.off()


p_nova <- adonis2(d_1 ~ soil_data$binned_sp, method = 'bray')

p_nova

p_nova <- adonis2(d_1 ~ soil_data$genetic_spec, method = 'bray')

p_nova

p_nova <- adonis2(d_1 ~ soil_data$site_id, method = 'bray')

p_nova

p_nova <- adonis2(d_1 ~ soil_data$BGR, method = 'bray')

p_nova


adonis <- adonis(ord ~ X * Y, data = your.metadata) #this is basically your #PERMANOVA analysis

mod <- betadisper(d_1, soil_data$BGR)
permutest(mod)

plot(mod)

boxplot(mod)

mod.HSD <- TukeyHSD(mod)
mod.HSD
plot(mod.HSD)

plot2 <- ggplot(data = data.frame(soil_PCoA$points),
                aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = soil_data$site_id))
plot2


plot3 <- ggplot(data = data.frame(soil_PCoA$points),
                aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = soil_data$BGR))
plot3


plot4 <- ggplot(data = data.frame(soil_PCoA$points),
                aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = soil_data$V2))
plot4