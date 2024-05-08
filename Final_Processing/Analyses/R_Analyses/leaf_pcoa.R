library('tidyverse')
library('vegan')
library('tabula')
library('gridExtra')
library('devtools')
library('pairwiseAdonis')

setwd("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/")

getwd()

leaf_data <- read.csv("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/Dataset_Creation/Datasets/leafGenus_plantDNA.csv", header = TRUE, sep = ",")

colnames(leaf_data[0:25])

leaf_taxa <- leaf_data[7:length(leaf_data)]


v2_list <- leaf_data$V2

v2_list <- round(v2_list)

leaf_data$genetic_spec <- v2_list

leaf_data['genetic_spec'][leaf_data['genetic_spec'] == 0] <- 'Exaltata'
leaf_data['genetic_spec'][leaf_data['genetic_spec'] == 1] <- 'Syriaca'

leaf_data$binned_sp <- leaf_data$V2

leaf_data['binned_sp'][leaf_data['binned_sp'] < 0.33] <- 0
leaf_data['binned_sp'][leaf_data['binned_sp'] > 0.66] <- 1
leaf_data['binned_sp'][leaf_data['binned_sp'] >= 0.33 & leaf_data['binned_sp'] <= 0.66] <- 2

leaf_data['binned_sp'][leaf_data['binned_sp'] == 0] <- 'Exaltata'
leaf_data['binned_sp'][leaf_data['binned_sp'] == 1] <- 'Syriaca'
leaf_data['binned_sp'][leaf_data['binned_sp'] == 2] <- 'Hybrid'

leaf_data['BGR'][leaf_data['BGR'] == 'HR'] <- 'Cole Mountain'



d_1 = vegdist(leaf_taxa, method="bray", na.rm = TRUE)

leaf_PCoA <- wcmdscale(d_1, eig = TRUE)

pdf(file="./Figures/leaf_micro_rel_abundance_binned_sp_PCoA.pdf", width = 7, height = 4)

p <- ggplot(data = data.frame(leaf_PCoA$points),
            aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = leaf_data$binned_sp)) + 
  labs(color='Plant Host Species') +
  labs(title="Relative Abundance of Leaf Microbes by Plant Species PCoA") +
  scale_color_manual(values = c('darkgreen', 'green', 'yellow'))
p

dev.off()

pdf(file="./Figures/leaf_micro_rel_abundance_plant_spec_PCoA.pdf", width = 7, height = 4)

p <- ggplot(data = data.frame(leaf_PCoA$points),
                aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = leaf_data$genetic_spec)) + 
               labs(color='Plant Genotype') +
               labs(title="Relative Abundance of Phyllosphere Microbes by Plant Genotype PCoA") +
  scale_color_manual(values = c('darkgreen', 'yellow'))
p


dev.off()

pdf(file="./Figures/leaf_micro_rel_abundance_site_PCoA.pdf", width = 7, height = 4)

p <- ggplot(data = data.frame(leaf_PCoA$points),
            aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = leaf_data$site_id)) + 
  labs(color='Plant Host Species') +
  labs(title="Relative Abundance of Leaf Microbes by Site PCoA")
p


dev.off()

pdf(file="./Figures/leaf_micro_rel_abundance_bgr_PCoA.pdf", width = 7, height = 4)

p <- ggplot(data = data.frame(leaf_PCoA$points),
            aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = leaf_data$BGR)) + 
  labs(color='Broad Geographic Location') +
  labs(title="Relative Abundance of Leaf Microbes by Broad Geographic Region PCoA")
p


dev.off()


p_nova <- adonis2(d_1 ~ leaf_data$binned_sp, method = 'bray')

p_nova

p_nova <- adonis2(d_1 ~ leaf_data$genetic_spec, method = 'bray')

p_nova

p_nova <- adonis2(d_1 ~ leaf_data$site_id, method = 'bray')

p_nova

p_nova <- adonis2(d_1 ~ leaf_data$BGR, method = 'bray')

p_nova


adonis <- adonis(ord ~ X * Y, data = your.metadata) #this is basically your #PERMANOVA analysis

mod <- betadisper(d_1, leaf_data$BGR)
permutest(mod)

plot(mod)

boxplot(mod)

mod.HSD <- TukeyHSD(mod)
mod.HSD
plot(mod.HSD)

plot2 <- ggplot(data = data.frame(leaf_PCoA$points),
                aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = leaf_data$site_id))
plot2


plot3 <- ggplot(data = data.frame(leaf_PCoA$points),
                aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = leaf_data$BGR))
plot3


plot4 <- ggplot(data = data.frame(leaf_PCoA$points),
                aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = leaf_data$V2))
plot4

#######

leaf_data <- read.csv("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/Dataset_Creation/Datasets/leafGenus_plantDNA.csv", header = TRUE, sep = ",")



colnames(leaf_data[0:25])

max(leaf_data[7:length(leaf_data)])


leaf_data <- leaf_data %>% replace(is.na(.), 0)

leaf_data[7:length(leaf_data)] <- (leaf_data[7:length(leaf_data)])/100


leaf_data[7:length(leaf_data)] <- ceiling(leaf_data[7:length(leaf_data)])

leaf_taxa <- leaf_data[7:length(leaf_data)]



v2_list <- leaf_data$V2

v2_list <- round(v2_list)

leaf_data$genetic_spec <- v2_list


leaf_data['genetic_spec'][leaf_data['genetic_spec'] == 0] <- 'Exaltata'
leaf_data['genetic_spec'][leaf_data['genetic_spec'] == 1] <- 'Syriaca'



d_2 = vegdist(leaf_taxa, method="jaccard", na.rm = TRUE)

leaf_PCoA_2 <- wcmdscale(d_2, eig = TRUE)

pdf(file="./Figures/leaf_microbiome_pres_absence_PCoA.pdf", width = 7, height = 4)

p <- ggplot(data = data.frame(leaf_PCoA_2$points),
            aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = leaf_data$genetic_spec)) + 
  labs(color='Plant Host Species') +
  labs(title="Presence/Absence of Leaf Microbes PCoA") +
  scale_color_manual(values = c('darkgreen', 'yellow'))
p


p_nova <- adonis2(d_2 ~ leaf_data$genetic_spec, method = 'jaccard')

p_nova


dev.off()

plot2 <- ggplot(data = data.frame(leaf_PCoA_2$points),
                aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = leaf_data$site_id))
plot2

plot3 <- ggplot(data = data.frame(leaf_PCoA_2$points),
                aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = leaf_data$BGR))
plot3

######

total_sums <- rowSums(leaf_taxa)
total_sums
max(total_sums)
min(total_sums)

leaf_data$total_sp <- total_sums

boxplot(leaf_data$total_sp ~ leaf_data$genetic_spec)

leaf_data$genetic_spec
leaf_data$total_sp

taxa_totals <- colSums(leaf_taxa)

hist(taxa_totals)

max(taxa_totals)
min(taxa_totals)