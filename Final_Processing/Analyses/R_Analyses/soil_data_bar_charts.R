library('tidyverse')
library('vegan')
library('tabula')
library('ggpubr')
library('CoreMicrobiomeR')
library('reshape2')

setwd("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/")

getwd()

soil_data <- read.csv("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/Dataset_Creation/Datasets/soilPhylum_plantDNA.csv", header = TRUE, sep = ",")

dim(soil_data)

taxa <- colnames(soil_data)

taxa <- taxa[7:length(taxa)]

soil_data <- soil_data %>% replace(is.na(.), 0)

head(taxa)
tail(taxa)

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

################

tot_sum = colSums(soil_data[, c(taxa)])

tot_sum <- tot_sum/sum(tot_sum)

val_list <- tot_sum

val_list <- sort(val_list, decreasing = TRUE)

sort(tot_sum, decreasing = TRUE)[0:5]

tot_sum <- sort(tot_sum, decreasing=TRUE)

df <- as.data.frame(tot_sum)

df$taxa_num <- seq(1, length(tot_sum))


pdf(file="./Figures/soil_micro_rel_abund_phyla.pdf", width = 15, height = 10)
ggplot(df, aes(x = taxa_num, y = tot_sum)) + 
  geom_point() +
  xlab("Phyla Number") +
  ylab("Relative Abundance") +
  ggtitle('Relative abundance Proportion by Phyla') +
  theme(plot.title.position = 'plot', 
        plot.title = element_text(hjust = 0.5))
dev.off()
sum(tot_sum[0:5])

tot_sum[0:6]

################

test_col <- subset(soil_data, site_id == unique(soil_data$site_id)[1])
col_sums = colSums(test_col[, c(taxa)])
col_sums <- col_sums/(sum(col_sums))
site_df <- as.data.frame(col_sums)

for (site_x in unique(soil_data$site_id)) {
  sub_df <- subset(soil_data, site_id == site_x)
  
  sub_df_sum = colSums(sub_df[, c(taxa)])
  sub_df_sum <- sub_df_sum/(sum(sub_df_sum))
  
  site_df[site_x] <- sub_df_sum
  dim(site_df)
}

colnames(site_df)

site_df <- subset(site_df, select = -c(col_sums))

testdf <- t(site_df)

df_long <- melt(testdf, id.vars = "Sample", variable.name = "Phyla")

site_plot <- ggplot(df_long, aes(x = Var1, y = value, fill = Var2)) + 
  geom_bar(stat = "identity")+
  xlab("Site") +
  ylab("Relative Abundance") +
  theme(axis.title.x=element_text(size=14,face="bold"),
        axis.title.y=element_text(size=14,face="bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.text = element_text(size=6))
site_plot <- site_plot +  guides(fill=guide_legend(title="Phyla"))
site_plot



edf <- subset(soil_data, genetic_spec == 'Exaltata')
sdf <- subset(soil_data, genetic_spec == 'Syriaca')

edf_sum = colSums(edf[, c(taxa)])
length(edf_sum)

sdf_sum = colSums(sdf[, c(taxa)])
length(sdf_sum)

edf_sum <- edf_sum/(sum(edf_sum))
sdf_sum <- sdf_sum/(sum(sdf_sum))

fin_df <-as.data.frame((edf_sum))


colnames(fin_df)
fin_df$syriaca <- (sdf_sum)

colnames(fin_df) <- c("Exaltata", "Syriaca")

fin_df <- as.data.frame(fin_df)

testdf <- t(fin_df)

df_long <- melt(testdf, id.vars = "Sample", variable.name = "Phyla")


sp_plot <- ggplot(df_long, aes(x = Var1, y = value, fill = Var2)) + 
  geom_bar(stat = "identity")+
  xlab("Genotype") +
  ylab("Relative Abundance") +
  theme(axis.title.x=element_text(size=14,face="bold"),
        axis.title.y=element_text(size=14,face="bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.text = element_text(size=6))
sp_plot <- sp_plot +  guides(fill=guide_legend(title="Phyla"))


################



edf <- subset(soil_data, binned_sp == 'Exaltata')
hdf <- subset(soil_data, binned_sp == 'Hybrid')
sdf <- subset(soil_data, binned_sp == 'Syriaca')


edf_sum = colSums(edf[, c(taxa)])
length(edf_sum)

hdf_sum = colSums(hdf[, c(taxa)])
length(hdf_sum)

sdf_sum = colSums(sdf[, c(taxa)])
length(sdf_sum)

edf_sum <- edf_sum/(sum(edf_sum))
hdf_sum <- hdf_sum/(sum(hdf_sum))
sdf_sum <- sdf_sum/(sum(sdf_sum))

fin_df <-as.data.frame((edf_sum))

fin_df$hybrid <- (hdf_sum)

fin_df$syriaca <- (sdf_sum)

colnames(fin_df) <- c("Exaltata", 'Hybrid', "Syriaca")

fin_df <- as.data.frame(fin_df)

testdf <- t(fin_df)

df_long <- melt(testdf, id.vars = "Sample", variable.name = "Phyla")


binned_plot <- ggplot(df_long, aes(x = Var1, y = value, fill = Var2)) + 
  geom_bar(stat = "identity")+
  xlab("Genotype") +
  ylab("Relative Abundance") +
  theme(axis.title.x=element_text(size=14,face="bold"),
        axis.title.y=element_text(size=14,face="bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.text = element_text(size=6))
binned_plot <- binned_plot +  guides(fill=guide_legend(title="Phyla"))



#################

bdf <- subset(soil_data, BGR == 'Blacksburg')
cdf <- subset(soil_data, BGR == 'Cole Mountain')
wdf <- subset(soil_data, BGR == 'Wintergreen')

dim(bdf)
dim(cdf)
dim(wdf)

bdf_sum = colSums(bdf[, c(taxa)])

cdf_sum = colSums(cdf[, c(taxa)])

wdf_sum = colSums(wdf[, c(taxa)])

bdf_sum <- bdf_sum/(sum(bdf_sum))
cdf_sum <- cdf_sum/(sum(cdf_sum))
wdf_sum <- wdf_sum/(sum(wdf_sum))

fin_df <-as.data.frame((bdf_sum))

fin_df$cole_mountain <- (cdf_sum)

fin_df$wintegreen <- (wdf_sum)

colnames(fin_df) <- c("Blacksburg", 'Cole Mountain', "Wintergreen")

fin_df <- as.data.frame(fin_df)

testdf <- t(fin_df)

df_long <- melt(testdf, id.vars = "Sample", variable.name = "Phyla")


bgr_plot <- ggplot(df_long, aes(x = Var1, y = value, fill = Var2)) + 
  geom_bar(stat = "identity")+
  xlab("Broad Geographic Location") +
  ylab("Relative Abundance") +
  theme(axis.title.x=element_text(size=14,face="bold"),
        axis.title.y=element_text(size=14,face="bold"),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.text = element_text(size=6))
bgr_plot <- bgr_plot +  guides(fill=guide_legend(title="Phyla"))

bgr_plot
pdf(file="./Figures/soil_micro_bgr_rel_abund_BC.pdf", width = 15, height = 10)
bgr_plot
dev.off()


pdf(file="./Figures/soil_micro_sp_rel_abund_BC.pdf", width = 15, height = 10)
sp_plot
dev.off()


pdf(file="./Figures/soil_micro_binned_rel_abund_BC.pdf", width = 15, height = 10)
binned_plot
dev.off()


pdf(file="./Figures/soil_micro_site_rel_abund_BC.pdf", width = 15, height = 10)
site_plot
dev.off()


