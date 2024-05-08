library('tidyverse')
library('vegan')
library('tabula')


setwd("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/")

getwd()

leaf_data <- read.csv("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/Dataset_Creation/Datasets/leafGenus_plantDNA.csv", header = TRUE, sep = ",")


taxa <- colnames(leaf_data)

taxa[0:20]

taxa <- taxa[7:length(taxa)]

leaf_data[taxa] <- mutate_if(leaf_data[taxa], is.numeric, ~ . * 10000)

leaf_data[taxa] <- mutate_if(leaf_data[taxa], is.numeric, round)

#df <- data.matrix(leaf_data[][taxa])
df <- leaf_data[][taxa]

df <- df %>% replace(is.na(.), 0)


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

#S
#J


###################

leaf_data$alpha <- alpha

leaf_data$simp <- simp

leaf_data$shannon <- shannon



one.way <- aov(alpha ~ binned_sp, data = leaf_data)

summary(one.way)

TukeyHSD(one.way, conf.level=.95) 

###################################


jpeg(file="./Figures/leaf_alpha_div_BGR.jpeg")
boxplot(alpha ~ genetic_spec + BGR, data = leaf_data, main = 'Leaf Microbiome Alpha Diversity by BGR',
        xlab = 'BGR', ylab = 'Alpha Diversity', col = c('darkgreen', 'yellow'))
dev.off()


jpeg(file="./Figures/leaf_shannon_div_BGR.jpeg")
boxplot(shannon ~ genetic_spec + BGR, data = leaf_data, main = 'Leaf Microbiome Shannon Index by BGR',
        xlab = 'BGR', ylab = 'Shannon Index', col = rainbow(10))
dev.off()


jpeg(file="./Figures/leaf_simp_div_BGR.jpeg")
boxplot(simp ~ BGR, data = leaf_data, main = 'Leaf Microbiome Simpson Index by BGR',
        xlab = 'BGR', ylab = 'Simpson Index', col = rainbow(10))
dev.off()



jpeg(file="./Figures/leaf_alpha_div_site.jpeg")
boxplot(alpha ~ site_id, data = leaf_data, main = 'Leaf Microbiome Alpha Diversity by Site',
        xlab = 'Site', ylab = 'Alpha Diversity', col = rainbow(10))
dev.off()

jpeg(file="./Figures/leaf_shannon_div_site.jpeg")
boxplot(shannon ~ site_id, data = leaf_data, main = 'Leaf Microbiome Shannon Index by Site',
        xlab = 'Site', ylab = 'Shannon Index', col = rainbow(10))
dev.off()

jpeg(file="./Figures/leaf_simp_div_site.jpeg")
boxplot(simp ~ site_id, data = leaf_data, xlab = 'Site', ylab = 'Simpson Index', col = rainbow(10))
dev.off()

v2_list <- leaf_data$V2

v2_list <- round(v2_list)

leaf_data$genetic_spec <- v2_list

leaf_data['genetic_spec'][leaf_data['genetic_spec'] == 0] <- 'Exaltata'
leaf_data['genetic_spec'][leaf_data['genetic_spec'] == 1] <- 'Syriaca'

dev.off()
plot(leaf_data$V2, leaf_data$alpha)

jpeg(file="./Figures/leaf_alpha_div_spec_hybrid.jpeg")
boxplot(alpha ~ binned_sp, data = leaf_data, main = 'Leaf Microbiome Alpha Diversity by Plant Species',
        xlab = 'Species', ylab = 'Alpha Diversity', col = c('darkgreen', 'green', 'yellow'))
dev.off()

jpeg(file="./Figures/leaf_alpha_div_species.jpeg")
boxplot(alpha ~ genetic_spec, data = leaf_data, main = 'Leaf Microbiome Alpha Diversity by Plant Host Species',
        xlab = 'Plant Species', ylab = 'Alpha Diversity', col = c('darkgreen', 'yellow'))
dev.off()

jpeg(file="./Figures/leaf_shannon_div_species.jpeg")
boxplot(shannon ~ genetic_spec, data = leaf_data, main = 'Leaf Microbiome Shannon Index by Plant Host Species',
        xlab = 'Plant Species', ylab = 'Shannon Index', col = c('darkgreen', 'yellow'))
dev.off()

jpeg(file="./Figures/leaf_simpson_div_species.jpeg")
boxplot(simp ~ genetic_spec, data = leaf_data, xlab = 'Plant Species', ylab = 'Simpson Index', col = c('darkgreen', 'yellow'))
dev.off()


t.test(alpha ~ genetic_spec, data = leaf_data)

t.test(simp ~ genetic_spec, data = leaf_data)

t.test(shannon ~ genetic_spec, data = leaf_data)


one.way <- aov(alpha ~ site_id, data = leaf_data)

summary(one.way)

one.way <- aov(simp ~ site_id, data = leaf_data)

summary(one.way)

one.way <- aov(shannon ~ site_id, data = leaf_data)

summary(one.way)



one.way <- aov(alpha ~ BGR, data = leaf_data)

summary(one.way)

one.way <- aov(simp ~ BGR, data = leaf_data)

summary(one.way)

one.way <- aov(shannon ~ BGR, data = leaf_data)

summary(one.way)


# 3 figures arranged in 3 rows and 1 column
dev.off()

one.way <- aov(alpha ~ binned_sp, data = leaf_data)

summary(one.way)

TukeyHSD(one.way, conf.level=.95) 


p <- ggplot(data = leaf_data, aes(x=binned_sp, y=alpha)) 
p <- p + geom_boxplot(aes(fill=binned_sp), outlier.shape = NA) + scale_fill_manual(values = c("Exaltata" = 'darkgreen',"Hybrid" = 'green', "Syriaca" = "yellow"))
p <- p + geom_point(aes(color=binned_sp), show.legend = F, 
                    position=position_jitterdodge(jitter.width = 0.25)) +
  scale_color_manual(values = c('black', 'black', 'black'))
p <- p + facet_wrap( ~ binned_sp, scales='free_x')
p <- p + xlab("Site") + ylab("Alpha Diversity") + ggtitle("Plant Microbiome Alpha Diversity")
p <- p + guides(fill=guide_legend(title="Plant Species"))

p


p <- ggplot(data = leaf_data, aes(x=site_id, y=alpha)) 
p <- p + geom_boxplot(aes(fill=BGR), outlier.shape = NA)# + scale_fill_manual(values = c("Exaltata" = 'darkgreen', "Syriaca" = "yellow"))
p <- p + geom_point(aes(color=BGR), show.legend = F, 
                    position=position_jitterdodge(jitter.width = 0.25)) +
                  scale_color_manual(values = c('black', 'black', 'black'))
p <- p + facet_wrap( ~ BGR, scales='free_x')
p <- p + xlab("Site") + ylab("Alpha Diversity") + ggtitle("Plant Microbiome Alpha Diversity")
p <- p + guides(fill=guide_legend(title="Broad Geographic Location"))

p


p <- ggplot(data = leaf_data, aes(x=BGR, y=alpha)) 
p <- p + geom_boxplot(aes(fill=genetic_spec), outlier.shape = NA) + scale_fill_manual(values = c("Exaltata" = 'darkgreen', "Syriaca" = "yellow"))
p <- p + geom_point(aes(color=genetic_spec), show.legend = F, 
                    position=position_jitterdodge(jitter.width = 0.25)) +
                    scale_color_manual(values = c("Exaltata" = "black", "Syriaca" = "black"),)
p <- p + facet_wrap( ~ BGR, scales='free_x')
p <- p + xlab("Geographic Location") + ylab("Alpha Diversity") + ggtitle("Plant Microbiome Alpha Diversity")
p <- p + guides(fill=guide_legend(title="Plant Species"))

p


p <- ggplot(data = leaf_data, aes(x=BGR, y=alpha)) 
p <- p + geom_boxplot(aes(fill=binned_sp)) + scale_fill_manual(values = c("Exaltata" = 'darkgreen', "Hybrid" = 'green', "Syriaca" = "yellow"))
p <- p + geom_point(aes(color=binned_sp), show.legend = F, 
                    position=position_jitterdodge(jitter.width = 0.25)) +
  scale_color_manual(values = c("Exaltata" = "black","Hybrid" = "black",  "Syriaca" = "black"))
p <- p + facet_wrap( ~ BGR, scales="free_x")
p <- p + xlab("Geographic Location") + ylab("Alpha Diversity") + ggtitle("")
p <- p + guides(fill=guide_legend(title="Plant Species"))

p






pdf(file="./Figures/leaf_microbiome_alpha_div_1.pdf", width = 7, height = 5)

par(mfrow=c(1,3), mar = c(6, 5, 7, 1), cex.axis = 0.4, cex.lab = 1.5)
boxplot(alpha ~ site_id, data = leaf_data, #main = 'leaf Microbiome Alpha Diversity by Site',
        xlab = 'Site', ylab = '', col = rainbow(10), las =2)
boxplot(alpha ~ binned_sp + BGR, data = leaf_data, #main = 'leaf Microbiome Alpha Diversity by BGR',
        xlab = 'BGR', ylab = '', col = c('darkgreen', 'green', 'yellow'))#,
        #names=c('Team A', 'Team B', 'Team C', 'Team A', 'Team B', 'Team C', 'Team A', 'Team B', 'Team C'))
  
boxplot(alpha ~ binned_sp, data = leaf_data, xlab = 'Plant Species', #main = 'leaf Microbiome Alpha Diversity by Plant Host Species',
        ylab = '', col = c('darkgreen', 'green', 'yellow'))
mtext("Leaf Microbiome Alpha Diversity is Influenced by Plant Host Species, Site, and BGR",   # Add main title
      side = 3,
      line = -5,
      outer = TRUE)

mtext("Alpha Diversity",   # Add main title
      side = 2,
      line = -2.25,
      outer = TRUE)
dev.off()




pdf(file="./Figures/leaf_microbiome_alpha_div_hybrid.pdf", width = 7, height = 5)

par(mfrow=c(1,3), mar = c(6, 5, 7, 1), cex.axis = 0.4, cex.lab = 1.5)
boxplot(alpha ~ site_id, data = leaf_data, #main = 'leaf Microbiome Alpha Diversity by Site',
        xlab = 'Site', ylab = '', col = rainbow(10), las =2)
boxplot(alpha ~ BGR, data = leaf_data, #main = 'leaf Microbiome Alpha Diversity by BGR',
        xlab = 'BGR', ylab = '', col = rainbow(10))
boxplot(alpha ~ binned_sp, data = leaf_data, xlab = 'Plant Species', #main = 'leaf Microbiome Alpha Diversity by Plant Host Species',
        ylab = '', col = c('darkgreen', 'green', 'yellow'))
mtext("Leaf Microbiome Alpha Diversity is Influenced by Plant Host Species, Site, and BGR",   # Add main title
      side = 3,
      line = -5,
      outer = TRUE)

mtext("Alpha Diversity",   # Add main title
      side = 2,
      line = -2.25,
      outer = TRUE)
dev.off()





pdf(file="./Figures/leaf_microbiome_shannon_index.pdf", width = 7, height = 4)

par(mfrow=c(1,3), mar = c(6, 3.5, 7, 0), cex.axis = 0.5, cex.lab = 1.5)
boxplot(shannon ~ site_id, data = leaf_data, #main = 'leaf Microbiome Alpha Diversity by Site',
        xlab = 'Site', ylab = '', col = rainbow(10), las=2)
boxplot(shannon ~ BGR, data = leaf_data, #main = 'leaf Microbiome Alpha Diversity by BGR',
        xlab = 'BGR', ylab = '', col = rainbow(10))
boxplot(shannon ~ genetic_spec, data = leaf_data, xlab = 'Plant Species', #main = 'leaf Microbiome Alpha Diversity by Plant Host Species',
        ylab = '', col = c('darkgreen', 'yellow'))
mtext("Leaf Microbiome Shannon Index is Influenced by Site and BGR",   # Add main title
      side = 3,
      line = -5,
      outer = TRUE)

mtext("Shannon Index",   # Add main title
      side = 2,
      line = -1.25,
      outer = TRUE)
dev.off()



