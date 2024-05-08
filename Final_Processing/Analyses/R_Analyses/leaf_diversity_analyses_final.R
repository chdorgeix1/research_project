library('tidyverse')
library('vegan')
library('tabula')
library('ggpubr')

setwd("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/")

getwd()

leaf_data <- read.csv("C:/Users/15404/Documents/GitHub/research_project/Final_Processing/Dataset_Creation/Datasets/leafGenus_plantDNA.csv", header = TRUE, sep = ",")


taxa <- colnames(leaf_data)

taxa[0:20]

taxa <- taxa[7:length(taxa)]

leaf_data[taxa] <- mutate_if(leaf_data[taxa], is.numeric, ~ . * 10000)

leaf_data[taxa] <- mutate_if(leaf_data[taxa], is.numeric, round)

df <- leaf_data[][taxa]

df <- df %>% replace(is.na(.), 0)

df[taxa][0:5]

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

######################

# Diversity Calculations

H <- diversity(df)
simp <- diversity(df, index = "simpson")
shannon <- diversity(df, index = 'shannon')
invsimp <- diversity(df, "inv")
unbias.simp <- simpson.unb(df)
alpha <- fisher.alpha(df)

#pairs(cbind(H, simp, invsimp, unbias.simp, alpha), col="blue")

S <- specnumber(df) ## rowSums(BCI > 0) does the same...
J <- H/log(S)

S

leaf_data$alpha <- alpha
leaf_data$simp <- simp
leaf_data$shannon <- shannon
leaf_data$spec_num <- S

#####################
# Alpha Diversity Graphics

give.n <- function(x){
  return(c(y = min(x)-3, label = length(x)))
  # experiment with the multiplier to find the perfect position
}

# Plant Species(E, H, + S) x Alpha Diversity Boxplot
pdf(file="./Figures/leaf_microbiome_plant_species.pdf", width = 7, height = 5)
p <- ggplot(data = leaf_data, aes(x=binned_sp, y=alpha)) + theme(plot.title.position = 'plot', 
                                                                 plot.title = element_text(hjust = 0.5)) 
p <- p + geom_boxplot(aes(fill=binned_sp), outlier.shape = NA) + 
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.75), size = 4)+
  scale_fill_manual(values = c("Exaltata" = 'darkgreen',"Hybrid" = 'green', "Syriaca" = "yellow"))
p <- p + geom_point(aes(color=binned_sp), show.legend = F, 
                    position=position_jitterdodge(jitter.width = 0.25)) +
  scale_color_manual(values = c('black', 'black', 'black'))
p <- p + facet_wrap( ~ binned_sp, scales='free_x')
p <- p + xlab("Plant Host Genotype") + ylab("Alpha Diversity") + ggtitle("Phyllosphere Microbiome Alpha Diversity")
p <- p + guides(fill=guide_legend(title="Plant Species"))
p
dev.off()

# Plant Species x Alpha Diversity ANOVA + TUKEY
one.way <- aov(alpha ~ binned_sp, data = leaf_data)
summary(one.way)
TukeyHSD(one.way, conf.level=.95) 

# Plant Species(Grouped by BGL) x Alpha Diversity Boxplot
pdf(file="./Figures/leaf_microbiome_species_bgl.pdf", width = 7, height = 5)
p <- ggplot(data = leaf_data, aes(x=BGR, y=alpha)) 
p <- p + geom_boxplot(aes(fill=genetic_spec), outlier.shape = NA) + scale_fill_manual(values = c("Exaltata" = 'darkgreen', "Syriaca" = "yellow"))
p <- p + geom_point(aes(color=genetic_spec), show.legend = F, 
                    position=position_jitterdodge(jitter.width = 0.25)) + 
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.75), size = 5) +
  scale_color_manual(values = c("Exaltata" = "black", "Syriaca" = "black"),)
p <- p + facet_wrap( ~ BGR, scales='free_x')
p <- p + xlab("Geographic Location") + ylab("Alpha Diversity") + ggtitle("Plant Microbiome Alpha Diversity")
p <- p + guides(fill=guide_legend(title="Plant Species"))
p

dev.off()

# Plant Speices within geographic locations Alpha Div ANOVA + TUKEY
df1 <- subset(leaf_data, BGR == 'Wintergreen')
one.way <- aov(alpha ~ genetic_spec, data = df1)
summary(one.way)
TukeyHSD(one.way, conf.level=.95)
df2 <- subset(leaf_data, BGR == 'Cole Mountain')
one.way <- aov(alpha ~ genetic_spec, data = df2)
summary(one.way)
TukeyHSD(one.way, conf.level=.95)


# Plant Species(E,H,S)(Grouped by BGL) x Alpha Diversity Boxplot
pdf(file="./Figures/leaf_microbiome_species_hybrids_bgl.pdf", width = 7, height = 5)
p <- ggplot(data = leaf_data, aes(x=BGR, y=alpha)) 
p <- p + geom_boxplot(aes(fill=binned_sp)) + 
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.75), size = 5) + 
  scale_fill_manual(values = c("Exaltata" = 'darkgreen', "Hybrid" = 'green', "Syriaca" = "yellow"))
p <- p + geom_point(aes(color=binned_sp), show.legend = F, 
                    position=position_jitterdodge(jitter.width = 0.25)) +
  scale_color_manual(values = c("Exaltata" = "black","Hybrid" = "black",  "Syriaca" = "black"))
p <- p + facet_wrap( ~ BGR, scales="free_x")
p <- p + xlab("Geographic Location") + ylab("Alpha Diversity") + ggtitle("")
p <- p + guides(fill=guide_legend(title="Plant Species"))
p

# Binned Speices within geographic locations Alpha Div ANOVA + TUKEY
df1 <- subset(leaf_data, BGR == 'Wintergreen')
one.way <- aov(alpha ~ binned_sp, data = df1)
summary(one.way)
TukeyHSD(one.way, conf.level=.95)

df2 <- subset(leaf_data, BGR == 'Cole Mountain')
one.way <- aov(alpha ~ binned_sp, data = df2)
summary(one.way)
TukeyHSD(one.way, conf.level=.95)


# Site(grouped by BGL) x Alpha Diversity Boxplot
pdf(file="./Figures/leaf_microbiome_alpha_div_site.pdf", width = 9, height = 5)
p <- ggplot(data = leaf_data, aes(x=site_id, y=alpha)) + theme(plot.title.position = 'plot', 
                                                               plot.title = element_text(hjust = 0.5)) 
p <- p + geom_boxplot(aes(fill=BGR), outlier.shape = NA) +
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.5), size = 4)

# + scale_fill_manual(values = c("Exaltata" = 'darkgreen', "Syriaca" = "yellow"))
p <- p + geom_point(aes(color=BGR), show.legend = F, 
                    position=position_jitterdodge(jitter.width = 0.5)) +
  scale_color_manual(values = c('black', 'black', 'black'))
p <- p + facet_wrap( ~ BGR, scales='free_x')
p <- p + xlab("Site") + ylab("Alpha Diversity") + ggtitle("Phyllosphere Microbiome Alpha Diversity")
p <- p + guides(fill=guide_legend(title="Broad Geographic Location"))
p
dev.off()

# Site x Alpha Diversity ANOVA + TUKEY
df1 <- subset(leaf_data, BGR == 'Blacksburg')
one.way <- aov(alpha ~ site_id, data = df1)
summary(one.way)
TukeyHSD(one.way, conf.level=.95)

dfdf2 <- subset(leaf_data, BGR == 'Cole Mountain')
one.way <- aov(alpha ~ site_id, data = df2)
summary(one.way)
TukeyHSD(one.way, conf.level=.95)

df3 <- subset(leaf_data, BGR == 'Wintergreen')
one.way <- aov(alpha ~ site_id, data = df3)
summary(one.way)
TukeyHSD(one.way, conf.level=.95)

one.way <- aov(alpha ~ site_id, data = leaf_data)
summary(one.way)
TukeyHSD(one.way, conf.level=.95) 

# BGR x Alpha Diversity Boxplot
pdf(file="./Figures/leaf_microbiome_alpha_div_bgl.pdf", width = 7, height = 5)
p <- ggplot(data = leaf_data, aes(x=BGR, y=alpha))  + theme(plot.title.position = 'plot', 
                                                            plot.title = element_text(hjust = 0.5)) 
p <- p + geom_boxplot(aes(fill=BGR), outlier.shape = NA) +
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.75), size = 5)
p <- p + geom_point(aes(color=BGR), show.legend = F, 
                    position=position_jitterdodge(jitter.width = 0.25)) +
  scale_color_manual(values = c('black', 'black', 'black'))
p <- p + facet_wrap( ~ BGR, scales='free_x')
p <- p + xlab("Broad Geographic Location") + ylab("Alpha Diversity") + ggtitle("Phyllosphere Microbiome Alpha Diversity")
#p <- p + guides(fill=guide_legend(title="Broad Geographic Location"))
p <- p + guides(fill="none")
p
dev.off()

# BGR x Alpha Diversity ANOVA + TUKEY
one.way <- aov(alpha ~ BGR, data = leaf_data)
summary(one.way)
TukeyHSD(one.way, conf.level=.95)











# Shannon Index Graphics

give.n <- function(x){
  return(c(y = min(x)-0.2, label = length(x)))
  # experiment with the multiplier to find the perfect position
}

# Plant Species(E, H, + S) x Shannon Diversity Boxplot
pdf(file="./Figures/leaf_micro_binned_sp_shannon.pdf", width = 7, height = 5)
p <- ggplot(data = leaf_data, aes(x=binned_sp, y=shannon)) + theme(plot.title.position = 'plot', 
                                                                   plot.title = element_text(hjust = 0.5)) 
p <- p + geom_boxplot(aes(fill=binned_sp), outlier.shape = NA) + 
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.5), size = 4)+
  scale_fill_manual(values = c("Exaltata" = 'darkgreen',"Hybrid" = 'green', "Syriaca" = "yellow"))
p <- p + geom_point(aes(color=binned_sp), show.legend = F, 
                    position=position_jitterdodge(jitter.width = 0.5)) +
  scale_color_manual(values = c('black', 'black', 'black'))
p <- p + facet_wrap( ~ binned_sp, scales='free_x')
p <- p + xlab("Plant Genotype") + ylab("Shannon Diversity") + ggtitle("Phyllosphere Microbiome Shannon Diversity")
p <- p + guides(fill=guide_legend(title="Plant Genotype"))
p
dev.off()

# Plant Species x shannon Diversity ANOVA + TUKEY
one.way <- aov(shannon ~ binned_sp, data = leaf_data)
summary(one.way)
TukeyHSD(one.way, conf.level=.95) 

# Plant Species(Grouped by BGL) x shannon Diversity Boxplot
pdf(file="./Figures/leaf_microbiome_species_bgl.pdf", width = 7, height = 5)
p <- ggplot(data = leaf_data, aes(x=BGR, y=shannon)) 
p <- p + geom_boxplot(aes(fill=genetic_spec), outlier.shape = NA) + scale_fill_manual(values = c("Exaltata" = 'darkgreen', "Syriaca" = "yellow"))
p <- p + geom_point(aes(color=genetic_spec), show.legend = F, 
                    position=position_jitterdodge(jitter.width = 0.25)) + 
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.75), size = 5) +
  scale_color_manual(values = c("Exaltata" = "black", "Syriaca" = "black"),)
p <- p + facet_wrap( ~ BGR, scales='free_x')
p <- p + xlab("Geographic Location") + ylab("Shannon Diversity") + ggtitle("Plant Microbiome Shannon Diversity")
p <- p + guides(fill=guide_legend(title="Plant Species"))
p

dev.off()

# Plant Speices within geographic locations shannon Div ANOVA + TUKEY
df1 <- subset(leaf_data, BGR == 'Wintergreen')
one.way <- aov(shannon ~ genetic_spec, data = df1)
summary(one.way)
TukeyHSD(one.way, conf.level=.95)
df2 <- subset(leaf_data, BGR == 'Cole Mountain')
one.way <- aov(shannon ~ genetic_spec, data = df2)
summary(one.way)
TukeyHSD(one.way, conf.level=.95)


# Plant Species(E,H,S)(Grouped by BGL) x shannon Diversity Boxplot
pdf(file="./Figures/leaf_microbiome_species_hybrids_bgl.pdf", width = 7, height = 5)
p <- ggplot(data = leaf_data, aes(x=BGR, y=shannon)) 
p <- p + geom_boxplot(aes(fill=binned_sp)) + 
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.75), size = 5) + 
  scale_fill_manual(values = c("Exaltata" = 'darkgreen', "Hybrid" = 'green', "Syriaca" = "yellow"))
p <- p + geom_point(aes(color=binned_sp), show.legend = F, 
                    position=position_jitterdodge(jitter.width = 0.25)) +
  scale_color_manual(values = c("Exaltata" = "black","Hybrid" = "black",  "Syriaca" = "black"))
p <- p + facet_wrap( ~ BGR, scales="free_x")
p <- p + xlab("Geographic Location") + ylab("Shannon Diversity") + ggtitle("Plant Microbiome Shannon Diversity")
p <- p + guides(fill=guide_legend(title="Plant Species"))
p

# Binned Speices within geographic locations shannon Div ANOVA + TUKEY
df1 <- subset(leaf_data, BGR == 'Wintergreen')
one.way <- aov(shannon ~ binned_sp, data = df1)
summary(one.way)
TukeyHSD(one.way, conf.level=.95)

df2 <- subset(leaf_data, BGR == 'Cole Mountain')
one.way <- aov(shannon ~ binned_sp, data = df2)
summary(one.way)
TukeyHSD(one.way, conf.level=.95)


# Site(grouped by BGL) x shannon Diversity Boxplot
pdf(file="./Figures/leaf_micro_shannon_site.pdf", width = 9, height = 5)
p <- ggplot(data = leaf_data, aes(x=site_id, y=shannon)) + theme(plot.title.position = 'plot', 
                                                                 plot.title = element_text(hjust = 0.5))  
p <- p + geom_boxplot(aes(fill=BGR), outlier.shape = NA) +
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.5), size = 4)

# + scale_fill_manual(values = c("Exaltata" = 'darkgreen', "Syriaca" = "yellow"))
p <- p + geom_point(aes(color=BGR), show.legend = F, 
                    position=position_jitterdodge(jitter.width = 0.5)) +
  scale_color_manual(values = c('black', 'black', 'black'))
p <- p + facet_wrap( ~ BGR, scales='free_x')
p <- p + xlab("Site") + ylab("Shannon Diversity") + ggtitle("Phyllosphere Microbiome Shannon Diversity")
p <- p + guides(fill=guide_legend(title="Broad Geographic Location"))
p
dev.off()

# Site x shannon Diversity ANOVA + TUKEY
one.way <- aov(shannon ~ site_id, data = leaf_data)
summary(one.way)
TukeyHSD(one.way, conf.level=.95) 

# Site x Alpha Diversity ANOVA + TUKEY
df1 <- subset(leaf_data, BGR == 'Blacksburg')
one.way <- aov(shannon ~ site_id, data = df1)
summary(one.way)
TukeyHSD(one.way, conf.level=.95)

df2 <- subset(leaf_data, BGR == 'Cole Mountain')
one.way <- aov(shannon ~ site_id, data = df2)
summary(one.way)
TukeyHSD(one.way, conf.level=.95)

df3 <- subset(leaf_data, BGR == 'Wintergreen')
one.way <- aov(shannon ~ site_id, data = df3)
summary(one.way)
TukeyHSD(one.way, conf.level=.95)

one.way <- aov(shannon ~ site_id, data = leaf_data)
summary(one.way)
TukeyHSD(one.way, conf.level=.95) 


# BGR x shannon Diversity Boxplot
p <- ggplot(data = leaf_data, aes(x=BGR, y=shannon)) 
p <- p + geom_boxplot(aes(fill=BGR), outlier.shape = NA) +
  stat_summary(fun.data = give.n, geom = "text", fun = median,
               position = position_dodge(width = 0.75), size = 5)
p <- p + geom_point(aes(color=BGR), show.legend = F, 
                    position=position_jitterdodge(jitter.width = 0.25)) +
  scale_color_manual(values = c('black', 'black', 'black'))
p <- p + facet_wrap( ~ BGR, scales='free_x')
p <- p + xlab("Site") + ylab("Shannon Diversity") + ggtitle("Plant Microbiome Shannon Diversity")
p <- p + guides(fill=guide_legend(title="Broad Geographic Location"))
p

# BGR x shannon Diversity ANOVA + TUKEY
one.way <- aov(shannon ~ BGR, data = leaf_data)
summary(one.way)
TukeyHSD(one.way, conf.level=.95)







