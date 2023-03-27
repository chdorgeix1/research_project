setwd("C:/Users/15505/Documents/GitHub/research_project/rarefaction_curve")

df <- read.csv("202109_phylum_with_site.csv", header = TRUE, sep = ",")

df_sub <- subset(df, site == 'KCKT')

nrow(df_sub)

for (i in 1:nrow(df_sub))
{
  print(sum(df_sub[i,] != 0))
}
