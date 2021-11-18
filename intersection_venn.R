rm(list = ls())
library(stringr)
library(qpcR)
temp_nm1 <- "cecum_otu"
temp_nm2 <- "fecal_otu"
temp_nm3 <- "ileum_otu"

df1 <- read.csv(paste("origin/00result_na_rm/cecum/",temp_nm1, "_asso_metabo_merge_origin", "_NA_rm_trimmed.csv",
                     sep = ""), header=T, check.names=F)
df2 <- read.csv(paste("origin/00result_na_rm/fecal/",temp_nm2, "_asso_metabo_merge_origin", "_NA_rm_trimmed.csv",
                      sep = ""), header=T, check.names=F)
df3 <- read.csv(paste("origin/00result_na_rm/ileum/",temp_nm3, "_asso_metabo_merge_origin", "_NA_rm_trimmed.csv",
                      sep = ""), header=T, check.names=F)
dim(df1); dim(df2); dim(df3)

# intersection based on metabolites
# among three gut position
fecal.cecum.ileum.metabo.incl <- intersect(intersect(unique(df1$metabo_mz), unique(df2$metabo_mz)), unique(df3$metabo_mz))
# between fecal and cecum
fecal.cecum.metabo.incl <- intersect(unique(df1$metabo_mz), unique(df2$metabo_mz))
# betweem fecal and ileum
fecal.ileum.metabo.incl <- intersect(unique(df2$metabo_mz), unique(df3$metabo_mz))
# between cecum and ileum
cecum.ileum.metabo.incl <- intersect(unique(df1$metabo_mz), unique(df3$metabo_mz))


# intersection based on microbial taxa or otu
fecal.cecum.ileum.micro.incl <- intersect(intersect(unique(df1$taxa), unique(df2$taxa)), unique(df3$taxa))
# between fecal and cecum
fecal.cecum.micro.incl <- intersect(unique(df1$taxa), unique(df2$taxa))
# betweem fecal and ileum
fecal.ileum.micro.incl <- intersect(unique(df2$taxa), unique(df3$taxa))
# between cecum and ileum
cecum.ileum.micro.incl <- intersect(unique(df1$taxa), unique(df3$taxa))

# output
metabo.micro.merge.incl <- qpcR:::cbind.na(fecal.cecum.ileum.metabo.incl, 
                                           fecal.cecum.metabo.incl,fecal.ileum.metabo.incl,
                                           cecum.ileum.metabo.incl,fecal.cecum.ileum.micro.incl,
                                           fecal.cecum.micro.incl,fecal.ileum.micro.incl,cecum.ileum.micro.incl)
write.csv(metabo.micro.merge.incl, paste("origin/00result_na_rm/intersection_sum/", "metabo_micro_merge_incl_otu.csv",
                                         sep = ""), row.names = F)






## second part: draw venn plot
# load library

rm(list = ls())
library(stringr)
library(qpcR)
library(VennDiagram)
temp_nm1 <- "cecum_species"
temp_nm2 <- "fecal_species"
temp_nm3 <- "ileum_species"

df1 <- read.csv(paste("origin/00result_na_rm/cecum/",temp_nm1, "_asso_metabo_merge_origin", "_NA_rm_trimmed.csv",
                      sep = ""), header=T, check.names=F)
df2 <- read.csv(paste("origin/00result_na_rm/fecal/",temp_nm2, "_asso_metabo_merge_origin", "_NA_rm_trimmed.csv",
                      sep = ""), header=T, check.names=F)
df3 <- read.csv(paste("origin/00result_na_rm/ileum/",temp_nm3, "_asso_metabo_merge_origin", "_NA_rm_trimmed.csv",
                      sep = ""), header=T, check.names=F)

metabo_set1 <- unique(df1$metabo_mz)
metabo_set2 <- unique(df2$metabo_mz)
metabo_set3 <- unique(df3$metabo_mz)

micro_set1 <- unique(df1$taxa)
micro_set2 <- unique(df2$taxa)
micro_set3 <- unique(df3$taxa)

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
# Chart
venn.diagram(
  x = list(micro_set1, micro_set2, micro_set3),
  category.names = c("cecum" , "fecal" , "ileum"),
  filename = paste("origin/00result_na_rm/intersection_sum/", "metabo_micro_merge_incl_species_micro.png",
                                         sep = ""),
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


## Application on rap Lyrics
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)
venn.diagram(
  x = list(
    metabo_set1, metabo_set2, metabo_set3
  ),
  category.names = c("cecum" , "fecal" , "ileum"),
  filename = paste("origin/00result_na_rm/intersection_sum/", "metabo_micro_merge_incl_genus2.png",
                   sep = ""),
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 0.6,
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  rotation = 1
)










### 科属种不同分类水平下分别得到了三个部位交集的代谢物，接下来我们
# 这些代谢物在三个分类水平交集的代谢物还剩下哪些。
rm(list = ls())
path.w <- "origin/00result_na_rm/intersection_sum/"
family <- read.csv(paste(path.w,"metabo_micro_merge_incl_family.csv",
                         sep = ""), header=T, check.names=F)
genus <- read.csv(paste(path.w,"metabo_micro_merge_incl_genus.csv",
                        sep = ""), header=T, check.names=F)
species <- read.csv(paste(path.w,"metabo_micro_merge_incl_species.csv",
                          sep = ""), header=T, check.names=F)
otu <- read.csv(paste(path.w, "metabo_micro_merge_incl_otu.csv",
                      sep = ""), header=T, check.names=F)
# 去除第一列中的NA，留下代谢物即可
f.metabo.noNA <- family$fecal.cecum.ileum.metabo.incl[!is.na(family$fecal.cecum.ileum.metabo.incl)]
g.metabo.noNA <- genus$fecal.cecum.ileum.metabo.incl[!is.na(genus$fecal.cecum.ileum.metabo.incl)]
s.metabo.noNA <- species$fecal.cecum.ileum.metabo.incl[!is.na(species$fecal.cecum.ileum.metabo.incl)]
otu.metabo.noNA <- otu$fecal.cecum.ileum.metabo.incl[!is.na(otu$fecal.cecum.ileum.metabo.incl)]
metabo.overlap.f.g.s <- intersect(intersect(f.metabo.noNA, g.metabo.noNA), s.metabo.noNA)
metabo.overlap.f.g.s.otu <- intersect(metabo.overlap.f.g.s, otu.metabo.noNA)
# 只输出科属种
write.csv(metabo.overlap.f.g.s, file = paste(path.w, "micro_asso_metabo_overlap_f_g_s.csv", sep = ""))
# 输出科属种、otu
write.csv(metabo.overlap.f.g.s.otu, file = paste(path.w, "micro_asso_metabo_overlap_f_g_s_otu.csv", sep = ""))







