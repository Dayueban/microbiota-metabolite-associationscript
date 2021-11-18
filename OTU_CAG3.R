## version 3版本的代码用于原始ileum, cecum, feca的OTU数据来构建cag
####基于R语言的SparCC分析
rm(list=ls())
#install.packages("robCompositions")
#library(robCompositions)
library(vegan)
temp_nm1 <- "../otu/ileum/"
temp_nm2 <- "../otu/cecum/"
temp_nm3 <- "../otu/fecal/"
temp_nm4 <- "fecal" # 不同部位时该参数需要改动
# load working datasheet
options (stringsAsFactors = FALSE)
#work_df <- read.csv("Ileum_OTU.csv",header = T,row.names = 1,check.names = F)
a <- read.csv(paste(temp_nm1, "ileum_otu_80.csv", sep = ""),
              header = T,row.names = 1,check.names = F)
b <- read.csv(paste(temp_nm2, "cecum_otu_80.csv", sep = ""),
              header = T,row.names = 1,check.names = F)
c <- read.csv(paste(temp_nm3, "fecal_otu_80.csv", sep = ""),
              header = T,row.names = 1,check.names = F)
#overlap_index <- intersect(intersect(rownames(a), rownames(b)), rownames(c))
meta_df <- c #a,b,c三个数据框根据所需要的肠道部位来改动
otu <- meta_df
# centered log ratio transformation
# otu.clr <- cenLR(otu, exp(10))
otu <- decostand(otu,"log")
#otu <- log10(otu)
# load library


# sparcc.otu<-sparcc(otu)##Sparcc计算
# cor函数计算相关矩阵
#cor.otu <- cor(otu, method = "spearman")
source("normalization.R")
source("sparcc.R")
sparcc.otu <- sparcc(otu)
cor.otu <- sparcc.otu$Cor
rownames(cor.otu) <- colnames(otu)
colnames(cor.otu) <- colnames(otu)
#r_result <- data.frame(row=rownames(s)[row(s)[upper.tri(s)]],col=colnames(s)[col(s)[upper.tri(s)]],corr=s[upper.tri(s)])
#write.table(r_result,"otu_cor.txt",sep="\t")
s.matrix <- cor.otu
dis.matrix <- as.dist(1 - s.matrix)
otu.ward <- hclust(dis.matrix,method="ward.D2")
plot(otu.ward,hang=-1)
k.temp = 8
rect.hclust(otu.ward, k = k.temp)

# 可以看到k = k.temp 为比较合理的结果，根据k.temp参数来选择
otu.cluster <- cutree(otu.ward, k = k.temp) ##k=clust数量
# avoid change the names of idv. when transposing the dataframe:otu
idv.name <- rownames(otu)
otu.tr <- as.data.frame(t(otu))
names(otu.tr) <- idv.name
otu.cluster.8 <- cbind(otu.cluster, otu.tr)
write.csv(otu.cluster.8, paste("../otu/CAG2/", temp_nm4, "otu_cluster", k.temp, "cag.csv", sep = "_"),
          row.names = TRUE)


### method1：配对比较分析,进行模块之间是否显著差异的检测
pairwise.adonis <-function(x,factors, sim.method, p.adjust.m){
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),] ~factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method);
    pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted =p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

adonis.1 <- pairwise.adonis(otu.cluster.8, otu.cluster.8$otu.cluster,sim.method="bray",p.adjust.m="bonferroni")
adonis.1

# With 8 CAGs
# 1. The OTUs that belong to each CAG are selected and their abundance is combined.
# In this way, the abundance of each CAG in a given sample is equal to the sum of the abundance of the OTUs that compose it.
# 2. A table is generated with the abundances of each CAG in each sample
CAG_n8 <- factor(otu.cluster)
CAG_n8_a = otu[,CAG_n8 == 1]
CAG_n8_b = otu[,CAG_n8 == 2]
CAG_n8_c = otu[,CAG_n8 == 3]
CAG_n8_d = otu[,CAG_n8 == 4]
CAG_n8_e = otu[,CAG_n8 == 5]
CAG_n8_f = otu[,CAG_n8 == 6]
CAG_n8_g = otu[,CAG_n8 == 7]
CAG_n8_h = otu[,CAG_n8 == 8]

CAG_n8_matrix = data.frame(a = rowSums(CAG_n8_a), b = rowSums(CAG_n8_b),
                           c = rowSums(CAG_n8_c), d = rowSums(CAG_n8_d),
                           e = rowSums(CAG_n8_e), f = rowSums(CAG_n8_f),
                           g = rowSums(CAG_n8_g), h = rowSums(CAG_n8_h))
#                           e = rowSums(CAG_n3_e), f = rowSums(CAG_n3_f),
#                           g = rowSums(CAG_n3_g))
# rename the CAG from 1 to 7
cag.names <- paste("CAG", seq(1,8), sep = "_")
names(CAG_n8_matrix) <- cag.names
#temp_nm5 <- "cecum"
# output the CAG table in ileum
write.csv(CAG_n8_matrix, paste("../otu/CAG2/", temp_nm4, "otu_CAG.csv", sep = "_"),
          row.names = TRUE)


############# 第二部分，基于OTU的CAG和基于代谢物的WGCNA关联分析 ###############

# load dataframe
#rm(list = ls())
# load libraries
library(psych)
library(dplyr)
library(reshape2)
library(tidyr)
temp_nm1 <- "../otu/CAG2/"
temp_nm2 <- "../metabolites_orig/WGCNA2/"
temp_nm3 <- "_fecal"

OTU.cag <- read.csv(paste(temp_nm1, temp_nm3 ,"_otu_CAG.csv", sep = ""),
                    header = T, row.names = 1, check.names = F)
metabo.WGCNA <- read.csv(paste(temp_nm2, "fecal_metab_meta_module_54.csv", sep = ""),
                         header = T, row.names = 1, check.names = F)

corr_df <- corr.test(OTU.cag, metabo.WGCNA, method = "spearman", adjust = "fdr")

# extract correlation coefficients 
corr_df_cor <- corr_df$r
corr_df_p <- corr_df$p

corr_df_cor <- data.frame(row=rownames(corr_df_cor),corr_df_cor,check.names = F) # create a column called "row" 
rownames(corr_df_cor) <- NULL
corr_df_p <- data.frame(row=rownames(corr_df_p),corr_df_p,check.names = F) # create a column called "row" 
rownames(corr_df_p) <- NULL
# Melt
nbar.m <- melt(corr_df_cor)
nbap.m <- melt(corr_df_p)

nbar.m$value2<-cut(nbar.m$value,breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),include.lowest=TRUE, label=c("(-0.75,-1)","(-0.5,-0.75)","(-0.25,-0.5)","(0,-0.25)","(0,0.25)","(0.25,0.5)","(0.5,0.75)","(0.75,1)")) # the label for the legend
nbap.m$value2<-cut(nbap.m$value,breaks=c(-Inf, 0.001, 0.01, 0.05),label=c("***", "**", "*")) 
nbar.m<-cbind.data.frame(nbar.m,nbap.m$value,nbap.m$value2) # adding the p value and its cut to the first dataset of R coefficients
names(nbar.m)[5]<-paste("valuep") # change the column names of the dataframe 
names(nbar.m)[6]<-paste("signif.")
names(nbar.m)[1:3] <- c("CAG", "WGCNA", "Coef.")
#dir.create("./module_accoc")
#temp_nm6 <- "./module_assoc/"

# remove the rows with NA in column sifnif.
nbar.m.nona <- nbar.m %>% drop_na(signif.)
write.csv(nbar.m.nona, paste("./module_assoc2/", temp_nm3, "_CAG_assoc_WGCNA.csv", sep = ""),row.names = F)


