# for manipulation the results of association between otu cags and metabolite WGCNA
# Version 2.0
# sweep out everything
rm(list = ls())
# load libraries
library(psych)
library(dplyr)
library(reshape2)
library(tidyr)
temp_nm1 <- "./module_assoc2/"
temp_nm2 <- "../metabolites_orig/WGCNA2/"
temp_nm3 <- "_cecum"
temp_nm4 <- "../otu/CAG2/"
temp_nm5 <- "./module_assoc2/manipulation/"

cag.WGCNA <- read.csv(paste(temp_nm1, temp_nm3 ,"_CAG_assoc_WGCNA.csv", sep = ""),
                      header = T, check.names = F)
WGCNA.assign <- read.csv(paste(temp_nm2, "cecum_metab_modulelog+p20bicor-132.csv", sep = ""),
                         header = T, check.names = F)
cag.assign <- read.csv(paste(temp_nm4, temp_nm3, "_otu_cluster_8_cag.csv", sep = ""),
                       header = T, row.names = 1, check.names = F)
#names(cag.assign)[1] <- "CAG"

# 增加一列，给cluster赋值加上CAG字样
cag.assign$CAG <- paste(rep("CAG", nrow(cag.assign)), cag.assign$otu.cluster, sep = "_")

# unique的CAG和module分别统计如下：
CAG_uniq <- unique(cag.WGCNA$CAG)
WGCNA_uniq <- unique(cag.WGCNA$WGCNA)
# 由于WGCNA.assign数据框中的mz value之间有空格，影响后面取列，因此重命名
names(WGCNA.assign)[2] <- "mz_value"

# 1)找到归属以某个CAG旗下的所有OTU，以及其注释分类是什么
cag2otu <- NULL
table2 <- NULL
table3 <- NULL
for (i in 1:length(CAG_uniq)) {
  cag2otu <- rownames(cag.assign)[which(CAG_uniq[i] == cag.assign$CAG)]
  cag2otu <- paste(cag2otu, collapse = ";") # 向量内部的元素连接在一起，用分号分割
  table2 <- cbind(CAG_uniq[i], cag2otu)
  table3 <- as.data.frame(rbind(table3, table2))
  
}
names(table3)[1] <- "CAG"

# 2)找到归属于某个module旗下的所有mz值，以及其annotation结果是什么
WGCNA2mz <- NULL
WGCNA2anno <- NULL
table4 <- NULL
table5 <- NULL

for (i in 1:length(WGCNA_uniq)) {
  WGCNA2mz <- WGCNA.assign$mz_value[which(WGCNA_uniq[i] == WGCNA.assign$module)]
  WGCNA2anno <- WGCNA.assign$annotation[which(WGCNA_uniq[i] == WGCNA.assign$module)]
  #WGCNA2anno.nona <- na.omit(WGCNA2anno) # 20200916新增加代码
  #WGCNA2anno.nona.f <- as.vector(WGCNA2anno.nona) # 20200916新增加代码
  #WGCNA2mz <- paste(WGCNA2mz, collapse = ";") # 得到module对应的mz值
  #WGCNA2anno <- paste(WGCNA2anno.nona.f, collapse = ";") # 得到每个模块对应mz值的annotation
  table4 <- cbind(WGCNA_uniq[i], WGCNA2mz, WGCNA2anno)
  table5 <- as.data.frame(rbind(table5, table4))
  
}
names(table5)[1] <- "WGCNA"


# 3)返回到cag.WGCNA表格，将cag和WGCNA的分类结果分别写到后面
# 3.1)先做cag
cag.assin2 <- NULL
cag.otu <- NULL
for (j in 1:nrow(cag.WGCNA)) {
  cag.assin2 <- table3$cag2otu[table3$CAG == cag.WGCNA$CAG[j]] 
  cag.otu <- c(cag.otu, cag.assin2)
}
cag.WGCNA$cag2otu <- cag.otu # cag.WGCNA表格增加一列

# 3.2)再做WGCNA
WGCNA.assign2 <- NULL
WGCNA.mz <- NULL
WGCNA.assign3 <- NULL
WGCNA.anno <- NULL
for (j in 1:nrow(cag.WGCNA)) {
  WGCNA.assign2 <- table5$WGCNA2mz[table5$WGCNA == cag.WGCNA$WGCNA[j]] 
  WGCNA.assign3 <- table5$WGCNA2anno[table5$WGCNA == cag.WGCNA$WGCNA[j]]
  WGCNA.mz <- c(WGCNA.mz, WGCNA.assign2)
  WGCNA.anno <- c(WGCNA.anno, WGCNA.assign3)
  
}
cag.WGCNA$WGCNA2mz <- WGCNA.mz
cag.WGCNA$WGCNA2anno <- WGCNA.anno

# output the final file
write.csv(cag.WGCNA, paste(temp_nm5, temp_nm3, "_cag_WGCNA_assoc_assign_anno2.csv",
                           sep = ""), row.names = F)












