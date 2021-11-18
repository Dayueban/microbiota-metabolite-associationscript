## part Ⅰ 过滤掉没有关联的结果，也就是Threshold为NA的

# This script was used for remove the metabolites which final_p_value presented by NA
rm(list=ls())
# 该代码需要加载stringr包
library(stringr)
options(stringsAsFactors = FALSE)

## 可更改变量 ##########
temp_nm1 <- "ileum_species"
temp_nm2 <- "filtered"         # origin或者filtered
########################

# 读取数据
df <- read.csv(paste("species/ileum/",temp_nm1, "_asso_metabo_merge_", temp_nm2, ".csv",
                     sep = ""), header=T, row.names = 1, check.names=F)
dim(df)
#rownames_df <- ileum_T[,1]
#ileum_df2 <- gut_df[,1:100]
total_table <- NULL
table2 <- data.frame(matrix(nrow = dim(df)[1])) # to construct one empty dataframe with number of rows as the dataframe df

#change the table2 colname with the first colname in original dataframe
colnames(table2) <- colnames(df)[1]
j=0
for(i in seq(10,dim(df)[2],10)){
  if(!is.na(df[,i][1])){
    j=j+1
    table2[,(j*10-9):(j*10)] <- df[,(i-9):i]
    #print(i)
  }
  #total_table <- cbind(total_table,table2)
}
rownames(table2) <- rownames(df)
# to correct the first colnames of table2
# table2数据框的第一列和第二列最后的Estimate_q和Estimate_b都是10个字符，
# 反向选择其位置获得相应的mz值
names(table2)[1] <- paste(str_sub(names(table2)[2],1,-11),
                          str_sub(names(table2)[1],-10,nchar(names(table2)[1]))
                          ,sep = "")
dim(table2)
write.csv(table2, paste("species/ileum/", temp_nm1, "_asso_metabo_merge_", 
                        temp_nm2, "_NA_rm.csv", sep = ""))

#######################################################################
#######################################################################

## part Ⅱ 计算关联结果&与每个代谢物显著关联的菌是什么

# 计算菌和代谢物mz值之间总共有多少个关联结果
rm(list = ls())
temp_nm1 <- "fecal_phylum"
temp_nm2 <- "origin"
df <- read.csv(paste("origin/00result_na_rm/fecal/",temp_nm1, "_asso_metabo_merge_", temp_nm2, "_NA_rm.csv",
                     sep = ""), header=T, row.names = 1, check.names=F)
asso_count = 0
table2 <- NULL
for(j in seq(10,dim(df)[2],10)){
  asso_count = asso_count + sum(df[,j-1] <= df[j][1])
}
print(asso_count)  # 打印出总关联结果

# calculating the frequency of genera associated with metabolites
library(qpcR)
library(stringr)
all_table <- NULL
colnames_index_all <- c()
#table2 <- data.frame(matrix(nrow = nrow(df)))
#i = 0
for(j in seq(10,dim(df)[2],10)){
  #if(!is.na(df[,j][1])){
  #i=i+1  
  table3 <- which(df[,j-1] <= df[,j][1])
  rownames_index <- rownames(df)[table3]
  colnames_index <- str_sub(names(df)[j-1], 1, -15)
  colnames_index_all <- c(colnames_index_all, colnames_index)
  all_table <- qpcR:::cbind.na(all_table, rownames_index)
}
# 先删除第一列,无用信息
all_table <- as.data.frame(all_table[,2:ncol(all_table)])
# 将有显著性结果的mz值作为列名赋值给all_table表
names(all_table) <- colnames_index_all
write.csv(all_table,paste("origin/00result_na_rm/fecal/", temp_nm1, "_asso_metabo_merge_", 
                          temp_nm2, "_NA_rm_each_micro_metabo.csv", sep = ""))

