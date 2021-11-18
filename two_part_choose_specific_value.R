rm(list = ls())
library(stringr)
temp_nm1 <- "cecum_family"
temp_nm2 <- "filtered"
df <- read.csv(paste("family/cecum/",temp_nm1, "_asso_metabo_merge_", temp_nm2, "_NA_rm.csv",
                     sep = ""), header=T, row.names = 1, check.names=F)
dim(df)
## 主诉：如何将显著性结果给挑选出来，以及根据P值来源确定“相关系数”为binary model源、quantitive model源、还是meta Z score

# calculating the frequency of genera associated with metabolites
library(qpcR)
library(stringr)
all_table <- NULL
taxa_all <- c()
metabo_all <- c()
Final_pvalue_all <- c()
Coref_score_all <- c()
Coref_type_all <- c()
Coref_type <- NULL
Coref_score <- NULL
for(j in seq(10,dim(df)[2],10)){
  table3 <- which(df[,j-1] <= df[,j][1])
  #metabo <- str_sub(names(df)[j-1], 1, -15)
  for(i in 1:length(table3)) {
    taxa <- rownames(df)[table3[i]]
    Final_pvalue <- df[table3[i], j-1]
    metabo <- str_sub(names(df)[j-1], 1, -15)
    if(df[table3[i], j-1] %in% df[table3[i],j-2]) {
      Coref_score <- df[table3[i], j-3]
      Coref_type <- "Meta"
    } else {
      if(df[table3[i], j-1] %in% df[table3[i],j-4]) {
        Coref_score <- df[table3[i], j-5]
        Coref_type <- "binary"
        
      } else {
        if(df[table3[i], j-1] %in% df[table3[i],j-6]) {
          Coref_score <- df[table3[i], j-7]
          Coref_type <- "quantitive"
        }
      }
    }
    taxa_all <- c(taxa_all, taxa)
    metabo_all <- c(metabo_all, metabo)
    Final_pvalue_all <- c(Final_pvalue_all, Final_pvalue)
    Coref_score_all <- c(Coref_score_all, Coref_score)
    Coref_type_all <- c(Coref_type_all, Coref_type)
    #all_table <- qpcR:::cbind.na(all_table, taxa_all, metabo_all, Final_pvalue_all,
    #                             Coref_score_all, Coref_type_all)
  }
  #cat(j, "\n") # 查看循环是否正常
}
all_table <- cbind(taxa_all, metabo_all, Final_pvalue_all,
                   Coref_score_all, Coref_type_all)
all_table <- as.data.frame(all_table)
names(all_table) <- c("taxa", "metabo_mz", "Final_Pvalue", "Corr_score", "Corr_type")
write.csv(all_table, paste("family/cecum/",temp_nm1, "_asso_metabo_merge_", temp_nm2, "_NA_rm_trimmed.csv",
                           sep = ""), row.names = F)
