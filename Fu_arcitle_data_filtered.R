rm(list = ls())
setwd("E:\\F6代血清代谢组学相关性分析结果\\20200810-microbiota-metabolites-assoc\\taxa\\species\\ileum")
a <- read.csv("ileum_species_80.csv", header = T, check.names = F,
              row.names = 1)
b <- read.csv("ileum_serum_metabolites_merge_filtered.csv",check.names = F,
              row.names = 1)
## 统计几个统计量
a1 <- as.data.frame(t(a))
sum1 <- apply(a1, 1, sum)
sum2 <- apply(a1, 1, function(x)length(x[x == 0]))
sumall <- sum(sum1)
num_sample <- ncol(a1)

## 

a1$ratio1 <- (sum1 / sumall) * 100
a1$ratio1 <- paste(sprintf("%.2f", a1$ratio1), "%", sep = "")
a1$ratio2 <- (sum2 / num_sample)
a1$ratio2 <- paste(sprintf("%.2f", a1$ratio2), "%", sep = "")

## to exclude family which with abundance less than 0.05% and present 
## in at least 20% of the study cohort
keep <- rownames(a1)[intersect(which(a1$ratio1 > 0.05),which(a1$ratio2 < 0.8))]

## 保留keep变量里面的那些family
a2 <- a1[keep, ]
a2$ratio1 <- NULL
a2$ratio2 <- NULL

a3 <- as.data.frame(t(a2))

# 菌群数据和代谢物数据ID取交集
sample_name <- intersect(rownames(a), rownames(b))
length(sample_name)
b1 <- b[sample_name,]
a1 <- a[sample_name,]

#dir.create("cecum")
write.csv(a1,"ileum_species_80.csv")
write.csv(b1,"ileum_serum_metabolites_merge_filtered.csv")

