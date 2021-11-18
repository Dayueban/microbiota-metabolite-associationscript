# sweep out everything
rm(list = ls())
# load libraries
library(forestplot)
library(plyr)
library(psych)
library(ggcorrplot)
library(pheatmap)
library(reshape2)
library(forestplot)
library(plyr)
library(psych)
#install.packages("ggcorrplot")
library(ggcorrplot)
library(pheatmap)
library(RNOmni) # 该包可以进行rank-based inverse transformation
library(Hmisc)
library(glmnet)
library(car)
library(glmnet)
#install.packages("caret")
library(caret)

# temp names of some variables
temp_nm1 <- "../otu/fecal/"
temp_nm2 <- "fecal"

# require the dataset
alpha.df <- read.csv(paste(temp_nm1, temp_nm2,"_alpha_diversity.csv",sep=""),
                     check.names=FALSE, header = T, row.names = 1)
metabo <- read.csv(paste(temp_nm1, temp_nm2,"_serum_metabolites_merge_origin.csv", sep = ""),
                   check.names=FALSE, header = T, row.names = 1)

# metabolites after screening
metabo2 <- read.csv("G:\\20200810-microbiota-metabolites-assoc\\taxa\\variance_explained\\fecal_otu_identified_metabolites.csv",
                    check.names = FALSE, header = T, row.names = 1)

# 表型校正
residual_all <- NULL
for(j in 3:length(metabo2)){               ####从第3列开始才是表型前5列为样品ID和固定效应
  fix<-lm(metabo2[,j] ~ as.factor(metabo2$slaughter_batch) + metabo2$sex)
  residual <- fix$residual
  residual_all <- cbind(residual_all, residual)
}
residual_all <- as.data.frame(residual_all)
names(residual_all) <- names(metabo2)[3:ncol(metabo2)]
rownames(residual_all) <- rownames(metabo2)



# match the id of two data.frame
ID.match <- intersect(rownames(alpha.df), rownames(residual_all))
# re-assignment the id to alpha.df and metabo
alpha.df <- alpha.df[ID.match, ]
residual_all <- residual_all[ID.match, ]

# combine alpha.df and metabo
fecal.residual_all.alpha <- cbind(residual_all, alpha.df)
#write.csv(fecal.residual_all.alpha, paste(temp_nm1, temp_nm2, "_serum_metabolites2_merge_origin-Mshannon.csv", sep = ""), row.names = TRUE)
dim(fecal.residual_all.alpha)

# perform Lasso regression
a <- fecal.residual_all.alpha
#a$slaughter_batch <- NULL
#a$sex <- NULL
set.seed(123)
beta_table <- NULL
R2_table <- NULL
R22_table <- NULL
pre_table <- NULL
for (i in 1:10) {
  idx <- createDataPartition(a$chao, times = 1, p = 0.9, list = FALSE)
  # 522 113
  X <- a[,1:(ncol(a)-3)]
  Y <- a[,ncol(a)-1] # shannon first, then chao, sobs
  # train dataset
  # scale and standardize metabolites
  X <- scale(X)
  X.train <- X[idx, ]
  X.test <- X[-idx,]
  # test dataset
  Y.train <- Y[idx]
  Y.test <- Y[-idx]
  cv.out <- cv.glmnet(X.train, Y.train, alpha = 1, family = "gaussian")
  Lasso_mod <- glmnet(X.train, Y.train, family = "gaussian", lambda = cv.out$lambda.min, alpha = 1) # the data of X should be matrix
  #best.lambda <- cv.out$lambda.min
  co <- coef(Lasso_mod)
  #co <- as.data.frame(as.matrix(co))
  #co 
  beta_table <- cbind(beta_table, co)
  pre <- Lasso_mod %>% predict(newx = X.test) # 预测值，也就是mShannon
  # 根据文章“A heritable  subset of the core rumen microbiome dictates dairy cow productivity and emissions”
  # 的方法
  R2 <- 1 - cv.out$cvm[which(cv.out$glmnet.fit$lambda == cv.out$lambda.min)] / var(Y.train)
  # 参考：https://www.r-bloggers.com/how-and-when-ridge-regression-with-glmnet/
  # 另外的方法
  R22 <- Lasso_mod$dev.ratio[which(Lasso_mod$lambda == cv.out$lambda.min)]
  R2_table <- cbind(R2_table, R2)
  R22_table <- cbind(R22_table, R22)
  pre_table <- cbind(pre_table, pre)
}
beta_table <- as.data.frame(as.matrix(beta_table)) #将dgCMatrix转为data.frame
# 删除Intercept行
beta_table1 <- beta_table[-1,]
incl <- apply(beta_table1, 1, function(x) length(x[x == 0]))
keep <- names(incl [incl < 10])
# keep the metabolites that with non-zero β-coefficient in at leat one of the ten LASSO models
beta_table2 <- beta_table1[keep,]
dim(beta_table2)
# 116 10

# keep the metabolites that with non-zero β-coefficient in ten LASSO model
keep2 <- names(incl [incl == 0])
beta_table3 <- beta_table1[keep2,]
dim(beta_table3)
# 5 10

# mean the R2 score
mean_R2 <- mean(R2_table)
# 0.066
mean_R22 <- mean(R22_table)
# 0.1036286
std_R2 <- sd(R2_table)
# 0.016


# 计算实际shannon和预测shannon之间的相关系数，作散点图
library(ggpubr)
shannon_mshannon <- as.data.frame(cbind(Y.test, apply(pre_table, 1, median)))
names(shannon_mshannon) <- c("Observed Shannon diversity", "mShannon diversity")
p <- ggscatter(shannon_mshannon, x = "mShannon diversity", y = "Observed Shannon diversity", add = "reg.line", 
          add.params = list(color = "blue")
)
ggsave("fecal_alpha_Mshannon.png", width = 8, height = 6, dpi = 300)

#p <- p + stat_cor(method = "pearson")

# 计算实际chao和预测chao之间的相关系数，作散点图
library(ggpubr)
shannon_mshannon <- as.data.frame(cbind(Y.test, apply(pre_table, 1, median)))
names(shannon_mshannon) <- c("Observed chao diversity", "mchao diversity")
p <- ggscatter(shannon_mshannon, x = "mchao diversity", y = "Observed chao diversity", add = "reg.line", 
               add.params = list(color = "blue")
)
ggsave("fecal_alpha_Mshannon.png", width = 8, height = 6, dpi = 300)










#' 不校正表型
# require the dataset
alpha.df <- read.csv(paste(temp_nm1, temp_nm2,"_alpha_diversity.csv",sep=""),
                     check.names=FALSE, header = T, row.names = 1)
metabo <- read.csv(paste(temp_nm1, temp_nm2,"_serum_metabolites_merge_origin.csv", sep = ""),
                   check.names=FALSE, header = T, row.names = 1)

# metabolites after screening
metabo2 <- read.csv("G:\\20200810-microbiota-metabolites-assoc\\taxa\\variance_explained\\fecal_otu_identified_metabolites.csv",
                    check.names = FALSE, header = T, row.names = 1)

# match the id of two data.frame
ID.match <- intersect(rownames(alpha.df), rownames(residual_all))
# re-assignment the id to alpha.df and metabo
alpha.df <- alpha.df[ID.match, ]
metabo2 <- metabo2[ID.match, ]

# combine alpha.df and metabo
fecal.metabo2.alpha <- cbind(metabo2, alpha.df)
#write.csv(fecal.metabo.alpha, paste(temp_nm1, temp_nm2, "_serum_metabolites_merge_origin-Mshannon.csv", sep = ""),
#          row.names = TRUE)
dim(fecal.metabo2.alpha)

# perform Lasso regression
a <- fecal.metabo2.alpha
a$slaughter_batch <- NULL
a$sex <- NULL
set.seed(123)
beta_table <- NULL
R2_table <- NULL
R22_table <- NULL
pre_table <- NULL
for (i in 1:10) {
  idx <- createDataPartition(a$shannon, times = 1, p = 0.9, list = FALSE)
  # 522 113
  X <- a[,1:(ncol(a)-3)]
  Y <- a[,ncol(a)] # shannon first, then chao, sobs
  # train dataset
  # scale and standardize metabolites
  X <- scale(X)
  X.train <- X[idx, ]
  X.test <- X[-idx,]
  # test dataset
  Y.train <- Y[idx]
  Y.test <- Y[-idx]
  cv.out <- cv.glmnet(X.train, Y.train, alpha = 1, family = "gaussian")
  Lasso_mod <- glmnet(X.train, Y.train, family = "gaussian", lambda = cv.out$lambda.min, alpha = 1) # the data of X should be matrix
  #best.lambda <- cv.out$lambda.min
  co <- coef(Lasso_mod)
  #co <- as.data.frame(as.matrix(co))
  #co 
  beta_table <- cbind(beta_table, co)
  pre <- Lasso_mod %>% predict(newx = X.test) # 预测值，也就是mShannon
  # 根据文章“A heritable  subset of the core rumen microbiome dictates dairy cow productivity and emissions”
  # 的方法
  R2 <- 1 - cv.out$cvm[which(cv.out$glmnet.fit$lambda == cv.out$lambda.min)] / var(Y.train)
  # 参考：https://www.r-bloggers.com/how-and-when-ridge-regression-with-glmnet/
  # 另外的方法
  R22 <- Lasso_mod$dev.ratio[which(Lasso_mod$lambda == cv.out$lambda.min)]
  R2_table <- cbind(R2_table, R2)
  R22_table <- cbind(R22_table, R22)
  pre_table <- cbind(pre_table, pre)
}
beta_table <- as.data.frame(as.matrix(beta_table)) #将dgCMatrix转为data.frame
# 删除Intercept行
beta_table1 <- beta_table[-1,]
incl <- apply(beta_table1, 1, function(x) length(x[x == 0]))
keep <- names(incl [incl < 10])
# keep the metabolites that with non-zero β-coefficient in at leat one of the ten LASSO models
beta_table2 <- beta_table1[keep,]
dim(beta_table2)
# 116 10

# keep the metabolites that with non-zero β-coefficient in ten LASSO model
keep2 <- names(incl [incl == 0])
beta_table3 <- beta_table1[keep2,]
dim(beta_table3)
# 5 10

# mean the R2 score
mean_R2 <- mean(R2_table)
# 0.066
mean_R22 <- mean(R22_table)
# 0.1036286
std_R2 <- sd(R2_table)
# 0.016


# 计算实际shannon和预测shannon之间的相关系数，作散点图
library(ggpubr)
shannon_mshannon <- as.data.frame(cbind(Y.test, apply(pre_table, 1, median)))
names(shannon_mshannon) <- c("Observed Shannon diversity", "mShannon diversity")
p <- ggscatter(shannon_mshannon, x = "mShannon diversity", y = "Observed Shannon diversity", add = "reg.line", 
               add.params = list(color = "blue")
)
ggsave("fecal_alpha_Mshannon.png", width = 8, height = 6, dpi = 300)

#p <- p + stat_cor(method = "pearson")

# 计算实际chao和预测chao之间的相关系数，作散点图
library(ggpubr)
shannon_mshannon <- as.data.frame(cbind(Y.test, apply(pre_table, 1, median)))
names(shannon_mshannon) <- c("Observed chao diversity", "mchao diversity")
p <- ggscatter(shannon_mshannon, x = "mchao diversity", y = "Observed chao diversity", add = "reg.line", 
               add.params = list(color = "blue")
)
ggsave("fecal_alpha_Mshannon.png", width = 8, height = 6, dpi = 300)

