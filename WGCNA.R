rm(list = ls())
temp_nm1 <- "species/ileum/"
temp_nm2 <- "species/cecum/"
temp_nm3 <- "species/fecal/"
#temp_nm2 <- "filtered"
#library (xlsx) ### Saving to spreadsheet
library (data.table) ### Fast read of large files into R
library (WGCNA) ### -	Clustering software. Previously reported work done using v1.34
library (flashClust) ### Clustering software
library (ppcor) ### Partial Spearman correlations, for confounder analysis. Previously reported work done using v1.0
library (gplots) ### Plotting
library (cowplot) ### Plotting; to arrange several plots on the same page
library (ggplot2) ### Plotting
library (plyr) ### Data transformations
library(stringr) ### Mainly Dataframe/Character manipulation

# load working datasheet
options (stringsAsFactors = FALSE)
#work_df <- read.csv("Ileum_OTU.csv",header = T,row.names = 1,check.names = F)
a <- read.csv(paste(temp_nm1, "ileum_serum_metabolites_merge_origin.csv", sep = ""),
              header = T,row.names = 1,check.names = F)
b <- read.csv(paste(temp_nm2, "cecum_serum_metabolites_merge_origin.csv", sep = ""),
              header = T,row.names = 1,check.names = F)
c <- read.csv(paste(temp_nm3, "fecal_serum_metabolites_merge_origin.csv", sep = ""),
              header = T,row.names = 1,check.names = F)
overlap_index <- intersect(intersect(rownames(a), rownames(b)), rownames(c))
meta_df <- a[overlap_index,]
dat <- meta_df[,3:ncol(meta_df)]


### Settings for WGCNA generally
cor_method          = "spearman" ### for association with clinical parameters
corFun_tmp          = "bicor"
cluster_method      = "average"
corOptions_list     = list (use = 'pairwise.complete.obs') 
corOptions_str      = "use = 'pairwise.complete.obs'"
BH_pval_asso_cutoff = 0.05
NetworkType         = "signed" ### Signed-network (as the PC1 and median profile does not make sense as a summary measure of a cluster with anticorrelated metabolites.)



# substitute the zero by imputation with the method of knn
library(impute)
# replacing all zeroes into NA
dat[dat == 0] <- NA
sum1 <- apply(dat, 2, function(x)length(x[x == 0]))
dat.imp <- impute.knn(as.matrix(dat),k=73)
dat.imp <- as.data.frame(dat.imp$data)
sum2 <- apply(dat.imp, 2, function(x)length(x[x == 0]))
table(sum2)

# dat_tmp = log2 (dat.imp) ### work in logarithmic space
dat_tmp <- log2(dat + 0.00001)
### Settings for WGCNA on polar metabolite measurements
### Once these are established the steps below can be run
RsquareCut_val      = 0.89 ### usually ranges 0.80-0.95 but requires inspecting curves
mergingThresh       = 0.20 ### Maximum dissimilarity of module eigengenes (i.e. 1-correlation) for merging modules.
minModuleSize       = 8 ### minimum number of metabolites constituting a cluster
SoftPower           = 8 ### beta-value, main parameter to optimize

### Calculate weighted adjacency matrix
# Choose a set of soft-thresholding powers
set.seed(123)
powers = c(c(1:10), seq(from = 12, to=30, by=1))
# Call the network topology analysis function
#Analysis of scale free topology for multiple soft thresholding powers. 
#The aim is to help the user pick an appropriate soft-thresholding power for network construction
sft = pickSoftThreshold(dat_tmp,powerVector=powers,verbose = 5)



A = adjacency (dat_tmp, power = SoftPower, type = NetworkType, corFnc = corFun_tmp, corOptions = corOptions_str)
colnames (A) = rownames (A) = colnames (dat_tmp)
### Define dissimilarity based on topological overlap
dissTOM = TOMdist (A, TOMType = NetworkType)
colnames (dissTOM) = rownames (dissTOM) = colnames (dat_tmp)

### Hierarchical clustering
metaTree = flashClust (as.dist (dissTOM), method = cluster_method)
### Define modules by cutting branches
moduleLabels1 = cutreeDynamic (dendro = metaTree, distM = dissTOM, method = "hybrid", deepSplit = 4, pamRespectsDendro = T, minClusterSize = minModuleSize)
moduleLabels1 = labels2colors (moduleLabels1)
### Automatically merge highly correlated modules
merge = mergeCloseModules (dat_tmp, moduleLabels1, corFnc = corFun_tmp, corOptions = corOptions_list, cutHeight = mergingThresh)
### Determine resulting merged module colors
moduleLabels2 = merge$colors
table(moduleLabels2)
length(table(moduleLabels2))
### Establish eigengenes of the newly merged modules, used for cluster overall abundances
MEs = merge$newMEs
### Choose final module assignments
moduleColorsMeta = moduleLabels2
names (moduleColorsMeta) = colnames (dat_tmp)
MEsMeta = orderMEs (MEs)
rownames (MEsMeta) = rownames (dat_tmp)

### Determine relevant descriptive statistics of established clusters
### kIN: within-module connectivity, determined by summing connectivity with all
###      other metabolites in the given cluster.
### kME: bicor-correlation between the metabolite profile and module eigenvector; 
### both measures of intramodular hub-metabolite status.
kIN <-      vector (length = ncol (dat_tmp)); names (kIN) = colnames (dat_tmp)
kME <-      vector (length = ncol (dat_tmp)); names (kME) = colnames (dat_tmp)
modules <-  vector (length = ncol (dat_tmp)); names (modules) = colnames (dat_tmp)

# cluster_mapping_file = read.table (file = "cluster_mapping_file.tab", header = T, sep = "\t", row.names = 1)
# cluster_mapping_file$label =  sapply (rownames (cluster_mapping_file),  function (x) paste (cluster_mapping_file [x, "New_Name"], cluster_mapping_file [x, "Description"], sep = ": "))
cluster_mapping_file <- cbind(names (table (moduleColorsMeta)), paste("M",formatC(seq("01",length(
  names (table (moduleColorsMeta))),by = 1),width=2, flag=0), sep = "_"))
cluster_mapping_file <- as.data.frame(cluster_mapping_file)
names(cluster_mapping_file) <- c("Module_Name", "New_Name")
grey_index <- which(cluster_mapping_file$Module_Name == "grey")
cluster_mapping_file <- cluster_mapping_file[-grey_index,]
cluster_mapping_file$New_Name <- paste("M",formatC(seq("01",nrow(cluster_mapping_file),by = 1),width=2, flag=0), sep = "_")
grey_module <- c("grey", "M_remaining")
cluster_mapping_file <- rbind(grey_module, cluster_mapping_file)
cluster_mapping_file$V2[which(cluster_mapping_file$V1 == "grey")] <- "M_remaining"


for (module in names (table (moduleColorsMeta))) {   
  
  all.metabolites = names (dat_tmp)
  inModule = (moduleColorsMeta == module)
  module.metabolites = names (moduleColorsMeta [inModule])
  modules [module.metabolites] = module 
  kIN [module.metabolites] = sapply (module.metabolites, function (x) sum (A [x, module.metabolites]) - 1)
  datKME = signedKME (dat_tmp, MEsMeta, corFnc = corFun_tmp, corOptions = corOptions_str)
  rownames (datKME) = colnames (dat_tmp)
  kME [module.metabolites] = datKME [module.metabolites, paste ("kME", module, sep = "")]   
  
}
# write.table(names (table (moduleColorsMeta)), "cfi_commonidv_metabo_module.tab",quote = FALSE,row.names = FALSE,
#            col.names = FALSE)
output = data.frame ("module" = modules, "kME" = kME, "kIN" = kIN, "cluster_name" = sapply (modules, function (m) cluster_mapping_file [paste0 ("M_ME", m), "New_Name"]))
















#ID_common_all_sample <- intersect(rownames(MEsMeta),rownames(meta_df))

# 为了将构建好的聚类颜色分别赋值给每个OTU
cluster_metabo <- paste(rep("ME",length(moduleColorsMeta)),as.vector(moduleColorsMeta),sep = "")
metabo_cluster <- rbind(meta_df,cluster_metabo)
rownames(metabo_cluster)[nrow(metabo_cluster)] <- "ME_cluster"
write.csv(metabo_cluster,paste(temp_nm1,"species_metabo_cluster_sof10_V1.csv"))
## Link OTU clusters to metabolites of interest
cor_OTU.metabolite <- list () ### Data structure for storing results of correlation tests under different setups

tmpMat = array (NA, c (ncol (MEsMeta), 26, 52))
dimnames (tmpMat) [[1]] = names (MEsMeta)
dimnames (tmpMat) [[2]] = names (meta_df[,-c(1,2)]) 
dimnames (tmpMat) [[3]] = paste(c ("estimate", "p.value"),rep(names (meta_df[,-c(1,2)]),each=2),sep="_")

### Associating metabolite clusters with HOMA-IR without de-confounding for BMI    
for (i in 3:ncol(meta_df)){
  tmpMat [, names(meta_df)[i], c (paste(c("estimate", "p.value"),c(names(meta_df)[i],names(meta_df)[i]),sep = "_"))] =
    t (apply (MEsMeta [ID_common_all_sample,], MARGIN = 2, FUN = function (x) 
      unlist (cor.test (meta_df [ID_common_all_sample, names(meta_df)[i]], x,
                        method = cor_method) [c ("estimate", "p.value")])))       
}

# transform the array to dataframe
data_df <- as.data.frame.table(tmpMat)

#data_df[which(data_df$Freq == NA),] <- NULL
# remove the rows with NA
data_df <- data_df[complete.cases(data_df), ]
write.csv(data_df,"meta_micro_bioassoc_sof6_V2.csv")

### Associating metabolite clusters with HOMA-IR while de-confounding for sex and batch

tmpMat_adj = array (NA, c (ncol (MEsMeta), 26, 52))
dimnames (tmpMat_adj) [[1]] = names (MEsMeta)
dimnames (tmpMat_adj) [[2]] = names (meta_df[,-c(1,2)]) 
dimnames (tmpMat_adj) [[3]] = paste(c ("estimate", "p.value"),rep(names (meta_df[,-c(1,2)]),each=2),sep="_")

for (i in 3:ncol(meta_df)){
  tmpMat_adj [, names(meta_df)[i], c (paste(c("estimate", "p.value"),c(names(meta_df)[i],names(meta_df)[i]),sep = "_"))] =
    t (apply (MEsMeta [ID_common_all_sample,], MARGIN = 2, FUN = function (x) 
      unlist (pcor.test (x = x, y = meta_df [ID_common_all_sample, names(meta_df)[i]], z = meta_df [ID_common_all_sample, c("batch","sex")],
                         method = cor_method) [c ("estimate", "p.value")])))       
}

# transform the array to dataframe
data_df_adj <- as.data.frame.table(tmpMat_adj)

#data_df[which(data_df$Freq == NA),] <- NULL
# remove the rows with NA
data_df_adj <- data_df_adj[complete.cases(data_df_adj), ]
write.csv(data_df_adj,"meta_micro_bioassoc_sof6_adj_V2.csv")




##--------------- part two plot the WGCNA for associations between ------------##
##---------------      OTUs cluster and metabolites                ------------##
rm(list = ls())
vari_name = "Fecal"
path = "E:/F6代血清代谢组学相关性分析结果/20180327重新分析/merge/OTUs"
setwd(paste(path,"/",vari_name,sep = ""))
# load the dataframe
asso_df <- read.csv("meta_micro_bioassoc_sof16_adj_V2.csv",header = T)
asso_df <- asso_df[order(asso_df$Estimate,decreasing = T),]
asso_df$value1<-cut(asso_df$Estimate,breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),include.lowest=TRUE, label=c("(-0.75,-1)","(-0.5,-0.75)","(-0.25,-0.5)","(0,-0.25)","(0,0.25)","(0.25,0.5)","(0.5,0.75)","(0.75,1)")) # the label for the legend
asso_df$value2<-cut(asso_df$pvalue,breaks=c(-Inf, 0.001, 0.01, 0.05),label=c("***", " ** ", "*")) 
names(asso_df)[5]<-paste("valuep") # change the column names of the dataframe 
names(asso_df)[6]<-paste("signif.")
pa <- ggplot(asso_df, aes(ME_cluster, Metabolite)) +
  geom_tile(aes(fill=Estimate),colour="white") +
  #scale_fill_brewer(palette = "RdYlGn",name="Correlation")# RColorBrewer package
  scale_fill_gradient2(low="blue3", high="red3", guide="colorbar",name="Correlation coefficient") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2,color="black",size=20),
        axis.text.y=element_text(size = 20, color = "black"),
        axis.ticks = element_blank(),
        axis.title.y=element_text(size = 25,color = "black"),axis.title.x=element_text(size = 25,colour = "black"))+
  #labs(title="Spearman correlation") + 
  xlab("OTU cluster") +
  ylab("Serum metabolites") +
  theme(plot.title = element_text(size = 20,color = "black")) +
  theme(axis.line=element_line(colour = "white"),
        #panel.border = element_rect(colour = "black", fill=NA, size=1,linetype = 1),
        legend.key.height=unit(0.8,"cm"),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(lineheight=0.8,size=16),
        legend.title=element_text(size=20))
# Adding the significance level stars using geom_text 
pp<- pa +
  geom_text(aes(label=signif.),size=6,na.rm=TRUE,nudge_y = -0.15) # nudge_y调节星号在方框的位置
ggsave(paste("WGCNA_",vari_name,"_adj_V2.tiff",sep = ""),width = 14,height = 30)




