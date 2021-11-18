rm(list = ls())
library(vegan)
library(regress)
library(gap)

temp_nm = "variance_explained/cecum"
pheno <- read.csv(paste(temp_nm,"otu_identified_metabolites.csv", sep = "_"), header=T,
                  check.names = F)
# ind <- sample(seq_len(nrow(pheno)), size=200)
otu <- read.csv(paste("../otu/cecum/cecum","otu_80.csv",sep="_"),header=T,sep = ",",
                check.names = F)

#pheno <- pheno[1:50,1:8] #just for test this code with small sample
#otu <- otu[1:50,1:30]



### 第二种方法


# 下面采用的是Fu的文章中的方法

for(i in 4:ncol(pheno)){
  mean <- mean(pheno[,i],na.rm=T)
  pheno[,i][which(is.na(pheno[,i]))] <- mean
}

Explained_variance<-function(pheno,otu,gradient){
  discovery_n<-round(nrow(pheno)*0.7,0)
  pheno_rank<-sample(1:nrow(pheno),discovery_n)
  discovery_set<-pheno[pheno_rank,]
  subset_otu<-otu[pheno_rank,]
  validation_set<-pheno[-(pheno_rank),]
  validation_otu<-otu[-(pheno_rank),]
  fix<-lm(discovery_set$traits~+discovery_set$slaughter_batch+discovery_set$sex)
  residual<-fix$residual
  table<-assoc(residual,subset_otu)
  if(any(table[,7]<=gradient)){
    mm<-which(table[,7]<=gradient)
    sub_table<-table[mm,]
    sig_otu<-validation_otu[c(names(mm))]
    r<-NULL
    for(i in 1:nrow(sig_otu)){
      regi<-sum(sig_otu[i,]*sub_table[1],na.rm=T)+sum(sub_table[2],na.rm=T)+sum(as.numeric(sig_otu[i,]>0),na.rm=T)
      r<-c(r,regi)
    }
    fix2<-lm(validation_set$traits~validation_set$slaughter_batch+validation_set$sex)
    residua2<-fix2$residual
    temp<-lm(residua2~r)
    temp1<-summary(temp)
    R2<-temp1$r.squared
    return(R2)
  }
}
assoc<-function(residual,OTU){
  all=NULL
  for(i in 2:length(OTU)){
    binary<-as.numeric(OTU[,i]>0)
    temp2<-lm(residual~binary)
    seat<-which(binary>0)
    residual2<-residual[c(seat)]
    present_otu<-OTU[,i][c(seat)]
    temp1<-lm(residual2~log10(present_otu)) 
    temp<-summary(temp1)
    summar<-summary(temp2)
    Estimate_q<-temp$coefficients[2]    
    Estimate_b<-summar$coefficients[2]  
    t_value_q<-temp$coefficients[6]     
    p_value_q<-temp$coefficients[8] 
    t_value_b<-summar$coefficients[6]
    p_value_b<-summar$coefficients[8]
    if(is.na(p_value_b)){
      p_meta<-NA
      Estimate_b<-NA
      meta_z_score<-abs(qnorm(p_value_q/2))*Estimate_q/abs(Estimate_q)
    }
    if(is.na(p_value_q)){
      p_meta<-NA
      Estimate_q<-NA
      meta_z_score<-abs(qnorm(p_value_b/2))*Estimate_b/abs(Estimate_b)
    }
    if(!(is.na(p_value_b))&&!(is.na(p_value_q))&& p_value_q<=p_value_b){
      meta_z_score<-(abs(qnorm(p_value_q/2))*Estimate_q/abs(Estimate_q)+abs(qnorm(p_value_b/2))*Estimate_b/abs(Estimate_b))/sqrt(2)            
      meta_z_score<-Estimate_b/abs(Estimate_b)*abs(meta_z_score)
      p_meta<-pnorm(-abs(meta_z_score))*2
    }
    if(!(is.na(p_value_b))&& !(is.na(p_value_q))&& p_value_q>p_value_b){
      meta_z_score<-(abs(qnorm(p_value_q/2))*Estimate_q/abs(Estimate_q)+abs(qnorm(p_value_b/2))*Estimate_b/abs(Estimate_b))/sqrt(2) 
      meta_z_score<-Estimate_q/abs(Estimate_q)*abs(meta_z_score)
      p_meta<-pnorm(-abs(meta_z_score))*2
    }
    final_p_value<-min(c(p_value_q,p_value_b,p_meta),na.rm=T)
    table<-cbind(Estimate_q,Estimate_b,p_value_q,p_value_b,meta_z_score,p_meta,final_p_value)
    all<-rbind(all,table)
  }
  rownames(all)<-colnames(OTU)[2:length(OTU)]
  return(all)
}

de<-c(0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1)
box_table <- data.frame(matrix(ncol = 9, nrow = 0))
#colnames(box_table) <- de

#随机挑取248个样本
ind <- sample(seq_len(nrow(pheno)), size=248)  ## sample函数，用于随机抽样
pheno_sub <- pheno[ind,]
otu_sub <- otu[ind,]

for(x in 4:ncol(pheno_sub)){
    mm<-cbind(pheno_sub[1:3],pheno_sub[,x])
    colnames(mm)<-c('ID','slaughter_batch','sex','traits')
    box<-NULL
    for(n in de){R1<-NULL
    for(j in 1:50){
      x1<-Explained_variance(mm,otu_sub,gradient=n)
      R1<-cbind(R1,x1)
    }
    exp1<-mean(R1,na.rm=T)
    box<-cbind(box,exp1)
    }
    rownames(box)<-colnames(pheno)[x]
    #colnames(box)<-de
    box_table <- rbind(box_table, box)
}

colnames(box_table) <- de
write.csv(box_table, paste(temp_nm, "/cecum_variance_explained20210402.csv",sep = ""))
