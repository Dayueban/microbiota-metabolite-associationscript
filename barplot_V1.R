rm(list = ls())
setwd("E:\\F6代血清代谢组学相关性分析结果\\20180327重新分析\\Paper_graph\\2_barplot")
meta_df <- read.table("three_gutregion_share_genus19.csv",header = T,sep = ',',check.names = F,row.names = 1)
meta_df<- sweep(meta_df,2,colSums(meta_df),'/')*100
meta_df <- as.matrix(meta_df) # barplot要求是矩阵或向量
library("RColorBrewer")
#col = brewer.pal(nrow(data),"Paired")
#colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628","#F781BF", "#999999", "#1B9E77", "#D95F02","#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#1A1A1A", "#67001F")
colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628","#F781BF", "#999999", "#1B9E77", "#D95F02","#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#66C2A5", "#FC8D62")
tiff("Genus_taxonomy_barplot0324.tiff",width=1280,height=1280,res=300)
par(mai=c(0.6,1,0.5,1.2))
bp = barplot(meta_df, beside=F, xlab="", col=colors, border=NA, xaxs="i", width = c(15,15,15),
             mgp=c(3,0.2,0.5), las=2, cex.names=1.2, axes=F, main = "Genus taxonomy",
             ylab="Relative abundance (%)",cex.lab=1.2,xaxt="n")
axis(2, las=2, cex.axis=1, mgp=c(3,0.4,0.5), tcl=-0.2)
text(bp, par("usr")[3]-0.4,labels = colnames(meta_df),xpd=TRUE,adj=1,srt=45)
abline(h=axTicks(2),lty=2,col=rgb(0,0,0,0.2))
legend(par("usr")[2],par("usr")[4]+1,cex=0.5,bty="n",x.intersp=0.4,pt.cex=1,rownames(meta_df),fill=colors,text.font=1,ncol=1,xpd=TRUE)
dev.off()

