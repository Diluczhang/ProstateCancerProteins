
# 对处理后的数据进行NMF分型

rm(list = ls())
options(stringsAsFactors = FALSE)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))



NMF = function(path){
  
  setwd(path)
  
  
  # nmf分型的代码
  library(NMF)
  library(doParallel)
  library(foreach)
  library(dplyr)
  load("1.rawdata/7.Protein.50.log2.Combat.T.Top25.RData")
 
  # 启用parallel作为foreach并行计算的后端
  cl <- makeCluster(8)
  registerDoParallel(cl)
  
  
  protein.nmf <-function(x){
    result<-nmf(x,rank = 3,nrun = 200)
    return(result)
  }
  
  
  
  
  Protein.50.log2.Combat.T.Top25.nmf <- list()
  for(i in c(1:100)){
    print(paste("Protein.50.log2.Combat.T.Top25.nmf",i,sep="_"))
    Protein.50.log2.Combat.T.Top25.nmf[[i]] <-nmf(Protein.50.log2.Combat.T.Top25,rank = 3,nrun = 200)
  }
  print("Protein.50.log2.Combat.T.Top25 finish")
  save(Protein.50.log2.Combat.T.Top25.nmf,file = "Protein.50.log2.Combat.T.Top25.nmf.RData")

  stopCluster(cl)   
  
  
}

path <- getwd()
NMF(path)




Group50Twoway <- function(NMFResult,filename){
  
  
  Group100 <-predict(NMFResult[[100]])
  Group100<-data.frame(Sample=names(Group100),Group=Group100)
  write.xlsx(Group100,file=paste0(filename,"第100次分组结果.xlsx"))
  
  allGroup<-names(Group100)
  for (i in 1:100) {
    group.tmp<- predict(NMFResult[[i]])
    allGroup<- cbind(allGroup,group.tmp)
  }
  colnames(allGroup)<-c("sample",seq(1:100))
  write.xlsx(allGroup,file=paste0(filename,"NMF100次重复结果.xlsx"))
  
  
  
  for(j in c(1:100)){
    protein.T.nmf.100 = matrix(0,dim(NMFResult[[j]]@consensus)[1],dim(NMFResult[[j]]@consensus)[2])
    mx <- NMFResult[[j]]@consensus
    protein.T.nmf.100 <- protein.T.nmf.100+mx
  }
  protein.T.nmf.100 <- protein.T.nmf.100/100

  #kmean分类预测
  km <- kmeans(protein.T.nmf.100,3,nstart=24)
  Group100_dup <- data.frame(Sample = names(km$cluster),Group = km$cluster)
  write.xlsx(Group100_dup,file=paste0(filename,"K均值合并100次分组结果.xlsx"))

} 


load("1.rawdata/Protein.50.log2.Combat.T.Top25.nmf.RData")

Group50Twoway(NMFResult = Protein.50.log2.Combat.T.Top25.nmf,
              filename = "Protein.50.log2.Combat.T.Top25.nmf_"
)




#### 分组比较
g1<- read.xlsx("1.rawdata/分组信息-Protein.50.log2.Combat.T.Top25.6.3.xlsx",rowNames = T)
g2<- read.xlsx("Protein.50.log2.Combat.T.Top25.nmf_第100次分组结果.xlsx")
g3<-  read.xlsx("Protein.50.log2.Combat.T.Top25.nmf_K均值合并100次分组结果.xlsx")


Group <- data.frame(Sample = g1$sample,
                    "subtype" = g1[,2],
                    "the hundredth time subtype" = g2[,2],
                    "combined 100 time subtype" = g3[,2])
row.names(Group) <- Group$Sample
Group <- Group[order(Group$subtype),]


Group$subtype<-paste0("Group",Group$subtype)
Group[,3]<-paste0("Group",Group[,3])
Group[,4]<-paste0("Group",Group[,4])

# 画图
library(ComplexHeatmap)
library(RColorBrewer)
unique(unlist(Group[,2:4]))

col_sub = c("#66C2A5","#FC8D62","#8DA0CB")
names(col_sub) <- c("Group1","Group2","Group3")

filename="NMF.subtype"
p <- Heatmap(t(Group[,2:4]), name = "Cluster",col=col_sub,
             width = unit(18, "cm"),
             height = unit(5, "cm"),
             cluster_columns = FALSE, cluster_rows = FALSE,
             #rect_gp = gpar(col = "white", lwd = 1),
             border = F,show_column_names = T,
             column_title_gp = gpar(fontsize = 15, family="serif"),
             column_names_gp = gpar(fontsize=4, family="serif"),
             row_names_gp = gpar(fontsize=12, family="serif"),
             row_names_side = "right",show_heatmap_legend = T,
             heatmap_legend_param = list(labels_gp = gpar(family="serif",fontsize = 12),
                                         title_gp = gpar(family="serif",fontsize = 12)))


png(filename = paste0(filename,".png"),width = 12,height = 5,units ="in",bg="white",res=300)
draw(p,heatmap_legend_side = "right")
dev.off()
pdf(file = paste0(filename,".pdf"),width = 12,height = 5)
draw(p,heatmap_legend_side = "right")
dev.off()
