#Seurat4Cluster 请按需修改参数
options(stringAsFactors=FALSE)
#rds结果文件路径
rdsFile = "/inputfile.rds"
#输出结果前缀
prefix = "samplename"
#输出结果路径
outpath = "Clusteroutpath"
#选择一个计算marker方法: MAST,wilcox,None(JustCluster)
MarkerGeneMethod = "wilcox"
#聚类方法选择kmean时的kmean值
Kmean = 10
#Logfc.Threshold卡值选择
threshold = 0.25
#最小pct卡值
minpct = 0.1
#聚类系数
resolution = 0.8
#输入指定聚类结果（Manual）,没有则为NULL
clusterFile = NULL
#算logfc的底数 2或e
base = "2"
OnlyPos = TRUE
if(base=="e"){
	base=exp(1)
}else if(base=="2"){
	base=2
}

library(SingleCellExperiment)
library(scater)
library(plyr)
library(reshape2)
library(Seurat)
library(mclust)
library(dplyr)
print("Start")
setwd(outpath)
seuset = readRDS(rdsFile)
assay = DefaultAssay(seuset)
dir.create("GraphClust")
if(!is.null(seuset@commands$RunUMAP)){
    dims = seuset@commands$RunUMAP@params$dims
}else if(!is.null(seuset@commands$RunTSNE)){
    dims = seuset@commands$RunTSNE@params$dims
}
seuset <- FindNeighbors(seuset, reduction = "pca",dims = dims)
seuset <- FindClusters(seuset, resolution = resolution)
table = table(Idents(seuset),seuset@meta.data$orig.ident)
print(table)
table = dcast(data.frame(table),Var1~Var2)
colnames(table)[1]="Cluster"
write.table(table,file=paste0(prefix,"_GraphClust.Statistics.txt"),sep="\t",row.names=F,quote=F)
data = data.frame(Cell = colnames(seuset),Cluster=Idents(seuset))
write.table(data,file=paste0(prefix,"_GraphClust.Summary_Cell.txt"),sep="\t",row.names=F,quot=F)
allmarkers <- FindAllMarkers(object = seuset, only.pos = onlypos, min.pct = minpct,logfc.threshold = threshold, test.use = MarkerGeneMethod,base=base)
write.table(data.frame(allmarkers$gene,allmarkers[,1:6]),file=paste0(prefix,"_GraphClust.AllMarkerGenes.txt"),sep="\t",row.names=F,quote=F)
saveRDS(seuset, file = paste0("../",prefix,"_GraphClust.seuset.rds"))

