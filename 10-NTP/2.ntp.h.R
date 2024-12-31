rm(list = ls())
options(stringsAsFactors = FALSE)
path=getwd()

library(openxlsx)
library(dplyr)
library(cmapR)

# 设定的参数
distance.selection<-"cosine"
weight.genes<-"F" # 写成FALSE会报错
num.resamplings<-"1000"
GenePattern.output<-"T"
random.seed<-6

# 
source("/ntp_function.R")

NTPez(
  input.exp.filename <- "./3.NTP分型top39/1.输入数据/Protein.50.log2.Gene.gct",
  input.features.filename <- "./3.NTP分型top39/1.输入数据/protein.DEP.top39.intersect.verifiy.txt",
  output.name  <- "Protein.50.log2.NPT.top39.human",
  dist.selection  <- distance.selection,
  temp.nn.wt  <- weight.genes,
  nresmpl  <- num.resamplings,
  GenePattern.output  <- GenePattern.output,
  rnd.seed  <- random.seed
)

save(list = ls(),file = "8.NTP分型top39.RData")



rm(list = ls())
options(stringsAsFactors = FALSE)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(openxlsx)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
path="./human/3.NTP分型top39/"
# NTP结果
NTP <- read.table(file = "./human/3.NTP分型top39/2.原始结果/Protein.50.log2.NPT.top39.human_prediction_result.xls",
                  header = TRUE,row.names = 1)
load("./human/2.toGENE/Protein.50.log2.Gene.RData")

GeneSignature <- read.table(paste0(path,"1.输入数据/protein.DEP.top39.intersect.verifiy.txt"),header = TRUE)

re<- data.frame(table(NTP$predict.label))
re$Var1<- paste0("Group",re$Var1)
write.xlsx(re,"Protein.50.log2.NPT.top39.human.table.xlsx")


###NTP热图
ParseNMF <- function(NMFResult,NTP,SampleInfo,filename,Protein,Feature){
  
  #### 获得分组信息并输出
  library(NMF)
  Group <- predict(NMFResult)
  Group <- data.frame(Sample = names(Group),
                      NMF.Group = paste0("Group",Group),
                      NTP.Group = paste0("PCS",NTP[names(Group),"predict.label"]))
  row.names(Group) <- Group$Sample
  library(openxlsx)
  write.xlsx(Group,file = paste0(filename,".xlsx"),
             asTable = TRUE,overwrite = TRUE,
             colNames = TRUE,rowNames = FALSE)
  
  # 输出NMF和NTP分型的对应表
  library(dplyr)
  Group.table <- table(Group[,c("NMF.Group","NTP.Group")])
  Group.table.rownames <- row.names(Group.table)
  Group.table.colnames <- colnames(Group.table)
  Group.table <- matrix(data = Group.table,nrow = nrow(Group.table),byrow = FALSE)
  row.names(Group.table) <- Group.table.rownames
  colnames(Group.table) <- Group.table.colnames
  write.xlsx(as.data.frame(Group.table),file = paste0(filename,".table.xlsx"),
             asTable = TRUE,overwrite = TRUE,
             colNames = TRUE,rowNames = TRUE)
  
  # 生存分析

    col_map<-c("#e75545", "#55bed6", "#25a38a","#f6b130","#70ae46","#4472c5","#9FB6CD")[1:GroupNum]

  
  # 输出热图
  Protein <- Protein[,Group$Sample]
  Protein.mean <- apply(Protein,1,function(x){
    tapply(x,INDEX = Group$NTP.Group,function(y){
      mean(y,na.rm = TRUE)
    })
  }) %>% t()
  
  Gene.label <- Feature[row.names(Protein.mean),"Group"]
  Protein.mean.scale <- apply(Protein.mean,1,scale) %>% t()
  colnames(Protein.mean.scale) <- colnames(Protein.mean)
  
  library(ComplexHeatmap)
  Heatmap.color <- c("#473C8B", "#EE4000", "#FFA500","#00A08A","#F98400", "#85D4E3","#F4B5BD", "#9C964A")
  Heatmap.color.list <- Heatmap.color[1:length(unique(Gene.label))]
  names(Heatmap.color.list) <- unique(Gene.label)
  
  ## 设置注释
  Heatmap.anno.Gene = rowAnnotation("NTP" = Gene.label,
                                    simple_anno_size = unit(0.5, "cm"),
                                    annotation_name_side = c("top"),
                                    annotation_name_gp = gpar(fontsize = 14,font = 1),
                                    col = list("NTP" = Heatmap.color.list),
                                    annotation_legend_param = list(grid_height = unit(1, "cm"),
                                                                   grid_width = unit(8, "mm"),
                                                                   labels_gp = gpar(fontsize = 14),
                                                                   title_gp = gpar(fontsize = 14)))
  library(circlize)
  col_fun = colorRamp2(c(-1,0,1), c("#233574","white", "#C01822")) # 刘博
  
  p <- Heatmap(Protein.mean.scale,name = "z-score",col = col_fun,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_row_dend = FALSE,show_row_names = TRUE,
               show_column_names = TRUE,show_column_dend = FALSE,column_names_side = "top",
               right_annotation = Heatmap.anno.Gene,
               heatmap_legend_param = list(legend_height = unit(3, "cm"),
                                           grid_width = unit(0.6, "cm"),
                                           labels_gp = gpar(fontsize = 14),
                                           title_gp = gpar(fontsize = 14)))
  
  png(filename = paste0(filename,".png"),
      width = (0.5/2.54*ncol(Protein.mean) + 2.5),height = (0.5/2.54*nrow(Protein.mean) + 1),
      units ="in",bg="white",res=300)
  print(p)
  dev.off()
  pdf(file = paste0(filename,".pdf"),
      width = (0.5/2.54*ncol(Protein.mean) + 2.5),height = (0.5/2.54*nrow(Protein.mean) + 1))
  print(p)
  dev.off()
}



###signature热图
NTP_downstrain<-function(Data,NTP,GeneSignature,filename){
  NTP <- NTP[order(NTP$predict.label),]
  NTP$predict.label<- factor(NTP$predict.label)
  
  Data <- Data[intersect(GeneSignature,rownames(Data)),rownames(NTP)]
  
  Data.zscore <- t(apply(log2(Data),1,scale))
  colnames(Data.zscore) <- colnames(Data)
  
  Heatmap.color <- c("#473C8B", "#EE4000", "#FFA500","#00A08A","#F98400", "#85D4E3","#F4B5BD", "#9C964A")
  Heatmap.color.list <- Heatmap.color[1:length(levels(NTP$predict.label))]
  names(Heatmap.color.list) <- levels(NTP$predict.label)
  
  ## 设置注释
  Heatmap.anno.Group = HeatmapAnnotation("Protein Subgroup" = NTP$predict.label,
                                         simple_anno_size = unit(1.5, "cm"),
                                         annotation_name_side = c("right"),
                                         annotation_name_gp = gpar(fontsize = 20,font = 2),
                                         col = list("Protein Subgroup" = Heatmap.color.list),
                                         annotation_legend_param = list(grid_height = unit(1, "cm"),
                                                                        grid_width = unit(8, "mm"),
                                                                        labels_gp = gpar(fontsize = 20),
                                                                        title_gp = gpar(fontsize = 20))
  )
  col_fun = colorRamp2(c(-1,0,1), c("#233574","white", "#C01822")) # 刘博
  
  p <- Heatmap(Data.zscore,name = "z-score",col = col_fun,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_row_dend = FALSE,show_row_names = FALSE,
               show_column_names = FALSE,show_column_dend = FALSE,
               top_annotation = Heatmap.anno.Group,
               column_split = NTP$predict.label,
               heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                           grid_width = unit(1, "cm"),
                                           labels_gp = gpar(fontsize = 20),
                                           title_gp = gpar(fontsize = 20)))
  
  png(filename = paste0(filename,".png"),width = 15,height = 24,units ="in",bg="white",res=300)
  print(p)
  dev.off()
  pdf(file = paste0(filename,".pdf"),width = 15,height = 24)
  print(p)
  dev.off()
  
}

NTP_downstrain(Protein.50.log2.T.Gene,NTP,GeneSignature$GeneName,"Protein.50.log2.T.verifiy.NPT.top39")






###2.gene panel热图
library(pheatmap)
#基因
gene_group<- rep(c("SEG1","SEG2","SEG3"),c(9,9,15))


load("./human/2.toGENE/Protein2Gene.RData")

Protein.50.log2.SEG<-Protein.50.log2.Gene[GeneSignature$GeneName,]
out<- cbind(names(Protein2Gene)[match(GeneSignature$GeneName,Protein2Gene)],data.frame(Protein.50.log2.SEG))
colnames(out)<- c("Protein ID",colnames(Protein.50.log2.SEG))
write.xlsx(out,paste0(path,"3.分析结果/Protein.50.log2.top39.human.xlsx"),
           overwrite = TRUE,colNames = TRUE,rowNames = T)

##分组
group<- data.frame(sample = rownames(NTP),group = paste0("Group",NTP$predict.label))



Data_mean<-c()
for(i in levels(factor(group$group))) {
  d <- Protein.50.log2.SEG[,group$group==i]
  row_mean<-apply(d,1,mean)
  Data_mean<-cbind(Data_mean,row_mean)
  
}
colnames(Data_mean)<-c("Subgroup1","Subgroup2","Subgroup3")

Geneorder<-read.table(paste0(path,"top39_33.txt"),header = TRUE)
plotData<- Data_mean[Geneorder$GeneName,]

Heatmap.color <- c("#473C8B", "#EE4000", "#FFA500","#00A08A","#F98400", "#85D4E3","#F4B5BD", "#9C964A")
color = colorRampPalette(c("navy", "white", "firebrick3"))


scale_dat<- t(apply(log2(plotData),1,scale))
colnames(scale_dat)<-colnames(plotData)

#color
ann_colors<-list(Group=c(Heatmap.color[1:3]))
names(ann_colors[[1]]) <- levels(factor(gene_group))

ann_row<-data.frame(Group=gene_group)
rownames(ann_row) <- GeneSignature$GeneName


#breaks_list <- seq(-2,2,length.out= 244)
p<-pheatmap(scale_dat,cluster_rows = F, cluster_cols = F,
            show_colnames = T,show_rownames = T,
            #annotation_col = ann_col,
            annotation_row = ann_row,
            annotation_colors = ann_colors,
            annotation_names_row = F,
            #annotation_legend = F,
            angle_col = "45",
            color = colorRampPalette(c("navy", "white", "firebrick3"))(n=244),
            #breaks = breaks_list,
            border = F)
pdf(paste0(path,"3.分析结果/Protein.RF.top39.human.heatmap.pdf"),width = 5,height = 6)
print(p)
dev.off()
png(paste0(path,"3.分析结果/Protein.RF.top39.human.heatmap.png"),width = 5,height = 6,units ="in",bg="white",res=300)
print(p)
dev.off()






###分组结果GSEA
ParseNMF <- function(Group,filename,data,database){
  
  Group$group<- factor(Group$group)
  #### 从获得的分组信息进行下游的整理
  NumGroup <- length(unique(Group$group))
  for(i in 1:NumGroup){
    assign(paste0("Group.",i),Group$Sample[Group$group == i])
  }
  NewColumn <- c() # 排序后的样本名（按照亚型）
  GroupNum <- c() # 各亚型的样本数目，用于t检验
  GroupLabel <- c() # 排序后每个样本对应的分组标签，用于anova和热图
  GroupLabel.uniq <- c() # 排序后的分组标签，用于绘制热图时分组颜色的设定
  
  for(i in 1:NumGroup){
    NewColumn <- c(NewColumn,get(paste0("Group.",i)))
    GroupNum <- c(GroupNum,length(get(paste0("Group.",i))))
    GroupLabel <- c(GroupLabel,rep(paste0("Subroup",i),length(get(paste0("Group.",i)))))
    GroupLabel.uniq <- c(GroupLabel.uniq,paste0("Subroup",i))
  }
  
  # 对每一组进行1vsRest的GSEA富集分析
  library(clusterProfiler)
  pathwayDatabase <- read.gmt(database)
  for(i in 1:NumGroup){
    
    SampleInGroup <- Group[Group$group==i,1]
    SampleNotInGroup <- colnames(data)[!(colnames(data) %in% SampleInGroup)]
    tem <- apply(data,1,function(x){
      log2(mean(x[SampleInGroup])/mean(x[SampleNotInGroup]))
    })
    tem <- sort(tem,decreasing = TRUE)
    save(tem,file = paste0(filename,"GSEA.Group",i,".fc.RData"))
    
    gsegmt <- GSEA(tem, TERM2GENE=pathwayDatabase, verbose=F,pvalueCutoff = 1,seed = TRUE,minGSSize = 1)
    gsegmt.sum <- summary(gsegmt)
    
    # 输出结果
    save(gsegmt,file = paste0(filename,"GSEA.Group",i,".RData"))
    
    library(openxlsx)
    write.xlsx(as.data.frame(gsegmt.sum),file = paste0(filename,"GSEA.Group",i,".xlsx"),
               asTable = TRUE,overwrite = TRUE)
    
    library(enrichplot)
    library(ggplot2)
    for(j in 1:nrow(gsegmt)){
      p <- gseaplot2(gsegmt, gsegmt.sum$ID[j],title = gsegmt.sum$ID[j])
      filename2 <- paste0(filename,".GSEA.Group",i,".",gsegmt.sum$ID[j])
      if(nchar(filename2) > 100){
        filename2 <- substr(filename2,start = 1,stop = 100)
      }
      ggsave(p,filename = paste0(filename2,".png"),width = 5, height = 5, units = "in")
      ggsave(p,filename = paste0(filename2,".pdf"),width = 5, height = 5, units = "in")
    }
    
    # 清理环境，避免问题
    rm(SampleInGroup,SampleNotInGroup,tem,gsegmt,gsegmt.sum,p)
  }
  
}

ParseNMFGSEAResult <- function(path,filename,filename.abbr,path2){
  
  library(clusterProfiler)
  n = 3
  for(i in 1:n){
    load(paste0(path,"/",filename,"Group",i,".RData"))
    assign(paste0("Table",i),summary(gsegmt))
  }
  
  CombineTable3 <- function(Table1,Table2,Table3,filename){
    library(openxlsx)
    
    ID <- unique(c(Table1$ID,Table2$ID,Table3$ID))
    
    Result.NES <- data.frame(NES.Group1 = rep(NA,length(ID)),
                             NES.Group2 = rep(NA,length(ID)),
                             NES.Group3 = rep(NA,length(ID)))
    row.names(Result.NES) <- ID
    Result.NES[Table1$ID,"NES.Group1"] <- Table1[,"NES"]
    Result.NES[Table2$ID,"NES.Group2"] <- Table2[,"NES"]
    Result.NES[Table3$ID,"NES.Group3"] <- Table3[,"NES"]
    write.xlsx(Result.NES,file = paste0(filename,".NES.xlsx"),
               asTable = TRUE,overwrite = TRUE,rowNames = TRUE,sheetName = "NES")
    
    Result.FDR <- data.frame(FDR.Group1 = rep(NA,length(ID)),
                             FDR.Group2 = rep(NA,length(ID)),
                             FDR.Group3 = rep(NA,length(ID)))
    row.names(Result.FDR) <- ID
    Result.FDR[Table1$ID,"FDR.Group1"] <- Table1[,"p.adjust"]
    Result.FDR[Table2$ID,"FDR.Group2"] <- Table2[,"p.adjust"]
    Result.FDR[Table3$ID,"FDR.Group3"] <- Table3[,"p.adjust"]
    write.xlsx(Result.FDR,file = paste0(filename,".FDR.xlsx"),
               asTable = TRUE,overwrite = TRUE,rowNames = TRUE,sheetName = "FDR")
    
    re <- list(NES = Result.NES,
               FDR = Result.FDR)
    return(re)
  }
  CombineTable4 <- function(Table1,Table2,Table3,Table4,filename){
    library(openxlsx)
    
    ID <- unique(c(Table1$ID,Table2$ID,Table3$ID,Table4$ID))
    
    Result.NES <- data.frame(NES.Group1 = rep(NA,length(ID)),
                             NES.Group2 = rep(NA,length(ID)),
                             NES.Group3 = rep(NA,length(ID)),
                             NES.Group4 = rep(NA,length(ID)))
    row.names(Result.NES) <- ID
    Result.NES[Table1$ID,"NES.Group1"] <- Table1[,"NES"]
    Result.NES[Table2$ID,"NES.Group2"] <- Table2[,"NES"]
    Result.NES[Table3$ID,"NES.Group3"] <- Table3[,"NES"]
    Result.NES[Table4$ID,"NES.Group4"] <- Table4[,"NES"]
    write.xlsx(Result.NES,file = paste0(filename,".NES.xlsx"),
               asTable = TRUE,overwrite = TRUE,rowNames = TRUE,sheetName = "NES")
    
    Result.FDR <- data.frame(FDR.Group1 = rep(NA,length(ID)),
                             FDR.Group2 = rep(NA,length(ID)),
                             FDR.Group3 = rep(NA,length(ID)),
                             FDR.Group4 = rep(NA,length(ID)))
    row.names(Result.FDR) <- ID
    Result.FDR[Table1$ID,"FDR.Group1"] <- Table1[,"p.adjust"]
    Result.FDR[Table2$ID,"FDR.Group2"] <- Table2[,"p.adjust"]
    Result.FDR[Table3$ID,"FDR.Group3"] <- Table3[,"p.adjust"]
    Result.FDR[Table4$ID,"FDR.Group4"] <- Table4[,"p.adjust"]
    write.xlsx(Result.FDR,file = paste0(filename,".FDR.xlsx"),
               asTable = TRUE,overwrite = TRUE,rowNames = TRUE,sheetName = "FDR")
    
    re <- list(NES = Result.NES,
               FDR = Result.FDR)
    return(re)
  }
  
  if(file.exists(path2)){
    # 在函数中改变路径，会改变外部工作路径
    current_wd <- getwd()
    setwd(path2)
  }else{
    dir.create(path2)
    current_wd <- getwd()
    setwd(path2)
  }
  if(n == 3){
    Re = CombineTable3(Table1,Table2,Table3,filename = paste0(filename.abbr,".Group","n"))
  }
  if(n == 4){
    Re = CombineTable4(Table1,Table2,Table3,Table4,filename = paste0(filename.abbr,".Group","n"))
  }
  save(Re,file = paste0(substring(filename,first = 1,last = nchar(filename)-5),".RData"))
  
  PlotResult <- function(NES,FDR,name,filename){
    
    # 对数据进行过滤
    # FDR.Flag <- apply(FDR,1,function(x){
    #   sum(x > 0.05,na.rm = TRUE) < sum(!is.na(x))
    # })
    # NES <- NES[FDR.Flag,]
    # FDR <- FDR[FDR.Flag,]
    
    # 对数据进行排序
    library(dplyr)
    colnames(NES) <- paste0("Group",1:ncol(NES))
    colnames(FDR) <- paste0("Group",1:ncol(FDR))
    if(ncol(NES) == 2){
      NES = arrange(NES,Group1,Group2)
      FDR = FDR[row.names(NES),]
    }
    if(ncol(NES) == 4){
      NES = arrange(NES,Group1,Group2,Group3,Group4)
      FDR = FDR[row.names(NES),]
    }
    
    # name参数有两个取值，NES或者FDR
    library(ComplexHeatmap)
    
    ## 亚型的颜色设置
    Heatmap.color <- c("#473C8B", "#EE4000", "#FFA500","#00A08A","#F98400", "#85D4E3","#F4B5BD", "#9C964A")
    Heatmap.color.list <- Heatmap.color[1:ncol(NES)]
    names(Heatmap.color.list) <- paste0("Group",1:ncol(NES))
    
    Heatmap.anno.Group = HeatmapAnnotation("Protein Subgroup" = paste0("Group",1:ncol(NES)),
                                           
                                           simple_anno_size = unit(0.8, "cm"),
                                           annotation_name_side = c("right"),
                                           annotation_name_gp = gpar(fontsize = 20,font = 2),
                                           
                                           col = list("Protein Subgroup" = Heatmap.color.list),
                                           annotation_legend_param = list(grid_height = unit(0.8, "cm"),
                                                                          grid_width = unit(8, "mm"),
                                                                          labels_gp = gpar(fontsize = 20),
                                                                          title_gp = gpar(fontsize = 20)))
    # 缺失值的图例
    lgd_list = list(
      Legend(labels = "NA", legend_gp = gpar(fill = "#5D5D5D"),
             grid_height = unit(0.8, "cm"), grid_width = unit(8, "mm"),
             labels_gp = gpar(fontsize = 20))
    )
    
    FDR2 <- FDR
    FDR2[FDR >= 0.05] = NA
    FDR2[FDR < 0.05] <- "*"
    
    library(circlize)
    MyCol = colorRamp2(c(-2,0,2), c("#233574","white", "#C01822")) # 刘博
    nameLabel = "NES"
    rect_gp = "white"
    
    p <- Heatmap(NES, name = nameLabel, 
                 col = MyCol,na_col = "#5D5D5D",rect_gp = gpar(col = rect_gp, lwd = 1),
                 cluster_rows = FALSE,show_row_names = TRUE,row_names_side = "left",
                 cluster_columns = FALSE,show_column_names = FALSE,
                 width = unit(1.8*ncol(NES), "cm"),height = unit((1.5*nrow(NES)),"cm"),
                 top_annotation = Heatmap.anno.Group,
                 heatmap_legend_param = list(legend_height = unit(3, "cm"),
                                             grid_width = unit(0.5, "cm"),
                                             labels_gp = gpar(fontsize = 20),
                                             title_gp = gpar(fontsize = 20)),
                 cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                   if(!is.na(FDR2[i, j])){
                     grid.text(FDR2[i, j], x, y)
                   }
                 })
    
    png(filename = paste0(filename,".heatmap.png"),
        width = (0.8/2.54*ncol(NES) + 15),height = (0.5/2.54*nrow(NES) + 3),
        units ="in",bg="white",res=300)
    draw(p)
    dev.off()
    pdf(file = paste0(filename,".heatmap.pdf"),
        width = (0.8/2.54*ncol(NES) + 15),height = (0.5/2.54*nrow(NES) + 3))
    draw(p)
    dev.off()
    
  }
  PlotResult(NES = Re$NES,FDR = Re$FDR,name = "NES",filename = substring(filename,first = 1,last = nchar(filename)-5))
  
  setwd(current_wd)
}

#ref1
group<-data.frame(Sample=rownames(NTP),
                  group=NTP$predict.label)
ParseNMF(group,
         filename = "Protein.50.log2.NPT.top39.GSEA.human",
         data = Protein.50.log2.Gene,
         database = "./select.pathways.gmt")

ParseNMFGSEAResult(path = "./human/4.GSEA/单独的富集结果",
                   filename = "Protein.50.log2.NPT.top39.GSEA.humanGSEA.",
                   filename.abbr = "NPT.huaman.Hallmark",
                   path2 = "./human/4.GSEA/合并后的结果")



