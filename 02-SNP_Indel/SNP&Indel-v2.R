
# 20220824
# 调整做图：
# 1. SNP参考maftools的绘图风格排序
# 2. 加入突变频率的注释
# 3. 字体放大
# 4. 行按照突变频率排序
# 对比不同蛋白组学亚型（top25的3型）患者特定突变（突变频率≥3%的基因）的差异（SNV和Indel加在一起画图），计算P值。

rm(list = ls())
options(stringsAsFactors = FALSE)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(dplyr)
library(tidyr)
library(reshape2)



maf = read.table(file = "1.rawdata/all.T.138.maf",
                 sep = "\t",quote = "",fill = TRUE,header = TRUE)
colnames(maf)

D <- maf[,c("Hugo_Symbol","Variant_Classification","Variant_Type","Sample")]
D2 <- filter(D, Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation",
                                              "Frame_Shift_Del","Splice_Site",
                                              "In_Frame_Del","Frame_Shift_Ins",
                                              "Translation_Start_Site","In_Frame_Ins",
                                              "Nonstop_Mutation"))
# 这里是9种，图上少了三种，原因是图上对应的基因没有注释到这样的突变
# 选这9种是因为文献都是这九种

table(D2$Variant_Type)
length(unique(D2$Sample))
#mut <- mut %>% distinct()  两种结果：1，一个基因多次突变。2，一个基因不同突变注释类型

# 用于核对
# "A0629_T"
D3 <- data.frame(Hugo_Symbol= character(0), Variant_Classification= character(0), Sample = character(0))
for(i in unique(D2$Sample)){
  print(i)
  a <- D2[,c("Hugo_Symbol","Variant_Classification","Sample")][D2$Sample == i,]
  
  # 需要考虑是否是同一种类型
  # 这是我修改的脚本
  b <- row.names(table(a))[apply(table(a),1,function(x){sum(x != 0)}) > 1]
  if(length(b) == 0){
    a <- a %>% distinct()
    D3 <- rbind(D3,a)
  }else{
    a[a$Hugo_Symbol %in% b,"Variant_Classification"] <- "Multi_Hit"
    a <- a %>% distinct()
    D3 <- rbind(D3,a)
  }
}
rm(i,a,b)
D3.data <- dcast(data=D3,Hugo_Symbol ~ Sample,value.var = "Variant_Classification")

row.names(D3.data) <- D3.data$Hugo_Symbol
D3.data <- D3.data[,2:ncol(D3.data)]
row.names(D3.data)[1] <- "SEPT9"

colnames(D3.data) <- sapply(colnames(D3.data),function(x){
  if(nchar(x) == 4){
    paste0("A0",substring(x,first = 2),"_T")
  }else{
    paste0(x,"_T")
  }
})

dim(D3.data)
library(openxlsx)
# write.xlsx(as.data.frame(D3.data),file = "4.SNP&Indelv2-用于核对的数据.xlsx",
#            asTable = TRUE,overwrite = TRUE,
#            colNames = TRUE,rowNames = TRUE)

table(D3$Variant_Classification)
table(as.matrix(D3.data[,2:ncol(D3)]))

D3.data["GOLGA8K",] %>% t() %>% table()

load("1.rawdata/Protein.50.log2.Combat.T.Top25.6.RData")
ParseNMF <- function(NMFResult,Data,filename){
  
  library(NMF)
  Group <- predict(NMFResult)
  Group <- data.frame(Sample = names(Group),Group = Group)
  row.names(Group) <- Group$Sample
  
  NumGroup <- length(unique(Group$Group))
  for(i in 1:NumGroup){
    tem <- intersect(Group$Sample[Group$Group == i],colnames(Data))
    assign(paste0("Group.",i),tem)
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
  
  Data <- Data[,NewColumn]
  Data2 <- Data
  Data2[!is.na(Data2)] <- 1
  Data2[is.na(Data2)] <- 0
  pvalue <- apply(Data2,1,function(x){
    # 筛选后，有些行可能没有突变，会导致报错
    if(length(unique(x)) == 1){
      1
    }else{
      # 这里是离散值，需要用fisher检验
      df.table <- table(data.frame(Protein.Group = GroupLabel,
                                   Data = x))
      fisher.test(df.table,alternative = "two.sided",
                  simulate.p.value=F,workspace = 2e7)$p.value
    }
  })
  
  Data.mut <- apply(Data,1,function(x){sum(!is.na(x))/length(x)})
  
  Data.out <- cbind(Data,pvalue,data.frame("Mutation percent" = Data.mut))
  library(openxlsx)
  write.xlsx(x = as.data.frame(Data.out),file = paste0(filename,".xlsx"),
             asTable = TRUE,overwrite = TRUE,
             rowNames = TRUE,colNames = TRUE)
  
  Data <- Data[order(Data.mut,decreasing = TRUE),]
  pvalue <- pvalue[order(Data.mut,decreasing = TRUE)]
  Data.mut <- Data.mut[order(Data.mut,decreasing = TRUE)]
  
  Data <- Data[Data.mut > 0.03,]
  pvalue <- pvalue[Data.mut > 0.03]
  Data.mut <- Data.mut[Data.mut > 0.03]

  # 画图
  library(ComplexHeatmap)
  # Heatmap.color <- c("#473C8B", "#EE4000", "#FFA500","#00A08A","#F98400", "#85D4E3","#F4B5BD", "#9C964A")
  # 更改top annotation颜色
  Heatmap.color <- c("#f68c80", "#92cbc2", "#7aa8d0")
  Heatmap.color.list <- Heatmap.color[1:length(GroupLabel.uniq)]
  names(Heatmap.color.list) <- GroupLabel.uniq
  
  pvalue.label <- sapply(pvalue,function(x){
    if(as.numeric(x) < 0.001){
      format(x,scientific = TRUE,digits = 2)
    }else{
      format(round(x,digits = 3),nsmall = 3)
    }
  })
  Freq.label <- sapply(Data.mut,function(x){
    if(as.numeric(x) < 0.005){
      "0%"
    }else{
      paste0(format(round(x*100,digits = 0),nsmall = 0),"%")
    }
  })
  
  hr<-rowAnnotation(text = anno_text(pvalue.label,
                                     gp = gpar(fontsize = 24, font = 2,
                                               col = ifelse(pvalue < 0.05,"red","black")), # 大小
                                     rot = 0,location = 0,just = "left",width = unit(25,"mm")),
                    Frequency = anno_text(Freq.label,
                                          gp = gpar(fontsize = 24, font = 2,
                                                    col = ifelse(pvalue < 0.05,"red","black")), # 大小
                                          rot = 0,location = 0,just = "left",width = unit(25,"mm")))
  hc = HeatmapAnnotation("Protein Subgroup" = GroupLabel,
                         simple_anno_size = unit(1, "cm"),
                         annotation_name_side = c("right"),
                         annotation_name_gp = gpar(family="serif",fontsize = 24,font = 2),
                         col = list("Protein Subgroup" = Heatmap.color.list),
                         annotation_legend_param = list(labels_gp = gpar(family="serif",fontsize = 24),
                                                        title_gp = gpar(family="serif",fontsize = 24),
                                                        grid_width = unit(1, "cm"),
                                                        grid_height = unit(1, "cm")))
  
  Mycol <- c("#2371A9", "#F5BA6E", "#ACD387","#F09694","#A3C9DC", 
             "#379838","#EE7C1A", "#9C964A","#AEDCDF")
  names(Mycol) <- c("Missense_Mutation","Splice_Site","Nonsense_Mutation","Frame_Shift_Ins","Frame_Shift_Del",
                    "Multi_Hit","In_Frame_Del","Nonstop_Mutation","In_Frame_Ins")
  
  p <- Heatmap(Data, name = "?",na_col = "#CBCCCC",col = Mycol,
               cluster_columns = FALSE, cluster_rows = FALSE,
               rect_gp = gpar(col = "white", lwd = 2),
               # column_split = GroupLabel,
               border = F,show_column_names = T,column_title_gp = gpar(fontsize = 24, family="serif"),
               column_names_gp = gpar(fontsize=24, family="serif"),
               row_names_gp = gpar(fontsize=24, family="serif"),
               row_names_side = "left",show_heatmap_legend = T,
               right_annotation = hr,
               top_annotation = hc,
               heatmap_legend_param = list(labels_gp = gpar(family="serif",fontsize = 24),
                                           title_gp = gpar(family="serif",fontsize = 24),
                                           title = " ",
                                           grid_width = unit(1, "cm"),
                                           grid_height = unit(1, "cm")))
  
  png(filename = paste0(filename,".png"),
      width = (0.6/2.54*ncol(Data) + 22),
      height = (1.5/2.54*nrow(Data) + 1),
      units ="in",bg="white",res=300)
  print(p)
  dev.off()
  pdf(file = paste0(filename,".pdf"),
      width = (0.6/2.54*ncol(Data) + 22),
      height = (1.5/2.54*nrow(Data) + 1))
  print(p)
  dev.off()
  
  save(list = ls(),file = paste0(filename,".RData"))
  
}

ParseNMF(NMFResult = Protein.50.log2.Combat.T.Top25.6$fit$`3`,
         Data = D3.data,
         filename = "Protein.50.log2.Combat.T.Top25.6.Group3-SNV_Indel")

