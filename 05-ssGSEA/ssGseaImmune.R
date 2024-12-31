library(cmapR)
library(openxlsx)
library(pacman)
library(NMF)
library(ComplexHeatmap)
library(circlize)
source("ssGSEA_function.R")
ParseNMFssGSEA <- function(NMFResult,filename,GCTfile,Immune.Group,RNAData){
  ssGSEA.GCT <- parse_gctx(GCTfile)
  ssGSEA.NES <- ssGSEA.GCT@mat
  ssGSEA.rdesc <- ssGSEA.GCT@rdesc
  ssGSEA.pvalue <- ssGSEA.rdesc[,grepl(pattern = "^pvalue",colnames(ssGSEA.rdesc),ignore.case = TRUE)]
  ssGSEA.FDR <- ssGSEA.rdesc[,grepl(pattern = "fdr",colnames(ssGSEA.rdesc),ignore.case = TRUE)]
  colnames(ssGSEA.pvalue) <- sapply(colnames(ssGSEA.pvalue),function(x){
    gsub(pattern = "pvalue.",replacement = "",x)
  })
  colnames(ssGSEA.FDR) <- sapply(colnames(ssGSEA.FDR),function(x){
    gsub(pattern = "fdr.pvalue.",replacement = "",x)
  })
  ssGSEA.FDR.label <- ssGSEA.FDR
  ssGSEA.FDR.label[ssGSEA.FDR.label >= 0.05] = NA
  ssGSEA.FDR.label[ssGSEA.FDR.label < 0.05] <- "*"
  
  # 列按照组别排列并对组间NES进行比较,两组用wilcox,两组以上用kruskal
  Group <- predict(NMFResult)
  Group <- data.frame(Sample = names(Group),Group = Group)
  row.names(Group) <- Group$Sample
  NumGroup <- length(unique(Group$Group))
  for(i in 1:NumGroup){
    assign(paste0("Group.",i),intersect(Group$Sample[Group$Group == i],colnames(RNAData)))
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
  ssGSEA.NES <- ssGSEA.NES[,NewColumn]
  ssGSEA.FDR.label <- ssGSEA.FDR.label[,NewColumn]
  # 计算pvalue
  if(NumGroup == 2){
    pvalue <- apply(ssGSEA.NES,1,function(x){
      wilcox.test(x[1:GroupNum[1]],x[(GroupNum[1] + 1):sum(GroupNum)],
                  alternative = "two.sided",paired = FALSE)$p.value
    })
  }else{
    pvalue <- apply(ssGSEA.NES,1,function(x){
      kruskal.test(x~GroupLabel)$p.value
    })
  }
  # 计算FDR
  FDR <- p.adjust(pvalue,method = "BH")
  
  # 行按组别排列
  Gene.Group1 <- intersect(Immune.Group[Immune.Group$Group == 0,"Pathway"],row.names(ssGSEA.NES))
  Gene.Group2 <- intersect(Immune.Group[Immune.Group$Group == 1,"Pathway"],row.names(ssGSEA.NES))
  Gene.Group3 <- intersect(Immune.Group[Immune.Group$Group == 2,"Pathway"],row.names(ssGSEA.NES))
  ssGSEA.NES <- ssGSEA.NES[c(Gene.Group1,Gene.Group2,Gene.Group3),]
  ssGSEA.FDR.label <- ssGSEA.FDR2[c(Gene.Group1,Gene.Group2,Gene.Group3),]
  pvalue <- pvalue[c(Gene.Group1,Gene.Group2,Gene.Group3)]
  FDR<- FDR[c(Gene.Group1,Gene.Group2,Gene.Group3)]
  
  # 输出结果
  write.xlsx(as.data.frame(cbind(ssGSEA.NES,pvalue,FDR)),
             file = paste0(filename,".total.xlsx"),
             asTable = TRUE,overwrite = TRUE,rowNames = TRUE)
  
  # FDR < 0.05进行筛选,并输出结果
  ssGSEA.NES.sig <- ssGSEA.NES[FDR < 0.05,]
  ssGSEA.FDR.sig <- ssGSEA.FDR.label[FDR < 0.05,]
  pvalue.sig <- pvalue[FDR < 0.05]
  FDR.sig <- FDR[FDR < 0.05]
  write.xlsx(as.data.frame(cbind(ssGSEA.NES.sig,pvalue.sig,FDR.sig)),
             file = paste0(filename,".sig.xlsx"),
             asTable = TRUE,overwrite = TRUE,rowNames = TRUE)

  ## 设置列注释
  Heatmap.color <- c("#473C8B", "#EE4000", "#FFA500","#00A08A","#F98400", "#85D4E3","#F4B5BD", "#9C964A")
  Heatmap.color.list <- Heatmap.color[1:NumGroup]
  names(Heatmap.color.list) <- GroupLabel.uniq
  hc = HeatmapAnnotation("Protein Subgroup" = GroupLabel,
                         simple_anno_size = unit(0.8, "cm"),
                         annotation_name_side = c("right"),
                         annotation_name_gp = gpar(fontsize = 20,font = 2),
                         col = list("Protein Subgroup" = Heatmap.color.list),
                         annotation_legend_param = list(grid_height = unit(0.8, "cm"),
                                                        grid_width = unit(8, "mm"),
                                                        labels_gp = gpar(fontsize = 20),
                                                        title_gp = gpar(fontsize = 20)))
  # 设置行注释
  Pathway.color <- c("#DFDFDF","#CB3127","#3281B2")
  names(Pathway.color) <- c("Pathway","ACTIVE","INHIBIT")
  FDR.label <- sapply(FDR,function(x){
    if(as.numeric(x) < 0.001){
      format(x,scientific = TRUE,digits = 2)
    }else{
      format(round(x,digits = 3),nsmall = 3)
    }
  })
  hr <- rowAnnotation(text = anno_text(FDR.label,
                                       gp = gpar(fontsize = 18, font = 2,col = ifelse(FDR < 0.05,"red","black")), # 大小
                                       rot = 0, # 角度
                                       location = 0, # 位置
                                       just = "left",
                                       width = unit(25,"mm")),
                      Group = c(rep("Pathway",length(Gene.Group1)),
                                rep("ACTIVE",length(Gene.Group2)),
                                rep("INHIBIT",length(Gene.Group3))),
                      simple_anno_size = unit(0.8, "cm"),
                      show_annotation_name = FALSE,
                      # annotation_name_side = c("bottome"),
                      # annotation_name_gp = gpar(fontsize = 20,font = 2),
                      col = list("Group" = Pathway.color),
                      annotation_legend_param = list(grid_height = unit(0.8, "cm"),
                                                     grid_width = unit(8, "mm"),
                                                     labels_gp = gpar(fontsize = 20),
                                                     title_gp = gpar(fontsize = 20)))
  FDR.label.sig <- sapply(FDR.sig,function(x){
    if(as.numeric(x) < 0.001){
      format(x,scientific = TRUE,digits = 2)
    }else{
      format(round(x,digits = 3),nsmall = 3)
    }
  })
  hr.sig <- rowAnnotation(text = anno_text(FDR.label.sig,
                                           gp = gpar(fontsize = 18, font = 2,col = ifelse(FDR.sig < 0.05,"red","black")), # 大小
                                           rot = 0, # 角度
                                           location = 0, # 位置
                                           just = "left",
                                           width = unit(25,"mm")),
                          Group = c(rep("Pathway",length(Gene.Group1[Gene.Group1 %in% Gene.Group.sig])),
                                    rep("ACTIVE",length(Gene.Group2[Gene.Group2 %in% Gene.Group.sig])),
                                    rep("INHIBIT",length(Gene.Group3[Gene.Group3 %in% Gene.Group.sig]))),
                          simple_anno_size = unit(0.8, "cm"),
                          show_annotation_name = FALSE,
                          # annotation_name_side = c("bottome"),
                          # annotation_name_gp = gpar(fontsize = 20,font = 2),
                          col = list("Group" = Pathway.color),
                          annotation_legend_param = list(grid_height = unit(0.8, "cm"),
                                                         grid_width = unit(8, "mm"),
                                                         labels_gp = gpar(fontsize = 20),
                                                         title_gp = gpar(fontsize = 20)))
  # NES scale
  ssGSEA.NES.scale <- t(apply(ssGSEA.NES, 1, scale))
  ssGSEA.NES.scale.sig <- t(apply(ssGSEA.NES.sig, 1, scale))
  # plot
  p <- Heatmap(ssGSEA.NES.scale,
               name = "NES z-score",
               col = colorRamp2(c(-2,0,2), c("#233574","white", "#C01822")),
               na_col = "#5D5D5D",rect_gp = gpar(col = "black", lwd = 1),
               cluster_rows = FALSE,show_row_dend = FALSE,
               show_row_names = TRUE,row_names_side = "left",
               cluster_columns = FALSE,show_column_names = FALSE,
               top_annotation = hc,right_annotation = hr,
               heatmap_legend_param = list(legend_height = unit(3, "cm"),
                                           grid_width = unit(0.5, "cm"),
                                           labels_gp = gpar(fontsize = 20),
                                           title_gp = gpar(fontsize = 20)))
  
  png(filename = paste0(filename,".heatmap.total.png"),
      width = (0.8/2.54*ncol(ssGSEA.NES) + 12),
      height = (0.5*nrow(ssGSEA.NES)/2.54 + 1),
      units ="in",bg="white",res=300)
  draw(p, padding = unit(c(1, 5, 1, 5), "cm"))
  dev.off()
  pdf(file = paste0(filename,".heatmap.total.pdf"),
      width = (0.8/2.54*ncol(ssGSEA.NES) + 12),
      height = (0.5*nrow(ssGSEA.NES)/2.54 + 1))
  draw(p, padding = unit(c(1, 5, 1, 5), "cm"))
  dev.off()
  p.sig <- Heatmap(ssGSEA.NES.scale.sig,
                   name = "NES z-score",
                   col = colorRamp2(c(-2,0,2), c("#233574","white", "#C01822")),
                   na_col = "#5D5D5D",rect_gp = gpar(col = "black", lwd = 1),
                   cluster_rows = FALSE,show_row_dend = FALSE,
                   show_row_names = TRUE,row_names_side = "left",
                   cluster_columns = FALSE,show_column_names = FALSE,
                   top_annotation = hc,right_annotation = hr.sig,
                   heatmap_legend_param = list(legend_height = unit(3, "cm"),
                                               grid_width = unit(0.5, "cm"),
                                               labels_gp = gpar(fontsize = 20),
                                               title_gp = gpar(fontsize = 20)))
  png(filename = paste0(filename,".heatmap.sig.png"),
      width = (0.8/2.54*ncol(ssGSEA.NES.sig) + 12),
      height = (0.5*nrow(ssGSEA.NES.sig)/2.54 + 1),
      units ="in",bg="white",res=300)
  draw(p.sig, padding = unit(c(1, 5, 1, 5), "cm"))
  dev.off()
  pdf(file = paste0(filename,".heatmap.sig.pdf"),
      width = (0.8/2.54*ncol(ssGSEA.NES.sig) + 12),
      height = (0.5*nrow(ssGSEA.NES.sig)/2.54 + 1))
  draw(p.sig, padding = unit(c(1, 5, 1, 5), "cm"))
  dev.off()
}
# Immune_new
load("1.rawdata/Immune_new.RData")
colnames(Immune)
Immune.Group <- data.frame(Pathway = colnames(Immune),Group = c(rep(0,19),rep(1,14),rep(2,12)))
load("1.rawdata/Protein.50.log2.Combat.T.Top25.6.RData")
##########################
dataTpm <- read.table("1.rawdata/PRAD_RNA_TPM.txt",header = T,row.names = 1,sep = "\t",quote = "")
gctTpm <- new("GCT",
              mat = data.matrix(dataTpm),
              rdesc = data.frame(Type="RNA",
                                 geneSymbol=row.names(dataTpm),
                                 stringsAsFactors = FALSE))
write_gct(gctTpm, ofile="RNA.T.gct", appenddim=FALSE)
ssGSEA_u(
  input.ds = "RNA.T.gct",
  output.prefix = "RNA.T.gct.ssGSEA.Immune_new",
  gene.set.databases = "1.rawdata/Immune_new.gmt",
  sample.norm.type = "rank",
  weight = 1,
  statistic = "area.under.RES",
  output.score.type = "NES",
  nperm = 1000,
  global.fdr = TRUE,
  min.overlap = 5,
  correl.type = "z.score",
  combine.mode = "combine.off"
)

ParseNMFssGSEA(NMFResult = Protein.50.log2.Combat.T.Top25.6$fit$`3`,
               filename = "RNA.T.ssGSEA.Immune-Protein.50.log2.Combat.T.Top25.6-Group3",
               GCTfile = "RNA.T.gct.ssGSEA.Immune_new-combined.gct",
               Immune.Group = Immune.Group,
               RNAData = dataTpm)