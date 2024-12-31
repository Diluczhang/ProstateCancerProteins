library(openxlsx)
library(dplyr)
library(NMF)
library(ComplexHeatmap)
library(circlize)
# input
Cytoband <- read.xlsx("1.rawdata/AmpDelFus.xlsx", rowNames = T)
allSam <- colnames(Cytoband)
# Protein Subgroup
load("1.rawdata/Protein.50.log2.Combat.T.Top25.6.RData")
NMFResult <- Protein.50.log2.Combat.T.Top25.6$fit$`3`
group <- predict(NMFResult) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "sample") %>% 
    rename(group = 2)
ParseGroup <- function(group,data = NULL){
    library(NMF)
    UniqueGroup <- unique(group$group)
    Group.Sample <- list()
    for(i in UniqueGroup){
        if(is.null(data)){
            tem <- group$sample[group$group == i]
        }else{
            tem <- intersect(group$sample[group$group == i],colnames(data))
        }
        Group.Sample <- c(Group.Sample,list(tem))
    }
    names(Group.Sample) <- UniqueGroup
    
    NewColumn <- c()
    GroupNum <- c()
    GroupLabel <- c()
    GroupLabel.uniq <- c()
    
    for(i in UniqueGroup){
        NewColumn <- c(NewColumn,Group.Sample[[i]])
        GroupNum <- c(GroupNum,length(Group.Sample[[i]]))
        GroupLabel <- c(GroupLabel,rep(i,length(Group.Sample[[i]])))
        GroupLabel.uniq <- c(GroupLabel.uniq,i)
    }
    
    Result <- list(NumGroup = length(Group.Sample),
                   NewColumn = NewColumn,
                   GroupNum = GroupNum,
                   GroupLabel = GroupLabel,
                   GroupLabel.uniq = GroupLabel.uniq,
                   Group.Sample = Group.Sample)
    return(Result)
}
groupInfo <- ParseGroup(group,data = NULL)
# pvalue : fisher test
Cytoband <- bind_rows(as.data.frame(t(data.frame(groupInfo$NewColumn, row.names = groupInfo$NewColumn))), Cytoband)
Cytoband <- Cytoband[-1, ]
pvalue <- apply(Cytoband,1,function(x){
    df.table <- table(na.omit(data.frame(Protein.Group = groupInfo$GroupLabel, Cytoband = x)))
    fisher.test(df.table,alternative = "two.sided",simulate.p.value=F,workspace = 2e7)$p.value
})
pvalueAll <- sapply(as.numeric(pvalue), function(x) ifelse(x > 0.001, format(round(x, 3), scientific = FALSE), format(signif(x, 2), scientific = TRUE)))
names(pvalueAll) <- names(pvalue)
################### plot
df <- Cytoband[, c(groupInfo$Group.Sample$`1`, groupInfo$Group.Sample$`2`, groupInfo$Group.Sample$`3`)]
Heatmap.color <- c("#f68c80", "#92cbc2", "#7aa8d0")
Heatmap.color.list <- Heatmap.color[1:length(groupInfo$GroupLabel.uniq)]
names(Heatmap.color.list) <- c("Subgroup1","Subgroup2","Subgroup3")

hc = HeatmapAnnotation("Protein Subgroup" = c(rep("Subgroup1", length(groupInfo$Group.Sample$`1`)),
                                              rep("Subgroup2", length(groupInfo$Group.Sample$`2`)),
                                              rep("Subgroup3", length(groupInfo$Group.Sample$`3`))),
                       simple_anno_size = unit(0.5, "cm"),
                       annotation_name_side = c("right"),
                       annotation_name_gp = gpar(family="serif",fontsize = 16,font = 2),
                       col = list("Protein Subgroup" = Heatmap.color.list),
                       annotation_legend_param = list(labels_gp = gpar(family="serif",fontsize = 16),
                                                      title_gp = gpar(family="serif",fontsize = 16)))


# amp/del
Mycol2 <- c("#236BA9", "#EEEDED", "#E13C2D")
names(Mycol2) <- c("-2", "0", "2")
hr2 <-rowAnnotation(text = anno_text(pvalueAll[1:7],
                                     gp = gpar(fontsize = 16, font = 2, col = "black"),
                                     rot = 0,location = 0,just = "left",width = unit(25,"mm")))
hm2 <- Heatmap(df[1:7,,drop = FALSE], name = "amp",na_col = "#CBCCCC",col = Mycol2,
               top_annotation = hc,
               right_annotation = hr2,
               row_names_gp = gpar(fontsize=16, family="serif"),
               row_names_side = "left",
               show_heatmap_legend = F)

# fusion
Mycol4 <- c("#EEEDED","orange")
names(Mycol4) <- c("0", "2")
hr4<-rowAnnotation(text = anno_text(pvalueAll[8],
                                    gp = gpar(fontsize = 16, font = 2,
                                              col = "black"), # 大小
                                    rot = 0,location = 0,just = "left",width = unit(25,"mm")))

hm4 <- Heatmap(df[8,,drop = FALSE], name = "fusion",na_col = "#CBCCCC",col = Mycol4,
               right_annotation = hr4,
               row_names_gp = gpar(fontsize=16, family="serif"),row_names_side = "left",
               show_heatmap_legend = F)

ht_list <- hm2 %v% hm4

# 图例

lgd1 = Legend(at = c(0,-2),
              direction = "vertical",
              title = "Del",
              title_position = "topleft",
              labels = c("No del","Del","NA"),
              grid_height = unit(10, "mm"), # 区别于前面的legend_height
              grid_width = unit(10, "mm"),
              legend_gp = gpar(fill = c("#EEEDED","#236BA9","#CBCCCC")),
              labels_gp = gpar(family="serif",fontsize = 16),
              title_gp = gpar(family="serif",fontsize = 16))
lgd2 = Legend(at = c(0,2),
              direction = "vertical",
              title = "Amp",
              title_position = "topleft",
              labels = c("No amp","Amp","NA"),
              grid_height = unit(10, "mm"), # 区别于前面的legend_height
              grid_width = unit(10, "mm"),
              legend_gp = gpar(fill = c("#EEEDED","#E13C2D","#CBCCCC")),
              labels_gp = gpar(family="serif",fontsize = 16),
              title_gp = gpar(family="serif",fontsize = 16))

lgd3 = Legend(at = c(0,2),
              direction = "vertical",
              title = "Fusion",
              title_position = "topleft",
              labels = c("No fusion","Fusion","NA"),
              grid_height = unit(10, "mm"), # 区别于前面的legend_height
              grid_width = unit(10, "mm"),
              legend_gp = gpar(fill = c("#EEEDED","orange","#CBCCCC")),
              labels_gp = gpar(family="serif",fontsize = 16),
              title_gp = gpar(family="serif",fontsize = 16))


png(filename = "AmpDelFusGroupPlot.png",
    width = 30,
    height = 5,
    units ="in",bg="white",res=300)
draw(ht_list, ht_gap = unit(0, "mm"), 
     heatmap_legend_list = list(lgd1,lgd2,lgd3))
dev.off()

pdf(file = "AmpDelFusGroupPlot.pdf",
    width = 30,
    height = 5)
draw(ht_list, ht_gap = unit(0, "mm"), 
     heatmap_legend_list = list(lgd1,lgd2,lgd3))
dev.off()

