library(Seurat)
library(ggplot2)

#GeneBubblePlot
gg = DotPlot(seuset,features=rev(reallist),cols=featureplot) + RotatedAxis()
gg = gg + scale_x_discrete(breaks=rev(reallist),labels=rev(plottitle))
ggsave(paste0("Genebubbleplot.pdf"),gg,width = bubbleWidth,height = bubbleHeight)
ggsave(paste0("Genebubbleplot.png"),gg,width = bubbleWidth,height = bubbleHeight)

#分面GeneBubblePlot
DotPlot(seuset,features=split(markers$gene,markers$cluster),cols=featureplot)+RotatedAxis()+ 
  theme(panel.border = element_rect(color = "black"),
        panel.spacing = unit(1, "mm"),
        axis.title = element_blank(),
        axis.text.y = element_blank())
ggsave(paste0("Genebubbleplot.pdf"),gg,width = bubbleWidth,height = bubbleHeight)

#UMAP
DimPlot(object = seuset,reduction = "umap",group.by = "ident")