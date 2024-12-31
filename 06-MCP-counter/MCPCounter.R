library(MCPcounter) # ?MCPcounter.estimate
library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggpubr)

# TPM
dataTpm <- read.table("1.rawdata/PRAD_RNA_TPM.txt",header = T,row.names = 1,sep = "\t",quote = "")
mcpTpm <- MCPcounter.estimate(dataTpm,featuresType="HUGO_symbols")
write.xlsx(as.data.frame(mcpTpm),"MCPcounter_result.xlsx",rowNames = T,colNames = T)
# 样本蛋白亚型
group <- read.xlsx("1.rawdata/分组信息-Protein.50.log2.Combat.T.Top25.6.3.sub.xlsx")
# plot
McpCounterBoxplot <- function(MCPCounter.Result = NULL,out = NULL,group = NULL){
  groupF <- group %>% filter(sample %in% colnames(MCPCounter.Result))
  groups <- split(groupF$sample, groupF$group)
  MCPCounter.Result <- MCPCounter.Result[,c(groups$Subgroup1,
                                            groups$Subgroup2,
                                            groups$Subgroup3)]
  df <- rbind(c(rep("Group1",length(groups$Subgroup1)),
                rep("Group2",length(groups$Subgroup2)),
                rep("Group3",length(groups$Subgroup3))),
              MCPCounter.Result)
  row.names(df)[1] <- "Group"
  df2 <- data.frame()
  for(i in 1:nrow(MCPCounter.Result)){
    tem <- data.frame(group = df[1,] %>% t() %>% as.vector(),
                      value = df[i + 1,] %>% t() %>% as.vector() %>% as.numeric() %>% log2(),
                      cell = rep(row.names(MCPCounter.Result)[i],ncol(MCPCounter.Result)))
    df2 <- rbind(df2,tem)
  }
  p <- ggplot(data = df2) +
    geom_boxplot(mapping = aes(x = cell, y = value, fill = group)) +
    theme_bw() +
    theme(plot.title = element_text(size = 20),
          axis.text.x = element_text(size = 13, angle = 60, hjust = 1),
          axis.text.y = element_text(size = 13),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20)) +
    labs(x = 'Cell Type', y = "log2(Estimated Propotion)") +
    stat_compare_means(aes(x = cell, y = value, group = group),
                       method = "kruskal.test",
                       label = "p.signif", 
                       label.y = max(df2$value) + 0.1)
  ggsave(filename = paste0(out,"/MCPCounter_boxplot.png"),plot = p,width = 8,height = 6,units = "in",dpi = 300,limitsize = FALSE)
  ggsave(filename = paste0(out,"/MCPCounter_boxplot.pdf"),plot = p,width = 8,height = 6,units = "in")                  
}
McpCounterBoxplot(MCPCounter.Result = mcpTpm,
                  out = "./",
                  group = group)
########################################################################
MCPCounter.Result <- read.xlsx("MCPcounter_result.xlsx",rowNames = T,colNames = T)
group <- read.xlsx("1.rawdata/分组信息-Protein.50.log2.Combat.T.Top25.6.3.sub.xlsx")
groupF <- group %>% filter(sample %in% colnames(MCPCounter.Result))
groups <- split(groupF$sample, groupF$group)
MCPCounter.Result <- MCPCounter.Result[,c(groups$Subgroup1,
                                          groups$Subgroup2,
                                          groups$Subgroup3)]
testRes <- data.frame()
for(i in row.names(MCPCounter.Result)){ # i <- 1
  kruskalRes <- kruskal.test(list(Subgroup1 = as.numeric(MCPCounter.Result[i,groups$Subgroup1]),
                                  Subgroup2 = as.numeric(MCPCounter.Result[i,groups$Subgroup2]),
                                  Subgroup3 = as.numeric(MCPCounter.Result[i,groups$Subgroup3])))
  testRes[i,"kruskal.test.pvalue"] <- kruskalRes$p.value
  testGro_1Gro_2 <- wilcox.test(as.numeric(MCPCounter.Result[i,groups$Subgroup1]),
                                as.numeric(MCPCounter.Result[i,groups$Subgroup2]))
  testGro_1Gro_3 <- wilcox.test(as.numeric(MCPCounter.Result[i,groups$Subgroup1]),
                                as.numeric(MCPCounter.Result[i,groups$Subgroup3]))
  testGro_2Gro_3 <- wilcox.test(as.numeric(MCPCounter.Result[i,groups$Subgroup2]),
                                as.numeric(MCPCounter.Result[i,groups$Subgroup3]))
  testRes[i,"wilcox.test.pvalue.Subgroup1_Subgroup2"] <- testGro_1Gro_2$p.value
  testRes[i,"wilcox.test.pvalue.Subgroup1_Subgroup3"] <- testGro_1Gro_3$p.value
  testRes[i,"wilcox.test.pvalue.Subgroup2_Subgroup3"] <- testGro_2Gro_3$p.value
}
testRes$kruskal.test.FDR <- p.adjust(testRes$kruskal.test.pvalue,method = "BH")
testRes <- testRes %>% select(kruskal.test.pvalue,kruskal.test.FDR,everything())
testRes <- testRes[order(rownames(testRes)), ]
write.xlsx(testRes,"MCPcounter_result_pvalue.xlsx",rowNames = T,colNames = T)
