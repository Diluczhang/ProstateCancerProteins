
# 20220829
# 对比3个蛋白亚型患者TMB和SCNA burden差异，给出P值
# 这个脚本解决SCNA burden

rm(list = ls())
options(stringsAsFactors = FALSE)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

load("1.rawdata/Protein.50.log2.Combat.T.Top25.6.RData")

library(dplyr)
library(openxlsx)

D <- read.table("1.rawdata/ALL.focal_data_by_genes.xls.txt",
                sep = "\t",header = TRUE,row.names = 1)
# 调整样本名，使样本名和蛋白组的样本名一致
colnames(D) <- sapply(colnames(D),function(x){
  ifelse(nchar(x) == 6,paste0("A0",substring(x,first = 2)),x)
})
D <- D[,3:ncol(D)]
SCNA.burden <- apply(D,2,function(x){
  tem <- x %>% as.numeric() %>% abs()
  sum(tem > 0.3)/length(tem)
})
SCNA.burden <- data.frame(SCNA.burden = SCNA.burden) %>% t() %>% as.data.frame()


library(NMF)
Group <- predict(Protein.50.log2.Combat.T.Top25.6$fit$`3`)
Group <- data.frame(Sample = names(Group),Group = Group)
row.names(Group) <- Group$Sample

ParseNMFGroup <- function(Group,Data){
  NumGroup <- length(unique(Group$Group))
  Group.Sample <- list()
  for(i in 1:NumGroup){
    assign(paste0("Group.",i),intersect(Group$Sample[Group$Group == i],colnames(Data)))
    Group.Sample <- c(Group.Sample,list(get(paste0("Group.",i))))
  }
  names(Group.Sample) <- paste0("Subroup",1:NumGroup)
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
  
  Result <- list(NumGroup = NumGroup,
                 NewColumn = NewColumn,
                 GroupNum = GroupNum,
                 GroupLabel = GroupLabel,
                 GroupLabel.uniq = GroupLabel.uniq,
                 Group.Sample = Group.Sample)
  return(Result)
}
Group.info <- ParseNMFGroup(Group,Data = SCNA.burden)

SCNA.burden <- SCNA.burden[,Group.info$NewColumn,drop = FALSE]

# 计算3组的pvalue
pvalue <- apply(SCNA.burden,1,function(x){
  kruskal.test(x~Group.info$GroupLabel)$p.value
})

# 计算两两相比的pvalue
# 计算两两相比的pvalue
comparisons <- combn(Group.info$GroupLabel.uniq,2, simplify = F)
pvalue.compare <- NULL
for(j in 1:length(comparisons)){
  comparisons.1 <- Group.info$Group.Sample[[which(names(Group.info$Group.Sample) == comparisons[[j]][1])]]
  comparisons.2 <- Group.info$Group.Sample[[which(names(Group.info$Group.Sample) == comparisons[[j]][2])]]
  
  SCNA.burden.tem <- SCNA.burden[,c(comparisons.1,comparisons.2),drop = FALSE]
  pvalue.tem <- apply(SCNA.burden.tem,1,function(x){
    wilcox.test(x[1:length(comparisons.1)],x[(length(comparisons.1) + 1):(length(comparisons.1) + length(comparisons.2))],
                alternative = "two.sided",paired = FALSE)$p.value
  })
  
  pvalue.compare <- c(pvalue.compare,pvalue.tem)
}
names(pvalue.compare) <- lapply(comparisons, function(x){
  paste(x,collapse = "vs")
}) %>% unlist()

library(openxlsx)
Data.out <- cbind(rbind(Group.info$GroupLabel,SCNA.burden),
                  c(NA,pvalue),
                  rbind(rep(NA,3),pvalue.compare))
colnames(Data.out)[ncol(Data.out)-3] <- "kruskal.test pvalue"
row.names(Data.out)[1] <- "Group"
row.names(Data.out)[2] <- "SCNA.burden"
write.xlsx(as.data.frame(Data.out),file = "Protein.50.log2.Combat.T.Top25.6.Group3-SCNA.burden.xlsx",
           asTable = TRUE,overwrite = TRUE,
           colNames = TRUE,rowNames = TRUE)


# 画图部分
library(ggsignif)
library(ggplot2)
library(dplyr)
library(reshape2)

PlotFun <- function(df,title){
  
  # 调整显著性曲线在y轴开始的位置
  y.position <- df$Gene %>% max()
  y.position.lower <- y.position %>% floor()
  if((y.position - y.position.lower) > 0.5){
    y.position <- df$Gene %>%  max() %>% ceiling()
  }else{
    y.position <- df$Gene %>% max() %>% ceiling() - 1
  }
  
  # 调整显著差异的步进间距
  # 当y轴跨度大时，步进间距调大
  interval <- df$Gene %>% log10() %>% max() %>% ceiling() - df$Gene %>% min() %>% floor()
  if(interval > 1.5){
    y.step_increase <- 0.1
  }else{
    y.step_increase <- 0.1
  }
  
  p <- ggplot(df,aes(x = GroupLabel,y = Gene))+
    geom_boxplot(width = 0.5,outlier.alpha = 0)+
    geom_jitter(aes(colour = GroupLabel),shape=16,size=2,position=position_jitter(0.2))+
    theme_bw()+
    theme(plot.title=element_text(size = 20),
          axis.text.x=element_text(size=13,angle=0),
          axis.text.y=element_text(size=13),
          axis.title.x=element_text(size = 20),
          axis.title.y=element_text(size = 20))+
    labs(title = "SCNA burden",x='Group', y= "SCNA burden (%)")+
    geom_signif(aes(x = GroupLabel,y = Gene),df,
                # y_position = y.position,
                # tip_length = c(0.05,0.05),
                comparisons = combn(Group.info$GroupLabel.uniq,2, simplify = F),
                step_increase = 0.1,
                map_signif_level = F,
                test = wilcox.test)
  return(p)
}

df <- cbind(t(SCNA.burden),Group.info$GroupLabel) %>% as.data.frame()
colnames(df) <- c("Gene","GroupLabel")
df$Gene <- as.numeric(df$Gene) * 100

p <- PlotFun(df)
ggsave(filename = "Protein.50.log2.Combat.T.Top25.6.Group3-SCNA.burden.png",
       plot = p,width = 6,height = 6,units = "in",dpi = 300,limitsize = FALSE)
ggsave(filename = "Protein.50.log2.Combat.T.Top25.6.Group3-SCNA.burden.pdf",
       plot = p,width = 6,height = 6,units = "in")





