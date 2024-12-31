setwd("./8-H650_SELECT_boxplot")

rm(list = ls())
options(stringsAsFactors = FALSE)
###loading data
library(openxlsx)
library(tidyverse)
library(ggplot2)
ParseGene <- function(Group,Data,Gene,filename){
    
    colnames(Group)<- c("Sample","Group")
    
    ParseNMFGroup <- function(Group){
        NumGroup <- length(unique(Group$Group))
        Group.Sample <- list()
        for(i in 1:NumGroup){
            assign(paste0("Group.",i),intersect(Group$Sample[Group$Group == i],colnames(Data)))
            Group.Sample <- c(Group.Sample,list(get(paste0("Group.",i))))
        }
        names(Group.Sample) <- paste0("Subgroup",1:NumGroup)
        NewColumn <- c() # 排序后的样本名（按照亚型）
        GroupNum <- c() # 各亚型的样本数目，用于t检验
        GroupLabel <- c() # 排序后每个样本对应的分组标签，用于anova和热图
        GroupLabel.uniq <- c() # 排序后的分组标签，用于绘制热图时分组颜色的设定
        
        for(i in 1:NumGroup){
            NewColumn <- c(NewColumn,get(paste0("Group.",i)))
            GroupNum <- c(GroupNum,length(get(paste0("Group.",i))))
            GroupLabel <- c(GroupLabel,rep(paste0("Subgroup",i),length(get(paste0("Group.",i)))))
            GroupLabel.uniq <- c(GroupLabel.uniq,paste0("Subgroup",i))
        }
        
        Result <- list(NumGroup = NumGroup,
                       NewColumn = NewColumn,
                       GroupNum = GroupNum,
                       GroupLabel = GroupLabel,
                       GroupLabel.uniq = GroupLabel.uniq,
                       Group.Sample = Group.Sample)
        return(Result)
    }
    Group.info <- ParseNMFGroup(Group)
    Data <- Data[,intersect(Group.info$NewColumn,colnames(Data))]
    
    if(sum(Gene %in% row.names(Data)) == 0){
        print(paste0(Gene," not quantified"))
        quit()
    }
    Data <- Data %>% filter(row.names(.) %in% Gene)
    print(Data)
    
    if(Group.info$NumGroup == 2){
        pvalue <- apply(Data,1,function(x){
            wilcox.test(x[1:Group.info$GroupNum[1]],x[(Group.info$GroupNum[1] + 1):sum(Group.info$GroupNum)],
                        alternative = "two.sided",paired = FALSE)$p.value
        })
    }else{
        pvalue <- apply(Data,1,function(x){
            kruskal.test(x~Group.info$GroupLabel)$p.value
        })
    }
    
    # 计算两两相比的pvalue
    comparisons <- combn(Group.info$GroupLabel.uniq,2, simplify = F)
    pvalue.compare <- NULL
    for(j in 1:length(comparisons)){
        comparisons.1 <- Group.info$Group.Sample[[which(names(Group.info$Group.Sample) == comparisons[[j]][1])]]
        comparisons.2 <- Group.info$Group.Sample[[which(names(Group.info$Group.Sample) == comparisons[[j]][2])]]
        
        Data.tem <- Data[,c(comparisons.1,comparisons.2)]
        pvalue.tem <- apply(Data.tem,1,function(x){
            wilcox.test(x[1:length(comparisons.1)],x[(length(comparisons.1) + 1):(length(comparisons.1) + length(comparisons.2))],
                        alternative = "two.sided",paired = FALSE)$p.value
        })
        
        pvalue.compare <- cbind(pvalue.compare,pvalue.tem)
    }
    colnames(pvalue.compare) <- lapply(comparisons, function(x){
        paste(x,collapse = "vs")
    }) %>% unlist()
    
    library(openxlsx)
    Data.out <- cbind(rbind(Group.info$GroupLabel,Data),
                      c(NA,pvalue),
                      rbind(rep(NA,3),pvalue.compare))
    colnames(Data.out)[ncol(Data.out)-3] <- "kruskal.test pvalue"
    row.names(Data.out)[1] <- "Group"
    write.xlsx(as.data.frame(Data.out),file = paste0(filename,".xlsx"),
               asTable = TRUE,overwrite = TRUE,
               colNames = TRUE,rowNames = TRUE)
    
    # 画图部分
    library(ggsignif)
    library(ggplot2)
    library(dplyr)
    library(reshape2)
    
    PlotFun <- function(df){
        
        # 调整显著性曲线在y轴开始的位置
        y.position <- df$Gene %>% log10() %>% max()
        y.position.lower <- y.position %>% floor()
        if((y.position - y.position.lower) > 0.5){
            y.position <- df$Gene %>% log10() %>% max() %>% ceiling()
        }else{
            y.position <- df$Gene %>% log10() %>% max() %>% ceiling() - 1
        }
        
        # 调整显著差异的步进间距
        # 当y轴跨度大时，步进间距调大
        interval <- df$Gene %>% log10() %>% max() %>% ceiling() - df$Gene %>% log10() %>% min() %>% floor()
        if(interval > 1.5){
            y.step_increase <- 0.02
        }else{
            y.step_increase <- 0.02
        }
        
        p <- ggplot(df,aes(x = GroupLabel,y = log10(Gene)))+
            geom_boxplot(width = 0.5,outlier.alpha = 0)+
            geom_jitter(aes(colour = GroupLabel),shape=16,size=2,
                        position=position_jitter(0.1))+
            theme_bw()+
            theme(plot.title=element_text(size = 20),
                  axis.text.x=element_text(size=13,angle=0),
                  axis.text.y=element_text(size=13),
                  axis.title.x=element_text(size = 20),
                  axis.title.y=element_text(size = 20))+
            labs(title = row.names(Data)[i],x='Group', y= "Expression Level")+
            geom_signif(aes(x = GroupLabel,y = log2(Gene)),
                        y_position = y.position,
                        tip_length = c(0.005,0.005),
                        comparisons = combn(Group.info$GroupLabel.uniq,2, simplify = F),
                        step_increase = y.step_increase,
                        map_signif_level = F,
                        test = wilcox.test)
        return(p)
    }
    
    if(sum(Gene %in% row.names(Data)) == 1){
        df <- as.data.frame(cbind(t(Data),Group.info$GroupLabel))
        colnames(df) <- c("Gene","GroupLabel")
        df$Gene <- as.numeric(df$Gene)
        i <- 1
        p <- PlotFun(df)
        ggsave(filename = paste0(filename,".png"),
               plot = p,width = 6,height = 6,units = "in",dpi = 300,limitsize = FALSE)
        ggsave(filename = paste0(filename,".pdf"),
               plot = p,width = 6,height = 6,units = "in")
        
    }else{
        for(i in 1:nrow(Data)){
            df <- rbind(Data[i,],Group.info$GroupLabel) %>% t() %>% as.data.frame()
            colnames(df) <- c("Gene","GroupLabel")
            df$Gene <- as.numeric(df$Gene)
            df<-na.omit(df)
            p <- PlotFun(df)
            ggsave(filename = paste0(filename,".",row.names(Data)[i],".png"),
                   plot = p,width = 6,height = 6,units = "in",dpi = 300,limitsize = FALSE)
            ggsave(filename = paste0(filename,".",row.names(Data)[i],".pdf"),
                   plot = p,width = 6,height = 6,limitsize = FALSE)
        }
    }
}

h650<-read.xlsx("1.rawdata/附件1_样本_定性.xlsx")
sample<- read.xlsx("1.rawdata/2.代谢数据分型.xlsx")
feature<-read.xlsx("1.rawdata/H650_SELECT.xlsx")
expr<- h650[,sample[,1]]
rownames(expr)<- h650$Name

ParseGene(Group=sample,
          Data = expr,
          Gene = feature$Name,
          filename = "H650_subGroup")



