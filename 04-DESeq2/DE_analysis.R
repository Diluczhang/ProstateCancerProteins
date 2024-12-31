#!/usr/bin/env Rscript
suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))

args <- commandArgs(T)

samplefile <- as.character(args[1]) # "1.rawdata\sample.list"
contrastfile <- as.character(args[2]) # "1.rawdata\contrast.list"
countfile <- as.character(args[3]) # "1.rawdata\gene_count.xls"
fpkmfile <- as.character(args[4]) # "1.rawdata\gene_fpkm.xls"
outdir <- as.character(args[5]) # "./"
desfile <- as.character(args[6])
fc <- as.numeric(args[7])
fdr <- as.numeric(args[8])
value <- as.numeric(args[9])
#=============================================================================================#

# 0代表没有description注释文件，默认是0
if (is.na(desfile)){
  desfile <- 0
}

# FC默认是1
if(is.na(fc)){
  fc <- 1
}

# fdr = 1 表明使用FDR校正；fdr = 0 表明不做FDR校正，默认是1
if(is.na(fdr)){
  fdr <- 1
}

# p.adj or pvalue阈值默认是0.05
if(is.na(value)){
  value <- 0.05
}

# samplefile <- "../sample.txt"
# contrastfile <- "../contrasts.txt"
# countfile <- "../gene_count.xls"

#=============================================================================================#

setwd(outdir)
sample_df <- read.table(file = samplefile, sep = "\t", header = F, stringsAsFactors = F)
contrast_df <- read.table(file = contrastfile, sep = "\t", header = F, stringsAsFactors = F, quote = "", colClasses = c("character", "character"))

count_matrix_all <- read.table(file = countfile, sep = "\t", header = T, row.names = 1, stringsAsFactors = F, check.names = F,quote = "")
fpkm_matrix_all <- read.table(file = fpkmfile, sep = "\t", header = T, row.names = 1, stringsAsFactors = F, check.names = F,quote = "")
count_matrix_all <- round(count_matrix_all)
# fpkm_matrix_all <- round(fpkm_matrix_all)

if (desfile != 0){
  filename <- basename(desfile)
  if (filename == "Annotation_summary.xls"){
    des_data <- read.table(file = desfile, sep = "\t", header = T, stringsAsFactors = F, quote = "",fill=TRUE)
    des_data <- des_data[,c(1,3,5,7,9,11,13,14,15)]
    names(des_data) <- c("gene_id", "NR_annotation", "Swissprot_annotation", "Pfam_annotation", "BP_annotation",
                         "MF_annotation", "CC_annotation", "KO_annotation", "Pathway_annotation")
  }else{
    des_data <- read.table(file = desfile, sep = "\t", header = T, stringsAsFactors = F, quote = "",fill=TRUE)
    des_data <- des_data[,1:3]
    names(des_data) <- c("gene_id", "symbol", "description")
  }
}
pwd <- getwd()
for (i in 1:nrow(contrast_df)) {
  info <- plyr::dlply(sample_df, "V1")
  sample_one <- contrast_df[i, "V1"]
  sample_two <- contrast_df[i, "V2"]
  sample_num <- table(sample_df$V1)
  
  group_folder <- paste0(sample_one, "_vs_", sample_two)
  DEGpath <- paste0(group_folder, "/DEG_list")
  if (dir.exists(DEGpath)) {
    setwd(DEGpath)
  }else{
    dir.create(DEGpath, recursive = T)
    setwd(DEGpath)
  }
  sample_index <- c(info[[sample_one]]$V2, info[[sample_two]]$V2)
  condition <- factor(c(rep(sample_one, as.numeric(sample_num[sample_one])), rep(sample_two, as.numeric(sample_num[sample_two]))),
                      levels = c(sample_one, sample_two))
  
  count_matrix <- count_matrix_all[,sample_index]
  fpkm_matrix <- fpkm_matrix_all[,sample_index]
  coldata <- data.frame(row.names = colnames(count_matrix), condition)
  
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design =~ condition)
  dds <- DESeq(dds)
  dds <- dds[rowSums(counts(dds)) > 1,]
  res <- results(dds, contrast = c("condition", sample_one, sample_two), independentFiltering = F)
  
  basemean_A <- round(rowMeans(counts(dds, normalized=T)[,coldata$condition == sample_one], na.rm = T), 2)
  basemean_B <- round(rowMeans(counts(dds, normalized=T)[,coldata$condition == sample_two], na.rm = T), 2)

  count_matrix_res <- count_matrix %>% filter(rownames(count_matrix) %in% rownames(res))
  fpkm_matrix_res <- fpkm_matrix %>% filter(rownames(fpkm_matrix) %in% rownames(res))
  # A_count_normalized <- counts(dds, normalized=T)[,coldata$condition == sample_one]

  resdata <- cbind(row.names(res), basemean_A, basemean_B,data.frame(res, stringsAsFactors = F)[,c(2,5,6)], count_matrix_res)
  resdata_fpkm <- cbind(row.names(res),data.frame(res,stringsAsFactors = F)[,2], fpkm_matrix_res)

  resdata <- dplyr::arrange(resdata, padj)
  resdata[,4] <- round(resdata[,4], 4)
  resdata[,5:6] <- signif(resdata[,5:6], 5)

  # resdata_fpkm <- dplyr::arrange(resdata_fpkm, log2FoldChange)
  # resdata_fpkm[,2] <- signif(resdata[,2],4)
  colnames(resdata_fpkm)[1:2] <- c("gene_id","log2FoldChange")

  colnames(resdata)[1:3] <- c("gene_id", paste0(sample_one, "_readcount"), paste0(sample_two, "_readcount"))
  
  if (desfile != 0){
    filename <- basename(desfile)
    if (filename == "Annotation_summary.xls"){
      resdata$gene_id <- as.character(resdata$gene_id)
      des_data$gene_id <- as.character(des_data$gene_id)
      resdata <- left_join(resdata, des_data, by = c("gene_id" = "gene_id"))
      resdata[,7:ncol(resdata)][is.na(resdata[,7:ncol(resdata)])] <- "-"
    }else{
      resdata$gene_id <- as.character(resdata$gene_id)
      des_data$gene_id <- as.character(des_data$gene_id)
      resdata <- left_join(resdata, des_data, by = c("gene_id" = "gene_id"))
      resdata$symbol <- ifelse(is.na(resdata$symbol) | resdata$symbol == "", "-", resdata$symbol)
      resdata$description <- ifelse(is.na(resdata$description) | resdata$description == "", "-", resdata$description)
    }
  }
  
  if (fdr == 1){
    resdata_DEG <- resdata[(resdata$padj < value & !is.na(resdata$padj)) & (resdata$log2FoldChange > fc | resdata$log2FoldChange < -fc), ]
    resdata_UP <- resdata_DEG[resdata_DEG$log2FoldChange > fc,]
    resdata_DOWN <- resdata_DEG[resdata_DEG$log2FoldChange < -fc,]
  }else{
    resdata_DEG <- resdata[(resdata$pvalue < value & !is.na(resdata$pvalue)) & (resdata$log2FoldChange > fc | resdata$log2FoldChange < -fc), ]
    resdata_UP <- resdata_DEG[resdata_DEG$log2FoldChange > fc,]
    resdata_DOWN <- resdata_DEG[resdata_DEG$log2FoldChange < -fc,]
  }
  
  # resdata_DEG <- resdata[(resdata$padj < 0.05 & !is.na(resdata$padj)) & (resdata$log2FoldChange > 1 | resdata$log2FoldChange < -1), ]
  # resdata_UP <- resdata_DEG[resdata_DEG$log2FoldChange > 1,]
  # resdata_DOWN <- resdata_DEG[resdata_DEG$log2FoldChange < -1,]
  
  output <- paste0(sample_one, "_vs_", sample_two)
  write.table(resdata_fpkm, file = paste0(output, "_FPKM.xls"), sep = "\t", row.names = F, col.names = T, quote = F) 
  write.table(resdata, file = paste0(output, ".result.xls"), sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(resdata_DEG, file = paste0(output, ".DEG.xls"), sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(resdata_UP, file = paste0(output, ".DEG_up.xls"), sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(resdata_DOWN, file = paste0(output, ".DEG_down.xls"), sep = "\t", row.names = F, col.names = T, quote = F)
  
  #=============================================================================================#
  #produce cluter.txt for every group
  
  clusterpath <- "../DEG_cluster"
  if (dir.exists(clusterpath)) {
    setwd(clusterpath)
  }else{
    dir.create(clusterpath, recursive = T)
    setwd(clusterpath)
  }
  
  cluster_matrix <- counts(dds, normalized=T)
  cluster_matrix <- data.frame(geneid = row.names(cluster_matrix), cluster_matrix, stringsAsFactors = F, check.names = F)
  cluster_matrix_DEG <- filter(cluster_matrix, geneid %in% resdata_DEG$gene_id)
  
  write.table(cluster_matrix_DEG, file = "cluster.txt", sep = "\t", row.names = F, col.names = T, quote = F)
  
  setwd(pwd)
}
