

options(stringsAsFactors=FALSE, scipen=100)

library(dplyr)
library(magrittr)
library(plyr)

library(stringr)

# exp<-"./web_app/data/disease/exp_data/TCGA-BRCA.txt"
# 表达谱构建
exp_data='C:/Users/Administrator/Desktop/XX/TCGA-LIHC.txt'
name=paste0('TCGA-LIHC','.csv')
save_path='C:/Users/Administrator/Desktop/XX/'
uigene='C:/Users/Administrator/Desktop/XX/UbiqGene.txt'
#Exp_Conversion <- function(exp_data,name,save_path) {
  library(org.Hs.eg.db)
  # 读取TCGA数据
  GeneID_EXP = read.table(exp_data, sep="\t", row.names=1, header=TRUE,check.names = F)
  
  ###--- 行列去重 ---###
  # TCGA基因ID的转换、去重
  rownames(GeneID_EXP) = gsub("\\..*", "", rownames(GeneID_EXP)) ##去掉版本号
  Symbol_ENSG = rownames(GeneID_EXP) %>% select(org.Hs.eg.db, ., columns="SYMBOL", keytype="ENSEMBL") %>% na.omit()

  #========================================================================================
  #管道连接符 分解如下：
  GeneSymbol_EXP <- merge(Symbol_ENSG, GeneID_EXP,by.x="ENSEMBL", by.y="row.names")#根据ID匹配,两个数据框生成一个新的数据框
  GeneSymbol_EXP$ENSEMBL = NULL#删除ID列
  #一个基因对应多个探针,取所有探针的平均值作为该基因的表达值
  GeneSymbol_EXP<-aggregate(.~SYMBOL,GeneSymbol_EXP,max)
  
  GeneSymbol_EXP$SYMBOL <-gsub("*.\\\ ///.*","", GeneSymbol_EXP$SYMBOL)#把"*.\\\ ///.*"替换成""
  GeneSymbol_EXP<- GeneSymbol_EXP[ GeneSymbol_EXP$SYMBOL!='',]#去除无genesymbol的空探针
  rownames(GeneSymbol_EXP)<-GeneSymbol_EXP$SYMBOL#将gene symbol作为行名
  GeneSymbol_EXP$SYMBOL=NULL
  #========================================================================================
  
  
  # 对于多次测试的样本（以“Rep”结尾），去重，取平均值
  unique_EXP = GeneSymbol_EXP %>% t() %>% as.data.frame() %>% mutate(Sample=str_sub(rownames(.), 1, 15))
  unique_EXP = aggregate(.~Sample, unique_EXP, max) %>% set_rownames(.$Sample) %>% t() %>% as.data.frame()
  unique_EXP = unique_EXP[-1,]
  ###--- 分组信息 ---###
  # 分组情况（11：normal   01：tumor）
  pheno = data.frame(Sample=colnames(unique_EXP)) %>% 
    mutate(Group=ifelse(substr(Sample, 14, 15)=='11', "normal", ifelse(substr(Sample, 14, 15)=='01', "tumor", NA))) %>%
    na.omit()
  pheno = pheno[order(pheno$Group,decreasing = T),]
  write.csv(pheno, paste0(save_path,"pheno.csv",sep = ''))
  
  # 写出count表达谱
  
  exprSet = unique_EXP[, pheno$Sample]
  write.csv(exprSet, paste0(save_path,name,"_exprSet.csv",sep = ''))
  
  
  # 用voom方法标准化表达谱
  library(limma)
  normalize_EXP = exprSet[, pheno$Sample] %>% apply(2, as.numeric) %>% voom(normalize="quantile") %>% .$E %>% 
    as.data.frame() %>% set_rownames(rownames(exprSet))
  write.csv(normalize_EXP, paste0(save_path,name,sep = ''))#标准化后的表达谱
  
  
  ###======================================== limma包差异分析 ========================================###
  library(limma)
  # 设置对比
  design = model.matrix(~0 + factor(pheno$Group)) %>% set_colnames(levels(factor(pheno$Group))) %>%
    set_rownames(pheno$Sample)
  # 构建差异比较矩阵
  contrast.matrix = makeContrasts(tumor-normal, levels=design)
  # 差异分析，b vs. a
  diff_all = lmFit(normalize_EXP, design) %>% contrasts.fit(contrast.matrix) %>% eBayes() %>% 
    topTable(coef=1, n=Inf, adjust.method="fdr", sort.by="P")
  diff_0.05 = diff_all[diff_all$adj.P.Val < 0.05, ]
  diff_0.01 = diff_all[diff_all$adj.P.Val < 0.01, ]
  
  write.csv(diff_all, paste0(save_path,"limma_DEG_all.csv" ,sep = ''))
  write.csv(diff_0.01,paste0(save_path,"limma_DEG_0.01.csv",sep = ''))
  write.csv(diff_0.05,paste0(save_path,"limma_DEG_0.05.csv",sep = ''))
  
  
  
  
  ###======================================== DESeq2包差异分析 ========================================###
  #if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
  #BiocManager::install("DESeq2")
  library(DESeq2)
  DESeq2_EXP = exprSet[, pheno$Sample] %>% apply(2, as.integer) %>% as.data.frame %>% set_rownames(rownames(exprSet))
  
  # 构建输入矩阵
  input_Matrix = DESeqDataSetFromMatrix(countData=DESeq2_EXP, colData=pheno, design=~Group)
  
  # DESeq2 进行差异分析
  dds = DESeq(input_Matrix, parallel=TRUE)
  
  # 输出差异结果
  diff_all = results(dds, contrast=c("Group", "tumor", "normal")) %>% as.data.frame() %>% na.omit() ##表示RSV-normal
  diff_0.01 = diff_all[diff_all$padj < 0.01, ]
  diff_0.05 = diff_all[diff_all$padj < 0.05, ]
  
  write.csv(diff_all,  paste0(save_path,"DESeq2_DEG_all.csv" ,sep = ''))
  write.csv(diff_0.01, paste0(save_path,"DESeq2_DEG_0.01.csv",sep = ''))
  write.csv(diff_0.05, paste0(save_path,"DESeq2_DEG_0.05.csv",sep = ''))
  
  # 输出标准化表达谱
  normalize_EXP = varianceStabilizingTransformation(as.matrix(DESeq2_EXP))
  write.csv(normalize_EXP, paste(save_path,"DESeq2_normalize_exprSet.csv",sep = ''))
  #=====
  limma_DEG_0.05 = read.csv("C:\\Users\\Administrator\\Desktop\\XX\\limma_DEG_0.05.csv", sep = ",", header = TRUE, row.names = 1)
  Ubiquitin_Gene <- read.table("C:\\Users\\Administrator\\Desktop\\XX\\UbiqGene.txt", sep="\t", header=TRUE,check.names = F)#202
  Ubiquitin_Gene$pathway = "Ubiquitination"
  intersect = intersect(Ubiquitin_Gene$Gene, rownames(limma_DEG_0.05))#154
  
  Ubiquitin_DEG = unique(merge(Ubiquitin_Gene,limma_DEG_0.05,by.x="Gene",by.y="row.names"))
  write.csv(Ubiquitin_DEG, paste0(save_path,"Ubiquitin_DEG.csv"))#154
  
  
