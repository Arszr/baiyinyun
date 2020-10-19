

options(stringsAsFactors=FALSE, scipen=100)

library(dplyr)
library(magrittr)
library(plyr)

library(stringr)

# exp<-"./web_app/data/disease/exp_data/TCGA-BRCA.txt"
# 表达谱构建
Exp_Conversion <- function(exp_data,name,save_path) {
    library(org.Hs.eg.db)
    # 读取TCGA数据
    GeneID_EXP = read.table(exp_data, sep="\t", row.names=1, header=TRUE,check.names = F)

    ###--- 行列去重 ---###
    # TCGA基因ID的转换、去重
    rownames(GeneID_EXP) = gsub("\\..*", "", rownames(GeneID_EXP)) ##去掉版本号
    Symbol_ENSG = rownames(GeneID_EXP) %>% select(org.Hs.eg.db, ., columns="SYMBOL", keytype="ENSEMBL") %>% na.omit()
    GeneSymbol_EXP = merge(Symbol_ENSG, GeneID_EXP, by.x="ENSEMBL", by.y="row.names") %>% dplyr::select(-ENSEMBL) %>% 
                    aggregate(.~SYMBOL, ., max) %>% .[grep("\\S", .$SYMBOL), ] %>% set_rownames(.$SYMBOL) %>% 
                    dplyr::select(-SYMBOL)

    #========================================================================================
    #管道连接符 分解如下：
    GeneSymbol_EXP <- merge(Symbol_ENSG, GeneID_EXP,by.x="ENSEMBL", by.y="row.names")#根据ID匹配,两个数据框生成一个新的数据框
    GeneSymbol_EXP$ENSEMBL = NULL#删除ID列
    #一个基因对应多个探针,取所有探针的平均值作为该基因的表达值
    GeneSymbol_EXP<-aggregate(.~SYMBOL,GeneSymbol_EXP,max)

    GeneSymbol_EXP<-GeneSymbol_EXP[GeneSymbol_EXP$SYMBOL!='',]#去除无genesymbol的空探针
    GeneSymbol_EXP = data.frame(GeneSymbol_EXP[-grep("///",gene_symbol_exprSet$"SYMBOL"),])#去除含"///"的多个基因
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
    write.csv(pheno, paste(save_path,"pheno.csv",sep = ''))

    # 写出count表达谱

    exprSet = unique_EXP[, pheno$Sample]
    write.csv(exprSet, paste(save_path,name,"_exprSet.csv",sep = ''))


    # 用voom方法标准化表达谱
    library(limma)
    normalize_EXP = exprSet[, pheno$Sample] %>% apply(2, as.numeric) %>% voom(normalize="quantile") %>% .$E %>% 
                    as.data.frame() %>% set_rownames(rownames(exprSet))
    write.csv(normalize_EXP, paste(save_path,"voom_normalize_exprSet.csv",sep = ''))#标准化后的表达谱


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

    write.csv(diff_all, paste(save_path,"limma_DEG_all.csv" ,sep = ''))
    write.csv(diff_0.01,paste(save_path,"limma_DEG_0.01.csv",sep = ''))
    write.csv(diff_0.05,paste(save_path,"limma_DEG_0.05.csv",sep = ''))




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

    write.csv(diff_all,  paste(save_path,"DESeq2_DEG_all.csv" ,sep = ''))
    write.csv(diff_0.01, paste(save_path,"DESeq2_DEG_0.01.csv",sep = ''))
    write.csv(diff_0.05, paste(save_path,"DESeq2_DEG_0.05.csv",sep = ''))

    # 输出标准化表达谱
    normalize_EXP = varianceStabilizingTransformation(as.matrix(DESeq2_EXP))
    write.csv(normalize_EXP, paste(save_path,"DESeq2_normalize_exprSet.csv",sep = ''))

    #提取差异基因与泛素化基因的交集
    DEG_0.05 <- paste(save_path,"limma_DEG_0.05.csv",sep = '')
    gene_set <- "./web_app/data/data/ubiq/UbiqGene.txt"

    DEG_0.05 = read.csv(DEG_0.05, sep = ",", header = TRUE, row.names = 1)
    gene_set <- read.table(gene_set, sep="\t", header=TRUE,check.names = F)#202
    intersect = intersect(gene_set$Gene, rownames(DEG_0.05))#154

    gene_set_DEG = unique(merge(gene_set,DEG_0.05,by.x="Gene",by.y="row.names"))
    write.csv(gene_set_DEG, paste(save_path,"gene_set_DEG.csv",sep = ''))#154
    return(1)
}

#提取差异基因与泛素化基因的交集
Difference <- function(difference_data_all,difference_data,aims_set,pheno,voom_normalize,save_path) {


    DEG_0.05 <- difference_data #"limma_DEG_0.05.csv"
    gene_set <- aims_set #"..\\..\\data\\泛素化Gene.txt"

    DEG_0.05 = read.csv(DEG_0.05, sep = ",", header = TRUE, row.names = 1)
    gene_set <- read.table(gene_set, sep="\t", header=TRUE,check.names = F)#202
    intersect = intersect(gene_set$Gene, rownames(DEG_0.05))#154

    gene_set_DEG = unique(merge(gene_set,DEG_0.05,by.x="Gene",by.y="row.names"))
    gene_dpath=paste(save_path,"gene_set_DEG.csv",sep = '')
    write.csv(gene_set_DEG, gene_dpath)#154


    #差异结果可视化：气泡火山图展示差异结果（将154个泛素化相关差异基因P值最显著的Top10标注基因名称）
    #=======================1、画气泡火山图====================================
    hot<-function(DEG_all,gene_set,path,logFCcut,p,pCut){

        library(ggplot2)
        library(ggrepel)
        library(ggthemes)
        library(gridExtra)

        Sys.setenv(LANGUAGE = "en") #显示英文报错信息
        options(stringsAsFactors = FALSE) #禁止chr转成factor
        # 输入文件
        # 全部基因差异表达分析结果

        #差异结果画火山图
        x <-read.csv(DEG_all, sep = ",", header = TRUE, row.names = 1)
        x$label<- rownames(x)
        # 突出展示感兴趣的基因:选出泛素化相关差异基因TOP10标注基因名称
        selectedGeneID <- read.csv(gene_set,row.names = 1)
        selectgenes <- selectedGeneID[order(selectedGeneID$adj.P.Val)[c(1:10)],]
        #selectgenes <- selectedGeneID[order(selectedGeneID$logFC)[c(1:10,(nrow(selectedGeneID)-9):nrow(selectedGeneID))],]
        selectgenes$label<- selectgenes$Gene
        colnames(selectgenes)[1] = "gsym"

        # 参数设置
        # 点的颜色和虚线的位置都由下面的阈值决定，根据具体需求调整。
        #logFCcut <- 0 #log2-foldchange
        #pCut <- 0.05
        #pvalCut <- 0.05 #P.value
        #adjPcut <- 0.05 #adj.P.Val
        #p <- "adj.P.Val" #选择P.value 或 adj.P.Val

        #置x，y軸的最大最小位置
        xmin <- 0
        xmax <- ceiling(max(-log10(x[,p])))
        ymin <- floor(min(x$logFC))
        ymax <- ceiling(max(x$logFC))

        # 基因名的颜色，需大于等于pathway的数量，这里自定义了足够多的颜色
        mycol <- c("blueviolet","#223D6C","darkgreen","chocolate4","#D20A13","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

        # 开始画图
        # 根据需要，可以把P.Value换成adj.P.Val

        # 簡單的setting for color
        x$color_transparent <- ifelse((x[,p] < pCut & x$logFC > logFCcut), "#ED4F4F", ifelse((x[,p] < pCut & x$logFC < -logFCcut), "#329E3F","grey"))
        # 簡單的setting for size
        x$size <- ifelse((x[,p] < pCut & abs(x$logFC) > logFCcut), 4, 2)

        # Construct the plot object
        p1 <- ggplot(data=x, aes( -log10(adj.P.Val), logFC,label = label, color = pathway)) +
        geom_point(alpha = 0.6, size = x$size, colour = x$color_transparent) +
        
        labs(x=bquote(~-Log[10]~("adj.P")), y=bquote(~Log[2]~"(fold change)"), title="Volcano picture of DEGs") + 
        ylim(c(ymin,ymax)) + 
        xlim(c(xmin, xmax)) + 
        
        #画阈值分界线
        geom_hline(yintercept = c(-logFCcut, logFCcut), color="grey40", 
                    linetype="longdash", lwd = 0.5) + #虚线的形状和粗细
        geom_vline(xintercept = -log10(pCut), color="grey40", 
                    linetype="longdash", lwd = 0.5) +
        
        theme_bw(base_size = 18) + #, base_family = "Times" #修改字体
        theme(legend.position="right",
                panel.grid=element_blank(),
                legend.title = element_blank(),
                legend.text= element_text(color="black",size=20),
                plot.title = element_text(hjust = 0.5,size=21),
                axis.text.x = element_text(color="black", size=20),
                axis.text.y = element_text(color="black", size=20),
                axis.title.x = element_text(color="black", size=20),
                axis.title.y = element_text(color="black", size=20))

        # 突出显示候选基因
        p2 <- p1 + 
        # 在感兴趣的基因外面画个黑色圈
        geom_point(data = selectgenes, alpha = 1, size = 4, shape = 1, 
                    stroke = 1, #圈粗细
                    color = "black") +
        
        # 显示感兴趣的基因的基因名
        scale_color_manual(values = mycol) + 
        geom_text_repel(data = selectgenes, 
                        show.legend = FALSE, #不显示图例
                        size = 5, box.padding = unit(0.35, "lines"), 
                        point.padding = unit(0.3, "lines")) +
        guides(color=guide_legend(title = NULL)) 

        # 显示pathway
        np <- length(unique(selectgenes$pathway))
        (labelsInfo <- data.frame(pathway = names(table(selectgenes$pathway)),
                                col = mycol[1:np]))
        p3 <-p2 + annotation_custom(tableGrob(labelsInfo$pathway, rows = c(rep("", np)), cols = "",
                                        theme = ttheme_minimal(base_colour = labelsInfo$col)),
                            ymin = ymax - 0.5, ymax = ymax, xmin = xmin + 0.5, xmax = xmax/4)
        #保存到pdf文件
        ggsave(paste(path,"Volcano_DEG.pdf",sep=''), p3, width = 8,height = 7)

    }
    #运行function
    hot(difference_data_all,gene_dpath,save_path,1,"adj.P.Val",0.05)




    #=========================2、画热图==================
    limma_DEG_0.05<-read.csv(difference_data,row.names = 1)
    normalize_exprSet<-read.csv(voom_normalize,row.names = 1,check.names = F)
    gene_set_DEG <- read.csv(gene_dpath,row.names = 1)

    DEG_expr<-normalize_exprSet[as.character(unique(gene_set_DEG$Gene)),]

    #GSVA分析
    library(msigdbr)
    library(dplyr)
    library(data.table)
    library(GSVA)
    library(limma)
    library(stringr)
    library(ggplot2)
    library(pheatmap)

    # 读入基因集的基因列表
    kegg = read.csv(gene_dpath,row.names = 1)
    kegg = data.frame(Genes="Genes_GSVA",t(unique(kegg[,1])))
    
    kegg_gene = list()
    for (i in 1:nrow(kegg)) {
    kegg_gene[[kegg[i, 1]]] = t(kegg[i, 2:ncol(kegg)]) %>% .[grep("\\S", .)]##等于空的不要
    }##转换列表，可以用(unlist(strsplit(test_1[,2],split="/")))
    gs = kegg_gene
    # 读入基因表达矩阵
    exprSet <- normalize_exprSet
    # 这一句就完成了GSVA分析
    gsva_es = gsva(as.matrix(exprSet), gs)
    # 保存到文件
    write.csv(gsva_es, paste(save_path,"gene_set_DEG_gsva.csv",sep = ''), quote = F)


    ##制作一个分组信息用于注释
    pheno<-read.csv(pheno,row.names = 1)
    row.names(pheno)=pheno[,1]

    GSVA = data.frame(t(gsva_es))
    annotation_col = merge(pheno,GSVA,by.x = "row.names",by.y = "row.names")
    row.names(annotation_col)=annotation_col[,1]
    annotation_col = annotation_col[,-c(1:2)]
    annotation_col = annotation_col[order(annotation_col$Genes),]
    annotation_col = annotation_col[order(annotation_col$Group),]
    DEG_expr = DEG_expr[,row.names(annotation_col)]

    annotation_col$Group<-factor(annotation_col$Group,levels = c("normal","tumor"))

    # 对应注释信息的颜色
    annColors <- list(Genes_GSVA = c("#283285","#DEDEDE","#651430"))

    #如果注释出界，可以通过调整格子比例和字体修正
    pheatmap(DEG_expr, #热图的数据
            cluster_rows = TRUE,#行聚类
            cluster_cols = F,#列聚类，可以看出样本之间的区分度
            #color = c("#283285","#DEDEDE","#651430"),
            annotation_col =annotation_col, #标注样本分类
            annotation_colors = annColors,
            annotation_legend=TRUE, # 显示注释
            show_rownames = F,# 显示行名
            show_colnames = F,# 显示列名
            fontsize_col = 16, # 列名大小
            fontsize_row = 16, # 行名大小（也可以不显示）
            fontsize = 16, # 图的基本字体大小，目的是调低图例的大小
            
            scale = "row", #以行来标准化，这个功能很不错
            main = "DEG_expression Heatmap",# 图的标题
            color =colorRampPalette(c("blue","white","red"))(200),
            legend_labels = c("down","none","up"), # 修改图例各位置的名称
            filename = paste(save_path,"DEG_expression Heatmap.pdf",sep = ''),
            width = 10,
            height = 7
            )

    return(c(gene_dpath,paste0(save_path,"gene_set_DEG_gsva.csv")))
}
# a <- Exp_Conversion('./web_app/data/disease/exp_data/TCGA-LIHC.txt','TCGA-LIHC','./web_app/temp/Arsz/Ubiquitination/TCGA-LIHC/')

# b <- Difference('./web_app/data/Difference/TCGA-BRCA/limma_DEG_all.csv','./web_app/data/Difference/TCGA-BRCA/limma_DEG_0.05.csv'
#     ,'./web_app/data/data/ubiq/UbiqGene.txt','./web_app/data/Difference/TCGA-BRCA/pheno.csv','./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/'
#     )

