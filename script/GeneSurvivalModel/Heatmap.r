Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)

Heatmap<-function(data.train,Cox_genes_OS_pValue,save_path){

	#=========================画热图==================
	#只展示预后显著的泛素化相关差异基因18个

	#读取训练集763个样本的RNA表达谱
	data.train <- read.csv(data.train,row.names = 1,check.names = F)
	exprSet <- data.frame(t(data.train[,-c(1:12)]),check.names = F)
	Cox_genes_OS_pValue <- read.csv(Cox_genes_OS_pValue,header = T)

	Cox_genes_expr<-exprSet[as.character(unique(Cox_genes_OS_pValue$Gene)),]



	# 读入基因集的基因列表
	genes = data.frame(Genes="Genes_GSVA",t(unique(Cox_genes_OS_pValue[,1])))
	
	gs = list()
	for (i in 1:nrow(genes)) {
	gs[[genes[i, 1]]] = t(genes[i, 2:ncol(genes)]) %>% .[grep("\\S", .)]##等于空的不要
	}##转换列表，可以用(unlist(strsplit(test_1[,2],split="/")))

	# 读入基因表达矩阵
	exprSet <- exprSet
	# 这一句就完成了GSVA分析
	gsva_es = gsva(as.matrix(exprSet), gs)
	# 保存到文件
	write.csv(gsva_es, paste(save_path,"Cox_gene_set_gsva.csv",sep=''), quote = F)


	##制作一个分组信息用于注释
	annotation_col = data.frame(t(gsva_es))
	annotation_col$sample = row.names(annotation_col)
	annotation_col = annotation_col[order(annotation_col$Genes_GSVA),]

	Cox_genes_expr = Cox_genes_expr[,row.names(annotation_col)]
	annotation_col$sample = NULL

	# 对应注释信息的颜色
	annColors <- list(Genes_GSVA = c("#FDE725", "#24868D","#440154"))

	#如果注释出界，可以通过调整格子比例和字体修正
	library(pheatmap)
	pheatmap(Cox_genes_expr, #热图的数据
			cluster_rows = TRUE,#行聚类
			cluster_cols = F,#列聚类，可以看出样本之间的区分度
			#color = c("#283285","#DEDEDE","#651430"),
			annotation_col =annotation_col, #标注样本分类
			annotation_colors = annColors,
			annotation_legend=TRUE, # 显示注释
			show_rownames = T,# 显示行名
			show_colnames = F,# 显示列名
			fontsize_col = 16, # 列名大小
			fontsize_row = 16, # 行名大小（也可以不显示）
			fontsize = 16, # 图的基本字体大小，目的是调低图例的大小
			
			scale = "row", #以行来标准化，这个功能很不错
			main = "COX_DEG_expression Heatmap",# 图的标题
			color =colorRampPalette(c("blue","white","red"))(200),
			#legend_labels = c("down","none","up")) # 修改图例各位置的名称
			filename = paste(save_path,"COX_DEG_expression Heatmap.pdf",sep=''),width = 11,height = 6)
}	

