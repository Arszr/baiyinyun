options(stringsAsFactors = F)

RelatedBubbles <- function(cox_case_exp,pdata,save_path) {

	#读取表达谱
	Cox_case_exp <- read.csv(cox_case_exp,row.names=1,header = T, check.names = F)
	pdata = read.csv(pdata,header = T, fill=TRUE,check.names = F)
	case_exp1 = data.frame(t(Cox_case_exp),check.names = F)
	case_exp1$sample = substr(row.names(case_exp1),1,12)

	#相关性分析
	library(Hmisc)
	Cor_result_Final = data.frame()
	for(j in 7:ncol(pdata)){
	pdata_EXP = na.omit(merge(pdata[,c(1,j)],case_exp1,by.x="A0_Samples",by.y="sample"))
	row.names(pdata_EXP) = pdata_EXP[,1]
	
	for(i in 3:ncol(pdata_EXP)){
	cordata = pdata_EXP[,c(2,i)]
	cor_result = rcorr(as.matrix(cordata))
	flattenCorrMatrix <- function(cormat, pmat) {
		ut <- upper.tri(cormat) 
		data.frame( row = rownames(cormat)[row(cormat)[ut]],column = rownames(cormat)[col(cormat)[ut]], cor =(cormat)[ut], p = pmat[ut])
	}
	result = flattenCorrMatrix(cor_result$r, cor_result$P)
	Cor_result_Final = rbind(Cor_result_Final,result)
	}
	}
	#Cor_result_Final = na.omit(Cor_result_Final)
	colnames(Cor_result_Final)[1:2] = c("phenotype","Gene")
	bubble <- paste(save_path,"phenotype_CoxGene_cor.txt",sep = '')
	write.table(Cor_result_Final,bubble,row.names = F,quote = F,sep="\t")#生成相关性系数文件#126



	#==============================================
	# bubbles气泡热图
	#数据整理:第一列是Y轴(phenotype)，第二列是X轴(Gene)，后面几列依次是cor,pvalue，logFC, count值等需要展示的数据
	pm <- read.table(bubble,sep = "\t",comment.char = "#", stringsAsFactors = F,header = T,check.names=F,fill=TRUE, quote="")

	#开始画图
	##melt dataframe and draw figues
	library(reshape2)
	library(ggplot2)

	x=ggplot(pm,aes(phenotype,Gene,size=-log10(p)))+
	geom_point(shape=21,aes(fill=cor),position =position_dodge(0))+
	theme_minimal()+
	scale_size_continuous(range=c(1,12))+
	scale_fill_gradientn(colours=c("#2381B3","white","#F0E366"),guide="legend")+
	theme(legend.text= element_text(size=18),
			legend.title= element_text(color="black", size=20),
			axis.text.x = element_text(color="black",angle=35,hjust=1, vjust=1, size=18),
			axis.text.y = element_text(color="black",size=18),
			axis.title.x = element_text(color="black", size=20),
			axis.title.y = element_text(color="black", size=20),
			#legend.position = "bottom",
			legend.box = "vertical",
			panel.grid =element_blank(),#不显示网格线
			legend.margin=margin(t= 0, unit='cm'),
			legend.spacing = unit(0,"in"))
	ggsave(file=paste(save_path,"Cor_bubbles.pdf",sep = ''),width = 10, height = 11)

	return(bubble)
}
# RelatedBubbles(
#             './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Cox_case_exp.csv',
#             './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/pdata.csv',
#             './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/'
#         )
