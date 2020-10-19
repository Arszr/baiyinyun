library(ranger)
library(randomForest)
library(survival)
library(survminer)
library(randomForestSRC)
options(stringsAsFactors = FALSE) #禁止chr转成factor
# 显示进程
display.progress = function (index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}  
Rsf <- function(train_data,signature_gene_path,save_path) {

	#输入文件:表达矩阵,临床信息
	trait <- na.omit(read.csv(train_data,check.names = F,stringsAsFactors = F,header = T,row.names = 1)[,-c(4:12)])
	#trait = trait[trait$A1_OS!=0,]#生存时间有0值会报错，去掉0的生存时间就好使了 
	# 过滤出8个最优模型的特征基因
	expr <- data.frame(t(trait[,-c(1:3)]),check.names = F)
	signature_gene<- read.table(signature_gene_path, sep = "\t", stringsAsFactors = F,header = T, check.names=F,fill=TRUE)
	exp <- expr[signature_gene$x,]#8
	cli <- trait[,c(1,3)]

	gene.sel <- signature_gene$x
	#加权随机森林进一步降维
	dt.rf <- cbind.data.frame(cli,t(exp))
	dt.rf <- dt.rf[,setdiff(colnames(dt.rf),"A1_OS")] # 二分类分类树，结局考虑overall survival

	ntree <- 1000
	mtry <- floor(sqrt(length(gene.sel)))
	weight <- 0.999999 # 算法要求无法取1但无限接近1，该参数表示变量被选择的概率，原文为100%。
	seed <- 10

	set.seed(seed) # 设置外部种子保证结果可重复
	surv.rf <- ranger(formula = event ~ ., 
					data = dt.rf,
					num.trees = ntree,
					mtry = mtry,
					importance = "impurity",
					split.select.weights = rep(weight,length(gene.sel)))

	# 变量重要性按照降序排列
	var.imp <- sort(ranger::importance(surv.rf),decreasing = T)

	#逐步回归筛选变量:根据重要性逐步纳入基因，每一次做当前基因组合下的随机森林并计算oob，取oob达到最小时的组合。
	ntree <- 1000 # 树的数目
	nPerm <- 50 # 扰动次数，一般为50
	var.now <- oob <- c()
	for (var in names(var.imp)) {
	cat(var,"\n")
	set.seed(seed)
	var.now <- c(var.now,var)
	swsfs.dt <- dt.rf[,c("event",var.now)]
	swsfs.dt$event <- factor(swsfs.dt$event)
	model_RF <- randomForest(event ~ ., 
							data = swsfs.dt,
							ntree = ntree, 
							nPerm = nPerm, 
							mtry = floor(sqrt(ncol(swsfs.dt))), 
							proximity = T,
							importance = F)
	oob <- c(oob,model_RF$err.rate[ntree,1])
	}
	names(oob) <- names(var.imp)
	signature <- var.imp[1:which.min(oob)] # 取oob达到最小时候的gene组合作为最终signature
	write.table(data.frame(signature = names(signature),importance = as.numeric(signature),stringsAsFactors = F),
				paste0(save_path,"signature with importance.txt"),sep = "\t",row.names = F,quote = F)

	#开始画图
	pdf(paste0(save_path,"oob from swsfs.pdf"),width = 6,height = 6)
	par(bty = "o", mgp = c(1.5,.33,0), mar = c(3,4,1,2),las = 1, tcl = -.25)
	plot(1:length(gene.sel),oob,
		xlab = "Number of Genes",
		ylab = "",
		type = "l",
		lty = 4,
		col = "red", # 可修改线的颜色
		cex = 1.5)
	mtext("OOB error rate",side = 2,line = 2.5,las = 3) # 添加y标签
	points(1:length(gene.sel),oob,# 加圆圈
		col = "red", # 可修改圆圈颜色
		pch = 19)
	abline(v = which.min(oob),lty = 2, col = "red") # 找到oob最小的位置添加垂直虚线
	dev.off()





	#条形图
	#输入文件：表达矩阵、临床信息
	trait <- na.omit(read.csv(train_data,check.names = F,stringsAsFactors = F,header = T,row.names = 1)[,-c(4:12)])
	#trait = trait[trait$A1_OS!=0,]#生存时间有0值会报错，去掉0的生存时间就好使了 
	# 加载数据（来自FigureYa128Prognostic 迭代Lasso分析后的特征基因）
	# 过滤出8个最优模型的特征基因
	expr <- data.frame(t(trait[,-c(1:3)]),check.names = F)
	signature_gene<- read.table(signature_gene_path, sep = "\t", stringsAsFactors = F,header = T, check.names=F,fill=TRUE)
	exp <- expr[signature_gene$x,]#8
	cli <- trait[,c(1,3)]

	#随机森林进一步降维
	gene.sel <- signature_gene$x
	tmp <- exp[gene.sel,]; 
	dt.rf <- cbind.data.frame(cli,t(tmp))

	surv.rf <- rfsrc(Surv(A1_OS, event) ~ ., 
					data = dt.rf,
					ntree = 1000,
					nodesize = 10,##该值建议多调整，使得到的重要性基因数与加权随机森林得到的重要性基因数一致
					splitrule = 'logrank',
					importance = T,
					proximity = T,
					forest = T,
					seed = 10)

	#开始画图：用包里自带的函数直接出图
	pdf(paste0(save_path,"error rate10.pdf"),width = 12,height = 7)
	plot(surv.rf)
	dev.off()

	# 相对重要性（relative importance）；其实就是把重要性划分到0-1区间内
	range01 <- function(x){(x-min(x))/(max(x)-min(x))}
	raw.imp <- surv.rf$importance
	rel.imp <- range01(raw.imp) # calculate relative importance
	# 输出重要性矩阵
	imp.res <- data.frame(gene = names(raw.imp),
						raw.importance = raw.imp,
						rel.importance = rel.imp,
						stringsAsFactors = F)
	write.csv(imp.res[order(imp.res$rel.importance,decreasing = T),],paste0(save_path,"importance result.csv"),row.names = F,quote = F)




	## 获得每个样本的 riskscore，进一步进行下游分析

	score_t <- data.frame(dt.rf[,c(1,2)],SRC_Score=surv.rf$predicted)
	score_t$time = as.numeric(score_t$A1_OS)/365
	write.csv(score_t, paste0(save_path,"randomForestSRC_riskscore.csv"), quote = F)

	cut <- surv_cutpoint(score_t,'A1_OS','event','SRC_Score')


	## 生存分析
	cat <- surv_categorize(cut)
	cat$time = as.numeric(cat$A1_OS)/365
	fit <- survfit(Surv(time,event)~SRC_Score,cat)
	mytheme <- theme_survminer(font.legend = c(16,"plain", "black"),
							font.x = c(16,"plain", "black"),
							font.y = c(16,"plain", "black")) ## 自定义主题
	ggsurvplot(fit,cat,
			palette = 'jco',
			size=1.5,
			pval=T,
			legend.labs=c("High","Low"), 
			legend.title='randomForestSRC_Score',
			xlab="Time(years)",
			ylab='Overall survival',
			ggtheme = mytheme)
	ggsave(file=paste0(save_path,"Survival Analysis of randomForestSRC.pdf"),width=8, height=6)

}

# Rsf('./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/train.csv','./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/signature_gene.txt','./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/')