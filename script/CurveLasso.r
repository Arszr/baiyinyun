#发散曲线

library(survival)
library(glmnet)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
#输入文件:表达矩阵和带followup的临床数�?
# 筛选既有表达矩阵又有followup的样�?
CurveLasso <- function(train_data,signature_gene,save_path) {
	trait <- na.omit(read.csv(train_data,check.names = F,stringsAsFactors = F,header = T,row.names = 1)[,-c(4:12)])
	trait = trait[trait$A1_OS!=0,]#生存时间�?0值会报错，去�?0的生存时间就好使�? 

	expr <- data.frame(t(trait[,-c(1:3)]),check.names = F)
	mysurv <- trait[,c(1:3)]
	# 过滤�?8个最优模型的特征基因
	signature_gene<- read.table(signature_gene, sep = "\t", stringsAsFactors = F,header = T, check.names=F,fill=TRUE)
	surv.expr <- expr[signature_gene$x,]#8
	set.seed(10)
	cvfit = cv.glmnet(t(surv.expr), Surv(mysurv$A1_OS,mysurv$event), 
					#10倍交叉验证，非必须限定条件，这篇文献有，其他文献大多没提
					nfold=10,
					family = "cox") 
	pdf(file = paste0(save_path,"Lasso_cox.pdf"), width = 7, height = 6)
	plot(cvfit)
	dev.off()
	#两个lambda值均可采用，具体lambda选值要根据自己实验设计而定�?
	#此处使用`lambda min`
	cvfit$lambda.min #最佳lambda�?
	# 0.003229972
	cvfit$lambda.1se #一倍SE内的更简洁的模型
	# 0.04370307
	fit <- glmnet(t(surv.expr), Surv(mysurv$A1_OS,mysurv$event), 
				family = "cox") 
	#用包自带的函数画�?
	plot(fit, label = TRUE)#发散曲线�?


	#修改画图函数
	#自定义颜�?
	mycol <- rep(c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767",
				"#BD6263","#8EA325","#A9D179","#84CAC0","#F5AE6B","#BCB8D3","#4387B5","#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA",
				"#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767"),2)

	#设置x轴最大�?
	xmax <- 2

	plotCoef_plus <- function (beta, norm, lambda, df, dev, label = FALSE, legend = FALSE, xvar = c("norm","lambda", "dev"), xlab = iname, ylab = "Coefficients", ...) {
		which = sort(unique(beta@i)+1)
		nwhich = length(which)
		switch(nwhich + 1, `0` = {
			warning("No plot produced since all coefficients zero")
			return()
		}, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
		beta = as.matrix(beta[which, , drop = FALSE])
		xvar = match.arg(xvar)
		switch(xvar, norm = {
			index = if (missing(norm)) apply(abs(beta), 2, sum) else norm
			iname = "L1 Norm"
			approx.f = 1
		}, lambda = {
			index = log(lambda)
			iname = "Log Lambda"
			approx.f = 0
		}, dev = {
			index = dev
			iname = "Fraction Deviance Explained"
			approx.f = 1
		})
		dotlist = list(...)
		type = dotlist$type
		
		if (legend){
			#在右侧留出画图例的地�?
			par(xpd = T, mar = par()$mar + c(0,0,0,6))
		}
		
		#修改bty，换个更好看的边框，还可以改成，o / n / 7 / l / c / u / ]
		if (is.null(type)) 
			matplot(index, t(beta), lty = 1, lwd = 2,
					xlab = xlab, ylab = ylab, 
					xlim = c(0, xmax), #设置x轴最大�?
					col = mycol,#线的颜色
					type = "l", cex.lab=1.2, cex.axis=1,
					bty="n", ...)#不画右边�?
		else matplot(index, t(beta), lty = 1, lwd = 2,
					xlab = xlab, ylab = ylab, 
					xlim = c(0, xmax), 
					col = mycol,
					type = "l", cex.lab=1.2, cex.axis=1,
					bty="n", ...)
		atdf = pretty(index)
		prettydf = approx(x = index, y = df, xout = atdf, rule = 2, 
							method = "constant", f = approx.f)$y
		axis(3, at = atdf, labels = prettydf, tcl = NA)
		
		if (label) {
			nnz = length(which)
			xpos = max(index)
			pos = 4
			if (xvar == "lambda") {
			xpos = min(index)
			pos = 2
			}
			xpos = rep(xpos, nnz)
			ypos = beta[, ncol(beta)]
			
			#原函数打印序号，修改为打印基因名
			text(xpos, ypos, paste(rownames(surv.expr)[which]),
				cex = 0.8, #基因名字体大�?
				#基因名的颜色跟线一�?
				col = mycol,
				#如果你不想要彩色的字，就用下面这�?
				#col = "black",
				pos = pos)
		}
		if (legend) {
			#画图�?
			legend("topright",
				inset=c(-0.12,0),#图例画到图外�?
				legend = rownames(surv.expr), #图例文字
				col = mycol, #图例线的颜色，与文字对应
				lwd = 3, #图例中线的粗�?
				cex = 1, #图例字体大小
				bty = "n") #不显示图例边�?
		}
		par(xpd=FALSE)
	}

	plot.glmnet_plus <- function (x, xvar = c("norm", "lambda", "dev"), label = FALSE, legend = FALSE,...) {
		xvar = match.arg(xvar)
		plotCoef_plus(x$beta, lambda = x$lambda, df = x$df, dev = x$dev.ratio, 
						label = label, legend = legend, xvar = xvar, ...)
	}


	#用修改后的函数画�?:在线的旁边显示基因名
	pdf(paste0(save_path,"lasso_name.pdf"),width = 7,height = 6)
	plot.glmnet_plus(fit, label = TRUE, #打印基因�?
					legend = FALSE) #不显示图�?

	#在图上画虚线
	#你想用哪个cutoff，就在“v = ”写上相应的数字
	#此处以lambda.min作为cutoff
	abline(v = cvfit$lambda.min, lty = 3, #线的类型，可以改�?0, 1, 2, 3, 4, 5, 6
		lwd = 2, #线的粗细
		col = "black") #线的颜色
	dev.off()


	coef.min = coef(cvfit, s = "lambda.min") 
	coef.min
	#提取选中的基因名
	active.min = which(coef.min != 0)
	geneids <- rownames(surv.expr)[active.min]
	geneids
	#提取选中的基因对应的coefficient
	index.min = coef.min[active.min]
	index.min
	#输出到文�?
	combine <- cbind(geneids, index.min)
	write.csv(combine,paste0(save_path,"gene_index.csv"))

	#输出用于nomogram作图的文�?:将纳入signature的变量拟合成一个变量，作为nomogram的输�?
	signature <- as.matrix(t(surv.expr[geneids,])) %*% as.matrix(index.min) 
	summary(signature)
	colnames(signature)[1] <- "lasso_Score"
	row.names = row.names(surv.expr)
	write.table(signature,paste0(save_path,"lasso_Score.txt"),row.names = T, quote = F,sep="\t")

}
# CurveLasso('./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/train.csv','./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/signature_gene.txt','./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/')
