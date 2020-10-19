Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
library(glmnet)
library(plyr)
library(dplyr)
library(tableone)
library(tidyverse)
library(stringr) 
library("forestplot")
library(caret)
library(survival)
library(magrittr)

# clinical='./web_app/data/disease/clinical/TCGA-BRCA.csv'
# gene_set='./web_app/data/disease/exp_data/TCGA-BRCA.csv'
# set_deg='./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Fig1/gene_set_DEG.csv'
# save_path='./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Fig1/'

# 输入clinical临床信息，表达谱，差异表达基因，保存路径
SingleFactor <-  function(clinical,gene_set,set_deg,save_path){

	Clinical = read.table(clinical,sep = ",",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE,check.names = F,quote ="")
	Clinical_col=colnames(Clinical)
    Gene_Exp <- read.csv(gene_set,row.names = 1,check.names=F)
    gene_set_DEG <- read.csv(set_deg,row.names = 1,check.names=F)

    Gene_Exp1 = Gene_Exp[as.character(unique(gene_set_DEG$Gene)),]
    Gene_Exp1 =data.frame(t(Gene_Exp1),check.names = F)
    Gene_Exp1$sample = row.names(Gene_Exp1)
    case<-Gene_Exp1[!grepl("\\-11",Gene_Exp1$sample),]
    case$sample=str_split(case$sample,'\\-0',simplify = T)[,1]

    trait<-merge(Clinical,case,by.x=Clinical_col[1],by.y = "sample")
    
	
	trait <- trait[!(trait[Clinical_col[2]] == 0), ]
	trait <- na.omit(trait)
	row.names(trait)=trait[,1]
	trait_col=colnames(trait)
    data_all <- paste(save_path,'data_all.csv',sep = '')
    write.csv(trait,file = data_all,row.names = F)#1090


    ##生成随机按7:3比例分组

    in_train = createDataPartition(trait[,trait_col[3]], p=0.7, list=FALSE)
    data.train = trait[in_train, ]
    data.test = trait[-in_train, ]
    train_path=paste(save_path,'train.csv',sep = '')
    test_path=paste(save_path,'test.csv',sep = '')
    write.csv(data.train, file = train_path, row.names = F)#764
    write.csv(data.test, file = test_path, row.names = F)#326
#     #单因素
#         # ===========================================================================
#     #1.data.train单因素Cox回归
    s <- length(Clinical_col)+1
    pdata1<-data.train[,c(2,3,s:ncol(data.train))]
	pdata1_col <- colnames(pdata1)
	pdata1[pdata1_col[2]]  <-  ifelse(pdata1[pdata1_col[2]]=="Alive",0,ifelse(pdata1[pdata1_col[2]]=="Dead",1,NA))
    #for(i in colnames(pdata1[,3:ncol(pdata1)])){
    #  pdata1[,i]<-ifelse(pdata1[,i]<median(pdata1[,i]),0,1)
    #} 
    #用以下循环进行单因素Cox回归分析，死亡风险预测
    Coxoutput_train=data.frame()
    for(i in colnames(pdata1)[3:ncol(pdata1)]){
    pdata2 <- pdata1[,c(pdata1_col[1], pdata1_col[2],i)]
    cox <- coxph(Surv(pdata2[,1], pdata2[,2]) ~ pdata2[,3], data = pdata2)
    coxSummary = summary(cox)
    Coxoutput_train=rbind(Coxoutput_train,cbind(factors=i,
                                    beta=round(coxSummary$coefficients[,1],3),
                                    z=coxSummary$coefficients[,"z"],
                                    pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                    CI95=paste0(round(coxSummary$conf.int[,3:4],3),collapse = '-'),
                                    HR=coxSummary$coefficients[,"exp(coef)"],
                                    lower=coxSummary$conf.int[,3],
                                    upper=coxSummary$conf.int[,4]))
    }
    Coxoutput_train$pvalue <- as.numeric(as.vector(Coxoutput_train$pvalue))
    Coxoutput_train <- Coxoutput_train[order(Coxoutput_train$pvalue),]%>% filter(pvalue < 0.05) #筛选P值显著的性状和基因
    #保存到文件
    Single_train_path=paste(save_path,'CoxSingle_train.csv',sep = '')
    write.csv(Coxoutput_train, Single_train_path, row.names = F)#18

    #=====================================
    #开始画图:森林图：
    
    #options(digits = 3)
    #读取输入文件，输入数据一般包括分类、样本数、风险比及置信区间（上限及下限）等，需要注意的是输入文件的布局(包括文字缩进)将展示在最后结果中，所见即所得。
    data <- read.csv(Single_train_path, stringsAsFactors=FALSE)

    #将要在图中展示的文本
    tabletext <- cbind(c("\nfactors",NA,NA,data$factors,NA),
                    c("\nBeta",NA,NA,data$beta,NA),
                    c("\npvalue",NA,NA,format(round(data$pvalue,3),nsmall=3),NA),
                    c("Hazard Ratio\n(95% CI)",NA,NA,paste(format(round(data$HR,3),nsmall=3)," (", data$CI95, ")",sep=""), NA),NA)
    #用默认参数画图
#     forestplot(labeltext=tabletext,mean=c(NA,NA,1,data$HR,NA),
#             lower=c(NA,NA,1,data$lower,NA),upper=c(NA,NA,1,data$upper,NA))

    #仔细调整参数
    x1 <- as.character(nrow(data)+5)
    #第三行顶部加黑线，引号内数字标记行位置
    lines <- paste('"3" = gpar(lwd=2, col="black"),', '"',x1,'"', '= gpar(lwd=2, col="black")',sep = '')
    lines1 <- paste0("list(",lines,")")
    linex <- eval(parse(text=lines1))

    pdf(file = paste(save_path,"SigCox_OS_forestplot.pdf",sep = ''), width = 10, height = 10)
    forestplot(labeltext=tabletext, #图中的文本
            mean=c(NA,NA,1,data$HR,NA),#HR
            lower=c(NA,NA,1,data$lower,NA), #95%置信区间下限
            upper=c(NA,NA,1,data$upper,NA),#95%置信区间上限
            title="Univariate Cox regression analysis",
            graph.pos=4,#图在表中的列位置
            graphwidth = unit(.4,"npc"),#图在表中的宽度比例
            fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
            col=fpColors(box="steelblue", lines="black", zero = "black"),#box颜色
            #boxsize=c(NA,NA,NA,data$Percent,NA)/75,#box大小根据样本量设置
            lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型
            zero=1,#zero线横坐标
            lwd.zero=2,#zero线宽
            grid = structure(c(data[1,]$Point.Estimate), gp = gpar(col = "black", lty=2,lwd=2)),#虚线(可多条)及其横坐标、颜色、线宽
            #xticks = c(floor(min(data$lower)), 1, 1+(ceiling(max(data$upper))-1)/2, ceiling(max(data$upper))),#横坐标刻度根据需要可随意设置
            lwd.xaxis=2,#X轴线宽
            xlab="",#X轴标题
            hrzl_lines=linex,#最后一行底部加黑线,""中数字为nrow(data)+5
            txt_gp=fpTxtGp(label=gpar(cex=1.2),#各种字体大小设置
                            ticks=gpar(cex=1.2),
                            xlab=gpar(cex = 1.2),
                            title=gpar(cex = 1.5)),
            #is.summary = c(T,rep(F,27)),#首行字体类型设置
            lineheight = unit(1,"cm"),#固定行高
            #align=c("l","c","c"),#每列文字的对齐方式，偶尔会用到
            #cex=10,
            colgap = unit(1,"cm"),#列间隙
            #mar=unit(rep(0.5, times = 4), "cm"),#图形页边距
            new_page = T)#是否新页
    dev.off()
#     #保存图片"Cox_OS_forestplot.pdf",width = 10,height = 10

    return(c(train_path,Single_train_path,test_path,data_all))

}


# a <- SingleFactor('./web_app/data/disease/clinical/TCGA-BRCA.csv','./web_app/data/disease/exp_data/TCGA-BRCA.csv','E:/python_code/Baiyinyun/web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Fig1/gene_set_DEG.csv','./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Fig1/')
# MultipleFactor("./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/train.csv","./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/CoxSingle_train.csv",'./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/')
# print(a)