library(survival) # 生存分析
library(glmnet) # LASSO回归
library(forestplot)

options(stringsAsFactors = FALSE) #禁止chr转成factor
#自定义函数显示进程
display.progress = function (index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}  
#输入文件:表达矩阵和生存信息
# 加载训练集
CoxIndependent <- function(train_data,save_path) {

trait <- na.omit(read.csv(train_data,check.names = F,stringsAsFactors = F,header = T,row.names = 1)[,-c(4:12)])
trait = trait[trait$A1_OS!=0,]#生存时间有0值会报错，去掉0的生存时间就好使了 

# 多变量cox+LASSO惩罚：用18个单因素良好的基因做多因素分析
# 读取18个与OS预后良好的基因
Coxoutput_train<- read.csv(paste0(save_path,"CoxSingle_train.csv"), stringsAsFactors=FALSE)
trait = trait[,c("A1_OS","event",Coxoutput_train$factors)]
# 保留既有表达数据又有生存数据的sample
tcga.expr <- data.frame(t(trait[,-c(1:3)]),check.names = F)
tcga.surv <- trait

set.seed(10)
cvfit = cv.glmnet(x = t(as.matrix(tcga.expr)), 
                  y = Surv(tcga.surv$A1_OS,tcga.surv$event),
                  nfold = 10,
                  family = "cox") 

myCoefs <- coef(cvfit, s="lambda.min");
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]

# 计算训练集riskscore（rs）
rs.tcga <- apply(t(tcga.expr[lasso_fea,]), 1, function(x) {x %*% myCoefs@x})
tcga.surv$riskscore <- rs.tcga[rownames(tcga.surv)]
write.csv(tcga.surv[,c(1,2,21)],file = paste0(save_path,'Cox_riskscore.csv'),row.names = T)#745

#--------------------#
# 计算rs的预后独立性 #
# 训练集预后独立性
mulcox.tcga <- summary(coxph(Surv(A1_OS, event) ~ ., data = tcga.surv))
mulcox.tcga <- data.frame(variable = rownames(mulcox.tcga$conf.int),
                          beta=round(mulcox.tcga$coefficients[,1],3),
                          z = mulcox.tcga$coefficients[,"z"],
                          pvalue = mulcox.tcga$coefficients[,5],
                          HR = mulcox.tcga$conf.int[,1],
                          CI95=paste0(round(mulcox.tcga$conf.int[,3],3),'-',round(mulcox.tcga$conf.int[,4],3)),
                          lower = mulcox.tcga$conf.int[,3],
                          upper = mulcox.tcga$conf.int[,4],
                          stringsAsFactors = F)
mulcox.tcga <- mulcox.tcga[as.character(Coxoutput_train$factors),]#按照单因素顺序排序
write.csv(mulcox.tcga,file = paste0(save_path,'MulCox.csv'),row.names = F)#18


#绘制森林图:
library("forestplot")
options(digits = 3)
#读取输入文件，输入数据一般包括分类、样本数、风险比及置信区间（上限及下限）等，需要注意的是输入文件的布局(包括文字缩进)将展示在最后结果中，所见即所得。
data <- read.csv(paste0(save_path,"MulCox.csv"), stringsAsFactors=FALSE)

#将要在图中展示的文本
tabletext <- cbind(c("\nVariable",NA,NA,data$variable,NA),
                   c("\nBeta",NA,NA,data$beta,NA),
                   c("\npvalue",NA,NA,format(round(data$pvalue,3),nsmall=3),NA),
                   c("Hazard Ratio\n(95% CI)",NA,NA,paste(format(round(data$HR,3),nsmall=3)," (", data$CI95, ")",sep=""), NA),NA)

pdf(paste0(save_path,"MulCox_OS_forestplot1.pdf"), width = 10, height = 10)
forestplot(labeltext=tabletext, #图中的文本
           mean=c(NA,NA,1,data$HR,NA),#HR
           lower=c(NA,NA,1,data$lower,NA), #95%置信区间下限
           upper=c(NA,NA,1,data$upper,NA),#95%置信区间上限
           title="Multivariate Cox regression analysis",
           graph.pos=5,#图在表中的列位置
           graphwidth = unit(.25,"npc"),#图在表中的宽度比例
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),#box颜色
           boxsize=0.4,#box大小固定
           lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           grid = structure(c(data[1,]$Point.Estimate), gp = gpar(col = "black", lty=2,lwd=2)),#虚线(可多条)及其横坐标、颜色、线宽
           #xticks = c(floor(min(data$lower)), 1, 1+(ceiling(max(data$upper))-1)/2, ceiling(max(data$upper))),#横坐标刻度根据需要可随意设置
           lwd.xaxis=2,#X轴线宽
           xlab="",#X轴标题
           hrzl_lines=list("3" = gpar(lwd=2, col="black"),#第三行顶部加黑线，引号内数字标记行位置
                           #"4" = gpar(lwd=60,lineend="butt", columns=c(1:4), col="#99999922"),#加阴影，弱项不建议使用
                           "23" = gpar(lwd=2, col="black")),#最后一行底部加黑线,""中数字为nrow(data)+5
           txt_gp=fpTxtGp(label=gpar(cex=1.2),#各种字体大小设置
                          ticks=gpar(cex=1.2),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.5)),
           #is.summary = c(T,rep(F,27)),#首行字体类型设置
           lineheight = unit(1,"cm"),#固定行高
           #align=c("l","c","c"),#每列文字的对齐方式，偶尔会用到
           #cex=10,
           colgap = unit(1,"cm"),#列间隙
           #mar=unit(rep(1.25, times = 4), "cm"),#图形页边距
           new_page = T)#是否新页
dev.off()


}

# CoxIndependent('./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/train.csv','./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/')