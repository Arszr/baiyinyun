library(survival)
library(glmnet)
library(pbapply)
library(survivalROC)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
# 显示进程
display.progress = function (index, totalN, breakN=20) {
  
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}    
# lasso回归
surv_lasso <- function(iter.times = NULL, surv.obj = NULL, expr.obj = NULL, nfolds = 10, alpha = 1, family = "cox") {
  # iter.times: pblapply的传入参数，用于迭代次数
  # surv.obj: surv对象，由Surv()函数得到；
  # expr.obj: 表达谱对象，注意行为特征，列为样本
  # nfolds：筛选最优lambda时的交叉验证次数，默认为10
  # alpha： 默认为1表示LASSO回归
  # family： 默认为"cox"
  cvfit = cv.glmnet(x = t(as.matrix(expr.obj)), 
                    y = surv.obj, 
                    nfolds = nfolds, # 10-fold交叉验证选取最优lambda
                    alpha = alpha, # alpha = 1 意味着 lasso
                    family = family) # 依赖cox模型
  # 取出最优lambda
  myCoefs <- coef(cvfit, s="lambda.min");
  lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )] # 取出非0的特征
  
  return(lasso_fea)
}

#输入文件:表达矩阵和带followup的临床数据
# 筛选既有表达矩阵又有followup的样本
LineLasso <- function(train_data,coxsingle_train,save_path) {

    trait <- na.omit(read.csv(train_data,check.names = F,stringsAsFactors = F,header = T,row.names = 1)[,-c(4:12)])
    trait = trait[trait$A1_OS!=0,]#生存时间有0值会报错，去掉0的生存时间就好使了 

    expr <- data.frame(t(trait[,-c(1:3)]),check.names = F)
    Sinfo <- trait[,c(1:3)]
    # 过滤出方差>var.cutoff的基因
    # 参数设置
    var.cutoff <- 0.2 # 例文为GEO数据，阈值设为0.2
    var <- apply(expr, 1, sd)
    expr.filtered <- expr[var > var.cutoff,] # 因为这里是TCGA数据，所以我设置的大一些
    #保存到文件
    write.csv(round(expr.filtered, 1), paste(save_path,"easy_input_expr.csv",sep=''), row.names = T,quote = F)
    write.csv(Sinfo, paste0(save_path,"easy_input_cli.csv"), row.names = T)

    #筛选基因
    expr <- read.csv(paste0(save_path,"easy_input_expr.csv"), check.names = F,stringsAsFactors = F,header = T,row.names = 1)
    Sinfo <- read.csv(paste0(save_path,"easy_input_cli.csv"), check.names = F,stringsAsFactors = F,header = T,row.names = 1)

    # 过滤出18个OS预后良好的基因
    Coxoutput_train<- read.csv(coxsingle_train, stringsAsFactors=FALSE)
    surv.expr <- expr[Coxoutput_train$factors,]#18

    ### 迭代LASSO挑选高频特征 ###
    iter.times <- 500 # 设置迭代次数，速度非常慢请耐心，例文是1000次
    # 运算lasso回归
    lasso_fea <- list()
    surv <- Surv(Sinfo$A1_OS, Sinfo$event)

    set.seed(10) # 外部设置种子，使得迭代过程是可重复的
    #下面这步运行时间较长，我们把它保存到lasso_fea.rda里
    lasso_fea <- pblapply(1:iter.times,
                        surv_lasso, 
                        surv.obj = surv, 
                        expr.obj = surv.expr)
    save(lasso_fea,file = paste0(save_path,"lasso_fea.rda")) # 保存该结果
    #这里直接加载上一步运行的结果
    load(paste0(save_path,"lasso_fea.rda")) # 加载该结果

    #######################
    ### 根据AUC挑选基因 ###
    genes <- sort(table(unlist(lasso_fea)), decreasing = T) # 根据基因出现的频次排序
    # 如果觉得出现次数较少的基因是不鲁棒的，也可以仅选择top基因
    freq.cutoff <- 50 
    genes <- names(genes[genes > freq.cutoff]) # 这里选择出现频次大于50的基因，认为是多次lasso的共识基因

    pred.time <- 5 # 查看5年ROC
    roc <- list() # 初始化roc列表
    auc <- c() # 初始化auc向量
    result=data.frame()
    for (i in 1:length(genes)) {
    gene <- genes[i]
    tmp <- data.frame(gene = as.numeric(surv.expr[gene,]),row.names = colnames(surv.expr),stringsAsFactors = F); colnames(tmp) = gene
    if(i == 1) { # 如果为第一个基因就把生存信息纳入数据框
        surv.dat <- cbind.data.frame(Sinfo[rownames(tmp),c("A1_OS","event")],tmp)
    } else {
        surv.dat <- cbind.data.frame(surv.dat,tmp)
    }
    cox <- coxph(Surv(A1_OS, event) ~ ., data = surv.dat) # 多变量cox比例风险模型（i=1时为单变量）
    riskScore <- predict(cox,type="risk",newdata=surv.dat) # 计算风险
    roc[[i]] <- survivalROC(Stime=surv.dat$A1_OS, 
                            status=surv.dat$event, 
                            marker = riskScore[rownames(surv.dat)], 
                            predict.time =pred.time*365, # 计算pred.time时刻的ROC，一般是五年生存
                            method="KM")
    auc <- c(auc,roc[[i]]$AUC) # 保存auc
    a=data.frame(gene,auc)
    result=rbind(result,a)
    }
    write.table(result,paste0(save_path,"auc.txt"), row.names = F, quote = F)

    #############################
    ### 最终signature与KM曲线 ###
    prog.sig <- genes[1:which.max(auc)]
    #保存到文件
    write.table(prog.sig,paste0(save_path,"signature_gene.txt"), row.names = F, quote = F)

    surv.dat <- t(surv.expr[prog.sig,])
    surv.dat <- cbind.data.frame(Sinfo[rownames(surv.dat),c("A1_OS","event"),],surv.dat)
    surv.dat$A1_OS <- surv.dat$A1_OS/365
    cox <- coxph(Surv(A1_OS, event) ~ ., data = surv.dat) # 多变量cox比例风险模型
    surv.dat$riskScore <- predict(cox,type="risk",newdata=surv.dat) # 计算风险

    surv.dat$Risk <- ifelse(riskScore > median(riskScore),"High","Low") # 根据风险值划分高低风险组
    fitd <- survdiff(Surv(A1_OS, event) ~ Risk, data=surv.dat, na.action=na.exclude)
    p.val <- 1-pchisq(fitd$chisq, length(fitd$n)-1) # 计算KM曲线p值
    fit <- survfit(Surv(A1_OS, event)~ Risk, data=surv.dat, type="kaplan-meier", error="greenwood", conf.type="plain", na.action=na.exclude)


    #开始画图:用base plot画图
    ### 设置颜色 ###
    red  <- "#E94244"
    blue <- "#4084BB"

    pdf(paste0(save_path,"LASSO_auc.pdf"),width = 7,height = 6.5)
    par(mfrow = c(1,1))
    # 通过不断纳入基因（根据频次排序的），计算AUC的变化情况
    par(bty="o", mgp = c(2,0.5,0), mar = c(4.1,4.1,2.1,2.1),tcl=-.25, font.main=3)
    plot(1:length(genes), auc, # 画AUC随基因入组的变化情况
        type = "l", lwd = 2, col = blue,cex.lab=1.5, cex.axis=1.2,
        ylim = c(0.6,1), xlab = "Genes ordered by frequency",ylab = "Area under the curve")
    points(which.max(auc), auc[which.max(auc)], pch = 16, cex = 2, col = red) # 点出峰值
    arrows(x0 = which.max(auc), y0 = auc[which.max(auc)] - 0.05,
        x1 = which.max(auc), y1 = auc[which.max(auc)] - 0.01,
        length = 0.1)
    text(which.max(auc), auc[which.max(auc)] - 0.08,
        labels = paste0("Number of genes: ",which.max(auc),"/n",
                        "Current AUC: ",round(auc[which.max(auc)],3)),cex = 1.3, col = red,adj = 0.5)
    invisible(dev.off())

}

# LineLasso <- LineLasso('./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/train.csv','./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/CoxSingle_train.csv','./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/')


