Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
library(survival)
library(survminer)
library(randomForestSRC)
library(timeROC)
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(caret)

timeR <- function(train_data,importance_gsva,cox_risk,lasso_risk,rfsrc_risk,save_path) {

gsva_es <- read.csv(importance_gsva,row.names = 1,check.names = F)
data.train <- read.csv(train_data,row.names = 1,check.names = F)

## timeROC

GSVA = data.frame(t(gsva_es),check.names=F)
score_t <- merge(data.train[,c(1,3)],GSVA,by=0)
row.names(score_t)=score_t[,1]
score_t = score_t[,-1]
colnames(score_t)[3] = "GSVA_Score"
score_t$time = as.numeric(score_t$A1_OS)/365

col <- c("#0073C2FF","#E43889","#10A754","orange") ## 自定义颜色
tt <- timeROC(score_t$time,score_t$event,score_t$GSVA_Score,
              cause = 1,weighting = 'marginal',
              times = c(1,5,10,15,20),ROC = T,iid = T)
tt$AUC

pdf(file=paste0(save_path,"GSVA_timeROC.pdf"),width=6, height=6)
plot(tt,time=1,title=FALSE,lwd=2,col=col[1],cex.lab=2, cex.axis=2)
plot(tt,time=5,col=col[2],add=TRUE,title=FALSE,lwd=2,cex.lab=2, cex.axis=2)
plot(tt,time=10,col=col[3],add=TRUE,title=FALSE,lwd=2,cex.lab=2, cex.axis=2)
plot(tt,time=15,col=col[4],add=TRUE,title=FALSE,lwd=2,cex.lab=2, cex.axis=2)

id <- c(paste0("GSVA_Score  1-Year AUC = ",round(tt$AUC[1],3)),
        paste0("GSVA_Score  5-Year AUC = ",round(tt$AUC[2],3)),
        paste0("GSVA_Score  10-Year AUC = ",round(tt$AUC[3],3)),
        paste0("GSVA_Score  15-Year AUC = ",round(tt$AUC[4],3))
)
legend("bottomright",id,
       fill=col[1:4],
       bty="o",cex=1.2,
       border = NA)
abline(0,1,lty=2,lwd=0.5)
dev.off()






#================================ 多因素COX 的 timeROC==========================
## 读取COX 的 riskscore文件
score_t <- read.csv(cox_risk,row.names = 1,check.names = F)
score_t$A1_OS = as.numeric(score_t$A1_OS)/365
## timeROC
library(timeROC)
col <- c("#0073C2FF","#E43889","#10A754","orange") ## 自定义颜色
tt <- timeROC(score_t$A1_OS,score_t$event,score_t$riskscore,
              cause = 1,weighting = 'marginal',
              times = c(1,5,10,15,20),ROC = T,iid = T)
tt$AUC

pdf(file=paste0(save_path,"COX_timeROC.pdf"),width=6, height=6)
plot(tt,time=1,title=FALSE,lwd=2,col=col[1],cex.lab=2, cex.axis=2)
plot(tt,time=5,col=col[2],add=TRUE,title=FALSE,lwd=2,cex.lab=2, cex.axis=2)
plot(tt,time=10,col=col[3],add=TRUE,title=FALSE,lwd=2,cex.lab=2, cex.axis=2)
plot(tt,time=15,col=col[4],add=TRUE,title=FALSE,lwd=2,cex.lab=2, cex.axis=2)

id <- c(paste0("COX_riskScore  1-Year AUC = ",round(tt$AUC[1],3)),
        paste0("COX_riskScore  5-Year AUC = ",round(tt$AUC[2],3)),
        paste0("COX_riskScore  10-Year AUC = ",round(tt$AUC[3],3)),
        paste0("COX_riskScore  15-Year AUC = ",round(tt$AUC[4],3))
)
legend("bottomright",id,
       fill=col[1:4],
       bty="o",cex=1.2,
       border = NA)
abline(0,1,lty=2,lwd=0.5)
dev.off()







#================================ Lasso 的 timeROC==========================
## 读取Lasso 的 Score文件
lasso_Score <- read.table(lasso_risk,row.names = 1,check.names = F)
score_t = merge(data.train[,c(1,3)],lasso_Score,by=0)
row.names(score_t) = score_t[,1]
score_t = score_t[,-1]
score_t$A1_OS = as.numeric(score_t$A1_OS)/365
## timeROC
library(timeROC)
col <- c("#0073C2FF","#E43889","#10A754","orange") ## 自定义颜色
tt <- timeROC(score_t$A1_OS,score_t$event,score_t$lasso_Score,
              cause = 1,weighting = 'marginal',
              times = c(1,5,10,15,20),ROC = T,iid = T)
tt$AUC

pdf(file=paste0(save_path,"Lasso_timeROC.pdf"),width=6, height=6)
plot(tt,time=1,title=FALSE,lwd=2,col=col[1],cex.lab=2, cex.axis=2)
plot(tt,time=5,col=col[2],add=TRUE,title=FALSE,lwd=2,cex.lab=2, cex.axis=2)
plot(tt,time=10,col=col[3],add=TRUE,title=FALSE,lwd=2,cex.lab=2, cex.axis=2)
plot(tt,time=15,col=col[4],add=TRUE,title=FALSE,lwd=2,cex.lab=2, cex.axis=2)

id <- c(paste0("Lasso_Score  1-Year AUC = ",round(tt$AUC[1],3)),
        paste0("Lasso_Score   5-Year AUC = ",round(tt$AUC[2],3)),
        paste0("Lasso_Score   10-Year AUC = ",round(tt$AUC[3],3)),
        paste0("Lasso_Score   15-Year AUC = ",round(tt$AUC[4],3))
)
legend("bottomright",id,
       fill=col[1:4],
       bty="o",cex=1.2,
       border = NA)
abline(0,1,lty=2,lwd=0.5)
dev.off()





#================================ randomForestSRC 的 timeROC==========================
## 读取randomForestSRC 的 riskscore文件
score_t <- read.csv(rfsrc_risk,row.names = 1,check.names = F)
## timeROC
library(timeROC)
col <- c("#0073C2FF","#E43889","#10A754","orange") ## 自定义颜色
tt <- timeROC(score_t$time,score_t$event,score_t$SRC_Score,
              cause = 1,weighting = 'marginal',
              times = c(1,5,10,15,20),ROC = T,iid = T)
tt$AUC

pdf(file=paste0(save_path,"randomForestSRC_timeROC.pdf"),width=6, height=6)
plot(tt,time=1,title=FALSE,lwd=2,col=col[1],cex.lab=2, cex.axis=2)
plot(tt,time=5,col=col[2],add=TRUE,title=FALSE,lwd=2,cex.lab=2, cex.axis=2)
plot(tt,time=10,col=col[3],add=TRUE,title=FALSE,lwd=2,cex.lab=2, cex.axis=2)
plot(tt,time=15,col=col[4],add=TRUE,title=FALSE,lwd=2,cex.lab=2, cex.axis=2)

id <- c(paste0("randomForestSRC_Score  1-Year AUC = ",round(tt$AUC[1],3)),
        paste0("randomForestSRC_Score   5-Year AUC = ",round(tt$AUC[2],3)),
        paste0("randomForestSRC_Score   10-Year AUC = ",round(tt$AUC[3],3)),
        paste0("randomForestSRC_Score   15-Year AUC = ",round(tt$AUC[4],3))
)
legend("bottomright",id,
       fill=col[1:4],
       bty="o",cex=1.2,
       border = NA)
abline(0,1,lty=2,lwd=0.5)
dev.off()


}

timeR(
        './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/train.csv',
        './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/importance_gene_gsva.csv',
        './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Cox_riskscore.csv',
        './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/lasso_Score.txt',
        './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/randomForestSRC_riskscore.csv',
        './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/'
        )