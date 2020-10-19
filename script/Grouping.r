
library(survival)
library(survminer)
library(ggplot2)
options(stringsAsFactors = FALSE) #禁止chr转成factor


#================================ 6个基因集GSVA 的 timeROC==========================
#GSVA分析
library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)

#读取训练集763个样本的RNA表达谱
Grouping <- function(train_data,save_path) {


data.train <- read.csv(train_data,row.names = 1,check.names = F)
# 6个错误率最低的重要性基因作为一个基因集
importance_gene<- read.table(paste0(save_path,"signature with importance.txt"), sep = "\t", stringsAsFactors = F,header = T, check.names=F,fill=TRUE)
# 读入基因集的基因列表
genes = data.frame(Genes="Genes_GSVA",t(unique(importance_gene[,1])))

gs = list()
for (i in 1:nrow(genes)) {
  gs[[genes[i, 1]]] = t(genes[i, 2:ncol(genes)]) %>% .[grep("\\S", .)]##等于空的不要
}##转换列表，可以用(unlist(strsplit(test_1[,2],split="/")))

# 读入基因表达矩阵
exprSet <- t(data.train[,-c(1:12)])
# 这一句就完成了GSVA分析
gsva_es = gsva(as.matrix(exprSet), gs)
# 保存到文件
write.csv(gsva_es, paste0(save_path,"importance_gene_gsva.csv"), quote = F)


#1、GSVA评分
importance_gene_gsva <- read.csv(paste0(save_path,'importance_gene_gsva.csv'),row.names = 1,check.names = F)
GSVA = data.frame(t(importance_gene_gsva),check.names=F)
data <- merge(data.train[,c(1,3)],GSVA,by=0)
row.names(data)=data[,1]
svdata = na.omit(data[,-1])

#对行按照基因表达水平排序，默认从低到高
sortsv<-svdata[order(svdata$Genes_GSVA),]

#输出:中位值分组的生存曲线、最佳分组生存曲线、遍历所有分组情况下的P值和Hazard Ratio的分布情况
#Median separation
#先根据表达水平的中位值分组，画生存曲线，保存
ssdf<-cbind(sortsv,data.frame(gp=ifelse(sortsv$Genes_GSVA>median(sortsv$Genes_GSVA),"high","low")))
ssdf$time = as.numeric(ssdf$A1_OS)/30

fit<-survfit(Surv(time,event)~gp,data=ssdf)
pdf(file=paste0(save_path,"GSVA_medianSeparation_OS.pdf"),width=8, height=7)
sc_median<-ggsurvplot(fit,data = ssdf, linetype = "strata", conf.int = F, pval = TRUE,
                      palette = c("#D95F02","#1B9E77"),
                      fontsize=20,font.legend=20,pval.size=7,
                      font.main = c(21),
                      font.x = c(20),
                      font.y = c(20),
                      font.tickslab = c(18), 
                      legend.title="",legend=c(0.7,0.9),
                      legend.labs=c("GSVA_High","GSVA_low"))+  
  labs(title = "Survival Analysis of Genes_GSVA", x = "Time/Months", y = "Survival probability")
dev.off()
#遍历所有分组情况，计算P值和Hazard Ratio，p值用于判断分组之间差异是否显著，而Hazard Ratio用于衡量分组之间的差异程度
pvals<-c()
hrs<-c()
for(n in 1:(nrow(sortsv)-1)){
ssdf<-cbind(sortsv,data.frame(gp=rep(c("low","high"),c(n,nrow(sortsv)-n))))
ssdf$time = as.numeric(ssdf$A1_OS)/30
diff<-survdiff(Surv(time,event)~gp,data=ssdf,rho = 0)
pv<-pchisq(diff$chisq,length(diff$n)-1,lower.tail=FALSE)
pvals<-c(pvals,pv)
hr<-diff$obs[1]*diff$exp[2]/(diff$obs[2]*diff$exp[1])
hrs<-c(hrs,hr)
}
#展示所有分组情况下的P值和Hazard Ratio的分布情况，水平虚线标记位置的P值为0.05，两条竖直虚线标记的HR为0.5和2
fd<-data.frame(Tag=1:(nrow(sortsv)-1),HR=hrs,Pvalue=pvals)
ggplot(fd,aes(x=log2(HR),y=-log10(Pvalue)))+
   geom_point(shape=21,aes(fill=Tag))+
   scale_fill_gradientn(colours=c("#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444"),guide="legend")+
   geom_hline(yintercept=(-log10(0.05)),linetype="dashed")+
   geom_vline(xintercept=c(-1,1),linetype="dashed")+
   annotate("text",y=-log10(pvals[which.min(pvals)]),x=log2(hrs[which.min(pvals)]),label="min-Pvalue")
#ggsave(file="..\\img\\Pvalue_Hazard-Ratio.pdf")

#Best separation:虚线左上角区域的点p值最小，是最佳的分组方式，分组情况如下
ggplot(fd,aes(x=log2(HR),y=-log10(Pvalue)))+
    geom_point(shape=21,aes(fill=Tag))+
	scale_fill_gradientn(colours=c("#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444"),guide="legend")+
	geom_hline(yintercept=(-log10(0.05)),linetype="dashed")+
	geom_vline(xintercept=c(-1,1),linetype="dashed")+
	geom_text(aes(label=paste(Tag,nrow(sortsv)-Tag,sep=":")),vjust=-1)

#画出最佳分组的生存曲线
ssdf<-cbind(sortsv,data.frame(gp=rep(c("low","high"),c(which.min(pvals),nrow(sortsv)-which.min(pvals)))))
ssdf$time = as.numeric(ssdf$A1_OS)/30
write.table(ssdf,paste0(save_path,"GSVA_BestSeparation.txt"),row.names = T,quote = F,sep="\t")# 生成txt文件

fit<-survfit(Surv(time,event)~gp,data=ssdf)
pdf(file=paste0(save_path,"GSVA_BestSeparation_OS.pdf"),width=8, height=7)
sc_minp<-ggsurvplot(fit,data = ssdf, linetype = "strata", conf.int = F, pval = TRUE,
                    palette = c("#D95F02","#1B9E77"),
                    fontsize=20,font.legend=20,pval.size=7,
                    font.main = c(21),
                    font.x = c(20),
                    font.y = c(20),
                    font.tickslab = c(18), 
                    legend.title="",legend=c(0.7,0.9),
                    legend.labs=c("GSVA_High","GSVA_low"))+  
  labs(title = "Survival Analysis of Genes_GSVA", x = "Time/Months", y = "Survival probability")

dev.off()



#2、COX评分
COX <- read.csv(paste0(save_path,'Cox_riskscore.csv'),row.names = 1,check.names = F,stringsAsFactors=F)
svdata = COX

#对行按照基因表达水平排序，默认从低到高
sortsv<-svdata[order(svdata$riskscore),]

#输出:中位值分组的生存曲线、最佳分组生存曲线、遍历所有分组情况下的P值和Hazard Ratio的分布情况
#Median separation
#先根据表达水平的中位值分组，画生存曲线，保存
ssdf<-cbind(sortsv,data.frame(gp=ifelse(sortsv$riskscore>median(sortsv$riskscore),"high","low")))

ssdf$time1  <-  as.numeric(ssdf$A1_OS)/30

fit<-survfit(Surv(time1,event)~gp,data=ssdf)
pdf(file=paste0(save_path,"COX_medianSeparation_OS.pdf"),width=8, height=7)
sc_median<-ggsurvplot(fit,data = ssdf, linetype = "strata", conf.int = F, pval = TRUE,
                      palette = c("#D95F02","#1B9E77"),
                      fontsize=20,font.legend=20,pval.size=7,
                      font.main = c(21),
                      font.x = c(20),
                      font.y = c(20),
                      font.tickslab = c(18), 
                      legend.title="",legend=c(0.7,0.9),
                      legend.labs=c("High_riskscore","low_riskscore"))+  
  labs(title = "Survival Analysis of COX_riskscore", x = "Time/Months", y = "Survival probability")

dev.off()
#遍历所有分组情况，计算P值和Hazard Ratio，p值用于判断分组之间差异是否显著，而Hazard Ratio用于衡量分组之间的差异程度
pvals<-c()
hrs<-c()
for(n in 1:(nrow(sortsv)-1)){
  ssdf<-cbind(sortsv,data.frame(gp=rep(c("low","high"),c(n,nrow(sortsv)-n))))
  ssdf$time = as.numeric(ssdf$A1_OS)/30
  diff<-survdiff(Surv(time,event)~gp,data=ssdf,rho = 0)
  pv<-pchisq(diff$chisq,length(diff$n)-1,lower.tail=FALSE)
  pvals<-c(pvals,pv)
  hr<-diff$obs[1]*diff$exp[2]/(diff$obs[2]*diff$exp[1])
  hrs<-c(hrs,hr)
}
#展示所有分组情况下的P值和Hazard Ratio的分布情况，水平虚线标记位置的P值为0.05，两条竖直虚线标记的HR为0.5和2
fd<-data.frame(Tag=1:(nrow(sortsv)-1),HR=hrs,Pvalue=pvals)
ggplot(fd,aes(x=log2(HR),y=-log10(Pvalue)))+
  geom_point(shape=21,aes(fill=Tag))+
  scale_fill_gradientn(colours=c("#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444"),guide="legend")+
  geom_hline(yintercept=(-log10(0.05)),linetype="dashed")+
  geom_vline(xintercept=c(-1,1),linetype="dashed")+
  annotate("text",y=-log10(pvals[which.min(pvals)]),x=log2(hrs[which.min(pvals)]),label="min-Pvalue")
#ggsave(file="..\\img\\Pvalue_Hazard-Ratio.pdf")

#Best separation:虚线左上角区域的点p值最小，是最佳的分组方式，分组情况如下
ggplot(fd,aes(x=log2(HR),y=-log10(Pvalue)))+
  geom_point(shape=21,aes(fill=Tag))+
  scale_fill_gradientn(colours=c("#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444"),guide="legend")+
  geom_hline(yintercept=(-log10(0.05)),linetype="dashed")+
  geom_vline(xintercept=c(-1,1),linetype="dashed")+
  geom_text(aes(label=paste(Tag,nrow(sortsv)-Tag,sep=":")),vjust=-1)

#画出最佳分组的生存曲线
ssdf<-cbind(sortsv,data.frame(gp=rep(c("low","high"),c(which.min(pvals),nrow(sortsv)-which.min(pvals)))))
ssdf$time = as.numeric(ssdf$A1_OS)/30
write.table(ssdf,paste0(save_path,"COX_BestSeparation.txt"),row.names = T,quote = F,sep="\t")# 生成txt文件

fit<-survfit(Surv(time,event)~gp,data=ssdf)
pdf(file=paste0(save_path,"COX_BestSeparation_OS.pdf"),width=8, height=7)
sc_minp<-ggsurvplot(fit,data = ssdf, linetype = "strata", conf.int = F, pval = TRUE,
                    palette = c("#D95F02","#1B9E77"),
                    fontsize=20,font.legend=20,pval.size=7,
                    font.main = c(21),
                    font.x = c(20),
                    font.y = c(20),
                    font.tickslab = c(18), 
                    legend.title="",legend=c(0.7,0.9),
                    legend.labs=c("High_riskscore","low_riskscore"))+  
  labs(title = "Survival Analysis of COX_riskscore", x = "Time/Months", y = "Survival probability")
sc_minp
dev.off()

#3、Lasso评分
Lasso <- read.table(paste0(save_path,'lasso_Score.txt'),sep = "\t", header = T, check.names=F,stringsAsFactors=F)

data <- merge(Lasso,data.train[,c(1,3)],by=0)

row.names(data)=data[,1]
svdata = na.omit(data[,-1])

#对行按照基因表达水平排序，默认从低到高
sortsv<-svdata[order(svdata$lasso_Score),]

#输出:中位值分组的生存曲线、最佳分组生存曲线、遍历所有分组情况下的P值和Hazard Ratio的分布情况
#Median separation
#先根据表达水平的中位值分组，画生存曲线，保存
ssdf<-cbind(sortsv,data.frame(gp=ifelse(sortsv$lasso_Score>median(sortsv$lasso_Score),"high","low")))
ssdf$time = as.numeric(ssdf$A1_OS)/30

fit<-survfit(Surv(time,event)~gp,data=ssdf)
pdf(file=paste0(save_path,"Lasso_medianSeparation_OS.pdf"),width=8, height=7)
sc_median<-ggsurvplot(fit,data = ssdf, linetype = "strata", conf.int = F, pval = TRUE,
                      palette = c("#D95F02","#1B9E77"),
                      fontsize=20,font.legend=20,pval.size=7,
                      font.main = c(21),
                      font.x = c(20),
                      font.y = c(20),
                      font.tickslab = c(18), 
                      legend.title="",legend=c(0.7,0.9),
                      legend.labs=c("High_riskscore","low_riskscore"))+  
  labs(title = "Survival Analysis of Lasso_riskscore", x = "Time/Months", y = "Survival probability")
dev.off()
#遍历所有分组情况，计算P值和Hazard Ratio，p值用于判断分组之间差异是否显著，而Hazard Ratio用于衡量分组之间的差异程度
pvals<-c()
hrs<-c()
for(n in 1:(nrow(sortsv)-1)){
  ssdf<-cbind(sortsv,data.frame(gp=rep(c("low","high"),c(n,nrow(sortsv)-n))))
  ssdf$time = as.numeric(ssdf$A1_OS)/30
  diff<-survdiff(Surv(time,event)~gp,data=ssdf,rho = 0)
  pv<-pchisq(diff$chisq,length(diff$n)-1,lower.tail=FALSE)
  pvals<-c(pvals,pv)
  hr<-diff$obs[1]*diff$exp[2]/(diff$obs[2]*diff$exp[1])
  hrs<-c(hrs,hr)
}
#展示所有分组情况下的P值和Hazard Ratio的分布情况，水平虚线标记位置的P值为0.05，两条竖直虚线标记的HR为0.5和2
fd<-data.frame(Tag=1:(nrow(sortsv)-1),HR=hrs,Pvalue=pvals)
ggplot(fd,aes(x=log2(HR),y=-log10(Pvalue)))+
  geom_point(shape=21,aes(fill=Tag))+
  scale_fill_gradientn(colours=c("#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444"),guide="legend")+
  geom_hline(yintercept=(-log10(0.05)),linetype="dashed")+
  geom_vline(xintercept=c(-1,1),linetype="dashed")+
  annotate("text",y=-log10(pvals[which.min(pvals)]),x=log2(hrs[which.min(pvals)]),label="min-Pvalue")
#ggsave(file="..\\img\\Pvalue_Hazard-Ratio.pdf")

#Best separation:虚线左上角区域的点p值最小，是最佳的分组方式，分组情况如下
ggplot(fd,aes(x=log2(HR),y=-log10(Pvalue)))+
  geom_point(shape=21,aes(fill=Tag))+
  scale_fill_gradientn(colours=c("#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444"),guide="legend")+
  geom_hline(yintercept=(-log10(0.05)),linetype="dashed")+
  geom_vline(xintercept=c(-1,1),linetype="dashed")+
  geom_text(aes(label=paste(Tag,nrow(sortsv)-Tag,sep=":")),vjust=-1)

#画出最佳分组的生存曲线
ssdf<-cbind(sortsv,data.frame(gp=rep(c("low","high"),c(which.min(pvals),nrow(sortsv)-which.min(pvals)))))
ssdf$time = as.numeric(ssdf$A1_OS)/30
write.table(ssdf,paste0(save_path,"Lasso_BestSeparation.txt"),row.names = T,quote = F,sep="\t")# 生成txt文件

fit<-survfit(Surv(time,event)~gp,data=ssdf)
pdf(file=paste0(save_path,"Lasso_BestSeparation_OS.pdf"),width=8, height=7)
sc_minp<-ggsurvplot(fit,data = ssdf, linetype = "strata", conf.int = F, pval = TRUE,
                    palette = c("#D95F02","#1B9E77"),
                    fontsize=20,font.legend=20,pval.size=7,
                    font.main = c(21),
                    font.x = c(20),
                    font.y = c(20),
                    font.tickslab = c(18), 
                    legend.title="",legend=c(0.7,0.9),
                    legend.labs=c("High_riskscore","low_riskscore"))+  
  labs(title = "Survival Analysis of Lasso_riskscore", x = "Time/Months", y = "Survival probability")
sc_minp
dev.off()





#4、randomForestSRC评分
randomForestSRC <- read.csv(paste0(save_path,'randomForestSRC_riskscore.csv'),row.names = 1,check.names = F,stringsAsFactors=F)
svdata = na.omit(randomForestSRC)

#对行按照基因表达水平排序，默认从低到高
sortsv<-svdata[order(svdata$SRC_Score),]

#输出:中位值分组的生存曲线、最佳分组生存曲线、遍历所有分组情况下的P值和Hazard Ratio的分布情况
#Median separation
#先根据表达水平的中位值分组，画生存曲线，保存
ssdf<-cbind(sortsv,data.frame(gp=ifelse(sortsv$SRC_Score>median(sortsv$SRC_Score),"high","low")))
ssdf$time = as.numeric(ssdf$A1_OS)/30

fit<-survfit(Surv(time,event)~gp,data=ssdf)
pdf(file=paste0(save_path,"randomForestSRC_medianSeparation_OS.pdf"),width=8, height=7)
sc_median<-ggsurvplot(fit,data = ssdf, linetype = "strata", conf.int = F, pval = TRUE,
                      palette = c("#D95F02","#1B9E77"),
                      fontsize=20,font.legend=20,pval.size=7,
                      font.main = c(21),
                      font.x = c(20),
                      font.y = c(20),
                      font.tickslab = c(18), 
                      legend.title="",legend=c(0.7,0.9),
                      legend.labs=c("High_riskscore","low_riskscore"))+  
  labs(title = "Survival Analysis of randomForestSRC_riskscore", x = "Time/Months", y = "Survival probability")
sc_median
dev.off()
#遍历所有分组情况，计算P值和Hazard Ratio，p值用于判断分组之间差异是否显著，而Hazard Ratio用于衡量分组之间的差异程度
pvals<-c()
hrs<-c()
for(n in 1:(nrow(sortsv)-1)){
  ssdf<-cbind(sortsv,data.frame(gp=rep(c("low","high"),c(n,nrow(sortsv)-n))))
  ssdf$time = as.numeric(ssdf$A1_OS)/30
  diff<-survdiff(Surv(time,event)~gp,data=ssdf,rho = 0)
  pv<-pchisq(diff$chisq,length(diff$n)-1,lower.tail=FALSE)
  pvals<-c(pvals,pv)
  hr<-diff$obs[1]*diff$exp[2]/(diff$obs[2]*diff$exp[1])
  hrs<-c(hrs,hr)
}
#展示所有分组情况下的P值和Hazard Ratio的分布情况，水平虚线标记位置的P值为0.05，两条竖直虚线标记的HR为0.5和2
fd<-data.frame(Tag=1:(nrow(sortsv)-1),HR=hrs,Pvalue=pvals)
ggplot(fd,aes(x=log2(HR),y=-log10(Pvalue)))+
  geom_point(shape=21,aes(fill=Tag))+
  scale_fill_gradientn(colours=c("#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444"),guide="legend")+
  geom_hline(yintercept=(-log10(0.05)),linetype="dashed")+
  geom_vline(xintercept=c(-1,1),linetype="dashed")+
  annotate("text",y=-log10(pvals[which.min(pvals)]),x=log2(hrs[which.min(pvals)]),label="min-Pvalue")
#ggsave(file="..\\img\\Pvalue_Hazard-Ratio.pdf")

#Best separation:虚线左上角区域的点p值最小，是最佳的分组方式，分组情况如下
ggplot(fd,aes(x=log2(HR),y=-log10(Pvalue)))+
  geom_point(shape=21,aes(fill=Tag))+
  scale_fill_gradientn(colours=c("#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444"),guide="legend")+
  geom_hline(yintercept=(-log10(0.05)),linetype="dashed")+
  geom_vline(xintercept=c(-1,1),linetype="dashed")+
  geom_text(aes(label=paste(Tag,nrow(sortsv)-Tag,sep=":")),vjust=-1)

#画出最佳分组的生存曲线
ssdf<-cbind(sortsv,data.frame(gp=rep(c("low","high"),c(which.min(pvals),nrow(sortsv)-which.min(pvals)))))
ssdf$time = as.numeric(ssdf$A1_OS)/30
write.table(ssdf,paste0(save_path,"randomForestSRC_BestSeparation.txt"),row.names = T,quote = F,sep="\t")# 生成txt文件

fit<-survfit(Surv(time,event)~gp,data=ssdf)
pdf(file=paste0(save_path,"randomForestSRC_BestSeparation_OS.pdf"),width=8, height=7)
sc_minp<-ggsurvplot(fit,data = ssdf, linetype = "strata", conf.int = F, pval = TRUE,
                    palette = c("#D95F02","#1B9E77"),
                    fontsize=20,font.legend=20,pval.size=7,
                    font.main = c(21),
                    font.x = c(20),
                    font.y = c(20),
                    font.tickslab = c(18), 
                    legend.title="",legend=c(0.7,0.9),
                    legend.labs=c("High_riskscore","low_riskscore"))+  
  labs(title = "Survival Analysis of randomForestSRC_riskscore", x = "Time/Months", y = "Survival probability")
sc_minp
dev.off()

}
# Grouping('./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/train.csv','./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/')