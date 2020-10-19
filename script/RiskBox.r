#===============================================================
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

#读取4组数据
RiskBox <- function(GSVAdata,coxdata,lassodata,rfsrcdata,save_path) {
GSVA <- read.table  (GSVAdata ,row.names = 1)
COX <- read.table   (coxdata  ,row.names = 1)
Lasso <- read.table (lassodata,row.names = 1)
randomForestSRC <- read.table(rfsrcdata,row.names = 1)

GSVA = data.frame(Sample=row.names(GSVA),GSVA[,3:4])
GSVA$type = "GSVA Score"

COX = data.frame(Sample=row.names(COX),COX[,3:4])
COX$type = "COX riskScore"

Lasso = data.frame(Sample=row.names(Lasso),Lasso[,c(1,4)])
Lasso$type = "Lasso riskScore"

randomForestSRC = data.frame(Sample=row.names(randomForestSRC),randomForestSRC[,c(3,5)])
randomForestSRC$type = "randomForestSRC riskScore"

colnames(GSVA)<-c("Sample","value","Group","type")
colnames(COX)<-colnames(GSVA)
colnames(Lasso)<-colnames(GSVA)
colnames(randomForestSRC)<-colnames(GSVA)

GSVA$Group<-factor(GSVA$Group,levels = c("low","high"))
COX$Group<-factor(COX$Group,levels = c("low","high"))
Lasso$Group<-factor(Lasso$Group,levels = c("low","high"))
randomForestSRC$Group<-factor(randomForestSRC$Group,levels = c("low","high"))

pdata<-rbind(GSVA,COX,Lasso,randomForestSRC)
pdata$Group<-factor(pdata$Group,levels = c("low","high"))


#多个箱线图拼成一个大图
require(tidyr)
require(ggplot2)
require(cowplot)
box1 <- lapply(unique(pdata$type), function(i) {
  dd <- pdata[pdata$type==i,]
  ## 用one way anova计算 p value
  res <- aov(value ~ Group, data = dd)
  pv1 <- summary(res)[[1]]$'Pr(>F)'[1]
  pv1.lab <- paste("one-way ANOVA p =", round(pv1,3))
  ## t test for comparing stage 3 vs stage 1
  pv2 <- t.test(value ~ Group, data=dd)$p.value
  pv2.lab <- paste("high vs low p =", round(pv2, 3))
  
  lab <- paste(unique(dd$type), pv1.lab, pv2.lab, sep="\n")
  label <- ggdraw() + draw_label(lab)
  plot_grid(label, 
            ggplot(dd, aes(Group, value, fill=Group)) +
              geom_boxplot() +
              #geom_jitter(size=0.1) + 
              scale_fill_manual(values = c("#0073C2","#EFBF00"))+ #设置填充的颜色
              theme_bw()+ #背景变为白色
              xlab("") + ylab(unique(dd$type)) +
              theme(legend.title = element_text(color="black",size=15),
                    legend.text= element_text(color="black",size=15),
                    axis.text.x = element_text(color="black", size=15),
                    axis.text.y = element_text(color="black", size=15),
                    axis.title.y = element_text(size = 15),
                    panel.grid = element_blank()),   #不显示网格线, 
            ncol=1, rel_heights=c(.3, 1))
})
plot_grid(plotlist=box1, ncol=4)
ggsave(file=paste0(save_path,"BestSeparation_box.pdf"),  width =15, height = 4)


}

# RiskBox("./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/GSVA_BestSeparation.txt","./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/COX_BestSeparation.txt" ,"./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Lasso_BestSeparation.txt","./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/randomForestSRC_BestSeparation.txt",'./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/')