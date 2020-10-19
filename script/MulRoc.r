Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

library("pROC")  
options(digits=3)

#读取4组数据
MulRoc <- function(gsvad,coxd,lassod,rfd,save_path) {
  GSVA <- read.table(gsvad,row.names = 1)
  COX <- read.table(coxd,row.names = 1)
  Lasso <- read.table(lassod,row.names = 1)
  randomForestSRC <- read.table(rfd,row.names = 1)

  GSVA$Group = ifelse(GSVA$gp=="high",1,0)
  COX$Group = ifelse(COX$gp=="high",1,0)
  Lasso$Group = ifelse(Lasso$gp=="high",1,0)
  randomForestSRC$Group = ifelse(randomForestSRC$gp=="high",1,0)

  GSVA$median_Group = ifelse(GSVA$Genes_GSVA > median(GSVA$Genes_GSVA),1,0)
  COX$median_Group = ifelse(COX$riskscore > median(COX$riskscore),1,0)
  Lasso$median_Group = ifelse(Lasso$lasso_Score > median(Lasso$lasso_Score),1,0)
  randomForestSRC$median_Group = ifelse(randomForestSRC$SRC_Score > median(randomForestSRC$SRC_Score),1,0)


  #输出最优分组下4组数据的多重ROC曲线
  #mycol <- c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
  pdf(file=paste0(save_path,"BestSeparation_ROC.pdf"),width=6,height=6)
  auc.out <- c()  
  x <- plot.roc(GSVA$Group, GSVA$Genes_GSVA,ylim=c(0,1),xlim=c(1,0),
                  #smooth=T,
                  ci=TRUE,col="blue",
                  main= "BestSeparation_ROC",
                  lwd=2, legacy.axes=T) 
    ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
    ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
    auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out1 <- rbind(auc.out,auc.ci)
    
  x <- plot.roc(COX$Group, COX$riskscore,add=T,
                #smooth=T,
                ci=TRUE,col="magenta",lwd=2, legacy.axes=T)
    ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
    ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
    auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out2 <- rbind(auc.out,auc.ci)
    
  x <- plot.roc(Lasso$Group, Lasso$lasso_Score,add=T,
                #smooth=T,
                ci=TRUE,col="green",lwd=2, legacy.axes=T)
    ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
    ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
    auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out3 <- rbind(auc.out,auc.ci)
    
  x <- plot.roc(randomForestSRC$Group, randomForestSRC$SRC_Score,add=T,
                #smooth=T,
                ci=TRUE,col="orange2",lwd=2, legacy.axes=T)
    ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
    ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
    auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
    auc.out4 <- rbind(auc.out,auc.ci)
    
  legend("bottomright", legend=c(paste("GSVA_BestSeparation, AUC=", auc.out1[,1], sep=""), 
                                  paste("COX_BestSeparation, AUC=", auc.out2[,1],  sep=""),
                                  paste("Lasso_BestSeparation, AUC=", auc.out3[,1],  sep=""),
                                  paste("randomForestSRC_BestSeparation, AUC=", auc.out4[,1],  sep="")), 
          col=c( "blue","magenta","green","orange2"), lwd=2,cex=1,bty="n")
  dev.off()



  #输出中位数分组下4组数据的多重ROC曲线
  #mycol <- c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
  pdf(file=paste0(save_path,"medianSeparation_ROC.pdf"),width=6,height=6)
  auc.out <- c()  
  x <- plot.roc(GSVA$median_Group, GSVA$Genes_GSVA,ylim=c(0,1),xlim=c(1,0),
                smooth=F,
                ci=TRUE,col="blue",
                main= "medianSeparation_ROC",
                lwd=2, legacy.axes=T) 
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out1 <- rbind(auc.out,auc.ci)

  x <- plot.roc(COX$median_Group, COX$riskscore,add=T,
                smooth=F,
                ci=TRUE,col="magenta",lwd=2, legacy.axes=T)
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out2 <- rbind(auc.out,auc.ci)

  x <- plot.roc(Lasso$median_Group, Lasso$lasso_Score,add=T,
                smooth=F,
                ci=TRUE,col="green",lwd=2, legacy.axes=T)
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out3 <- rbind(auc.out,auc.ci)

  x <- plot.roc(randomForestSRC$median_Group, randomForestSRC$SRC_Score,add=T,
                smooth=F,
                ci=TRUE,col="orange2",lwd=2, legacy.axes=T)
  ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
  ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限
  auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out4 <- rbind(auc.out,auc.ci)

  legend("bottomright", legend=c(paste("GSVA_medianSeparation, AUC=", auc.out1[,1], sep=""), 
                                paste("COX_medianSeparation, AUC=", auc.out2[,1],  sep=""),
                                paste("Lasso_medianSeparation, AUC=", auc.out3[,1],  sep=""),
                                paste("randomForestSRC_medianSeparation, AUC=", auc.out4[,1],  sep="")), 
        col=c( "blue","magenta","green","orange2"), lwd=2,cex=1,bty="n")
  dev.off()






  ### 设置颜色 ###
  roc.GSVA <- plot.roc(GSVA$Group, GSVA$Genes_GSVA,ylim=c(0,1),xlim=c(1,0),
          #smooth=T,
          ci=TRUE,col="blue",
          main= "BestSeparation_ROC",
          lwd=2, legacy.axes=T) 
  roc.COX <- plot.roc(COX$Group, COX$riskscore, ylim=c(0,1),xlim=c(1,0),
                      smooth=F, #绘制平滑曲线
                      ci=TRUE, 
                      col="magenta",lwd=2, legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1

  legend.paste <- c(paste0("GSVA_BestSeparation AUC: ",round(roc.GSVA$auc,3), " (",paste0(round(roc.GSVA$ci[1],3),"-",round(roc.GSVA$ci[3],3)),")"),
                    paste0("COX_BestSeparation AUC: ",round(roc.COX$auc,3)," (",paste0(round(roc.COX$ci[1],3),"-",round(roc.COX$ci[3],3)),")"))


  pdf("multipanelROC.pdf",width = 4.5,height = 5)
  plot(1-roc.GSVA$specificities, roc.GSVA$sensitivities, #画训练集
      col="blue", xlab="1-Specificity (FPR)", main="", ylab="Sensitivity (TPR)",
      lwd=2, type="l",  xlim=c(0,1),ylim=c(0,1))
  lines(x=1-roc.COX$specificities,y=roc.COX$sensitivities, #补测试集
        lwd=2,type="l",col=ggplot2::alpha("magenta",0.7)) #曲线初始和终末可能重合，设置透明色

  lines(x=c(0,1),y=c(0,1),lwd=1.5,lty=2,col="grey40") #补中斜线
  legend("bottomright", bty="n", #添加图例
        fill=c("blue","magenta"), 
        legend.paste,
        cex=.8, border=NA, y.intersp=1, x.intersp=0.2 )
  invisible(dev.off())
}

MulRoc(
  "./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/GSVA_BestSeparation.txt" ,
  "./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/COX_BestSeparation.txt" , 
  "./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Lasso_BestSeparation.txt",
  "./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/randomForestSRC_BestSeparation.txt",
  './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/'
)
