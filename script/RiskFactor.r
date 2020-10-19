library(reshape2)
library(ggplot2)
library(scales)
library(cowplot)

options(stringsAsFactors = FALSE) #禁止chr转成factor
#输入文件:data包含生存信息、表达量、risk score;   bestvars特征基因名称
#读取训练集763个样本的RNA表达谱
RiskF <- function(train_data,save_path) {

  data.train <- read.csv(train_data,row.names = 1,check.names = F)

  #1、GSVA评分
  GSVA_BestSeparation <- read.table(paste0(save_path,'GSVA_BestSeparation.txt'),sep = "\t", stringsAsFactors = F,header = T, check.names=F,fill=TRUE)
  data <- merge(GSVA_BestSeparation,data.train[,-c(1:12)],by=0)
  row.names(data)=data[,1]
  data = na.omit(data[order(data$Genes_GSVA),-1])

  # 6个错误率最低的重要性基因作为bestvars特征基因
  bestvars <- read.table(paste0(save_path,"signature with importance.txt"), sep = "\t", stringsAsFactors = F,header = T, check.names=F,fill=TRUE)$signature

  #按risk score分组：这里用中值做cutoff分组。
  #还可以用X-tile找最佳分组cutoff；或用FigureYa4bestSeparation找到P value最小的分组；
  #或用FigureYa35batch_bestSeparation从一群基因中批量找最佳分组，并筛选出高低两组有显著差异的基因。
  # risk score，用于画顶部散点图
  rs <- data$Genes_GSVA
  names(rs) <- rownames(data)
  # 使用最佳分组
  rs_data <- data.frame(x=1:length(rs),rs=as.numeric(sort(rs)),Risk=data$gp)
  # follow-up，用于画中间B图
  surv_data <- data.frame(x=1:length(rs),
                      t=data[names(sort(rs)),'A1_OS']/365*12,
                      s=data[names(sort(rs)),'event']) 
  surv_data$Status <- as.factor(ifelse(surv_data$s==0,'Alive','Death'))
  # 提取signature对应的data，并按risk score排序，用于画底部热图
  exp_data <- data[names(sort(rs)),which(colnames(data) %in% bestvars)]

  #开始画图：分别画出最上方的risk score、中间的follow-up和最下面的signature，最后拼图。
  #A - risk score
  plot.A <- ggplot(rs_data, aes(x=x,y=rs))+
    geom_point(aes(col=Risk),size=0.5)+
    scale_color_manual(labels=c("GSVA_High","GSVA_Low"), 
                      #guide_legend(guide = NULL), #如果不想画图例就删掉#
                      name="GSVA Score", values =c("#FF9289", "#02A2BD")) + 
    
    # 画竖向虚线
    geom_segment(aes(x = sum(rs_data$Risk=="low"),
                    y = -1, 
                    xend = sum(rs_data$Risk=="low"), 
                    yend = max(rs_data$rs)), linetype="dashed", size = 0.6)+

    theme(axis.title.x=element_blank()) +
    scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
    labs(y="GSVA Score",x="",fill="Risk") +
    #scale_colour_discrete(name="Risk scores") +
    theme_classic() +
    theme(axis.ticks.x=element_blank(),
          axis.line = element_blank(), #如果想像example2那样画坐标轴，就删掉这行
          axis.text.x=element_blank())

  #B - follow-up
  plot.B <- ggplot(surv_data,aes(x=x,y=t))+
    geom_point(aes(col=Status),size=0.5)+
    geom_vline(aes(xintercept=sum(rs_data$Risk=="low")),size=0.6,linetype="dashed")+
    scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
    scale_color_manual(labels=c("Alive","Dead"),
                      values =c("#00A087FF","#DC0000FF"))+
    labs(y="OS(months)",x="")+
    theme_classic()+
    theme(axis.ticks.x=element_blank(),
          axis.line = element_blank(), #如果想像example2那样不画坐标轴，就删掉前面的#
          axis.text.x=element_blank())

  #C - signature
  tmp <- t(scale(exp_data))
  tmp[tmp > 1] = 1
  tmp[tmp < -1] = -1
  reorder_cormat <- function(cormat){
    dd <- dist(cormat)
    hc <- hclust(dd,method = "average")
    cormat <-cormat[hc$order,]
  }
  tmp1 <- reorder_cormat(tmp)
  tmp.m <- melt(tmp1)
  p2 <-ggplot(tmp.m, aes(Var2, Var1),size=0.5) + 
    geom_tile(aes(fill = value)) 

  plot.C <- p2 + scale_fill_gradient2(name="Genes\nexpression", low="#3C3CFF", high="#FF4A4A", mid="white") +
    labs(x = "", y = "")+
    theme_classic()+
    theme(legend.title = element_text(size = 12),legend.position = "right",
          axis.ticks=element_blank(), axis.text.x=element_blank(),
          axis.line = element_blank())+
    #geom_vline(aes(xintercept=sum(rs_data$Risk=="Low-risk")), linetype="dashed",size=0.6,col="black")+
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))

  annotation_col<- data.frame(Var1 = rs_data$Risk, 
                              Var2=c(1:length(colnames(tmp))),
                              value=1)
  plot.C
  #分类
  p4 <- ggplot(tmp.m, aes(Var2, Var1),size=0.5)+
    geom_raster(aes(Var2,value,fill = Var1), data= annotation_col)+ 
    labs(x = "", y = "")+
    theme_void()+
    theme(legend.position = "right", legend.title = element_blank())+
    scale_x_discrete(name="")+
    scale_fill_manual(labels=c("GSVA_High","GSVA_Low"),
                      values =c("#FF9289", "#02A2BD"))+
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))

  #拼图
  plot_grid(plot.A,plot.B,p4,NULL,plot.C,
            labels = c("A","B","C",""),
            rel_heights = c(1,1,0.3,-0.05,0.8),
            label_x=0,
            label_y=1,
            align = 'v',ncol = 1,greedy = F)
  # 保存到文件
  ggsave(paste0(save_path,"GSVA_BestSeparation_risk.pdf"), width = 7, height = 9)

  #2、COX评分
  COX_BestSeparation <- read.table(paste0(save_path,'COX_BestSeparation.txt'),sep = "\t", stringsAsFactors = F,header = T, check.names=F,fill=TRUE)
  data <- merge(COX_BestSeparation,data.train[,-c(1:12)],by=0)
  row.names(data)=data[,1]
  data = na.omit(data[order(data$riskscore),-1])

  # 6个错误率最低的重要性基因作为bestvars特征基因
  bestvars <- read.table(paste0(save_path,"signature with importance.txt"), sep = "\t", stringsAsFactors = F,header = T, check.names=F,fill=TRUE)$signature

  #按risk score分组：这里用中值做cutoff分组。
  #还可以用X-tile找最佳分组cutoff；或用FigureYa4bestSeparation找到P value最小的分组；
  #或用FigureYa35batch_bestSeparation从一群基因中批量找最佳分组，并筛选出高低两组有显著差异的基因。
  # risk score，用于画顶部散点图
  rs <- data$riskscore
  names(rs) <- rownames(data)
  # 使用最佳分组
  rs_data <- data.frame(x=1:length(rs),rs=as.numeric(sort(rs)),Risk=data$gp)
  # follow-up，用于画中间B图
  surv_data <- data.frame(x=1:length(rs),
                          t=data[names(sort(rs)),'A1_OS']/365*12,
                          s=data[names(sort(rs)),'event']) 
  surv_data$Status <- as.factor(ifelse(surv_data$s==0,'Alive','Death'))
  # 提取signature对应的data，并按risk score排序，用于画底部热图
  exp_data <- data[names(sort(rs)),which(colnames(data) %in% bestvars)]

  #开始画图：分别画出最上方的risk score、中间的follow-up和最下面的signature，最后拼图。
  #A - risk score
  plot.A <- ggplot(rs_data, aes(x=x,y=rs))+
    geom_point(aes(col=Risk),size=0.5)+
    scale_color_manual(labels=c("High-risk","Low-risk"), 
                      #guide_legend(guide = NULL), #如果不想画图例就删掉#
                      name="COX Risk score", values =c("#FF9289", "#02A2BD")) + 
    
    # 画竖向虚线
    geom_segment(aes(x = sum(rs_data$Risk=="low"),
                    y = 0, 
                    xend = sum(rs_data$Risk=="low"), 
                    yend = max(rs_data$rs)), linetype="dashed", size = 0.6)+

  theme(axis.title.x=element_blank()) +
    scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
    labs(y="COX Risk score",x="",fill="Risk") +
    #scale_colour_discrete(name="Risk scores") +
    theme_classic() +
    theme(axis.ticks.x=element_blank(),
          axis.line = element_blank(), #如果想像example2那样画坐标轴，就删掉这行
          axis.text.x=element_blank())

  #B - follow-up
  plot.B <- ggplot(surv_data,aes(x=x,y=t))+
    geom_point(aes(col=Status),size=0.5)+
    geom_vline(aes(xintercept=sum(rs_data$Risk=="low")),size=0.6,linetype="dashed")+
    scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
    scale_color_manual(labels=c("Alive","Dead"),
                      values =c("#00A087FF","#DC0000FF"))+
    labs(y="OS(months)",x="")+
    theme_classic()+
    theme(axis.ticks.x=element_blank(),
          axis.line = element_blank(), #如果想像example2那样不画坐标轴，就删掉前面的#
          axis.text.x=element_blank())

  #C - signature
  tmp <- t(scale(exp_data))
  tmp[tmp > 1] = 1
  tmp[tmp < -1] = -1
  reorder_cormat <- function(cormat){
    dd <- dist(cormat)
    hc <- hclust(dd,method = "average")
    cormat <-cormat[hc$order,]
  }
  tmp1 <- reorder_cormat(tmp)
  tmp.m <- melt(tmp1)
  p2 <-ggplot(tmp.m, aes(Var2, Var1),size=0.5) + 
    geom_tile(aes(fill = value)) 

  plot.C <- p2 + scale_fill_gradient2(name="Genes\nexpression", low="#3C3CFF", high="#FF4A4A", mid="white") +
    labs(x = "", y = "")+
    theme_classic()+
    theme(legend.title = element_text(size = 12),legend.position = "right",
          axis.ticks=element_blank(), axis.text.x=element_blank(),
          axis.line = element_blank())+
    #geom_vline(aes(xintercept=sum(rs_data$Risk=="Low-risk")), linetype="dashed",size=0.6,col="black")+
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))

  annotation_col<- data.frame(Var1 = rs_data$Risk, 
                              Var2=c(1:length(colnames(tmp))),
                              value=1)

  #分类
  p4 <- ggplot(tmp.m, aes(Var2, Var1),size=0.5)+
    geom_raster(aes(Var2,value,fill = Var1), data= annotation_col)+ 
    labs(x = "", y = "")+
    theme_void()+
    theme(legend.position = "right", legend.title = element_blank())+
    scale_x_discrete(name="")+
    scale_fill_manual(labels=c("High-risk","Low-risk"),
                      values =c("#FF9289", "#02A2BD"))+
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))

  #拼图
  plot_grid(plot.A,plot.B,p4,NULL,plot.C,
            labels = c("A","B","C",""),
            rel_heights = c(1,1,0.3,-0.05,0.8),
            label_x=0,
            label_y=1,
            align = 'v',ncol = 1,greedy = F)
  # 保存到文件
  ggsave(paste0(save_path,"COX_BestSeparation_risk.pdf"), width = 7, height = 9)

  #3、Lasso评分
  Lasso_BestSeparation <- read.table(paste0(save_path,'Lasso_BestSeparation.txt'),sep = "\t", stringsAsFactors = F,header = T, check.names=F,fill=TRUE)
  data <- merge(Lasso_BestSeparation,data.train[,-c(1:12)],by=0)
  row.names(data)=data[,1]
  data = na.omit(data[order(data$lasso_Score),-1])

  # 6个错误率最低的重要性基因作为bestvars特征基因
  bestvars <- read.table(paste0(save_path,"signature with importance.txt"), sep = "\t", stringsAsFactors = F,header = T, check.names=F,fill=TRUE)$signature

  #按risk score分组：这里用中值做cutoff分组。
  #还可以用X-tile找最佳分组cutoff；或用FigureYa4bestSeparation找到P value最小的分组；
  #或用FigureYa35batch_bestSeparation从一群基因中批量找最佳分组，并筛选出高低两组有显著差异的基因。
  # risk score，用于画顶部散点图
  rs <- data$lasso_Score
  names(rs) <- rownames(data)
  # 使用最佳分组
  rs_data <- data.frame(x=1:length(rs),rs=as.numeric(sort(rs)),Risk=data$gp)
  # follow-up，用于画中间B图
  surv_data <- data.frame(x=1:length(rs),
                          t=data[names(sort(rs)),'A1_OS']/365*12,
                          s=data[names(sort(rs)),'event']) 
  surv_data$Status <- as.factor(ifelse(surv_data$s==0,'Alive','Death'))
  # 提取signature对应的data，并按risk score排序，用于画底部热图
  exp_data <- data[names(sort(rs)),which(colnames(data) %in% bestvars)]

  #开始画图：分别画出最上方的risk score、中间的follow-up和最下面的signature，最后拼图。
  #A - risk score
  plot.A <- ggplot(rs_data, aes(x=x,y=rs))+
    geom_point(aes(col=Risk),size=0.5)+
    scale_color_manual(labels=c("High-risk","Low-risk"), 
                      #guide_legend(guide = NULL), #如果不想画图例就删掉#
                      name="Lasso Risk score", values =c("#FF9289", "#02A2BD")) + 
    
    # 画竖向虚线
    geom_segment(aes(x = sum(rs_data$Risk=="low"),
                    y = 0, 
                    xend = sum(rs_data$Risk=="low"), 
                    yend = max(rs_data$rs)), linetype="dashed", size = 0.6)+

  theme(axis.title.x=element_blank()) +
    scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
    labs(y="Lasso Risk score",x="",fill="Risk") +
    #scale_colour_discrete(name="Risk scores") +
    theme_classic() +
    theme(axis.ticks.x=element_blank(),
          axis.line = element_blank(), #如果想像example2那样画坐标轴，就删掉这行
          axis.text.x=element_blank())

  #B - follow-up
  plot.B <- ggplot(surv_data,aes(x=x,y=t))+
    geom_point(aes(col=Status),size=0.5)+
    geom_vline(aes(xintercept=sum(rs_data$Risk=="low")),size=0.6,linetype="dashed")+
    scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
    scale_color_manual(labels=c("Alive","Dead"),
                      values =c("#00A087FF","#DC0000FF"))+
    labs(y="OS(months)",x="")+
    theme_classic()+
    theme(axis.ticks.x=element_blank(),
          axis.line = element_blank(), #如果想像example2那样不画坐标轴，就删掉前面的#
          axis.text.x=element_blank())

  #C - signature
  tmp <- t(scale(exp_data))
  tmp[tmp > 1] = 1
  tmp[tmp < -1] = -1
  reorder_cormat <- function(cormat){
    dd <- dist(cormat)
    hc <- hclust(dd,method = "average")
    cormat <-cormat[hc$order,]
  }
  tmp1 <- reorder_cormat(tmp)
  tmp.m <- melt(tmp1)
  p2 <-ggplot(tmp.m, aes(Var2, Var1),size=0.5) + 
    geom_tile(aes(fill = value)) 

  plot.C <- p2 + scale_fill_gradient2(name="Genes\nexpression", low="#3C3CFF", high="#FF4A4A", mid="white") +
    labs(x = "", y = "")+
    theme_classic()+
    theme(legend.title = element_text(size = 12),legend.position = "right",
          axis.ticks=element_blank(), axis.text.x=element_blank(),
          axis.line = element_blank())+
    #geom_vline(aes(xintercept=sum(rs_data$Risk=="Low-risk")), linetype="dashed",size=0.6,col="black")+
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))

  annotation_col<- data.frame(Var1 = rs_data$Risk, 
                              Var2=c(1:length(colnames(tmp))),
                              value=1)

  #分类
  p4 <- ggplot(tmp.m, aes(Var2, Var1),size=0.5)+
    geom_raster(aes(Var2,value,fill = Var1), data= annotation_col)+ 
    labs(x = "", y = "")+
    theme_void()+
    theme(legend.position = "right", legend.title = element_blank())+
    scale_x_discrete(name="")+
    scale_fill_manual(labels=c("High-risk","Low-risk"),
                      values =c("#FF9289", "#02A2BD"))+
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))

  #拼图
  plot_grid(plot.A,plot.B,p4,NULL,plot.C,
            labels = c("A","B","C",""),
            rel_heights = c(1,1,0.3,-0.05,0.8),
            label_x=0,
            label_y=1,
            align = 'v',ncol = 1,greedy = F)
  # 保存到文件
  ggsave(paste0(save_path,"Lasso_BestSeparation_risk.pdf"), width = 7, height = 9)


  #4、randomForestSRC评分
  randomForestSRC_BestSeparation <- read.table(paste0(save_path,'randomForestSRC_BestSeparation.txt'),sep = "\t", stringsAsFactors = F,header = T, check.names=F,fill=TRUE)
  data <- merge(randomForestSRC_BestSeparation,data.train[,-c(1:12)],by=0)
  row.names(data)=data[,1]
  data = na.omit(data[order(data$SRC_Score),-1])

  # 6个错误率最低的重要性基因作为bestvars特征基因
  bestvars <- read.table(paste0(save_path,"signature with importance.txt"), sep = "\t", stringsAsFactors = F,header = T, check.names=F,fill=TRUE)$signature

  #按risk score分组：这里用中值做cutoff分组。
  #还可以用X-tile找最佳分组cutoff；或用FigureYa4bestSeparation找到P value最小的分组；
  #或用FigureYa35batch_bestSeparation从一群基因中批量找最佳分组，并筛选出高低两组有显著差异的基因。
  # risk score，用于画顶部散点图
  rs <- data$SRC_Score
  names(rs) <- rownames(data)
  # 使用最佳分组
  rs_data <- data.frame(x=1:length(rs),rs=as.numeric(sort(rs)),Risk=data$gp)
  # follow-up，用于画中间B图
  surv_data <- data.frame(x=1:length(rs),
                          t=data[names(sort(rs)),'A1_OS']/365*12,
                          s=data[names(sort(rs)),'event']) 
  surv_data$Status <- as.factor(ifelse(surv_data$s==0,'Alive','Death'))
  # 提取signature对应的data，并按risk score排序，用于画底部热图
  exp_data <- data[names(sort(rs)),which(colnames(data) %in% bestvars)]

  #开始画图：分别画出最上方的risk score、中间的follow-up和最下面的signature，最后拼图。
  #A - risk score
  plot.A <- ggplot(rs_data, aes(x=x,y=rs))+
    geom_point(aes(col=Risk),size=0.5)+
    scale_color_manual(labels=c("High-risk","Low-risk"), 
                      #guide_legend(guide = NULL), #如果不想画图例就删掉#
                      name="randomForestSRC\nRisk score", values =c("#FF9289", "#02A2BD")) + 
    
    # 画竖向虚线
    geom_segment(aes(x = sum(rs_data$Risk=="low"),
                    y = 0, 
                    xend = sum(rs_data$Risk=="low"), 
                    yend = max(rs_data$rs)), linetype="dashed", size = 0.6)+

  theme(axis.title.x=element_blank()) +
    scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
    labs(y="randomForestSRC Risk score",x="",fill="Risk") +
    #scale_colour_discrete(name="Risk scores") +
    theme_classic() +
    theme(axis.ticks.x=element_blank(),
          axis.line = element_blank(), #如果想像example2那样画坐标轴，就删掉这行
          axis.text.x=element_blank())
  #B - follow-up
  plot.B <- ggplot(surv_data,aes(x=x,y=t))+
    geom_point(aes(col=Status),size=0.5)+
    geom_vline(aes(xintercept=sum(rs_data$Risk=="low")),size=0.6,linetype="dashed")+
    scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
    scale_color_manual(labels=c("Alive","Dead"),
                      values =c("#00A087FF","#DC0000FF"))+
    labs(y="OS(months)",x="")+
    theme_classic()+
    theme(axis.ticks.x=element_blank(),
          axis.line = element_blank(), #如果想像example2那样不画坐标轴，就删掉前面的#
          axis.text.x=element_blank())

  #C - signature
  tmp <- t(scale(exp_data))
  tmp[tmp > 1] = 1
  tmp[tmp < -1] = -1
  reorder_cormat <- function(cormat){
    dd <- dist(cormat)
    hc <- hclust(dd,method = "average")
    cormat <-cormat[hc$order,]
  }
  tmp1 <- reorder_cormat(tmp)
  tmp.m <- melt(tmp1)
  p2 <-ggplot(tmp.m, aes(Var2, Var1),size=0.5) + 
    geom_tile(aes(fill = value)) 

  plot.C <- p2 + scale_fill_gradient2(name="Genes\nexpression", low="#3C3CFF", high="#FF4A4A", mid="white") +
    labs(x = "", y = "")+
    theme_classic()+
    theme(legend.title = element_text(size = 12),legend.position = "right",
          axis.ticks=element_blank(), axis.text.x=element_blank(),
          axis.line = element_blank())+
    #geom_vline(aes(xintercept=sum(rs_data$Risk=="Low-risk")), linetype="dashed",size=0.6,col="black")+
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))

  annotation_col<- data.frame(Var1 = rs_data$Risk, 
                              Var2=c(1:length(colnames(tmp))),
                              value=1)
  #分类
  p4 <- ggplot(tmp.m, aes(Var2, Var1),size=0.5)+
    geom_raster(aes(Var2,value,fill = Var1), data= annotation_col)+ 
    labs(x = "", y = "")+
    theme_void()+
    theme(legend.position = "right", legend.title = element_blank())+
    scale_x_discrete(name="")+
    scale_fill_manual(labels=c("High-risk","Low-risk"),
                      values =c("#FF9289", "#02A2BD"))+
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))

  #拼图
  plot_grid(plot.A,plot.B,p4,NULL,plot.C,
            labels = c("A","B","C",""),
            rel_heights = c(1,1,0.3,-0.05,0.8),
            label_x=0,
            label_y=1,
            align = 'v',ncol = 1,greedy = F)
  # 保存到文件
  ggsave(paste0(save_path,"randomForestSRC_BestSeparation_risk.pdf"), width = 7, height = 9)

}

# RiskF('./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/train.csv','./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/')