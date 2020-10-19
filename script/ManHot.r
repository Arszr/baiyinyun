#============================================================
#人体循环系统表达图

#devtools::install_github("jespermaag/gganatogram")
library(gganatogram)
library(stringr)
library(gridExtra)
library(ggplot2)
library(ggpolypath)
library(gganatogram)
library(dplyr)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
#输入文件预处理
#如果你的数据已经整理成very_easy_input.csv的格式，就可以跳过这步，进入“开始画图”。
#此处以FigureYa55panCancer输出的的基因在TCGA + GTEx的表达量TPM值（easy_input.csv文件）为例，需要把TCGA的癌症名称跟包里的organ对应上。


#1、GSVA评分
ManHot <- function(gavad,coxd,lassod,rfd,organd,save_path) {
GSVA_BestSeparation <- read.table(gavad,sep = "\t", stringsAsFactors = F,header = T, check.names=F,fill=TRUE)
TCGA_organ <- read.csv(organd, header = T, check.names=F)

## 将GSVA评分按高低组分别取平均值
df = GSVA_BestSeparation[,3:4]
df<-aggregate(.~gp,df,mean)
df$group = "Genes_GSVA"
colnames(df)=c("type","tpm","group")

organ = TCGA_organ[TCGA_organ$TCGA=="BRCA",]
df$organ=organ$organ
df = df[,c(4,1,2)]
df$type<-factor(df$type,levels = c("low","high"))

### 把TCGA癌症名称缩写换成gganatogram包里的organ
#根据背景知识，已整理成TCGA_organ
#Name列来源：https://cn.bing.com/translator
TCGA.organ.tpm <- df
# 男人
hgMale_key_tpm <- merge(hgMale_key, TCGA.organ.tpm, by = "organ")
hgMale_key_tpm$value <- hgMale_key_tpm$tpm #替换掉原来的value
hgMale_key_tpm$tpm <- NULL

#开始画图
# very_easy_input.csv，organ对应的数值。
# 第一列是组织器官名，必须跟包里的organ一致
# 第三列color，每种organ给一种颜色
# 第二列组织所在的系统
# 第四列数值，可以是基因表达量，或者其他临床指标。
# 第五列分组tumor和normal
# 其中第一列和第四列为必需,此处以人类为例，分别画男人和女人，对比tumor和normal。
# 还可以男女都画，然后ps成一半男一半女。其他物种按照very_easy_input_*.csv的格式整理好数据，就可以套用了

#画男人
hgMale_key_tpm$type <- hgMale_key_tpm$type.y
hgMale_key_tpm$type<-factor(hgMale_key_tpm$type,levels = c("low","high"))

hgMale <- gganatogram(data=hgMale_key_tpm, fillOutline='white', organism='human', sex='male', fill="value") + 
  facet_wrap(~type) +
  scale_fill_gradient(low = "#1E90FF", high = "#FF0000") +
  labs(fill = "GSVA Score",title="  Genes_GSVA Score") + 
  coord_cartesian(ylim = c(-120, 0)) +
  theme_void()+
  theme(legend.text= element_text(color="black", size=12),
        legend.title=element_text(colour="black", size=14),
        title  = element_text(color="black", size=15)) +
  theme(plot.title = element_text(hjust = 0.5))#标题居中
hgMale

#组图并保存到pdf文件:此处只保存男性图
pdf(paste0(save_path,"gganatogram_Genes_GSVA.pdf"),width=9,height=5.5)
grid.arrange(hgMale, ncol=1)
dev.off()

#2、COX评分
COX_BestSeparation <- read.table(coxd,sep = "\t", stringsAsFactors = F,header = T, check.names=F,fill=TRUE)
TCGA_organ <- read.csv(organd, header = T, check.names=F)

## 将GSVA评分按高低组分别取平均值
df = COX_BestSeparation[,3:4]
df<-aggregate(.~gp,df,mean)
df$group = "COX"
colnames(df)=c("type","tpm","group")

organ = TCGA_organ[TCGA_organ$TCGA=="BRCA",]
df$organ=organ$organ
df = df[,c(4,1,2)]
df$type<-factor(df$type,levels = c("low","high"))

### 把TCGA癌症名称缩写换成gganatogram包里的organ
#根据背景知识，已整理成TCGA_organ
#Name列来源：https://cn.bing.com/translator
TCGA.organ.tpm <- df
# 男人
hgMale_key_tpm <- merge(hgMale_key, TCGA.organ.tpm, by = "organ")
hgMale_key_tpm$value <- hgMale_key_tpm$tpm #替换掉原来的value
hgMale_key_tpm$tpm <- NULL

#开始画图
# very_easy_input.csv，organ对应的数值。
# 第一列是组织器官名，必须跟包里的organ一致
# 第三列color，每种organ给一种颜色
# 第二列组织所在的系统
# 第四列数值，可以是基因表达量，或者其他临床指标。
# 第五列分组tumor和normal
# 其中第一列和第四列为必需,此处以人类为例，分别画男人和女人，对比tumor和normal。
# 还可以男女都画，然后ps成一半男一半女。其他物种按照very_easy_input_*.csv的格式整理好数据，就可以套用了

#画男人
hgMale_key_tpm$type <- hgMale_key_tpm$type.y
hgMale_key_tpm$type<-factor(hgMale_key_tpm$type,levels = c("low","high"))

hgMale <- gganatogram(data=hgMale_key_tpm, fillOutline='white', organism='human', sex='male', fill="value") + 
  facet_wrap(~type) +
  scale_fill_gradient(low = "#1E90FF", high = "#FF0000") +
  labs(fill = "COX riskScore",title="  COX riskScore") + 
  coord_cartesian(ylim = c(-120, 0)) +
  theme_void()+
  theme(legend.text= element_text(color="black", size=12),
        legend.title=element_text(colour="black", size=14),
        title  = element_text(color="black", size=15)) +
  theme(plot.title = element_text(hjust = 0.5))#标题居中
hgMale

#组图并保存到pdf文件:此处只保存男性图
pdf(paste0(save_path,"gganatogram_COX riskScore.pdf"),width=9,height=5.5)
grid.arrange(hgMale, ncol=1)
dev.off()




#3、Lasso评分
Lasso_BestSeparation <- read.table(lassod,sep = "\t", stringsAsFactors = F,header = T, check.names=F,fill=TRUE)
TCGA_organ <- read.csv(organd, header = T, check.names=F)

## 将GSVA评分按高低组分别取平均值
df = Lasso_BestSeparation[,3:4]
df<-aggregate(.~gp,df,mean)
df$group = "Lasso"
colnames(df)=c("type","tpm","group")

organ = TCGA_organ[TCGA_organ$TCGA=="BRCA",]
df$organ=organ$organ
df = df[,c(4,1,2)]
df$type<-factor(df$type,levels = c("low","high"))

### 把TCGA癌症名称缩写换成gganatogram包里的organ
#根据背景知识，已整理成TCGA_organ
#Name列来源：https://cn.bing.com/translator
TCGA.organ.tpm <- df
# 男人
hgMale_key_tpm <- merge(hgMale_key, TCGA.organ.tpm, by = "organ")
hgMale_key_tpm$value <- hgMale_key_tpm$tpm #替换掉原来的value
hgMale_key_tpm$tpm <- NULL

#开始画图
# very_easy_input.csv，organ对应的数值。
# 第一列是组织器官名，必须跟包里的organ一致
# 第三列color，每种organ给一种颜色
# 第二列组织所在的系统
# 第四列数值，可以是基因表达量，或者其他临床指标。
# 第五列分组tumor和normal
# 其中第一列和第四列为必需,此处以人类为例，分别画男人和女人，对比tumor和normal。
# 还可以男女都画，然后ps成一半男一半女。其他物种按照very_easy_input_*.csv的格式整理好数据，就可以套用了

#画男人
hgMale_key_tpm$type <- hgMale_key_tpm$type.y
hgMale_key_tpm$type<-factor(hgMale_key_tpm$type,levels = c("low","high"))

hgMale <- gganatogram(data=hgMale_key_tpm, fillOutline='white', organism='human', sex='male', fill="value") + 
  facet_wrap(~type) +
  scale_fill_gradient(low = "#1E90FF", high = "#FF0000") +
  labs(fill = "Lasso riskScore",title="  Lasso riskScore") + 
  coord_cartesian(ylim = c(-120, 0)) +
  theme_void()+
  theme(legend.text= element_text(color="black", size=12),
        legend.title=element_text(colour="black", size=14),
        title  = element_text(color="black", size=15)) +
  theme(plot.title = element_text(hjust = 0.5))#标题居中
hgMale

#组图并保存到pdf文件:此处只保存男性图
pdf(paste0(save_path,"gganatogram_Lasso_riskScore.pdf"),width=9,height=5.5)
grid.arrange(hgMale, ncol=1)
dev.off()





#4、randomForestSRC评分
randomForestSRC_BestSeparation <- read.table(rfd,sep = "\t", stringsAsFactors = F,header = T, check.names=F,fill=TRUE)
TCGA_organ <- read.csv(organd, header = T, check.names=F)

## 将GSVA评分按高低组分别取平均值
df = randomForestSRC_BestSeparation[,c(3,5)]
df<-aggregate(.~gp,df,mean)
colnames(df)=c("type","tpm")

organ = TCGA_organ[TCGA_organ$TCGA=="BRCA",]
df$organ=organ$organ
df$type<-factor(df$type,levels = c("low","high"))

### 把TCGA癌症名称缩写换成gganatogram包里的organ
#根据背景知识，已整理成TCGA_organ
#Name列来源：https://cn.bing.com/translator
TCGA.organ.tpm <- df
# 男人
hgMale_key_tpm <- merge(hgMale_key, TCGA.organ.tpm, by = "organ")
hgMale_key_tpm$value <- hgMale_key_tpm$tpm #替换掉原来的value
hgMale_key_tpm$tpm <- NULL

#开始画图
# very_easy_input.csv，organ对应的数值。
# 第一列是组织器官名，必须跟包里的organ一致
# 第三列color，每种organ给一种颜色
# 第二列组织所在的系统
# 第四列数值，可以是基因表达量，或者其他临床指标。
# 第五列分组tumor和normal
# 其中第一列和第四列为必需,此处以人类为例，分别画男人和女人，对比tumor和normal。
# 还可以男女都画，然后ps成一半男一半女。其他物种按照very_easy_input_*.csv的格式整理好数据，就可以套用了

#画男人
hgMale_key_tpm$type <- hgMale_key_tpm$type.y
hgMale_key_tpm$type<-factor(hgMale_key_tpm$type,levels = c("low","high"))

hgMale <- gganatogram(data=hgMale_key_tpm, fillOutline='white', organism='human', sex='male', fill="value") + 
  facet_wrap(~type) +
  scale_fill_gradient(low = "#1E90FF", high = "#FF0000") +
  labs(fill = "randomForestSRC riskScore",title="randomForestSRC\nriskScore") + 
  coord_cartesian(ylim = c(-120, 0)) +
  theme_void()+
  theme(legend.text= element_text(color="black", size=12),
        legend.title=element_text(colour="black", size=14),
        title  = element_text(color="black", size=15)) +
  theme(plot.title = element_text(hjust = 0.5))#标题居中
hgMale

#组图并保存到pdf文件:此处只保存男性图
pdf(paste0(save_path,"gganatogram_randomForestSRC_riskScore.pdf"),width=9,height=5)
grid.arrange(hgMale, ncol=1)
dev.off()

}

# ManHot(
#   "./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/GSVA_BestSeparation.txt" ,
#   "./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/COX_BestSeparation.txt" , 
#   "./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Lasso_BestSeparation.txt",
#   "./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/randomForestSRC_BestSeparation.txt",
#   './web_app/data/Difference/TCGA-BRCA/TCGA_organ.csv',
#   './web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/'
# )
