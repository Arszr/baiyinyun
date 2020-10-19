library(survival)
library(survminer)
library(randomForestSRC)
library(timeROC)
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(caret)


clin <- read.csv('./web_app/temp/泛素化/2.单因素COX分析/all_dataset.csv')
genename <- read.csv('./web_app/temp/泛素化/2.单因素COX分析/Coxoutput_train.csv')

colnames(clin)[2] <- 'OS.time'
colnames(clin)[4] <- 'OS'
clin_col <- colnames(clin)
gene <- genename[,1]
s_data <- clin[,c(clin_col[1],clin_col[2],clin_col[4],gene)]

s_data <- na.omit(s_data)

s_data[,clin_col[2]] <- s_data[,clin_col[2]]/365 ##转化为年为单位的时间

s_data <- s_data[s_data$OS.time!=0,] ##删除生存时间为0的样本
s_data <- s_data[!duplicated(s_data[,1]),]
#length(unique(s_data[,clin_col[1]])) ## 查看有无重复样本
rownames(s_data)=s_data[,1]  #取出第一列
s_data=s_data[,-1]          #将第一列删除
head(s_data)
set.seed(2020911)
size <- createDataPartition(s_data$OS,p=0.7,list = F)  # 随机选择70%的数据作为Train data
train <- s_data[size,]
validation <- s_data[-size,]
s_data_col <- colnames(s_data)
window_block <- s_data_col[3:length(s_data_col)]
##构建随机生存森林
min_obbs <- data.frame(obb=0.0,gene='')
for (k in 3:length(s_data_col)-2) {
  
  
    print(k)
    t_min_obb <- data.frame(obb=0.0,gene='')
    
    for (num_ in 0:(length(s_data_col)-2-k)) {
      
        print(window_block[(1+num_):(k+num_)])
        x_col <- window_block[(1+num_):(k+num_)]
        
        Surv_fun <- as.formula(paste0('Surv(OS.time,OS)~',paste(x_col,sep = '',collapse = '+')))
        rsf_t <- rfsrc(Surv_fun,data = train,
                        ntree = 1000,
                        nodesize = 25,##该值建议多调整  
                        splitrule = 'logrank',
                        importance = T,
                        proximity = T,
                        forest = T,
                        seed = 2020911)
        min_obb <- min(rsf_t$err.block.rate)
        t_min_obb <- rbind(t_min_obb,c(min_obb,paste(x_col,collapse = '\t')))

    }
    t_min_obb <- t_min_obb[-1,]
    key <- which(t_min_obb$obb== min(t_min_obb$obb), arr.ind = TRUE)
    min_obbs<- rbind(min_obbs,t_min_obb[key,])
    
    #break
    #   plot(rsf_t)
}
min_obbs <- min_obbs[-1,]
write.csv(min_obbs,'obb_set.csv',row.names=F)

path=paste(save_path,'obb_set.pdf',sep='/')
# pdf(paste(save_path,'line.pdf',sep='/'))
pdf(dca_path,width = 6,height = 6)
par(bty = "o", mgp = c(1.5,.33,0), mar = c(3,4,1,2),las = 1, tcl = -.25)
plot(1:length(min_obbs$obb),min_obbs$obb,
     xlab = "Number of Genes",
     ylab = "",
     type = "l",
     lty = 4,
     col = "#00aef9", # 可修改线的颜色
     cex = 1.5)
mtext("OOB error rate",side = 2,line = 2.5,las = 3) # 添加y标签
points(1:length(min_obbs$obb),min_obbs$obb,# 加圆圈
       col = "#00aef9", # 可修改圆圈颜色
       pch = 19)
abline(v = which.min(min_obbs$obb),lty = 2, col = "#00aef9") 

dev.off()

v <- which.min(min_obbs$obb)
str_ <- paste('gene',min_obbs$gene[v],sep = '\n')
str_ <- gsub("\t","\n",str_)
str_ <- gsub('"','',str_)
write.table(str_,'min_obb_set.csv',quote = F,col.names = F,row.names=F)
