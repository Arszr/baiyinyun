options(stringsAsFactors = F)
library(tidyverse)
library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(stringr)
library(ggplot2)
library(survival)
library(survminer)
library(future.apply)

Survival_plot <- function(ossur, save_path) {
    ossur <- ossur
    gsva_km_path <- save_path

    for (i in length(colnames(ossur))) {
        g <- colnames(ossur)[i]
        ossur$GSVA <- ifelse(as.numeric(ossur[, g]) > median(as.numeric(ossur[, g])), "high", "low")

        sur_data <- survfit(Surv(time, A2_Event) ~ GSVA, data = ossur)

        ggsurvplot(sur_data,
            data = ossur,
            conf.int = F, pval = TRUE,
            fontsize = 20, font.legend = 20, pval.size = 7,
            font.main = c(20),
            font.x = c(20),
            font.y = c(20),
            font.tickslab = c(18)
        ) +
            labs(title = paste0("OS Survival Analysis of ", g), x = "Time/Months", y = "Survival probability")
        ggsave(paste0(gsva_km_path, colnames(ossur)[i], "_OS.pdf"), width = 8, height = 6)
    }
}

# exp_data="web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Fig1/train.csv"
# cox_train="./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Fig1/CoxSingle_train.csv"
# gene_deg="./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Fig1/gene_set_DEG.csv"
# pheno_path="./web_app/data/Difference/TCGA-BRCA/pheno.csv"
# c_pdata="./web_app/data/disease/clinical/TCGA-BRCA.csv"
# save_path="./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Fig1/"

Survival_ <- function(exp_data, cox_train, gene_deg, pheno_path, c_pdata, save_path) {
    Exp <- read.csv(exp_data, row.names = 1, header = T, check.names = F)[,-c(1:12)]
    Exp1 <- data.frame(t(Exp), check.names = F)
    Exp1_col <- colnames(Exp1)
    Coxoutput_train <- read.csv(cox_train, header = T, fill = TRUE, check.names = F)
    Coxoutput_train_col <- colnames(Coxoutput_train)
    gene_set_DEG <- read.csv(gene_deg, row.names = 1, check.names = F)
    gene_set_DEG_col <- colnames(gene_set_DEG)
    Cox_genes <- merge(Coxoutput_train, gene_set_DEG, by.x = Coxoutput_train_col[1], by.y = gene_set_DEG_col[1])

    Cox_genes_exp = merge(Cox_genes[,1:2],Exp1,by.x=Coxoutput_train_col[1],by.y="row.names")#13
    row.names(Cox_genes_exp) = Cox_genes_exp[,1]
    Cox_genes_exp = Cox_genes_exp[,-c(1:2)]

    coxexp <- paste0(save_path, "Cox_genes_exp.csv")
    write.csv(Cox_genes_exp, file = coxexp, row.names = T, quote = F) # 13

    Cox_genes_exp1 <- data.frame(t(Cox_genes_exp), check.names = F)
    Cox_genes_exp1$sample <- row.names(Cox_genes_exp1)
    
    meta <- read.csv(c_pdata, header = T, fill = TRUE, check.names = F)
    # print(1)
    OS <- na.omit(meta[, c(1:3)])
    os_col <- colnames(OS)
    OS[,os_col[2]] <- as.numeric(as.character(OS[,os_col[2]])) ### 鍥犲瓙鏁板€煎<U+FFFD>??
    OS$time <- abs(OS[os_col[2]]) / 30 ## abs缁濆<U+FFFD>?<U+FFFD>??

    # 棰勫悗鏄捐憲宸紓鍩哄洜涓庣敓瀛樹俊鎭悎<U+FFFD>??
    phe <- OS
    phe <- unique(merge(phe, Cox_genes_exp1, by.x = os_col[1], by.y = "sample"))
    phe_col <- colnames(phe)
    Cox_genes_exp1 <- Cox_genes_exp[, phe[,phe_col[1]]]


    # 鎵归噺杈撳嚭OS鐢熷瓨鏇<U+FFFD>?嚎閲岀殑P<U+FFFD>??
    # install.packages("future.apply")

    plan(multiprocess)

    Gene <- data.frame(unique(row.names(Cox_genes_exp1)))
    row.names(Gene) <- Gene[, 1]
    genes <- row.names(Gene)

    phe[phe_col[3]]  <-  ifelse(phe[phe_col[3]]=="Alive",0,ifelse(phe[phe_col[3]]=="Dead",1,NA))
    phe[phe_col[4]] <- unlist(phe[phe_col[4]])
    system.time(res1 <- future_lapply(1:length(genes), function(i) {
        group <- ifelse(as.numeric(phe[, genes[i]]) > median(as.numeric(phe[, genes[i]])), "high", "low")
        # if(length(table(group))==1) return(NULL)
        surv <- as.formula(paste("Surv(time, A2_Event)~", "group"))
        data <- cbind(phe[, 3:4], group)
        x <- survdiff(surv, data = data)
        pValue <- 1 - pchisq(x$chisq, df = 1)
        return(c(genes[i], pValue))
    }))
    res1 <- data.frame(do.call(rbind, res1))
    colnames(res1) <- c("Gene", "OS_pValue")
    res1 <- res1[order(res1$OS_pValue), ]
    pvalue_p <- paste(save_path, "Cox_genes_OS_pValue.csv", sep = "")
    write.csv(res1, pvalue_p, row.names = F, quote = F) # 18
    # 批量输出基因的生存曲线

    for (i in 5:ncol(phe)) {
        g <- colnames(phe)[i]
        phe$Exp <- ifelse(as.numeric(phe[, g]) > median(as.numeric(phe[, g])), "high", "low")
        ggsurvplot(survfit(Surv(time, A2_Event) ~ Exp, data = phe),
            conf.int = F, pval = TRUE,
            data = phe,
            legend.title = "",
            fontsize = 20, font.legend = 20, pval.size = 7,
            font.main = c(20),
            font.x = c(20),
            font.y = c(20),
            font.tickslab = c(18)
        ) +
            labs(title = paste0("OS Survival Analysis of ", g), x = "Time/Months", y = "Survival probability")
        ggsave(paste0(save_path, colnames(phe)[i], "_OS.pdf"), width = 8, height = 6)
    }

    # GSVA鍒嗘<U+FFFD>?

    # 璇<U+FFFD>?<U+FFFD>叆鍩哄洜闆嗙殑鍩哄洜鍒楄<U+FFFD>?
    gene_set <- data.frame(Genes = "Genes_GSVA", t(unique(row.names(Cox_genes_exp))))

    gs <- list()
    for (i in 1:nrow(gene_set)) {
        gs[[gene_set[i, 1]]] <- t(gene_set[i, 2:ncol(gene_set)]) %>% .[grep("\\S", .)] ## 绛<U+FFFD>?<U+FFFD>簬绌虹殑涓嶈
    } ## 杞崲鍒楄〃锛屽彲浠ョ<U+FFFD>?(unlist(strsplit(test_1[,2],split="/")))
    # 璇<U+FFFD>?<U+FFFD>叆鍩哄洜琛ㄨ揪<U+FFFD>?╅樀
    exprSet <- Exp1
    # 杩欎竴鍙ュ氨瀹屾垚浜咷SVA鍒嗘<U+FFFD>?
    gsva_es <- gsva(as.matrix(exprSet), gs)
    # 淇濆瓨鍒版枃<U+FFFD>??
    gsva_path <- paste(save_path, "Cox_gene_set_gsva.csv", sep = "")
    write.csv(gsva_es, gsva_path, quote = F)


    # 鐢熷瓨淇℃伅涓嶨SVA鍚堝<U+FFFD>?
    GSVA <- read.csv(gsva_path, row.names = 1, header = T, fill = TRUE, check.names = F)
    GSVA <- data.frame(t(GSVA), check.names = F)
    GSVA$sample <- substr(row.names(GSVA), 1, 12)
    ossur <- OS
    ossur_col <- colnames(ossur)
    # print(phe)
    ossur <- unique(merge(ossur, GSVA, by.x = ossur_col[1], by.y = "sample"))
    ossur_col <- colnames(ossur)
    ossur[ossur_col[3]]  <-  ifelse(ossur[ossur_col[3]]=="Alive",0,ifelse(ossur[ossur_col[3]]=="Dead",1,NA))
    ossur[ossur_col[4]] <- unlist(ossur[ossur_col[4]])
    Survival_plot(ossur, save_path)

    return(c(coxexp,pvalue_p, gsva_path))
}

# x <- Survival_(
#     "web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Fig1/train.csv",
#     "./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Fig1/CoxSingle_train.csv",
#     "./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Fig1/gene_set_DEG.csv",
#     "./web_app/data/Difference/TCGA-BRCA/pheno.csv",
#     "./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Fig1/pdata.csv",
#     "./web_app/temp/Arsz/Ubiquitination/TCGA-BRCA/Fig1/"
# )