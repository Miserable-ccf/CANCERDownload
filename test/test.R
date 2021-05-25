library(survivalROC)
library(survival)
library(survminer)
library(lattice)
library(DESeq2)
library(TCGAbiolinks)
library(tidyverse)
library(CANCERDownload)

project <-   "TCGA-COAD"

count <- getCount(project)
clinical <- getClinical(project)
miRNACount <- getmiRNACount(project)
clinical_count <- clinical_count(clinical,count@count)
deg <- DeSeq2Analysis(count,project = project)
lnRAN_mRNA_res <- lnRAN_mRNA(deg)
lnRAN_mRNA_res@lnRNA
lnRAN_mRNA_res@mRNA
lnRAN_mRNA_count_res <- lnRAN_mRNA(count)
lnRAN_mRNA_count_res@lnRNA
lnRAN_mRNA_count_res@mRNA


if(F){
  gene <- deg@deg%>%
    arrange(desc(log2FoldChange))%>%
    rownames_to_column("symbol")%>%
    dplyr::slice(1:5)%>%
    pull("symbol")
}else{
  gene <- "METTL3,METTL14,METTL16,WTAP,VIRMA,RBM15,RBM15B,ZC3H13,FTO,ALKBH5,YTHDC1,YTHDC2,IGF2BP1,IGF2B,P2,IGF2BP3,YTHDF1,YTHDF2,YTHDF3,HNRNPC,HNRNPA2B1,RBMX"
  gene <- strsplit(gene,split = ",")[[1]]
}


gene <- gene[gene %in% rownames(count@count)]
CoxSingle(clinical,count@count,gene)
coxMulti_res <- CoxMulti(clinical,count@count, gene,cutoff=Inf)
ROC(coxMulti_res@data,predict.time=4)
Survival(coxMulti_res)



## 完整的流程
mircode <- read_tsv("/home/wangyang/workspace/bioinfo_analysis/Rscript/data/lncRNA_miRCODE_miRNA.csv")%>%
  plyr::rename(c(lncRNA="gene_symbol"))








