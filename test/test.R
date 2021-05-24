library(survivalROC)
library(survival)
library(survminer)
library(lattice)
library(DESeq2)
library(TCGAbiolinks)
library(tidyverse)
library(CANCERDownload)
count <- getCount("TCGA-ESCA")
clinical <- getClinical("TCGA-ESCA")
miRNACount <- getmiRNACount("TCGA-ESCA")
clinical_count <- clinical_count(clinical,count@count)
deg <- DeSeq2Analysis(count,project = "TCGA-ESCA")
lnRAN_mRNA_res <- lnRAN_mRNA(deg)
lnRAN_mRNA_res@lnRNA
lnRAN_mRNA_count_res <- lnRAN_mRNA(count)
lnRAN_mRNA_count_res@lnRNA


gene <- deg@deg%>%
  arrange(desc(log2FoldChange))%>%
  rownames_to_column("symbol")%>%
  dplyr::slice(1:5)%>%
  pull("symbol")

CoxSingle(clinical,count@count,gene)
coxMulti_res <- CoxMulti(clinical,count@count,"ERBB2")
ROC(coxMulti_res@data)
Survival(coxMulti_res)









