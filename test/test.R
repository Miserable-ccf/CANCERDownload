library(CANCERDownload)
count <- getCount("TCGA-ESCA")
clinical <- getClinical("TCGA-ESCA")
miRNACount <- getmiRNACount("TCGA-ESCA")
clinical_count <- clinical_count(clinical,count@count)
deg <- DeSeq2Analysis(count)
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
coxMulti_res <- CoxMulti(clinical,count@count,gene)
ROC(coxMulti_res)
Survival(coxMulti_res)
