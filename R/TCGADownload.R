{
  TCGA_query <- function(project,data.type,workflow.type,filename){
    if(!file.exists(paste0("result/",filename,".rda"))){
      message("重新处理：",paste0("result/",filename,".rda"))
      query<- GDCquery(project = project,
                       data.category = "Transcriptome Profiling",
                       data.type = data.type,
                       workflow.type = workflow.type)
      saveRDS(query,file = paste0("result/",filename,".rda"))
      return(query)
    }else{
      message("加载本地数据: ",paste0("result/",filename,".rda"))
      return(readRDS(paste0("result/",filename,".rda")))
    }
  }
  #### 获得metadata信息

  TCAG_prepare <- function(query,filename){
    if(!file.exists(paste0("result/",filename,".rda"))){
      message("重新处理：",paste0("result/",filename,".rda"))
      data <- GDCprepare(query = query,
                         summarizedExperiment = F,
                         save = F)
      saveRDS(data,file =paste0("result/",filename,".rda"))
      return(data)
    }else{
      message("加载本地数据: ",paste0("result/",filename,".rda"))
      return(readRDS(paste0("result/",filename,".rda")))
    }

  }
  DESeq2_Analysis <- function(expr,filename,metadata){
    if(!file.exists(paste0("result/",filename,".rda"))){
      ## 差异基因表达分析
      library(DESeq2)
      if(!identical(colnames(expr),rownames(metadata))){
        message("expr的列名与metadata行名不一致！！")
        message("expr 列数:",dim(expr)[2])
        message("metadata行数:",dim(metadata)[1])
        expr <- expr[,rownames(metadata)]
      }

      metadata$group <- factor(metadata$group)
      dds <- DESeqDataSetFromMatrix(countData = expr,
                                    colData = metadata,
                                    design = ~ group)
      ## PCA图
      #vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
      #vsd_df <- assay(vsd)
      # pdf(file = paste0("result/",filename,".pdf"))
      #DESeq2::plotPCA(vsd, intgroup = "group")
      #dev.off()
      #message("写入图片:",paste0("result/",filename,".pdf"))
      keep <- rowSums(counts(dds) > 0) >= min(table(metadata$group))
      dds_filt <- dds[keep, ]
      dds2 <- DESeq(dds_filt, parallel = T)
      message(resultsNames(dds2))
      res <- results(dds2)
      deg <- as.data.frame(res)
      message("写入文件：",paste0("result/",filename,".rda")," ",paste0("result/",filename,".csv"))
      saveRDS(deg,file =paste0("result/",filename,".rda"))
      write.csv(deg,file=paste0("result/",filename,".csv"),quote = F,row.names = F)
      return(deg)
    }else{
      message("加载本地数据: ",paste0("result/",filename,".rda"))
      return(readRDS(paste0("result/",filename,".rda")))
    }

  }
  TCGA_metadata <- function(path,filename) {
    metadata <- jsonlite::read_json(path, simplifyVector = T)
    metadata <- tibble::tibble(
      file_name = metadata$file_name,
      md5sum = metadata$md5sum,
      cases = bind_rows(metadata$associated_entities)$entity_submitter_id,
      TCGA_id_full = cases,
      TCGA_id = stringr::str_sub(cases, 1, 16),
      patient_id = stringr::str_sub(TCGA_id, 1, 12),
      tissue_type_id = stringr::str_sub(TCGA_id, 14, 15),
      tissue_type = sapply(tissue_type_id, function(x) {
        switch(x,
               "01" = "Primary Solid Tumor",
               "02" = "Recurrent Solid Tumor",
               "03" = "Primary Blood Derived Cancer - Peripheral Blood",
               "05" = "Additional - New Primary",
               "06" = "Metastatic",
               "07" = "Additional Metastatic",
               "11" = "Solid Tissue Normal")}),
      group = ifelse(tissue_type_id == "11", "Normal", "Tumor"))%>%
      column_to_rownames("TCGA_id_full")
    saveRDS(metadata,file =paste0("result/",filename,".rda"))
    write.csv(metadata,file = paste0("result/",filename,".csv"),quote = F)
    return(metadata)
  }
  get_metadata <- function(metadata){
    metadata <- tibble::tibble(
      TCGA_id_full = metadata$cases,
      cases=metadata$cases,
      TCGA_id = stringr::str_sub(TCGA_id_full, 1, 16),
      patient_id = stringr::str_sub(TCGA_id, 1, 12),
      tissue_type_id = stringr::str_sub(TCGA_id, 14, 15),
      tissue_type = sapply(tissue_type_id, function(x) {
        switch(x,
               "01" = "Primary Solid Tumor",
               "02" = "Recurrent Solid Tumor",
               "03" = "Primary Blood Derived Cancer - Peripheral Blood",
               "05" = "Additional - New Primary",
               "06" = "Metastatic",
               "07" = "Additional Metastatic",
               "11" = "Solid Tissue Normal")}),
      group = ifelse(tissue_type_id == "11", "Normal", "Tumor"))%>%
      column_to_rownames("TCGA_id_full")
    return(metadata)
  }
}

#' TCGA Download



#' TCGA Download
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' TCGA()
TCGA <- function(){
  message("Hello TCGA")
  message( getwd())
}
TCGA()

#' Download clinical
#' @export
getClinical <- function(project,filename=NULL){
  if(!file.exists(project)){
    dir.create(project)
  }
  if(is.null(filename)){
    filename = paste0(project,"-clinical")
  }
  if(!file.exists(paste0(project,"/",filename,".rds"))){
    message("从网络下载 ",project," 数据")
    clinical <- GDCquery_clinic(project = project, type = "clinical")%>%
      mutate(futime=ifelse(vital_status=='Alive',days_to_last_follow_up,days_to_death),
             fustat = factor(vital_status,levels = c("Alive","Dead"),labels = c(0,1)),
             fustat=as.numeric(as.character(fustat)),
             futime_year=futime/365,
             patient_id=submitter_id,
             tumor_stage=case_when(tumor_stage=="stage i"~"i",
                                   grepl("stage i[^i ^v]+",tumor_stage)~"i",
                                   tumor_stage=="stage ii"~"ii",
                                   grepl("stage ii[^i]+",tumor_stage)~"ii",
                                   grepl("stage iii",tumor_stage)~"iii",
                                   grepl("stage iv",tumor_stage)~"iv",
                                   tumor_stage=="not reported"~"unknow"
             ))%>%
      dplyr::select(c("patient_id","futime","fustat","futime_year","tumor_stage","tumor_stage"))
    saveRDS(clinical,file = paste0(project,"/",filename,".rds"))
    return(clinical)
  }else{
    message("加载本地的 ",project," 数据")
    return(readRDS(paste0(project,"/",filename,".rds")))
  }
}

#' Download Count
#' @export
getCount <- function(project,filename=NULL,workflow.type = "HTSeq - Counts"){
  if(!file.exists(project)){
    dir.create(project)
  }
  if(is.null(filename)){
    filename = paste0(project,"-count-metadata")
  }
  if(!file.exists(paste0(project,"/",filename,".rds"))){
    message("从网络下载 ",project," 数据")
    query<- GDCquery(project = project,
                               data.category = "Transcriptome Profiling",
                               data.type = "Gene Expression Quantification",
                               workflow.type = workflow.type)
    message("TCGA Query Project: ",query$project[[1]])
    message("TCGA data.category: ",query$data.category)
    message("TCGA data.type: ",query$data.type)
    message("TCGA workflow.type: ",query$workflow.type)
    metadata <- query[[1]][[1]]%>%
      get_metadata()
    metadata_sample <- metadata%>%
      dplyr::select("cases")%>%
      mutate(sample= substr(cases,14,15))
    message("metadata中样本有 ",length(unique(rownames(metadata)))," 个, 矩阵大小为: ",
            paste0(dim(metadata),collapse = ","))
    message("TCGA 样本信息[01 原发实体瘤, 11正常组织]: ",
            paste0(names(table(metadata_sample$sample)),collapse = ","),"对应数量",
            paste0(table(metadata_sample$sample),collapse = ","))
    GDCdownload(query,directory=project)
    count <- GDCprepare(query = query,
                       summarizedExperiment = F,
                       save = F,
                       directory=project)
    #data <- data[,rownames(metadata)]
    message("表达矩阵中有样本 ",length(unique(colnames(count)))," 个, 矩阵大小为: ",paste0(dim(count),collapse = ","),
         "表达矩阵的列名与metadata行名: ",identical(rownames(metadata),colnames(count)))

    ### 基因注释
    message("开始使用http://wangyang-bucket.oss-cn-beijing.aliyuncs.com/gencode.gene.info.v22.tsv进行基因注释基因")
    gff_v22 <- read_tsv("http://wangyang-bucket.oss-cn-beijing.aliyuncs.com/gencode.gene.info.v22.tsv")
    id2symbol <- gff_v22 %>%
      dplyr::select(1, 2,"gene_type")
    count_type <-count%>%
      {.[1:(nrow(count)-5),]}%>%
      plyr::rename(c(X1="gene_id"))%>%
      inner_join(id2symbol,by = "gene_id")%>%
      {.[!duplicated(.$gene_name),]}%>%
      remove_rownames()%>%
      column_to_rownames("gene_name")%>%
      dplyr::select(-gene_id)
    count <- count_type %>% dplyr::select(-gene_type)
    cancerTable <- new("CANCERTable",count=count,metadata=metadata,gene_type=count_type$gene_type)
    saveRDS(cancerTable,file = paste0(project,"/",filename,".rds"))
    return(cancerTable)
  }else{
    message("加载本地的 ",project," 数据")
    return(readRDS(paste0(project,"/",filename,".rds")))
  }
}

if(F){

  project="TCGA-ESCA"
  a <- getCount("TCGA-ESCA",workflow.type="HTSeq - FPKM",filename = "fpkm-count")
  a <- getCount("TCGA-ESCA")

}



#' Download miRNACount
#' @export
getmiRNACount <- function(project,filename=NULL){
  if(!file.exists(project)){
    dir.create(project)
  }
  if(is.null(filename)){
    filename = project
  }
  if(!file.exists(paste0(project,"/",filename,"-miRNA-count-metadata.rds"))){
    message("从网络下载 ",project," 数据")
    query<- GDCquery(project = project,
                     data.category = "Transcriptome Profiling",
                     data.type = "miRNA Expression Quantification",
                     workflow.type = "BCGSC miRNA Profiling")
    message("TCGA Query Project: ",query$project[[1]])
    message("TCGA data.category: ",query$data.category)
    message("TCGA data.type: ",query$data.type)
    message("TCGA workflow.type: ",query$workflow.type)
    metadata <- query[[1]][[1]]%>%
      get_metadata()
    metadata_sample <- metadata%>%
      dplyr::select("cases")%>%
      mutate(sample= substr(cases,14,15))
    message("metadata中样本有 ",length(unique(rownames(metadata)))," 个, 矩阵大小为: ",
            paste0(dim(metadata),collapse = ","))
    message("TCGA 样本信息[01 原发实体瘤, 11正常组织]: ",
            paste0(names(table(metadata_sample$sample)),collapse = ","),"对应数量",
            paste0(table(metadata_sample$sample),collapse = ","))
    GDCdownload(query,directory=project)
    count <- GDCprepare(query = query,
                        summarizedExperiment = F,
                        save = F,
                        directory=project)%>%
      dplyr::select("miRNA_ID",starts_with("read_count"))%>%
      rename_at(vars(contains("read_count")), ~ substr(.,12,length(.)))%>%
      column_to_rownames("miRNA_ID")

    #data <- data[,rownames(metadata)]
    message("表达矩阵中有样本 ",length(unique(colnames(count)))," 个, 矩阵大小为: ",paste0(dim(count),collapse = ","),
            "表达矩阵的列名与metadata行名: ",identical(rownames(metadata),colnames(count)))


    cancerTable <- new("CANCERTable",count=count,metadata=metadata,gene_type="miRNA")
    saveRDS(cancerTable,file = paste0(project,"/",filename,"-miRNA-count-metadata.rds"))
    return(cancerTable)
  }else{
    message("加载本地的 ",project," 数据")
    return(readRDS(paste0(project,"/",filename,"-miRNA-count-metadata.rds")))
  }
}

#' clinical_count
#' @export
clinical_count <- function(clinical,count,symbol=NULL){
  if(is.null(symbol)){
    count <- count[1:5,]
  }else{
    count <- count[symbol,]
  }
  data <- count%>%
    dplyr::select(!matches("TCGA-([A-Z 0-9]*?)-([A-Z 0-9]*?)-(11.*?)-(.*?)"))%>%
    t()%>%
    as.data.frame()%>%
    rownames_to_column("TCGA_full_id")%>%
    mutate(patient_id = stringr::str_sub(TCGA_full_id, 1, 12),
           tissue_type_id = stringr::str_sub(TCGA_full_id, 14, 15))%>%
    left_join(clinical,by="patient_id")

  #message(paste0(table(d$tissue_type_id),collapse = ","))
  return(data)
}
if(F){
  clinical_count <- clinical_count(clinical = clinical,count = a@count)

}
#' DeSeq2Analysis
#' @export
DeSeq2Analysis <-function(count,project,filename=NULL){

  if(F){
    count <- a
    filename=NULL
  }
  if(!file.exists(project)){
    dir.create(project)
  }
  if(is.null(filename)){
    filename = paste0(project,"-deseq2")
  }
  if(!file.exists(paste0(project,"/",filename,".rds"))){
    message("重新进行差异基因处理")
    expr <- count@count
    metadata <- count@metadata
    if(!identical(colnames(expr),rownames(metadata))){
      message("expr的列名与metadata行名不一致！！")
      message("expr 列数:",dim(expr)[2])
      message("metadata行数:",dim(metadata)[1])
      expr <- expr[,rownames(metadata)]
    }

    metadata$group <- factor(metadata$group)
    dds <- DESeqDataSetFromMatrix(countData = expr,
                                  colData = metadata,
                                  design = ~ group)

    keep <- rowSums(counts(dds) > 0) >= min(table(metadata$group))
    dds_filt <- dds[keep, ]
    dds2 <- DESeq(dds_filt, parallel = T)
    message(resultsNames(dds2))
    res <- results(dds2)
    deg <- as.data.frame(res)

    count_gene_type <- data.frame(gene=rownames(count@count),gene_type=count@gene_type)
    deg_type <- deg%>%rownames_to_column("gene")%>%
      inner_join(count_gene_type,by="gene")%>%
      column_to_rownames("gene")
    deg <-deg_type%>%dplyr::select(-gene_type)
    DEGTable <-  new("DEGTable",deg=deg,gene_type=deg_type$gene_type)
    saveRDS(DEGTable,file = paste0(project,"/",filename,".rds"))
    return(DEGTable)
  }else{
    message("加载本地的数据")
    return(readRDS(paste0(project,"/",filename,".rds")))
  }
}

if(F){
  library(DESeq2)
  deg <- DeSeq2Analysis(a)
  identical(rownames(b@count),rownames(deg@deg))
}

setGeneric("lnRAN_mRNA",function(count_type) standardGeneric("lnRAN_mRNA"))

#' lnRAN_mRNA DEGTable
#' @export
setMethod("lnRAN_mRNA","DEGTable",function(count_type){
  if(F){
    count_type <- deg
  }
  lncRNA_types <- "3prime_overlapping_ncrna, antisense, lincRNA, macro_lncRNA, non_coding, processed_transcript, sense_intronic, sense_overlapping"
  lncRNA_types <- unlist(str_split(lncRNA_types, pattern = ", "))
  mRNA <- count_type@deg[count_type@gene_type=="protein_coding",]
  lncRNA <- count_type@deg[count_type@gene_type %in% lncRNA_types,]
  lnRAN_mRNA <- new("lnRAN_mRNA",lnRNA=lncRNA,mRNA=mRNA)
  return(lnRAN_mRNA)

})

#' lnRAN_mRNA CANCERTable
#' @export
setMethod("lnRAN_mRNA","CANCERTable",function(count_type){
  lncRNA_types <- "3prime_overlapping_ncrna, antisense, lincRNA, macro_lncRNA, non_coding, processed_transcript, sense_intronic, sense_overlapping"
  lncRNA_types <- unlist(str_split(lncRNA_types, pattern = ", "))
  mRNA <- count_type@count[count_type@gene_type=="protein_coding",]
  lncRNA <- count_type@count[count_type@gene_type %in% lncRNA_types,]
  lnRAN_mRNA <- new("lnRAN_mRNA",lnRNA=lncRNA,mRNA=mRNA,metadata=count_type@metadata)
  return(lnRAN_mRNA)
})
if(F){
  lnRAN_mRNA <- lnRAN_mRNA(deg)
  dim( lnRAN_mRNA@lnRNA)
  dim( lnRAN_mRNA@mRNA)
  dim( lnRAN_mRNA@metadata)
}

CoxSingle_ <- function(data,gene){
  outTab <- data.frame()
  for(i in gene){
    cox <- coxph(Surv(futime_year,fustat) ~ data[,i],data=data)
    coxSummary <- summary(cox)
    outTab <- rbind(outTab,cbind(id=i,
                                 HR=coxSummary$conf.int[,"exp(coef)"],
                                 HR.95L=coxSummary$conf.int[,"lower .95"],
                                 HR.95H=coxSummary$conf.int[,"upper .95"],
                                 pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
  }
  return(outTab)
}

#' CoxSingle
#' @export
CoxSingle <- function(clinical,count,gene){
  outTab <- data.frame()
  data <- clinical_count(clinical = clinical,count = count,symbol=gene)
  outTab<- CoxSingle_(data,gene)
  return(outTab)
}

#' CoxMulti
#' @export
CoxMulti <- function(clinical,count,gene){
  if(F){
    clinical <- clinical
    a <- getCount("TCGA-ESCA")
    count <- a@count
    gene <- deg@deg%>%
      arrange(desc(log2FoldChange))%>%
      rownames_to_column("symbol")%>%
      dplyr::slice(1:5)%>%
      pull("symbol")
  }
  data <- clinical_count(clinical = clinical,count = count,symbol=gene)
  coxSingle <- CoxSingle_(data,gene)
  gene <- coxSingle$id
  formula <-as.formula(paste0("Surv(futime_year,fustat)~",paste0(gene,collapse = "+")))
  cox <- coxph(formula,data=data)
  coxSummary <- summary(cox)
  outTab <- data.frame()
  outTab <- rbind(outTab,cbind(HR=coxSummary$conf.int[,"exp(coef)"],
                               HR.95L=coxSummary$conf.int[,"lower .95"],
                               HR.95H=coxSummary$conf.int[,"upper .95"],
                               pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
  riskScore <- predict(cox,type = "risk",newdata = data)

  res_data<- data%>%
    mutate(
      riskScoreNum=riskScore,
      riskScore = as.vector(ifelse(riskScore>median(riskScore),"high","low")))%>%
    dplyr::select(c("TCGA_full_id","futime","futime_year","fustat","riskScoreNum","riskScore"))
  if(F){
    ## 生存曲线
    fit <- surv_fit(Surv(futime_year,fustat) ~ riskScore, data = res_data)
    res <- ggsurvplot(fit, pval=T, risk.table=F,
                      risk.table.height = 0.3,data = res_data)
    print(res)
    ## ROC曲线
    roc <- survivalROC(Stime=res_data$futime_year ,status=res_data$fustat,marker=res_data$riskScoreNum,
                       predict.time=3,method="KM")
    library(lattice)
    res <- xyplot(TP~FP,roc,main=paste("ROC curve(","AUC=",round(roc$AUC,3),")"),
                  xlab="False postive rate",ylab="True positive rate",panel = function(x,y){
                    panel.xyplot(x,y,type="l",lwd=2,cex.main=1.3,cex.lab=1.2,cex.axis=1.2,font=1.2,col="red")
                    panel.abline(c(0,1))
                  })
    print(res)
  }
  coxMulti <- new("CoxMulti",single=coxSingle,multi=outTab,data=res_data)
  return(coxMulti)
}

#' ROC
#' @export
ROC <- function(data,predict.time=3){
  library(lattice)
  roc <- survivalROC(Stime=data$futime_year ,status=data$fustat,marker=data$riskScoreNum,
                     predict.time=predict.time,method="KM")
  res <- xyplot(TP~FP,roc,main=paste("ROC curve(","AUC=",round(roc$AUC,3),")"),
                xlab="False postive rate",ylab="True positive rate",panel = function(x,y){
                  panel.xyplot(x,y,type="l",lwd=2,cex.main=1.3,cex.lab=1.2,cex.axis=1.2,font=1.2,col="red")
                  panel.abline(c(0,1))
                })
  return(res)
}

#' Survival
#' @export
Survival <- function(data){
  data <- data@data
  fit <- surv_fit(Surv(futime_year,fustat) ~ riskScore, data = data)
  res <- ggsurvplot(fit, pval=T, risk.table=F,
                    risk.table.height = 0.3,data = data)
  return(res)
}


if(F){
  library(survivalROC)
  library(survival)
  library(survminer)
  deg <- DeSeq2Analysis(a)
  #  a <- getCount("TCGA-ESCA",workflow.type="HTSeq - FPKM",filename = "fpkm-count")
  a <- getCount("TCGA-ESCA")
  gene <- deg@deg%>%
    arrange(desc(log2FoldChange))%>%
    rownames_to_column("symbol")%>%
    dplyr::slice(1:5)%>%
    pull("symbol")

  clinical_count <- clinical_count(clinical = clinical,count = a@count,symbol=gene)
  #coxph(Surv(futime_year,fustat) ~ MAGEC2,data=clinical_count)
  CoxSingle(clinical,count,gene)
  boxplot(clinical_count[,"MAGEC2"],outline=F)
  coxMulti_res <- CoxMulti(clinical,count,gene)
  ROC(coxMulti_res)
  Survival(coxMulti_res)
}



#### 提取mRNA和lnRNA
