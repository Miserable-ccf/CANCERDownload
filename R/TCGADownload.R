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
  get_metadata <- function(metadata,filename){
    message("重新处理：",paste0("result/",filename,".rda"))
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
    saveRDS(metadata,file =paste0("result/",filename,".rda"))
    write.csv(metadata,file = paste0("result/",filename,".csv"),quote = F)
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
clinical <- function(project,filename=NULL){
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






