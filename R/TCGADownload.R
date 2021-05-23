library("optparse")
library(TCGAbiolinks)
library(tidyverse)
library(reshape2)





report <- function(...){
  if(!exists("IS_DEBUG")){
    IS_DEBUG <<-T
  }
  if(!exists("LOG_FILENAME")){
    LOG_FILENAME <<- "report.txt"
  }
  if(exists("IS_EMPTY")&&IS_EMPTY==T){
    cat("",file=LOG_FILENAME)
  }
  report <- function(...){
    if(IS_DEBUG){
      cat(...,"\n")
    }else{
      cat(...,"\n",file=LOG_FILENAME,append=T,sep="")
    }
  }
}



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

