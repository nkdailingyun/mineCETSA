#' ms_2D_caldiff
#'
#' Function to reformat dataframe or matrix data into an expressionset (eSet) format
#'
#' @param data Normalized dataset to calculate the relative protein abundance differences
#' @param treatmentlevel a vector of treatment labels, the control should be the first element,
#' such as c("DMSO","TNFa","AT26533")
#'
#' @import dplyr Biobase
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- ms_2D_caldiff(MOLM, treatmentlevel=c("DMSO","TNFa","AT26533"))
#' }
#'
#'

ms_2D_caldiff <- function(data, treatmentlevel=NULL) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (length(treatmentlevel)==0) {
    stop("Need to specify a vector of treatments")
  }
  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
    data$description <- NULL
  }
  if (length(grep("countNum", names(data)))) {
    countinfo <- unique(data[ ,c("id","sumUniPeps","sumPSMs","countNum")])
  }

  data <- gather(data[ ,-c((ncol(data)-2):ncol(data))], condition, reading, -id)
  if (length(unlist(strsplit(data$condition[1], "_")))==4) {
    data <- tidyr::separate(data, condition, into=c("set","temperature","replicate","treatment"), sep="_")
    if(!identical(sort(unique(data$treatment)),sort(treatmentlevel))) {
      stop("Need to specify the right treatments")
    }
    data$treatment <- factor(data$treatment, levels=treatmentlevel)
    data1 <- plyr::ddply(data, c("set","temperature","replicate","id"), function(data) {
      data<-data[order(data$treatment), ]
      base=data$reading[1]
      data<-mutate(data, reading=reading-base)
    })
    data1 <- tidyr::unite(data1, condition, set, temperature, replicate, treatment, sep="_")
  } else if (length(unlist(strsplit(data$condition[1], "_")))==3) {
    data <- tidyr::separate(data, condition, into=c("temperature","replicate","treatment"), sep="_")
    if(!identical(sort(unique(data$treatment)),sort(treatmentlevel))) {
      stop("Need to specify the right treatments")
    }
    data$treatment <- factor(data$treatment, levels=treatmentlevel)
    data1 <- plyr::ddply(data, c("temperature","replicate","id"), function(data) {
      data<-data[order(data$treatment), ]
      base=data$reading[1]
      data<-mutate(data, reading=reading-base)
    })
    data1 <- tidyr::unite(data1, condition, temperature, replicate, treatment, sep="_")
  } else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }

  data1 <- tidyr::spread(data1, condition, reading)
  data1 <- merge(data1, countinfo)
  data1 <- merge(proteininfo, data1)

  if (length(attr(data1,"outdir"))==0 & length(outdir)>0) {
    attr(data1,"outdir") <- outdir
  }
  return(data1)
}
