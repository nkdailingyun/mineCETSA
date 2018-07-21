#' ms_2D_caldiff
#'
#' Function to calculate the pair-wise (per replicate and temperature) protein abundance differences
#'
#' @param data Normalized dataset to calculate the relative protein abundance differences
#' @param treatmentlevel a vector of treatment labels, the control should be the first element,
#' such as c("DMSO","TNFa","AT26533")
#' @param withinrep whether the calculation of the relative protein abundance difference should
#' still within the same biorep, default set to TRUE, when the bioreps are balanced
#'
#' @import dplyr Biobase
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- ms_2D_caldiff(MOLM, treatmentlevel=c("DMSO","TNFa","AT26533"))
#' }
#'
#'

ms_2D_caldiff <- function(data, treatmentlevel=NULL, withinrep=TRUE) {

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
    data <- data[ ,!(names(data) %in% c("sumUniPeps","sumPSMs","countNum"))]
  }

  data <- gather(data, condition, reading, -id)
  if (length(unlist(strsplit(data$condition[1], "_")))==4) {
    data <- tidyr::separate(data, condition, into=c("set","temperature","replicate","treatment"), sep="_")
    if(!identical(sort(unique(data$treatment)),sort(treatmentlevel))) {
      stop("Need to specify the right treatments")
    }
    data$treatment <- factor(data$treatment, levels=treatmentlevel)
    if (withinrep) {
      data1 <- plyr::ddply(data, c("set","temperature","replicate","id"), function(data) {
        data<-data[order(data$treatment), ]
        base=data$reading[1]
        data<-mutate(data, reading=reading-base)
      })
    } else{
      data1 <- plyr::ddply(data1, c("set","temperature","id"), function(data) {
        data<-data[order(data$treatment), ]
        base=mean(subset(data, treatment==treatmentlevel[1])$reading,na.rm=T)
        data<-mutate(data, reading=ifelse(treatment==treatmentlevel[1], 0.0, reading-base))
      })
    }
    data1 <- tidyr::unite(data1, condition, set, temperature, replicate, treatment, sep="_")
  } else if (length(unlist(strsplit(data$condition[1], "_")))==3) {
    data <- tidyr::separate(data, condition, into=c("temperature","replicate","treatment"), sep="_")
    if(!identical(sort(unique(data$treatment)),sort(treatmentlevel))) {
      stop("Need to specify the right treatments")
    }
    data$treatment <- factor(data$treatment, levels=treatmentlevel)
    if (withinrep) {
      data1 <- plyr::ddply(data, c("temperature","replicate","id"), function(data) {
        data<-data[order(data$treatment), ]
        base=data$reading[1]
        data<-mutate(data, reading=reading-base)
      })
    } else {
      data1 <- plyr::ddply(data, c("temperature","id"), function(data) {
        data<-data[order(data$treatment), ]
        base=mean(subset(data, treatment==treatmentlevel[1])$reading,na.rm=T)
        data<-mutate(data, reading=ifelse(treatment==treatmentlevel[1], 0.0, reading-base))
      })
    }
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
