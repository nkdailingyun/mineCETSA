#' ms_clean
#'
#' Function to clean up read-in dataset, mainly to remove not-quantified entries
#'
#' @param data dataset to be cleaned up
#' @param nread number of reading channels or sample treatements, default value 10
#' @param remkeratin whether to remove Keratin protein, which is generally considered as a common contaminant, default set to TRUE
#' @param remsinglecondprot whether the orphan proteins that appear only in one of the datasets should be removed, default set to FALSE
#'
#' @import dplyr
#'
#' @export
#' @return a dataframe after clean-up
#' @examples \dontrun{
#'  ITDRdata_cleaned <- ms_clean(ITDRdata)
#' }
#'
#'

ms_clean <- function(data, nread=10, remkeratin=TRUE, remsinglecondprot=FALSE) {

  # add variable name to output
  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  nrowdata <- nrow(data)
  row_na <- apply(data, 1, function(x) {any(is.na(x[c(1:(nread+3))]))})
  if (sum(row_na)!=0) {
    data_na <- data[which(row_na), ]
    data <- data[-which(row_na), ]
    ms_filewrite(data_na, paste0(dataname, "_Missing Values.txt"), outdir=outdir)
  } else {
    print("Your data is complete in quantification, there is no NA values!")
  }

  # remove single condition proteins if desired
  if (remsinglecondprot) {
    counttable <- data %>% group_by(id) %>% summarize(count=n()) %>% filter(count > 1)
    fkeep <- which(data$id %in% counttable$id)
    data <- data[fkeep, ]
    if (nrow(data) == 0) {
      stop("Single conditioned dataset cannot set remsinglecondprot to TRUE!")
    }
  }
  # to remove contaminating keratin proteins
  if (remkeratin) {
    pattern <- grep("Keratin", data$description)
    # print(length(pattern))
    if (length(pattern) > 0) {
      data <- data[-pattern, ]
    }
    if (nrow(data) == 0) {
      stop("After removing Keratin the dataset is empty.")
    }
  }
  print(paste("The data composition under each experimental condition (after cleanup) is:"))
  print(table(data$condition))
  if (length(attr(data,"outdir"))==0  & length(outdir)>0) {
    attr(data,"outdir") <- outdir
  }
  return(data)
}


#' ms_conditionrename
#'
#' Function to customize the condition names in dataset,
#' make sure both order and length of \code{incondition} and \code{outcondition} should match exactly.
#'
#' @param data dataset of which the condition names to be customized
#' @param incondition a vector of current condition naming to be changed
#' @param outcondition a vector of new condition naming to be applied
#'
#'
#' @export
#' @return a dataframe after condition name customization
#' @examples \dontrun{
#'  ITDRdata <- ms_conditionrename(ITDRdata, incondition=c("staurosporine.3","staurosporine.4"), outcondition=c("ST.1","ST.2"))
#' }
#'
#'
ms_conditionrename <- function(data, incondition=NULL, outcondition=NULL) {
  #if(length(unique(data$condition)) != length(incondition)){
  #	stop("Condition lengths do not match!")
  #}
  if (length(incondition) != length(outcondition)) {
    stop("The number of in and out conditions provided don't match!")
  }
  for (i in 1: length(incondition)) {
    if (!(incondition[i] %in% (unique(data$condition)))) {
      stop("The to-be-renamed conditions don't match the presenting ones!")
    } else {
      data$condition <- gsub(paste0("^",incondition[i],"$"), as.character(outcondition[i]), data$condition)
    }
  }
  print("The data composition under each experimental condition (after renaming) is:")
  print(table(data$condition))
  return(data)
}
