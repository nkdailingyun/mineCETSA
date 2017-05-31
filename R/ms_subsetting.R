#' ms_subsetting
#'
#' Function to parse data against a list to subset the data
#'
#' @param data dataset to be subsetted
#' @param hitidlist a list of hit UniprotID to parse against the dataset
#' @param dose a vector of doses of compound applied to ITDR samples, in the same order as channels
#' @param isfile whether the provided hitidlist is in a txt file under column name "id", default set to TRUE
#' @param allisoform whether to retrieve other isoforms of the same parental Uniprot ID, default set to TRUE
#' @param revsel short for reverse selection, when set to TRUE, the output is the dataset after removing the targets, default is set to FALSE
#'
#'
#' @importFrom readr read_tsv
#' @export
#' @return a dataframe
#' @examples \dontrun{
#' IITDRdata_subset <- ms_subsetting(ITDRdata_f[[1]], hitidlist="hit_list.txt", isfile=TRUE)
#' ITDRdata_subset <- ms_subsetting(ITDRdata_f[[1]], hitidlist=c("P00000", "P12345-6"), isfile=FALSE)
#' }
#'
#'


ms_subsetting <- function(data, hitidlist, isfile=TRUE, allisoform=TRUE, revsel=FALSE) {

  #the hitlist contains the hits Uniprot id under a column name of "id"
  #listif <- read.table( file=hitlist, sep="\t", header=T, quote="", na.string="", as.is=T, check.names=F )$id
  dataname <- deparse(substitute(data))
  if (isfile) {
    listid <- readr::read_tsv(file=hitidlist)$id
  } else {
    listid <- hitidlist
  }
  #return(listid)
  listid <- gsub("_","-",listid) # It happens that sometimes "-" was read in as "_"
  if (allisoform) { listid <- gsub("-[0-9]*","",listid) }
  listid <- gsub("^\\s+|\\s+$","",listid) # To make sure no empty space was introduced in original input
  listid <- na.omit(listid) # remove NAs
  llength <- length(listid)
  print(paste0("Found ", llength, " valid hit ids to parse again the data ", dataname))

  fkeep <- NULL
  for (i in 1:llength) {
    if (!allisoform) {hits <- grep(paste0("^", listid[i], "$"), data$id, value=FALSE)}
    else {hits <- grep(paste0("^", listid[i]), data$id, value=FALSE)}
    fkeep <- c(fkeep, hits)
    # print(fkeep)
  }
  if (length(fkeep) > 0) {
    # remove double entries
    fkeep <- unique(fkeep)
    # keep listed proteins or do the reverse selection
    if (revsel) {
      data <- data[-fkeep, ]
    } else {
      data <- data[fkeep, ]
    }
    print(paste0("Retrieved ", nrow(data), " data entries from ", dataname))
  } else {
    stop("Opps, no matches were found, pls double check!")
  }
  return(data)
}
