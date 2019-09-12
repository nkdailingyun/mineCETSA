#' ms_isoform_resolve
#'
#' Function to perform automated isoform ambiguity cleanup using the principle of parsimony
#' when two or more isoforms are identified among the search result datasets
#' but for each individual dataset, only one isoform is identified.
#'
#' @param data dataset to be clead-up from isoform ambiguity
#'
#' @import dplyr
#' @export
#' @return a dataframe
#' @examples \dontrun{
#' data_cleaned1 <- ms_isoform_resolve(data_cleaned)
#' }
#'
#'

ms_isoform_resolve <- function(data) {

  # add variable name to output
  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  ncondition <- length(unique(data$condition))

  # To look for the ids found in same conditions
  uniqueid <- unique(gsub("-[0-9]+", "", data$id))

  questionid <- NULL
  for (uniid in uniqueid) {
    subdata <- data[grep(uniid, data$id), ]$condition
    if (anyDuplicated(subdata)) {
      questionid <- c(questionid, uniid)
    } else {
      next
    }
  }
  #return(questionid)
  uniqueid1 <- setdiff(uniqueid, questionid)

  counttable <- data %>% group_by(id) %>% summarize(count=n())
  counttable$uniid <- gsub("-[0-9]+", "", counttable$id)

  # to look for the isoforms could be automatically resolved based on parsimony principle
  counttable1 <- subset(counttable, uniid %in% uniqueid1)
  counttable2 <- counttable1 %>% group_by(uniid) %>% mutate(unicount=n()) %>%
    filter(unicount>1) %>% arrange(uniid, desc(count))
  #return(counttable2)
  resolvetable <- counttable2 %>% top_n(1, count) %>% filter(!duplicated(uniid))#top_n(-1, id)
  # note the problem with ranking of character column
  resolvetable <- merge(resolvetable, unique(data[ ,c(1,2)]), all.x=TRUE)
  ambitable <- counttable2[ ,c(1:2)]
  ambitable <- merge(ambitable, unique(data[ ,c(1,2)]), all.x=TRUE)
  names(ambitable)[c(1:2)] <- c("isoforms", "frequency")
  cat(nrow(ambitable), "Isoforms were identified from dataset", dataname, ".\n")
  cat("Check out the details in the current working directory.\n")
  write.table(ambitable, paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"), dataname, "_solved_isoforms.txt"),
              sep="\t", row.names=FALSE, quote=FALSE)

  #return(resolvetable)
  for (i in seq_len(nrow(resolvetable))) {
    names <- resolvetable[i,c(1,3,5)]
    id <- names[2]
    data[grep(paste0("^", id), data$id), c(1,2)] <- names[c(1,3)]
  }

  cat(nrow(ambitable), "Isoforms were solved.\n")
  message("To further solve isoform ambiguity, proceed to use ms_isoform_consolidate() function.")

  questionid_pos <- NULL
  for (id in questionid) {
    questionid_pos <- c(questionid_pos, grep(id, data$id))
    #questionid_full <- c(questionid_full, unique(grep(id, data$id, value=TRUE)))
  }
  questionid_table <- data[questionid_pos, ]
  questionid_table2 <- questionid_table %>% group_by(id) %>%
    summarize(frequency=n(), totalCountNum=sum(countNum), uniid=unique(id)) %>%
    arrange(id, desc(totalCountNum))
  questionid_table2$uniid <- gsub("-[0-9]+", "", questionid_table2$uniid)
  questionid_table3 <- questionid_table2 %>% group_by(uniid) %>%
    top_n(1, frequency) %>% top_n(1, totalCountNum) %>% filter(frequency != ncondition)
  questionid_table3 <- questionid_table3[ ,c(1,4)]
  names(questionid_table3) <- c("Tobe_id","groupid")
  questionid_table4 <- merge(questionid_table2, questionid_table3, by.x="uniid", by.y="groupid")
  questionid_table4$uniid <- NULL

  cat(nrow(questionid_table3), "base IDs were found with possible protein grouping ambuiguity.\n")
  cat("Carefully check and modify the suggested to_be_consoidated matching table in the working directory.\n")

  write.table(questionid_table4, paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"), dataname, "_tobe_consolidated.txt"),
              sep="\t", row.names=FALSE, quote=FALSE)
  #return(questionid_table)
  if (length(attr(data,"outdir"))==0 & length(outdir)>0) {
    attr(data,"outdir") <- outdir
  }
  ms_filewrite(data, paste0(dataname,"_isoform_resolved.txt"))

  return(data)
}


#' ms_isoform_consolidate
#'
#' Function to perform isoform ambiguity cleanup step 2 using the principle of parsimony
#' when two or more isoforms are identified within the same search result dataset
#'
#' Note: Data must contain "sumUniPeps","sumPSMs","countNum" columns,
#' otherwise the mock named columns should be provided.
#'
#' @param data dataset to be clead-up from isoform ambiguity
#' @param matchtable an isoform substitution matching table
#' @param nread number of reading channels or sample treatements,
#' default value 10
#' @param withabd whether the dataset contains abundance data as in ITDR/ITTR data,
#' by default set to FALSE
#' @param weightbycountnum by default to consolidate the isoforms
#'
#' @importFrom tibble as_tibble
#' @importFrom gtools mixedorder mixedsort
#' @importFrom plyr . ddply
#' @import dplyr
#' @import tidyr
#' @export
#' @return a dataframe
#' @examples \dontrun{
#' data_cleaned2 <- ms_isoform_consolidate(data_cleaned1, matchtable="./subfolder/tobe_consolidated.txt")
#' data_cleaned2 <- ms_isoform_consolidate(data_cleaned1, matchtable="./subfolder/tobe_consolidated.txt", withabd=TRUE)
#' }
#'
#'

ms_isoform_consolidate <- function(data, matchtable, nread=10, withabd=FALSE, weightbycountnum=TRUE) {

  # add variable name to output
  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  matchtable <- read.delim(file=matchtable, quote="", na.string="", as.is=T, check.names=F)
  matchtable$id <- gsub("_", "-", matchtable$id)
  matchtable$Tobe_id <- gsub("_", "-", matchtable$Tobe_id)
  stopifnot("id" %in% names(matchtable), "Tobe_id" %in% names(matchtable))
  #return(matchtable)

  proteininfo <- unique(data[ ,c(1:2)])
  data$description <- NULL

  if (withabd) {
    d1 <- tidyr::gather(data[ ,c(1:(2*nread+2))], treatment, reading, -id, -condition)
  } else {
    d1 <- tidyr::gather(data[ ,c(1:(nread+2))], treatment, reading, -id, -condition)
  }
  if (weightbycountnum) {
    d1 <- merge(d1, data[ ,c("id","condition","sumUniPeps","sumPSMs","countNum")])
  } else {
    stop("To be implemented")
  }

  for (i in 1:nrow(matchtable)) {
    unsolvedid <- matchtable[i, "id"]
    tobeid <- matchtable[i,"Tobe_id"]
    d1[grep(paste0("^", unsolvedid, "$"), d1$id), "id"] <- tobeid
  }

  d1_averaged <- plyr::ddply(d1, plyr::.(id,condition,treatment), .drop=TRUE,
                            .fun = function(xx) {
                              c(mean = weighted.mean(xx[["reading"]], xx[["countNum"]], na.rm=TRUE),
                                sumUniPeps_new = sum(xx[["sumUniPeps"]]),
                                sumPSMs_new = sum(xx[["sumPSMs"]]),
                                countNum_new = sum(xx[["countNum"]])
                              )
                            }
  )
  #return(d1_averaged)
  d1_averaged_w <- tidyr::spread(d1_averaged[ ,c(1:4)], treatment, mean)
  d1_averaged_w <- merge(d1_averaged_w, unique(d1_averaged[ ,c(1,2,5:7)]), all=FALSE)

  data <- merge(proteininfo, d1_averaged_w)
  #data <- data[ ,c(1, ncol(data), 2:(ncol(data)-1))]
  names(data) <- gsub("_new","", names(data))
  if (withabd) {
    data <- data[ ,c(1:3, gtools::mixedorder(names(data)[c(4:(2*nread+3))])+3, (2*nread+4):ncol(data))]
  } else {
    data <- data[ ,c(1:3, gtools::mixedorder(names(data)[c(4:(nread+3))])+3, (nread+4):ncol(data))]
  }
  if (length(attr(data,"outdir"))==0 & length(outdir)>0) {
    attr(data,"outdir") <- outdir
  }
  ms_filewrite(data, paste0(dataname,"_isoform_consolidated.txt"))
  return(data)

}


#' ms_isoform_match
#'
#' Function to perform an matching of isoform ambiguity between two datasets, ie when
#' only another isoform of the protein in one dataset is found in the reference dataset,
#' this isoform is changed to be same as the one present in reference dataset
#'
#' @param data dataset to be matched to reference dataset to relieve the problem of isoform ambiguity
#' @param refdata dataset used as the check reference, typically the IMPRINTS dataset
#'
#' @import dplyr
#' @export
#' @return a dataframe
#' @examples \dontrun{
#' data_meltcurve1 <- ms_isoform_match(data_meltcurve, data_IMPRINTS)
#' }
#'
#'
ms_isoform_match <- function(data, refdata) {

  # add variable name to output
  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  # To look for the ids found in same conditions
  commonid <- intersect(data$id,refdata$id)
  commonidunique <- gsub("-[0-9]+", "",commonid)
  commonuniqueid <- intersect(unique(gsub("-[0-9]+", "", data$id)),
                              unique(gsub("-[0-9]+", "", refdata$id)))
  ambiid <- setdiff(refdata$id,data$id) # ids only present in reference data
  ambiuniqueid <- gsub("-[0-9]+", "",ambiid) # unique ids only present in reference data
  #print(intersect(ambiuniqueid, commonidunique))
  ambiuniqueid1 <- setdiff(ambiuniqueid, commonidunique) # to remove the already matched isoform
  ambiuniqueid2 <- intersect(ambiuniqueid1, commonuniqueid) # to keep the ones with possible isoform
  cat(length(ambiuniqueid2), "isoforms were found with possible matches.\n")

  refdata$uniid <- gsub("-[0-9]+", "", refdata$id)
  refdata1 <- subset(refdata, uniid%in%ambiuniqueid2)
  #print(nrow(refdata1))
  data$uniid <- gsub("-[0-9]+", "", data$id)
  data1 <- subset(data, uniid%in%ambiuniqueid2)
  # when there is multiple matches, keep the one with maximal countNum
  data1 <- data1 %>% group_by(uniid) %>% top_n(1, countNum) %>% ungroup()
  #print(nrow(data1))

  matchtable <- merge(data1[,c("id","description","uniid")],refdata1[,c("id","description","uniid")],by="uniid")
  if (length(attr(matchtable,"outdir"))==0 & length(outdir)>0) {
    attr(matchtable,"outdir") <- outdir
  }
  ms_filewrite(matchtable, paste0(dataname,"_isoform_matched.txt"))

  # to update the data file with the match table information
  for (i in seq_len(nrow(matchtable))) {
    names <- matchtable[i,c(1:5)]
    id <- names[2]
    data[grep(paste0("^", id), data$id), c(1,2)] <- names[c(4,5)]
  }
  data$uniid <- NULL
  if (length(attr(data,"outdir"))==0 & length(outdir)>0) {
    attr(data,"outdir") <- outdir
  }
  return(data)
}
