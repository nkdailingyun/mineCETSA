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
  outdir <- data$outdir[1]

  # to prevent the change of sub-directory folder
  if (!length(outdir)) {
    outdir <- paste0(dataname,"_",format(Sys.time(), "%y%m%d_%H%M"))
    dir.create(outdir)
  } else if (dir.exists(outdir)==FALSE) {
    dir.create(outdir)
  }

  data$outdir <- NULL

  # To look for the ids found in same conditions
  uniqueid <- unique(gsub("-[0-9]+", "", data$id))

  questionid <- NULL
  for (uniid in uniqueid) {
    subdata <- data[grep(uniid, data$id), ]$condition
    if (anyDuplicated(subdata)) {
      questionid <- c(questionid, uniid)
    } else{
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
  resolvetable <- counttable2 %>% top_n(1, count) %>% top_n(-1, id)
  resolvetable <- merge(resolvetable, unique(data[ ,c(1,2)]), all.x=TRUE)
  ambitable <- counttable2[ ,c(1:2)]
  ambitable <- merge(ambitable, unique(data[ ,c(1,2)]), all.x=TRUE)
  names(ambitable)[c(1:2)] <- c("isoforms", "frequency")
  print(paste0(nrow(ambitable), " Isoforms were identified from dataset ", dataname, "."))
  print("Check out the details in the current working directory.")
  write.table(ambitable, paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"), dataname, "_solved_isoforms.txt"),
              sep="\t", row.names=FALSE, quote=FALSE)

  #return(resolvetable)
  for (i in 1:nrow(resolvetable)) {
    names <- resolvetable[i,c(1,3,5)]
    id <- names[2]
    data[grep(paste0("^", id), data$id), c(1,2)] <- names[c(1,3)]
  }

  print(paste0(nrow(ambitable), " Isoforms were solved."))

  print(paste0("Still, Pay attention to these ", length(questionid), " base IDs for possible protein grouping ambuiguity: "))
  print(questionid)

  print("To further solve isoform ambiguity issue, double check and modify the matching table accordingly in the current working directory.")
  print("Then proceed to use ms_isoform_consolidate() function.")

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
    top_n(1, frequency) %>% top_n(1, totalCountNum)
  questionid_table3 <- questionid_table3[ ,c(1,4)]
  names(questionid_table3) <- c("Tobe_id","groupid")
  questionid_table4 <- merge(questionid_table2, questionid_table3, by.x="uniid", by.y="groupid")
  questionid_table4$uniid <- NULL

  write.table(questionid_table4, paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"), dataname, "_tobe_consolidated.txt"),
              sep="\t", row.names=FALSE, quote=FALSE)
  #return(questionid_table)
  data$outdir <- outdir
  ms_filewrite(data, paste0(dataname,"_isoform_solved1.txt"))

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
  outdir <- data$outdir[1]

  # to prevent the change of sub-directory folder
  if (!length(outdir)) {
    outdir <- paste0(dataname,"_",format(Sys.time(), "%y%m%d_%H%M"))
    dir.create(outdir)
  } else if (dir.exists(outdir)==FALSE) {
    dir.create(outdir)
  }
  data$outdir <- NULL

  matchtable <- read.delim(file=matchtable, quote="", na.string="", as.is=T, check.names=F)
  matchtable$id <- gsub("_", "-", matchtable$id)
  matchtable$Tobe_id <- gsub("_", "-", matchtable$Tobe_id)
  stopifnot("id" %in% names(matchtable), "Tobe_id" %in% names(matchtable))
  #return(matchtable)

  originalname <- unique(data[ ,c(1:2)])
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
    d1[grep(paste0("^", unsolvedid, "$"), d1$id), 1] <- tobeid
  }

  d1_aveaged <- plyr::ddply(d1, plyr::.(id,condition,treatment), .drop=TRUE,
                            .fun = function(xx) {
                              c(mean = weighted.mean(xx[["reading"]], xx[["countNum"]], na.rm=TRUE),
                                sumUniPeps_new = sum(xx[["sumUniPeps"]]),
                                sumPSMs_new = sum(xx[["sumPSMs"]]),
                                countNum_new = sum(xx[["countNum"]])
                              )
                            }
  )
  #return(d1_aveaged)
  d1_aveaged_w <- tidyr::spread(d1_aveaged[ ,c(1:4)], treatment, mean)
  d1_aveaged_w <- merge(d1_aveaged_w, unique(d1_aveaged[ ,c(1,2,5:7)]), all=FALSE)

  data <- merge(d1_aveaged_w, originalname)
  data <- data[ ,c(1, ncol(data), 2:(ncol(data)-1))]
  names(data) <- gsub("_new","", names(data))
  data$outdir <- outdir
  return(data)

}
