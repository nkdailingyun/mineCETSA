#' ms_melt_curve_summary
#'
#' Function to summarize the melt curve fitting parameters and
#' perform statistical analysis for each protein
#'
#' @param data dataset after ms_fitting() function
#' @param nread number of reading channels or sample treatements,
#' default value 10
#' @param contrast a character to indicate the contrasting treatment conditions
#' typically in a format similar as "Drug-Ctrl"
#' @param writetofile whether to keep a local file copy of the summary,
#' default set to TRUE
#'
#' @import tidyr
#' @import dplyr
#' @importFrom tibble as_tibble
#' @export
#' @return a dataframe
#' @examples \dontrun{
#' LY_summary <- ms_melt_curve_summary(LY_fitted)
#' }
#'
#'

ms_melt_curve_summary <- function(data, nread=10, contrast=NULL, writetofile=TRUE) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  ncond <- length(unique(data$condition)) # unique condition names including replicate info
  ctrl <- strsplit(contrast,"-")[[1]][2]
  treatment <- strsplit(contrast,"-")[[1]][1]
  if (anyNA(c(ctrl,treatment))) { stop("Make sure you correctly specify the contrast") }

  if (length(grep("countNum", names(data)))) {
    countinfo <- unique(data[ ,c("id","sumUniPeps","sumPSMs","countNum")])
    data1 <- data[ ,!(names(data) %in% c("sumUniPeps","sumPSMs","countNum"))]
  }

  if (length(grep("description", names(data1)))) {
    proteininfo <- unique(data1[ ,c("id","description")])
    data1$description <- NULL
  }

  # data1 <- tidyr::unite(data1, "id", id, description, sep="_")
  # Use this to accommodate dataset which use gene instead of full protein description

  # subsetting complete replicates:
  keep1 <- which(table(data1$id) == ncond)
  if(length(keep1)==0) { stop("No proteins contain complete replicate") }
  nkeep1 <- names(keep1)
  fkeep1 <- which(data1$id %in% nkeep1)
  data_complete <- data1[fkeep1, ]
  cat("The number of proteins with complete replicates is:",
               length(unique(data_complete$id)), "\n")
  cat("Only these proteins will be summarized.\n")

  data2 <- tidyr::separate(data_complete, condition, into=c("condition", "replicate"), sep="\\.")
  ncondition <- length(unique(data2$condition))
  nreplicate <- length(unique(data2$replicate))
  uniquecond <- unique(data2[ ,c("condition", "replicate")])
  row.names(uniquecond) <- NULL
  uniquecondlength <- nrow(uniquecond)
  cat("Replicates information were extracted as follows:\n")
  print(as.data.frame(uniquecond))

  # Calculate euclidean distance metrics
  dism <- plyr::ddply(data_complete, "id", calEDscore, nread, nreplicate)
  dism <- dism[ ,c("id","variance","distance","EDscore")]

  # Calculate delta Tm
  Tmtable <- data2 %>% group_by(id, condition) %>% summarise(Tm=mean(Tm,na.rm=TRUE))
  Tmtable <- spread(Tmtable, condition, Tm)
  Tmtable <- mutate(Tmtable, diff=eval(parse(text=treatment))-eval(parse(text=ctrl)))
  names(Tmtable)[c(2:4)] <- paste0("Tm_",names(Tmtable)[c(2:4)])

  # merge into final reporting table
  dism <- merge(dism, Tmtable)
  dism <- merge(proteininfo, dism)
  # dism <- tidyr::separate(dism, id, into=c("id","description"), sep="_")

  # nametempvector <- names(data1)[3:(nread+2)]
  # numtempvector <- as.numeric(nametempvector)
  # data1 <- data1[ ,-c(3:(nread+2))]
  # print(head(data1))

  return(dism)

}
