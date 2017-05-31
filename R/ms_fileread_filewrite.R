#' ms_filewrite
#'
#' Function to write dataset into tab delimited files
#'
#' @param data dataset to save into txt file
#' @param filename name of the txt file
#' @param outdir the subfolder name to save the txt file
#' @param withdescription whether the dataset contains the description column, default set to TRUE
#'
#' @import tidyr
#' @importFrom readr write_tsv
#'
#' @export
#' @return NULL
#' @examples \dontrun{
#'  ms_filewrite(data, "file_to_mydata.txt")
#' }
#'
#'


ms_filewrite <- function(data, filename, outdir=NULL, withdescription=TRUE) {

  # add variable name to output
  dataname <- deparse(substitute(data))
  #print(outdir)
  if (length(grep("outdir", names(data)))) {
    outdir <- data$outdir[1]
    if (dir.exists(outdir)==FALSE) {
      dir.create(outdir)
    }
  } else if (!length(outdir)) {
    outdir <- paste0(dataname,"_",format(Sys.time(), "%y%m%d_%H%M"))
    dir.create(outdir)
  }

  if (!length(grep("description", names(data))) & withdescription) {
    #if(printout){print("ids in the original data are in composite format")}
    # list <- strsplit(as.character(data$id), "\n")
    # df <- ldply(list)
    # colnames(df) <- c("id", "description")
    # data <- cbind(df, data[ ,-1])
    data <- tidyr::separate(data, id, into=c("id", "description"), sep="\n")
  }
  data$outdir <- NULL

  if (length(outdir)) {
    readr::write_tsv(data, paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"), filename))
    #write.table(data, paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"), filename), sep="\t", row.names=FALSE, quote=FALSE);
  } else {
    readr::write_tsv(data, paste0(format(Sys.time(), "%y%m%d_%H%M_"), filename))
    #write.table(data, paste0(format(Sys.time(), "%y%m%d_%H%M_"), filename), sep="\t", row.names=FALSE, quote=FALSE);
  }
}

#' ms_fileread
#'
#' Function to read in dataset from tab delimited files.
#' The file should be within a subfolder under current working directory
#' The input such as in the format of "./subfolder/datafile.txt", where the first . indicated the current working directory
#'
#' @param file name of the txt file, in the format of "./subfolder/datafile.txt"
#'
#' @importFrom readr read_tsv locale
#'
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'  ms_fileread("./subfolder/datafile.txt")
#' }
#'
#'

ms_fileread <- function(file) {

  # file contains the subfolder information from working directory
  # such as "./subfolder/datafile.txt"
  # The first . indicated the current working directory
  outdir <- strsplit(file,"/")[[1]][2]
  data <- readr::read_tsv(file=file) #, sep="\t", header=T, quote="", na.string="", as.is=T, check.names=F );
  #for (i in 4:ncol(data)) { data[ ,i] <- as.numeric(data[ ,i]) }
  data$outdir <- outdir
  return(data)

}
