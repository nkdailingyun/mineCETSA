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
  if (!length(outdir)) {
    if (length(grep("outdir",names(data)))) {
      outdir <- data$outdir[1]
      data$outdir <- NULL
    } else {
      outdir <- attr(data,"outdir")
    }
  }

  # to prevent the change of sub-directory folder
  if (!length(outdir)) {
    outdir <- paste0(dataname,"_",format(Sys.time(), "%y%m%d_%H%M"))
    dir.create(outdir)
  } else if (dir.exists(outdir)==FALSE) {
    dir.create(outdir)
  }
  #print(outdir)

  if (!length(grep("description", names(data))) & withdescription) {
    data <- tidyr::separate(data, id, into=c("id", "description"), sep="\n")
  }

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
#' @param filename name of the txt file, in the format of "./subfolder/datafile.txt"
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

ms_fileread <- function(filename) {

  # file contains the subfolder information from working directory
  # such as "./subfolder/datafile.txt"
  # The first . indicated the current working directory
  outdir <- strsplit(filename,"/")[[1]][2]
  data <- readr::read_tsv(file=filename)
  if (length(attr(data,"outdir"))==0 & length(outdir)>0) {
    attr(data,"outdir") <- outdir
  }
  return(data)

}
