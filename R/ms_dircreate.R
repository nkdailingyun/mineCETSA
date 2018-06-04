#' ms_dircreate
#'
#' Function to create a sub-working directory to store relevant results
#'
#' @param dirname a file name to be used for sub-working directory naming
#' @param data the dataset to be appended with the sub-working directory name
#'
#' @keywords internal
#' @return a dataframe
#'
#'

ms_dircreate <- function(dirname, data) {

  # create outfolder & variable
  dirname <- gsub(":", "_", dirname)
  dirname <- gsub(" ", "_", dirname)
  dirname <- gsub("\"", "", dirname)
  dirname <- gsub("\\.[A-z]*$", "", dirname)
  dirname <- gsub("\\.", "", dirname)
  dirname <- paste0(dirname,"_",format(Sys.time(), "%y%m%d_%H%M"))
  # print(dirname)

  if (length(attr(data,"outdir"))==0 & length(dirname)>0) {
    attr(data,"outdir") <- dirname
  }
  dir.create(dirname)
  return(data)
}


#' ms_directory
#'
#' Function to extract the directory information from data or assign a new
#' directory to data
#'
#' @param data the dataset to be appended with the sub-working directory name
#'
#' @keywords internal
#' @return a working directory
#'

ms_directory <- function(data, dataname) {

  if (length(grep("outdir",names(data)))) {
    outdir <- data$outdir[1]
    data$outdir <- NULL
  } else {
    outdir <- attr(data,"outdir")
  }

  if (!length(outdir)) {
    dirname <- gsub(":", "_", dataname)
    dirname <- gsub(" ", "_", dirname)
    dirname <- gsub("\"", "", dirname)
    dirname <- gsub("\\.[A-z]*$", "", dirname)
    dirname <- gsub("\\.", "", dirname)
    outdir <- paste0(dirname,"_",format(Sys.time(), "%y%m%d_%H%M"))
    dir.create(outdir)
  } else if (dir.exists(outdir)==FALSE) {
    dir.create(outdir)
  }
  return(list(outdir=outdir, data=data))
}
