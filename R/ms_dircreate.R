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
  data$outdir <- dirname
  dir.create(dirname)
  return(data)
}
