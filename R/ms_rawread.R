#' ms_rawread
#'
#' Function to parse and read in CETSA melt curve data from tab delimited files
#' exported from Proteome Discoverer
#'
#' @param filevector a file name or a vector of filenems to import
#' @param fchoose whether to choose file interactively, default set to FALSE
#' @param temp a vector of heating temperature applied to CETSA melt samples,
#' in the same order as channels
#' @param nread number of reading channels, should match the number of channels
#' used, default value 10
#' @param abdread whether to read in protein abundance data, default set to TRUE
#' @param PD21 whether the data is searched using Proteome Discoverer 2.1,
#' default set to TRUE
#' @param refchannel names of reference channel used in Proteome Discoverer
#' search, default value 126
#' @param channels names of the read-in channels, default value
#' c("126","127N","127C","128N","128C","129N","129C","130N","130C","131")
#'
#' @seealso \code{\link{ms_ITTR_rawread}} for time response data
#' @seealso \code{\link{ms_ITDR_rawread}} for dose response data
#'
#' @import dplyr
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'  LY <- ms_rawread(c("LY_K562_Ctrl_rep1_PD21_Proteins copy.txt",
#'  "LY_K562_Ctrl_rep2_PD21_Proteins copy.txt",
#'  "LY_K562_Treatment_rep1_PD21_Proteins copy.txt",
#'  "LY_K562_Treatment_rep2_PD21_Proteins copy.txt"))
#' }
#'
#'
#'
ms_rawread <- function(filevector, fchoose=FALSE, temp=c(37,40,43,46,49,52,55,58,61,64), nread=10, abdread=FALSE, PD21=TRUE,
                       refchannel="126", channels=c("126","127N","127C","128N","128C","129N","129C","130N","130C","131")) {

  flength <- length(filevector)
  if (flength < 2) {
    dirname <- deparse(substitute(filevector))
    data <- ms_innerread(filevector, fchoose, treatment=temp, nread, abdread, PD21, refchannel, channels)
    data <- ms_dircreate(dirname, data)
    print("The data composition under each experimental condition is:")
    print(table(data$condition))
    return(data)
  } else {
    filename <- filevector[1]
    dirname <- deparse(substitute(filename))
    indata <- ms_innerread(filevector[1], fchoose, treatment=temp, nread, abdread, PD21, refchannel, channels)
    indata <- mutate(indata, condition = paste0(condition,".1"))
    outdata <- indata
    for (i in 2:flength) {
      indata <- ms_innerread(filevector[i], fchoose, treatment=temp, nread, abdread, PD21, refchannel, channels)
      indata <- mutate(indata, condition = paste0(condition, ".", i))
      outdata <- rbind(x=outdata, y=indata, by=NULL)
    }
    outdata <- ms_dircreate(paste0("merged_",dirname), outdata)
    print("The data composition under each experimental condition (read in) is:")
    print(table(outdata$condition))
    return(outdata)
  }
}
