#' ms_ITDR_rawread
#'
#' Function to parse and read in ITDR data from tab delimited files exported from Proteome Discoverer
#'
#' @param filevector a file name or a vector of filenems to import
#' @param fchoose whether to choose file interactively, default set to FALSE
#' @param dose a vector of doses of compound applied to ITDR samples, in the same order as channels
#' @param nread number of reading channels, should match the number of channels used, default value 10
#' @param abdread whether to read in protein abundance data, default set to TRUE
#' @param PDversion which version of Proteome Discoverer the data is searched, possible values 20,21,22,24
#' @param refchannel names of reference channel used in Proteome Discoverer search, default value 126
#' @param channels names of the read-in channels, default value NULL, it would automatically
#' match the provided channel number when it is 10 or 11
#'
#' @seealso \code{\link{ms_ITTR_rawread}} for time response data
#'
#' @import dplyr
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'  ITDRdata <- ms_ITDR_rawread(c("SN1.txt", "SN2.txt", "ST1.txt", "ST2.txt"), PDversion=21,
#'  dose=c(20,5,1.25,0.3123,0.078125,0.019531,0.004883,0.001221,0.000305,0))
#' }
#'
#'
#'
ms_ITDR_rawread <- function(filevector, fchoose=FALSE,
                            dose=c(4,1,0.250,0.0625,0.015625,0.00391,0.00098,0.000245,0.00006,0),
                            nread=10, abdread=TRUE, PDversion=21, refchannel="126",
                            channels=NULL) {

  if (nread==10 & length(channels)==0) {
    channels=c("126","127N","127C","128N","128C","129N","129C","130N","130C","131")
  } else if (nread==11 & length(channels)==0) {
    channels=c("126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C")
  } else if (nread==16 & length(channels)==0) {
    channels=c("126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C","132N","132C","133N","133C","134N")
  } else if (nread!=length(channels) | nread!=length(dose)) {
    stop("Please provide a vector of used TMT channels")
  }

  flength <- length(filevector)
  if (flength < 2) {
    dirname <- deparse(substitute(filevector))
    dirname_l <- unlist(strsplit(dirname, split="/"))
    dirname <- dirname_l[length(dirname_l)]
    data <- ms_innerread(filevector, fchoose, treatment=dose, nread, abdread, PDversion, refchannel, channels)
    data <- ms_dircreate(dirname, data)
    cat("The data composition under each experimental condition (read in) is:\n")
    print(table(data$condition))
    return(data)
  } else {
    filename <- filevector[1]
    dirname <- deparse(substitute(filename))
    dirname_l <- unlist(strsplit(dirname, split="/"))
    dirname <- dirname_l[length(dirname_l)]
    indata <- ms_innerread(filevector[1], fchoose, treatment=dose, nread, abdread, PDversion, refchannel, channels)
    indata <- mutate(indata, condition = paste0(condition,".1"))
    outdata <- indata
    for (i in 2:flength) {
      indata <- ms_innerread(filevector[i], fchoose, treatment=dose, nread, abdread, PDversion, refchannel, channels)
      indata <- mutate(indata, condition = paste0(condition, ".", i))
      outdata <- rbind(x=outdata, y=indata, by=NULL)
    }
    outdata <- ms_dircreate(paste0("merged_",dirname), outdata)
    cat("The data composition under each experimental condition (read in) is:\n")
    print(table(outdata$condition))
    return(outdata)
  }
}
