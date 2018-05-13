#' ms_2D_rawread
#'
#' Function to parse and read in CETSA data from tab delimited files exported from Proteome Discoverer
#'
#' @param filevector a file name or a vector of filenems to import
#' @param fchoose whether to choose file interactively, default set to FALSE
#' @param treatment a vector of treatment names applied to CETSA samples, in the same order as channels,
#' in a typical 2D setting, the treatment name should contain both replicate and treatment information,
#' preferably concatenated with an underline
#' @param nread number of reading channels, should match the number of channels used, default value 10
#' @param abdread whether to read in protein abundance data, default set to TRUE
#' @param PD21 whether the data is searched using Proteome Discoverer 2.1, default set to TRUE
#' @param refchannel names of reference channel used in Proteome Discoverer search, default value 126
#' @param channels names of the read-in channels, default value c("126","127N","127C","128N","128C","129N","129C","130N","130C","131")
#'
#' @seealso \code{\link{ms_ITTR_rawread}} for time response data
#' @seealso \code{\link{ms_IDTR_rawread}} for dose response data
#'
#' @import dplyr
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'  MOLM13 <- ms_2D_rawread(c("M13_37C_frac_Proteins.txt","M13_52C_frac_Proteins.txt","M13_58C_frac_Proteins.txt"), nread=10,
#'  treatment=("B1_DMSO","B1_TNFa","B1_AT26533","B2_DMSO","B2_TNFa","B2_AT26533","B3_DMSO","B3_TNFa","B3_AT26533","Mix"))
#' }
#'
#'
#'
ms_2D_rawread <- function(filevector, fchoose=FALSE, treatment=NULL, nread=10,
                          abdread=TRUE, PD21=TRUE, noratiodata=TRUE, refchannel="126",
                          channels=c("126","127N","127C","128N","128C","129N","129C","130N","130C","131")) {

  if (length(treatment)==0) {
    stop("Need to specify the treatment conditions in the same order as the TMT channel arrangement")
  }
  if (length(treatment)!=nread | length(treatment)!=length(channels)) {
    stop("The number of elements in treatment condition vector should be equal to the number of read-in channels")
  }
  if (length(filevector)==0) {
    stop("Need to specify the input data file names")
  }
  flength <- length(filevector)

  if (flength < 2) {
    dirname <- deparse(substitute(filevector))
    dirname_l <- unlist(strsplit(dirname, split="/"))
    dirname <- dirname_l[length(dirname_l)]
    data <- ms_innerread(filevector, fchoose, treatment, nread, abdread, PD21, refchannel, channels)
    data <- ms_dircreate(dirname, data)
    outdir <- attr(data,"outdir")
    if (noratiodata) {
      data <- data[ ,-c(4:(nread+3))]
      names(data) <- gsub("Abundance_", "", names(data))
    }
    if (length(attr(data,"outdir"))==0 & length(outdir)>0) {
      attr(data,"outdir") <- outdir
    }
    print("The data composition under each experimental condition (read in) is:")
    print(table(data$condition))
    return(data)
  } else {
    filename <- filevector[1]
    dirname <- deparse(substitute(filename))
    dirname_l <- unlist(strsplit(dirname, split="/"))
    dirname <- dirname_l[length(dirname_l)]
    indata <- ms_innerread(filevector[1], fchoose, treatment, nread, abdread, PD21, refchannel, channels)
    indata <- mutate(indata, condition = paste0(condition,".1"))
    outdata <- indata
    for (i in 2:flength) {
      indata <- ms_innerread(filevector[i], fchoose, treatment, nread, abdread, PD21, refchannel, channels)
      indata <- mutate(indata, condition = paste0(condition, ".", i))
      outdata <- rbind(x=outdata, y=indata, by=NULL)
    }
    outdata <- ms_dircreate(paste0("merged_",dirname), outdata)
    outdir <- attr(outdata,"outdir")
    if (noratiodata) {
      outdata <- outdata[ ,-c(4:(nread+3))]
      names(outdata) <- gsub("Abundance_", "", names(outdata))
    }
    if (length(attr(outdata,"outdir"))==0 & length(outdir)>0) {
      attr(outdata,"outdir") <- outdir
    }
    print("The data composition under each experimental condition (read in) is:")
    print(table(outdata$condition))
    return(outdata)
  }
}
