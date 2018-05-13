#' ms_ITTR_scaling
#'
#' Function to apply systematic scaling to the ITTR dataset
#'
#' @param data dataset to be scaled
#' @param nread number of reading channels or sample treatements, default value 10
#' @param abdnorm whether to apply protein abundance level normalization,
#' default set to TRUE
#' @param reftolowest whether to check and re-arrange the treatment dose,
#' (or time) in ascending order, using the readings from lowest dose (or time)
#' group as the reference to derive ratios, default set to TRUE
#' @param remloadc whether to remove loading control sample, default set to FALSE
#' @param loadcname the header name of loading control sample
#' @param numcharmix whether the treatment names contains both character and
#' numeric values
#' @param writefactortofile whether to save a copy of scaling factors,
#' default set to TRUE
#' @param bottomlabel textual label at the bottom of the plot
#' @param filename name for the file
#'
#' @export
#'
#' @return scaled dataset in dataframe format
#' @examples \dontrun{
#' ITTRdata_scaled <- ms_ITTR_scaling(ITTRdata_cleaned)
#' }
#'
#'
ms_ITTR_scaling <- function(data, nread=10, abdnorm=TRUE, reftolowest=TRUE,
                            remloadc=FALSE, loadcname="C", numcharmix=FALSE,
                            writefactortofile=TRUE, bottomlabel="Treatment time",
                            filename="ITTR_normalization_factors.txt") {

  dataname <- deparse(substitute(data))
  ms_isothermal_scaling(data, dataname, nread, abdnorm, reftolowest, remloadc, loadcname,
                        numcharmix,writefactortofile, bottomlabel, filename)
}
