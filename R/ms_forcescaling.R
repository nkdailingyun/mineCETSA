#' ms_force_scaling
#'
#' Function to apply forced systematic scaling to CETSA melt curve dataset, i.e.,
#' to force scale the median of dataset to be same as the median of a subset of dataset (reference dataset)
#'
#' @param data dataset to be scaled
#' @param refcondition condition names to specify the reference dataset
#' @param nread number of reading channels or sample treatements
#' @param reorder whether to check and re-arrange the treatment temperature
#' in ascending order, using the readings from lowest temperature
#' group as the reference to derive ratios
#' @param writefactortofile whether to save a copy of scaling factors
#' @param bottomlabel textual label at the bottom of the plot
#' @param filename name for the file
#'
#'
#' @importFrom tibble as_tibble
#' @import tidyr
#' @import dplyr
#' @import ggplot2
#'
#' @export
#'
#' @return a dataframe
#' @examples \dontrun{
#' LY_forcescaled <- ms_forcescaling(LY_cleaned)
#' }
#'
#'
ms_forcescaling <- function(data, refcondition=NULL, nread=10, reorder=FALSE,
                            writefactortofile=TRUE, bottomlabel="Temperature",
                            filename="CETSA_normalization_factors.txt") {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (reorder) {
    # make sure the temperature is in ascending trend
    int_data <- data[ ,c(4:(nread+3))]
    int_data <- int_data[ ,order(as.numeric(names(int_data)), decreasing=FALSE)]
    nametempvector <- names(int_data)
    # Normalize blank treatement to 1
    int_data <- as_tibble(t(apply(int_data, 1, function(x) x/x[1])))
    names(int_data) <- nametempvector
    data <- cbind(data[ ,c(1:3)], int_data, data[ ,c((nread+4):(ncol(data)))])
  } else {
    nametempvector <- names(data[c(4:(nread+3))])
  }
  numtempvector <- as.numeric(nametempvector)

  ms_innerplotbox(data[ ,c(1,3,4:(nread+3))],
                  filename=paste0(dataname,"_Pre_Normalization_wholeset_trend.pdf"),
                  xlabel=bottomlabel, ylabel="Pre_Normalization_Ratio",
                  isratio=TRUE, isothermalstyle=FALSE, outdir=outdir)

  # calculate median for protein ratio data
  mdata <- ms_innerscale(data, nread, sumf=FALSE)
  mdata$set <- "Pre-normalization Ratio"
  labels <- mdata$condition
  mrow <- nrow(mdata)

  force_fitted_y <- c()

  if (length(refcondition)) {
    if (toupper(refcondition[1])=="ALL") {
      refdata <- data
    } else {
      ref = NULL
      for (refcond in refcondition) {
        pattern <- grep(pattern=paste0("^",refcond,"$"), data$condition)
        if (length(pattern)==0) { stop("Make sure you specify the right conditions from the dataset") }
        ref <- c(ref, pattern)
      }
      refdata <- data[ref, ]
    }
    #print(nrow(refdata))
  } else {
    stop("You need to specify the selected reference conditions or ALL by specifying refcondition argument")
  }

  for (j in 1:nread) {
    force_fitted_y[j] <- as.numeric(median(refdata[[j+3]], na.rm=TRUE))
  }
  #print(force_fitted_y)
  fit.dat <-try(drm(formula = force_fitted_y ~ numtempvector, fct = LL.4()))

  if(class(fit.dat) != "try-error") { force_fitted_y <- fitted.values(fit.dat) }
  #print(force_fitted_y)

  force_fitted_y <- force_fitted_y*(1/force_fitted_y[1]) # y-shift to top =1
  force_fitted <- matrix(nrow = mrow, ncol = nread)
  for (i in 1:mrow) { force_fitted[i, ] <- force_fitted_y }

  # Calculation of Normalization factor (Fitted median / Raw median at each temperature point).
  Normfactor <- (force_fitted / mdata[ ,c(4:(nread+3))])

  # Apply Normalization factor to every value of each individual dataset
  uniqueCondition <- unique(data$condition)
  niter <- length(uniqueCondition)
  outdata <- data
  for (i in 1:niter) {
    pattern <- grep(pattern=paste0("^",uniqueCondition[i],"$"), data$condition, value=FALSE)
    tmpdata <- data[pattern, ]
    for (j in 1:nread) {
      tmpdata[ ,j+3] <- tmpdata[ ,j+3] * Normfactor[i,j]
      outdata[pattern, ] <- tmpdata
    }
  }

  ms_innerplotbox(outdata[ ,c(1,3,4:(nread+3))],
                  filename=paste0(dataname,"_Post_Normalization_wholeset_trend.pdf"),
                  xlabel=bottomlabel, ylabel="Post_Normalization_Ratio",
                  isratio=TRUE, isothermalstyle=FALSE, outdir=outdir)

  #return(outdata)
  # Print Normalization Factors into file:
  Normfactor$condition <- labels
  Normfactor <- Normfactor[ ,c((nread+1),1:nread)]

  if (writefactortofile) {
    filename <- paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),dataname,"_",filename)
    write("Normalization factors: \n", filename)
    for (i in 1:mrow) {
      write(unlist(Normfactor[i, ]), filename, sep=" ", append=TRUE)
      write("\n",filename, append=TRUE)
    }
  }

  # Print Normalization Factors:
  message("Normalization factors: \n")
  print(Normfactor, row.names=FALSE)

  # Generate median value plots pre- and post-scaling:
  scalemdata <- ms_innerscale(outdata, nread, sumf=FALSE)
  scalemdata$set <- "Post-normalization Ratio"

  mdata_l <- tidyr::gather(mdata[ ,-c(1,2)], treatment, reading, -set, -condition)
  scalemdata_l <- tidyr::gather(scalemdata[ ,-c(1,2)], treatment, reading, -set, -condition)
  mediandata <- rbind(mdata_l, scalemdata_l)
  mediandata$set <- factor(mediandata$set, levels=c("Pre-normalization Ratio", "Post-normalization Ratio"))
  mediandata$treatment <- factor(mediandata$treatment,
                                 levels=sort(as.numeric(unique(mediandata$treatment)), decreasing=FALSE))

  ms_innerplotmedian(mediandata, filename=paste0(dataname, "_forcescaling_median.pdf"),
                     xlabel=bottomlabel, ylabel="Median value of soluble proteins",
                     isothermalstyle=FALSE, outdir=outdir)

  if (length(attr(outdata,"outdir"))==0 & length(outdir)>0) {
    attr(outdata,"outdir") <- outdir
  }
  ms_filewrite(outdata, paste0(dataname,"_","Forcescaled_data.txt"), outdir=outdir)
  return(as_tibble(outdata))
}
