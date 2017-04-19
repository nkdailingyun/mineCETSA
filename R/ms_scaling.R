#' ms_scaling
#'
#' Function to apply systematic scaling to CETSA melt curve dataset
#'
#' @param data dataset to be scaled
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
#' LY_scaled <- ms_scaling(LY_cleaned)
#' }
#'
#'
ms_scaling <- function(data, nread=10, reorder=FALSE, writefactortofile=TRUE,
                       bottomlabel="Temperature",
                       filename="CETSA_normalization_factors.txt") {

  dataname <- deparse(substitute(data))
  outdir <- data$outdir[1]
  data$outdir <- NULL

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

  # Calculate Tm value, R2 and Slope, Fit Median to 4 point Sigmodal Curve of the median ratios
  Tmresult <- NULL
  R2result <- NULL
  Sloperesult <- NULL

  # capture fitted values
  fitted_y <- matrix(nrow = mrow, ncol = nread)
  for(i in 1:mrow){
    y <- as.numeric(mdata[i,c(4:(nread+3))])
    x <- numtempvector
    fit.dat <-try(drm(formula = y ~ x, fct = LL.4()))

    if(class(fit.dat) != "try-error")
    {
      coeffs= data.frame(coefficients(fit.dat))
      slope = coeffs[1,1]
      Tm = coeffs[4,1]
      fitted_y[i,] = fitted.values(fit.dat)
      R2 = 1-sum((residuals(fit.dat)^2))/sum((y-mean(y))^2)
    }else{
      Tm=NA; R2=NA; slope=NA;
    }
    Tmresult[i]= Tm
    R2result[i]= R2
    Sloperesult[i]= slope
  }

  # Calculation of Fitting factor (Fitted median / Raw median at each temperature point).
  Fittingfactor <- (fitted_y / mdata[ ,c(4:(nread+3))])

  # Calculation of Scaling factor (1 / Fitted median at the lowest temperature point (i.e., 37oC))
  Scalingfactor <- (1 / fitted_y[ ,1])

  # Apply correction factors (1/ overall median) to every value of each individual dataset
  uniqueCondition <- unique(data$condition)
  niter <- length(uniqueCondition)
  outdata <- data
  for (i in 1:niter) {
    pattern <- grep(pattern=paste0("^",uniqueCondition[i],"$"), data$condition, value=FALSE)
    tmpdata <- data[pattern, ]
    for (j in 1:nread) {
      tmpdata[ ,(j+3)] <- tmpdata[ ,(j+3)] * Scalingfactor[i] * Fittingfactor[i,j]
      outdata[pattern, ] <- tmpdata
    }
  }
  ms_innerplotbox(outdata[ ,c(1,3,4:(nread+3))],
                  filename=paste0(dataname,"_Post_Normalization_wholeset_trend.pdf"),
                  xlabel=bottomlabel, ylabel="Post_Normalization_Ratio",
                  isratio=TRUE, isothermalstyle=FALSE, outdir=outdir)
  #return(outdata)
  # Print Normalization Factors into file:
  Fittingfactor$condition <- labels
  Fittingfactor <- Fittingfactor[ ,c((nread+1),1:nread)]

  # Generate median value plots pre- and post-scaling:
  scalemdata <- ms_innerscale(outdata, nread, sumf=FALSE)
  scalemdata$set <- "Post-normalization Ratio"

  mdata_l <- tidyr::gather(mdata[ ,-c(1,2)], treatment, reading, -set, -condition)
  scalemdata_l <- tidyr::gather(scalemdata[ ,-c(1,2)], treatment, reading, -set, -condition)
  mediandata <- rbind(mdata_l, scalemdata_l)
  mediandata$set <- factor(mediandata$set, levels=c("Pre-normalization Ratio", "Post-normalization Ratio"))
  mediandata$treatment <- factor(mediandata$treatment,
                                 levels=sort(as.numeric(unique(mediandata$treatment)), decreasing=FALSE))

  ms_innerplotmedian(mediandata, filename=paste0(dataname, "_scaling_median.pdf"),
                     xlabel=bottomlabel, ylabel="Median value of soluble proteins",
                     isothermalstyle=FALSE, outdir=outdir)

  if (writefactortofile) {
    filename <- paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),dataname,"_",filename)
    write("Fittingfactors: \n", filename)
    for (i in 1:mrow) {
      write(unlist(Fittingfactor[i, ]), filename, sep=" ", append=TRUE)
      write("\n",filename, append=TRUE)
    }
    write("\nScalingfactors: \n", filename, sep="\n", append=TRUE)
    write(Scalingfactor, filename, sep=" ", append=TRUE)
    write("\n", filename, append=TRUE)
    write("Pre-Normalization medians: \n", filename,append=TRUE)
    for (i in 1:mrow) {
      write(unlist(mdata[i,c(3,4:(nread+3))]), filename, sep=" ", append=TRUE)
      write("\n",filename, append=TRUE)
    }
    write("Post-Normalization medians: \n", filename, append=TRUE)
    for (i in 1:mrow) {
      write(unlist(scalemdata[i,c(3,4:(nread+3))]), filename, sep=" ", append=TRUE)
      write("\n",filename, append=TRUE)
    }
  }

  # Print Fitting and Scaling Factors:
  message("Fittingfactors: \n")
  print(Fittingfactor, row.names = FALSE)

  message("\nScalingfactors: \n")
  cat(Scalingfactor, "\n")

  outdata$outdir <- outdir
  ms_filewrite(outdata, paste0(dataname,"_","Scaled_data.txt"), outdir=outdir)
  return(as_tibble(outdata))
}
