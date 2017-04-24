#' ms_isothermal_scaling
#'
#' Internal function to apply systematic scaling to the ITDR or ITTR dataset
#'
#' @param data dataset to be scaled
#' @param nread number of reading channels or sample treatements
#' @param abdnorm whether to apply protein abundance level normalization
#' @param reftolowest whether to check and re-arrange the treatment dose
#' (or time) in ascending order, using the readings from lowest dose (or time)
#' group as the reference to derive ratios
#' @param remloadc whether to remove loading control sample
#' @param loadcname the header name of loading control sample
#' @param numcharmix whether the treatment names contains both character and
#' numeric values
#' @param writefactortofile whether to save a copy of scaling factors
#' @param bottomlabel textual label at the bottom of the plot
#' @param filename name for the file
#'
#' @keywords internal
#'
#' @importFrom tibble as_tibble
#' @importFrom gtools mixedorder mixedsort
#' @import tidyr
#' @import dplyr
#' @import ggplot2
#'
#' @return a dataframe
#' @examples \dontrun{
#' }
#'
#'

ms_isothermal_scaling <- function(data, nread, abdnorm, reftolowest, remloadc,
                                  loadcname, numcharmix, writefactortofile,
                                  bottomlabel, filename) {

  dataname <- deparse(substitute(data))
  outdir <- data$outdir[1]
  data$outdir <- NULL

  abdcol <- grep("Abundance", names(data), value=FALSE)
  if(length(abdcol) > 0 & abdnorm) {
    #if(!length(abdcol)){stop("Make sure there is abundance values, or add argument: abdnorm=FALSE")}
    abddata <- data[ ,c(1:3,abdcol)]
    ms_innerplotbox(abddata[ ,-2],
                    filename=paste0(dataname,"_Raw_Abundance_wholeset_trend.pdf"),
                    xlabel=bottomlabel, ylabel="Raw protein abundance (10^)",
                    isratio=FALSE, isothermalstyle=TRUE, outdir=outdir)

    mabddata <- ms_innerscale(abddata, nread, sumf=TRUE) # total abundance in each channel
    #return(mabddata)
    abdmean <- mean(unlist(mabddata[ ,c(4:(nread+3))]), na.rm=TRUE) # mean abundance
    abdNormfactor <- abdmean/mabddata[ ,c(4:(nread+3))]
    #return(abdNormfactor)
    # Apply correction factors to adjust the protein abundance in each channel
    uniqueCondition <- unique(abddata$condition)
    niter <- length(uniqueCondition)
    outdata <- abddata
    for (i in 1:niter) {
      pattern <- grep(pattern=paste0("^",uniqueCondition[i],"$"), abddata$condition, value=FALSE)
      tmpdata <- abddata[pattern, ]
      for (j in 1:nread) {
        tmpdata[ ,j+3] <- tmpdata[ ,j+3] * abdNormfactor[i,j]
        outdata[pattern, ] <- tmpdata
      }
    }
    ms_innerplotbox(outdata[ ,-2],
                    filename=paste0(dataname,"_Normalized_Abundance_wholeset_trend.pdf"),
                    xlabel=bottomlabel, ylabel="Normalized protein abundance (10^)",
                    isratio=FALSE, isothermalstyle=TRUE, outdir=outdir)
    ms_filewrite(outdata, paste0(dataname,"_","Scaled_abundance_data.txt"), outdir=outdir)
    #return(outdata)

    # Remove protein adundance data
    rdata <- data[ ,-abdcol]
    # calculate median for protein ratio data
    mdata1 <- ms_innerscale(rdata, nread, sumf=FALSE)
    mdata1$set <- "Pre-Abundance adjustment"
    # Apply correction factors to adjust the raw protein ratio in each channel
    outdata <- rdata
    for (i in 1:niter) {
      pattern <- grep(pattern=paste0("^",uniqueCondition[i],"$"), rdata$condition, value=FALSE)
      tmpdata <- rdata[pattern, ]
      for (j in 1:nread) {
        tmpdata[ ,j+3] <- tmpdata[ ,j+3] * abdNormfactor[i,j]
        outdata[pattern, ] <- tmpdata
      }
    }
    rdata <- outdata # adundance adjusted ratio generated.
    # calculate median for protein ratio data
    mdata2 <- ms_innerscale(rdata, nread, sumf=FALSE)
    mdata2$set <- "Post-Abundance adjustment"
    mdata <- rbind(mdata1, mdata2)
    setorder <- c("Pre-Abundance adjustment", "Post-Abundance adjustment")
  } else {
    abdcol <- grep("Abundance", names(data), value=FALSE)
    if (length(abdcol)) {
      rdata <- data[ ,-abdcol]
    } else {
      rdata <- data
    }
    # calculate median for protein ratio data
    mdata <- ms_innerscale(rdata, nread, sumf=FALSE)
    mdata$set <- "Raw Ratio"
    setorder <- c("Raw Ratio")
  }

  # make sure the dose is in ascending trend
  int_data <- rdata[ ,c(4:(nread+3))]
  if (numcharmix) {
    int_data <- int_data[ ,gtools::mixedorder(names(int_data))]
  } else {
    int_data <- int_data[ ,order(as.numeric(names(int_data)), decreasing=FALSE)]
  }
  rdata <- cbind(rdata[ ,c(1:3)], int_data, rdata[ ,c((nread+4):(ncol(rdata)))])
  namedosevector <- names(rdata[c(4:(nread+3))])
  #print(namedosevector)
  #print(mdata)

  # Normalize blank treatement to 1
  if (reftolowest) {
    int_data <- tibble::as_tibble(t(apply(int_data, 1, function(x) x/x[1])))
    names(int_data) <- namedosevector
    rdata <- cbind(rdata[ ,c(1:3)], int_data, rdata[ ,c((nread+4):(ncol(rdata)))])
  }

  # calculate median for protein ratio data
  mdata_new <- ms_innerscale(rdata, nread, sumf=FALSE)

  labels <- mdata_new$condition
  mrow <- nrow(mdata_new)

  # Calculation of Normalization factors (1/ overall median at each dose point).
  Normfactor <- 1/ mdata_new[ ,4:(nread+3)]

  # Apply correction factors (1/ overall median) to every value of each individual dataset
  uniqueCondition <- unique(rdata$condition)
  niter <- length(uniqueCondition)
  outdata <- rdata
  for (i in 1:niter) {
    pattern <- grep(pattern=paste0("^",uniqueCondition[i],"$"), rdata$condition, value=FALSE)
    tmpdata <- rdata[pattern, ]
    for (j in 1:nread) {
      tmpdata[ ,j+3] <- tmpdata[ ,j+3] * Normfactor[i,j]
      outdata[pattern, ] <- tmpdata
    }
  }
  ms_innerplotbox(outdata[ ,c(1,3,4:(nread+3))],
                  filename=paste0(dataname,"_Normalized_ratio_wholeset_trend.pdf"),
                  xlabel=bottomlabel, ylabel="Normalized Ratio", isratio=TRUE,
                  isothermalstyle=TRUE, outdir=outdir)
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

  # Generate median value plots post-scaling:
  scalemdata <- ms_innerscale(outdata, nread, sumf=FALSE)
  scalemdata$set <- "Post-normalization Ratio"
  setorder <- c(setorder, "Post-normalization Ratio")

  mdata <- tidyr::gather(mdata[ ,-c(1,2)], treatment, reading, -set, -condition)
  scalemdata <- tidyr::gather(scalemdata[ ,-c(1,2)], treatment, reading, -set, -condition)
  mediandata <- rbind(mdata, scalemdata)
  mediandata$set <- factor(mediandata$set, levels=setorder)
  if (numcharmix) {
    #print(unique(mediandata$treatment))
    #print(gtools::mixedsort(unique(mediandata$treatment)))
    mediandata$treatment <- factor(mediandata$treatment,
                                   levels=gtools::mixedsort(unique(mediandata$treatment)))
  } else {
    mediandata$treatment <- factor(mediandata$treatment,
                                   levels=sort(as.numeric(unique(mediandata$treatment)), decreasing=FALSE))
  }
  ms_innerplotmedian(mediandata, filename=paste0(dataname, "_scaling_median.pdf"),
                     xlabel=bottomlabel, ylabel="Median value of soluble proteins",
                     isothermalstyle=TRUE, outdir=outdir)

  if (remloadc) {
    # to remove the loading control
    loadcpos <- grep(pattern=paste0("^",loadcname,"$"), names(outdata), value=FALSE)
    int_data <- outdata[ ,-loadcpos][ ,c(4:(nread+2))]
    int_data <- tibble::as_tibble(t(apply(int_data, 1, function(x) x/x[1])))
    names(int_data) <- setdiff(namedosevector, loadcname)
    outdata <- cbind(outdata[ ,c(1:3)], int_data, outdata[,c((nread+4):(ncol(outdata)))])
  }

  outdata$outdir <- outdir
  ms_filewrite(outdata, paste0(dataname,"_","Scaled_data.txt"), outdir=outdir)
  return(tibble::as_tibble(outdata))
}
