#' ms_fitting
#'
#' Function to perform sigmoidal melt curve fitting to generate fitting
#' parameters for each entry/curve
#'
#' @param data dataset to be fitted
#' @param nread number of reading channels or sample treatements,
#' default value 10
#' @param topasone whether the top plateau has to be fixed, i.e., 1.0
#' @param halforiginal whether calculate Tm based on the point with absolute
#' reading of 0.5, default is TRUE; otherwise, calculate Tm based on the
#' half way between top and bottom plateaus, which is essentially
#' the inflection point
#' @param writetofile whether to keep a local file copy of the original data
#' with fitting parameters, default set to TRUE
#' @param keepfittedvalue whether to keep the fitted data at each temperature
#' points with fitting parameters, default set to FALSE
#'
#' @import drc
#' @importFrom tibble as_tibble
#' @export
#' @return a dataframe
#' @examples \dontrun{
#' LY_fitted <- ms_fitting(LY_scaled)
#' }
#'
#'

ms_fitting <- function(data, nread=10, topasone=TRUE, halforiginal=TRUE,
                       writetofile=TRUE, keepfittedvalue=FALSE) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  nrowdata <- nrow(data)
  nametempvector <- names(data)[4:(nread+3)]
  numtempvector <- as.numeric(nametempvector)
  # Calculate Tm, R2, Slope, and Bottom plateau
  Tmresult <- NULL
  R2result <- NULL
  Sloperesult <- NULL
  Plateauresult <- NULL
  fitted_y <- matrix(data=NA, nrow=nrowdata, ncol=nread)

  if (topasone==TRUE) {
    ep <- 3
    top <- 1.0
  } else {
    ep <- 4
    top <- NA
  }

  pb <- txtProgressBar(min=0, max=nrowdata, style=3, initial="")
  message("Curve fitting in progress...")
  # seems the control in drc doesnot work, so write error message to a file
  zz <- file(paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),"curve_fitting_record.txt"), open = "wt")
  sink(zz, type = "message")
  on.exit(sink(type="message"))
  for (i in seq_len(nrowdata)) {
    message(paste0("row ",i," protein ",data$id[i]))
    y <- as.numeric(data[i,c(4:(nread+3))])
    #x <- numtempvector
    valueindex = which(!is.na(y))
    fit.dat <- try(drc::drm(formula=y ~ numtempvector,
                            fct=drc::LL.4(fixed=c(NA,NA,top,NA)),
                            na.action=na.omit,
                            control=drmc(noMessage=TRUE)), silent=TRUE, outFile=zz)
    if (class(fit.dat) != "try-error") {
      coeffs <- data.frame(coefficients(fit.dat))
      slope <- coeffs[1,1]#b
      plateau <- coeffs[2,1]#c
      if (is.na(top)) { top <- coeffs[3,1] }
      if (halforiginal) {
        Tm <- drc::ED(fit.dat, respLev=0.5, reference="upper",
                      type="absolute", display=FALSE)[1,1]
      } else {
        # print(round(0.5*(top+plateau),2))
        Tm <- drc::ED(fit.dat, respLev=0.5*(top+plateau), reference="upper",
                      type="absolute", display=FALSE)[1,1]
        # Tm <- coeffs[ep,1]
      }
      fitted_y[i, valueindex] <- fitted.values(fit.dat)
      y1 <- na.omit(y)
      R2 <- 1 - sum((residuals(fit.dat)^2))/sum((y1-mean(y1))^2)
      # RSE <- summary(fit.dat)$rseMat[[1]]
    } else {
      Tm = NA; R2 = NA; slope = NA; plateau = NA; #RSE = NA;
      #message("The function failed to fit data into a typical melt curve...")
    }
    Tmresult[i] = Tm
    R2result[i] = R2
    Sloperesult[i] = -slope
    Plateauresult[i] = plateau
    # print(Tm)
    setTxtProgressBar(pb, i)
  }
  close(pb) # to close the progress bar

  # Merge & Export fitting parameters file
  colnames(fitted_y) <- nametempvector
  Fitted <- tibble::as_tibble(fitted_y)
  Fitted <- cbind(data[ ,c(1:3)], Fitted)
  if (ncol(data) > (nread+3)) {
    #extra columns other than id, description, condition plus reading
    Fitted <- cbind(Fitted, data[ ,c((nread+4):(ncol(data)))])
  }
  Fitted["Tm"] = Tmresult
  Fitted["R2"] = R2result
  Fitted["Slope"] = Sloperesult
  Fitted["Plateau"] = Plateauresult

  data <- cbind(data, Fitted[ ,c((ncol(Fitted)-3):(ncol(Fitted)))])
  if (writetofile) {
    ms_filewrite(data, paste0(dataname, "_data plus fitting parameters.txt"),
                 outdir=outdir)
  }

  if (length(attr(data,"outdir"))==0 & length(outdir)>0) {
    attr(data,"outdir") <- outdir
  }
  if (length(attr(Fitted,"outdir"))==0 & length(outdir)>0) {
    attr(Fitted,"outdir") <- outdir
  }
  if (keepfittedvalue) {
    return(list(Rawdata=data, Fitteddata=Fitted))
  } else {
    return(data)
  }
}
