#' ms_fitting
#'
#' Function to perform sigmoidal melt curve fitting to generate fitting
#' parameters for each entry/curve
#'
#' @param data dataset to be fitted
#' @param nread number of reading channels or sample treatements,
#' default value 10
#' @param topasone whether the top plateau has to be fixed, i.e., 1.0
#' @param halfmelt whether calculate Tm based on the half way between
#' top and bottom plateaus, default is FALSE
#' @param halforiginal whether calculate Tm based on the point with absolute
#' reading of 0.5, default is FALSE
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

ms_fitting <- function(data, nread=10, topasone=TRUE,
                       halfmelt=FALSE, halforiginal=FALSE,
                       writetofile=TRUE, keepfittedvalue=FALSE) {

  dataname <- deparse(substitute(data))
  outdir <- data$outdir[1]
  data$outdir <- NULL

  # to prevent the change of sub-directory folder
  if (!length(outdir)) {
    outdir <- paste0(dataname,"_",format(Sys.time(), "%y%m%d_%H%M"))
    dir.create(outdir)
  } else if (dir.exists(outdir)==FALSE) {
    dir.create(outdir)
  }

  nrowdata <- nrow(data)
  nametempvector <- names(data[4:(nread+3)])
  numtempvector <- as.numeric(nametempvector)
  # Calculate EC value, R2, Slope
  Tmresult <- NULL
  R2result <- NULL
  Sloperesult <- NULL
  RSEresult <- NULL
  fitted_y <- matrix(nrow = nrowdata, ncol = nread)

  if (topasone==TRUE) {
    ep <- 3
    top <- 1.0
  } else {
    ep <- 4
    top <- NA
  }
  for (i in 1:nrowdata) {
    y <- as.numeric(data[i,c(4:(nread+3))])
    x <- numtempvector
    fit.dat <- try(drc::drm(formula = y ~ x, fct = drc::LL.4(fixed=c(NA,NA,top,NA))))
    if (class(fit.dat) != "try-error") {
      coeffs <- data.frame(coefficients(fit.dat))
      slope <- coeffs[1,1]
      Tm <- coeffs[ep,1]
      if (halforiginal) {
        Tm <- ED(fit.dat, respLev=0.5, reference="control", type="absolute", display=FALSE)
      }else if (halfmelt & topasone) {
        Tm <- ED(fit.dat, respLev=0.5*(top+coeffs[2,1]), interval="delta", reference="upper", type="absolute", uref=top, display=FALSE)
      }else if (halfmelt) {
        Tm <- ED(fit.dat, respLev=0.5*(coeffs[3,1]+coeffs[2,1]), interval="delta", reference="upper", type="absolute", uref=coeffs[3,1], display=FALSE)
      }
      fitted_y[i, ] <- fitted.values(fit.dat)
      R2 <- 1 - sum((residuals(fit.dat)^2))/sum((y-mean(y))^2)
      RSE <- summary(fit.dat)$rseMat[[1]]
    } else {
      Tm = NA; R2 = NA; slope = NA; RSE=NA;
    }
    Tmresult[i] = Tm
    R2result[i] = R2
    Sloperesult[i] = slope
    RSEresult[i] = RSE
  }

  # Merge & Export fitting parameters file
  colnames(fitted_y) <- nametempvector
  Fitted <- tibble::as_tibble(fitted_y)
  Fitted <- cbind(data[, c(1:3)], Fitted)
  if (ncol(data) > (nread+3)) {
    #extra columns other than id, description, condition plus reading
    Fitted <- cbind(Fitted, data[ ,c((nread+4):(ncol(data)))])
  }
  Fitted["Tm"] = Tmresult
  Fitted["R2"] = R2result
  Fitted["Slope"] = Sloperesult
  Fitted["RSE"] = RSEresult

  data <- cbind(data, Fitted[ ,c((ncol(Fitted)-3):(ncol(Fitted)))])
  if (writetofile) {
    ms_filewrite(data, paste0(dataname, "_data plus fitting parameters.txt"),
                 outdir=outdir)
  }

  if (length(outdir) > 0) {
    data$outdir <- outdir
    Fitted$outdir <- outdir
  }
  if (keepfittedvalue) {
    return(list(Rawdata=data, Fitteddata=Fitted))
  } else {
    return(data)
  }
}
