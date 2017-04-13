#' ms_ITDR_fitting
#'
#' Function to perform dose-response function fitting to generate fitting
#' parameters for each entry/curve
#'
#' @param data dataset to be fitted
#' @param nread number of reading channels or sample treatements,
#' default value 10
#' @param fc short for fold change, indicate the response level the fitting
#' function used to back-calculate the corresponding effective concentration
#' @param forcestart whether the base of dose-reponse curve start from 1.0,
#' default set to FALSE to allow free baseline
#' @param writetofile whether to keep a local file copy of the original data
#' with fitting parameters, default set to TRUE
#' @param keepfittedvalue whether to keep the fitted data at each dose points
#' with fitting parameters, default set to FALSE
#'
#' @import drc
#' @importFrom tibble as_tibble
#' @export
#' @return a dataframe
#' @examples \dontrun{
#' ITDRdata_fitted <- ms_ITDR_fitting(ITDRdata_scaled)
#' }
#'
#'


ms_ITDR_fitting <- function(data, nread=10, fc=0.3, forcestart=FALSE,
                            writetofile=TRUE, keepfittedvalue=FALSE) {

  dataname <- deparse(substitute(data))
  outdir <- data$outdir[1]
  data$outdir <- NULL

  nrowdata= nrow(data)
  nametempvector <- names(data[4:(nread+3)])
  # Calculate EC value, R2, Slope
  ECresult <- NULL
  R2result <- NULL
  Sloperesult <- NULL
  fitted_y <- matrix(nrow = nrowdata, ncol = nread)
  numtempvector <- as.numeric(nametempvector)
  for (i in 1:nrowdata) {
    y <- as.numeric(data[i,c(4:(nread+3))])
    x <- numtempvector
    if (forcestart & mean(y) > 1.0) {
      fit.dat <- try(drm(formula = y ~ x, fct = LL.4(fixed=c(NA,1,NA,NA))))
    } else if (forcestart & mean(y) <= 1.0){
      fit.dat <- try(drm(formula = y ~ x, fct = LL.4(fixed=c(NA,NA,1,NA))))
    } else {
      fit.dat <- try(drm(formula = y ~ x, fct = LL.4()))
    }
    #The forced fitting to 1, turned out not ideal,
    #maybe only appliable for "nice" curves

    if (class(fit.dat) != "try-error") {
      coeffs <- data.frame(coefficients(fit.dat))
      slope <- coeffs[1,1]
      fitted_y[i, ] <- fitted.values(fit.dat)
      if (mean(fitted.values(fit.dat)) > 1.0) {
        EC <- ED(fit.dat, respLev=1+fc, interval="delta", reference="control",
                 type="absolute", lref=1.0, display=FALSE)#coeffs[4,1]
      } else {
        EC <- ED(fit.dat, respLev=1/(1+fc), interval="delta", reference="upper",
                 type="absolute", uref=1.0, display=FALSE)#coeffs[4,1]
      }
      R2 <- 1 - sum((residuals(fit.dat)^2))/sum((y-mean(y))^2)
    } else {
      fit.dat <- try(lm(formula = y ~ x))
      if (class(fit.dat) != "try-error") {
        fitted_y[i, ] <- fitted.values(fit.dat)
        R2 <- summary(fit.dat)$r.squared
        slope <- fit.dat$coefficient[[2]]
        EC <- NA
      } else {
        EC = NA; R2 = NA; slope = NA;
      }
    }
    ECresult[i] = EC
    R2result[i] = R2
    Sloperesult[i] = slope
  }

  # Merge & Export fitting parameters file
  colnames(fitted_y) <- nametempvector
  Fitted <- as_tibble(fitted_y)
  Fitted <- cbind(data[, c(1:3)], Fitted)
  if (ncol(data) > (nread+3)) {
    #extra columns other than id, description, condition plus reading
    Fitted <- cbind(Fitted, data[ ,c((nread+4):(ncol(data)))])
  }
  #ECname <- gsub("\\.", "", past0(EC, as.character(fc)))
  Fitted["EC"] = ECresult
  Fitted["R2"] = R2result
  Fitted["Slope"] = Sloperesult

  data <- cbind(data, Fitted[ ,c((ncol(Fitted)-2):(ncol(Fitted)))])
  if (writetofile) {
    ms_filewrite(data, paste0(dataname, "_data plus fitting parameters.txt"),
                 outdir=outdir)
  }

  if (length(outdir) > 0) {
    data$outdir <- outdir
    Fitted$outdir <- outdir
  }
  if (keepfittedvalue){
    return(list(Rawdata=data, Fitteddata=Fitted))
  } else {
    return(data)
  }

}

#' ms_ITTR_fitting
#'
#' Function to perform time-response function fitting to generate fitting
#' parameters for each entry/curve
#'
#' @param data dataset to be fitted
#' @param nread number of reading channels or sample treatements,
#' default value 10
#' @param fc short for fold change, indicate the response level the fitting
#' function used to back-calculate the corresponding effective concentration
#' @param forcestart whether the base of time-reponse curve start from 1.0,
#' default set to FALSE to allow free baseline
#' @param writetofile whether to keep a local file copy of the original data
#' with fitting parameters, default set to TRUE
#' @param keepfittedvalue whether to keep the fitted data at each time points
#' with fitting parameters, default set to FALSE
#'
#' @import drc
#' @importFrom tibble as_tibble
#' @export
#' @return a dataframe
#' @examples \dontrun{
#' ITTRdata_fitted <- ms_ITTR_fitting(ITDRdata_scaled)
#' }
#'
#'

ms_ITTR_fitting <- function(data, nread=10, fc=0.3, forcestart=FALSE,
                            writetofile=TRUE, keepfittedvalue=FALSE) {

  dataname <- deparse(substitute(data))
  outdir <- data$outdir[1]
  data$outdir <- NULL

  nrowdata= nrow(data)
  nametempvector <- names(data[4:(nread+3)])
  # Calculate EC value, R2, Slope
  ETresult <- NULL
  R2result <- NULL
  Sloperesult <- NULL
  fitted_y <- matrix(nrow = nrowdata, ncol = nread)
  numtempvector <- as.numeric(nametempvector)
  for (i in 1:nrowdata) {
    y <- as.numeric(data[i,c(4:(nread+3))])
    x <- numtempvector
    if (forcestart & mean(y) > 1.0) {
      fit.dat <- try(drm(formula = y ~ x, fct = LL.4(fixed=c(NA,1,NA,NA))))
    } else if (forcestart & mean(y) <= 1.0){
      fit.dat <- try(drm(formula = y ~ x, fct = LL.4(fixed=c(NA,NA,1,NA))))
    } else {
      fit.dat <- try(drm(formula = y ~ x, fct = LL.4()))
    }
    #The forced fitting to 1, turned out not ideal,
    #maybe only appliable for "nice" curves

    if (class(fit.dat) != "try-error") {
      coeffs <- data.frame(coefficients(fit.dat))
      slope <- coeffs[1,1]
      fitted_y[i, ] <- fitted.values(fit.dat)
      if (mean(fitted.values(fit.dat)) > 1.0) {
        ET <- ED(fit.dat, respLev=1+fc, interval="delta", reference="control",
                 type="absolute", lref=1.0, display=FALSE)#coeffs[4,1]
      } else {
        ET <- ED(fit.dat, respLev=1/(1+fc), interval="delta", reference="upper",
                 type="absolute", uref=1.0, display=FALSE)#coeffs[4,1]
      }
      R2 <- 1 - sum((residuals(fit.dat)^2))/sum((y-mean(y))^2)
    } else {
      fit.dat <- try(lm(formula = y ~ x))
      if (class(fit.dat) != "try-error") {
        fitted_y[i, ] <- fitted.values(fit.dat)
        R2 <- summary(fit.dat)$r.squared
        slope <- fit.dat$coefficient[[2]]
        ET <- NA
      } else {
        ET = NA; R2 = NA; slope = NA;
      }
    }
    ETresult[i] = ET
    R2result[i] = R2
    Sloperesult[i] = slope
  }

  # Merge & Export fitting parameters file
  colnames(fitted_y) <- nametempvector
  Fitted <- as_tibble(fitted_y)
  Fitted <- cbind(data[, c(1:3)], Fitted)
  if (ncol(data) > (nread+3)) {
    #extra columns other than id, description, condition plus reading
    Fitted <- cbind(Fitted, data[,c((nread+4):(ncol(data)))])
  }
  #ECname <- gsub("\\.", "", past0(EC, as.character(fc)))
  Fitted["ET"] = ETresult
  Fitted["R2"] = R2result
  Fitted["Slope"] = Sloperesult

  data <- cbind(data, Fitted[ ,c((ncol(Fitted)-2):(ncol(Fitted)))])
  if (writetofile) {
    ms_filewrite(data, paste0(dataname, "_data plus fitting parameters.txt"),
                 outdir=outdir)
  }

  if (length(outdir) > 0) {
    data$outdir <- outdir
    Fitted$outdir <- outdir
  }
  if (keepfittedvalue){
    return(list(Rawdata=data, Fitteddata=Fitted))
  } else {
    return(data)
  }
}
