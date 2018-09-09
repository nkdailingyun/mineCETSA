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
#' @param calMTT whether to calculate the Mininal Time Threshold(MTT), this is
#' useful when do follow-up analysis such as R2-AUC plot
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

ms_ITTR_fitting <- function(data, nread=10, fc=0.3, calMTT=FALSE,
                            nbaseline=3, baselineMAD=NULL, nMAD=2.5,
                            writetofile=TRUE, keepfittedvalue=FALSE) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  nrowdata <- nrow(data)
  nametempvector <- names(data[4:(nread+3)])

  if (calMTT) {
    # to extract the baseline variance information from the dataset
    if(!length(baselineMAD)) {
      baselineMAD <- round(mad(unlist(data[ ,c(4:(3+nbaseline))]), na.rm=T), 4)
    }
    if (baselineMAD==0) {
      stop("It is impossible to calculate MTT value when baseline variance
           is 0, should possibly change to a larger nbaseline value")
    } else {
      print(paste0("The baseline variance (based on the first ", nbaseline,
                   " points) is ", baselineMAD, "."))
    }
    forcestart = TRUE
    #The forced fitting to 1 as shown below, turned out not ideal,
    #maybe only appliable for "nice" curves
  } else {
    if(!length(baselineMAD)) { baselineMAD = 0 }
    forcestart = FALSE
  }

  # Calculate EC value, R2, Slope
  ETresult <- NULL
  R2result <- NULL
  Sloperesult <- NULL
  fitted_y <- matrix(nrow = nrowdata, ncol = nread)
  numtempvector <- as.numeric(nametempvector)
  for (i in 1:nrowdata) {
    y <- as.numeric(data[i,c(4:(nread+3))])
    x <- numtempvector
    valueindex = which(!is.na(y))
    if (forcestart & mean(y,na.rm=T) > 1.0) {
      fit.dat <- try(drm(formula = y ~ x, fct = LL.4(fixed=c(NA,1,NA,NA)), na.action=na.omit))
    } else if (forcestart & mean(y,na.rm=T) <= 1.0){
      fit.dat <- try(drm(formula = y ~ x, fct = LL.4(fixed=c(NA,NA,1,NA)), na.action=na.omit))
    } else {
      fit.dat <- try(drm(formula = y ~ x, fct = LL.4()))
    }

    if (class(fit.dat) != "try-error") {
      coeffs <- data.frame(coefficients(fit.dat))
      slope <- coeffs[1,1]
      fitted_y[i, valueindex] <- fitted.values(fit.dat)
      if (mean(fitted.values(fit.dat),na.rm=T) > 1.0) {
        if (calMTT) {
          ET <- ED(fit.dat, respLev=(1+nMAD*baselineMAD), interval="delta", reference="control",
                   type="absolute", lref=1.0, display=FALSE)#coeffs[4,1]
        } else {
          ET <- ED(fit.dat, respLev=(1+nMAD*baselineMAD)*(1+fc), interval="delta", reference="control",
                   type="absolute", lref=1.0, display=FALSE)#coeffs[4,1]
        }
      } else {
        if (calMTT) {
          ET <- ED(fit.dat, respLev=(1-nMAD*baselineMAD), interval="delta", reference="control",
                   type="absolute", lref=1.0, display=FALSE)#coeffs[4,1]
        } else {
          ET <- ED(fit.dat, respLev=(1-nMAD*baselineMAD)/(1+fc), interval="delta", reference="upper",
                   type="absolute", uref=1.0, display=FALSE)#coeffs[4,1]
        }
      }
      R2 <- 1 - sum((residuals(fit.dat)^2))/sum((y-mean(y))^2)
    } else {
      fit.dat <- try(lm(formula = y ~ x, na.action=na.omit))
      if (class(fit.dat) != "try-error") {
        fitted_y[i, valueindex] <- fitted.values(fit.dat)
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
    if (calMTT) {
      ms_filewrite(data, paste0(dataname, "_plus_MTT.txt"),
                   outdir=outdir)
    } else {
      ms_filewrite(data, paste0(dataname, "_plus_fitting_parameters.txt"),
                   outdir=outdir)
    }
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
