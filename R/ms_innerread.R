#' ms_innerread
#'
#' Internal function to parse and read in data from tab delimited files exported from Proteome Discoverer
#'
#' @param file file name to import
#' @param fchoose whether to choose file interactively
#' @param treatment names of treatments (temperature, dose, time) applied to samples
#' @param nread number of reading channels, should match the number of channels used
#' @param abdread whether to read in protein abundance data
#' @param PD21 whether the data is searched using Proteome Discoverer 2.1
#' @param refchannel names of reference channel used in Proteome Discoverer search, such as 126
#' @param channels names of the read-in channels
#' @param ... other arguments ignored (for compatibility with generic)
#'
#' @keywords internal
#'
#' @importFrom readr read_tsv
#'
#' @return a dataframe
#' @examples \dontrun{
#'  ms_innerread("file.txt")
#' }
#'
#'
#'
ms_innerread <- function(file, fchoose, treatment, nread, abdread, PD21, refchannel, channels) {

  # Treatment and channel check
  if (length(treatment) != nread | length(channels) != nread) {
    stop("Make sure you specify the right number of channels!")
  }
  # File Check
  if (length(file) == 0) {
    stop("No valid file(s) loaded!")
  }

  # File Loading
  if (fchoose) {
    data <- readr::read_tsv(file.choose())
  } else {
    data <- readr::read_tsv(file=file)
  }

  # Clean up column names
  names(data) <- gsub(pattern="_", "", names(data))
  names(data) <- gsub(pattern=":", "", names(data))
  names(data) <- gsub(pattern="\\(", "", names(data))
  names(data) <- gsub(pattern="\\)", "", names(data))

  nrowdata <- nrow(data)
  colnames <- names(data)
  collength <- length(names(data))

  # remove standard error and variability

  pattern <- grep("Standard Error", names(data), value=FALSE)
  if (length(pattern) > 0) {
    data <- data[ ,-pattern]
    colnames <- names(data)
  }

  pattern <- grep("Variability", names(data), value=FALSE)
  if (length(pattern) > 0) {
    data <- data[ ,-pattern]
    colnames <- names(data)
  }

  # create condition vector
  #conditionchannel <- setdiff(channels, refchannel)[1]
  if (PD21) {
    conditions <- grep(pattern=paste0("^Abundance Ratio [A-z0-9,. -]+ / ", refchannel, "[A-z0-9,. -]+$"), colnames, value=TRUE)
    conditions <- gsub(pattern=paste0("^Abundance Ratio [A-z0-9,. -]+ / ", refchannel, ", "), "", conditions)
  } else {
    # create condition vector
    conditions <- grep(pattern=paste0("^[A-z0-9,. -]+ / ", refchannel, "[A-z0-9,. -]+$"), colnames, value=TRUE)
    conditions <- gsub(pattern=paste0("^[A-z0-9,. -]+ / ", refchannel,", "), "", conditions)
  }
  conditions <- unique(conditions)
  #print(conditions)

  # Move Accession
  tmppos <- grep(pattern="^Accession$", names(data), value=FALSE)
  collength <- length(names(data))
  data <- data[ ,c(tmppos,1:collength)]

  # Move Description
  tmppos <- grep(pattern="^Description$", names(data), value=FALSE)
  collength <- length(names(data))
  data <- data[ ,c(1,tmppos,2:collength)]

  # Condition row
  if (length(conditions) == 1) {
    data$condition <- rep(conditions, nrowdata)
    collength <- length(names(data))
    data <- data[ ,c(1:2,collength,3:(collength-1))]
  } else {
    stop("Make sure the condition was correctedly specified. Specify one unique condition for each input file.")
  }

  # Move reading values and rename them to correct position
  j = 1
  for (i in 1:nread) {
    if(PD21) {
      tmppos <- grep(pattern=paste0("^Abundance Ratio ",channels[i],", ",conditions," / ",refchannel,", ",conditions), names(data), value=FALSE)
    } else {
      tmppos <- grep(pattern=paste0(channels[i],", ",conditions," / ",refchannel,", ",conditions), names(data), value=FALSE)
    }
    if (length(tmppos)) {
      collength <- length(names(data))
      data <- data[ ,c(1:(j+2),tmppos,(j+3):collength)]
      j = j + 1
    }else if(channels[i] == refchannel) {
      # Set up start reference channel
      data$Reference <- rep(1.0, nrowdata)
      collength <- length(names(data))
      data <- data[ ,c(1:(j+2),(collength),(j+3):(collength-1))]
      j = j + 1
    }
  }

  # to read in protein abundance raw data
  if (PD21 & abdread) {
    for (i in 1:nread) {
      #get column for each channel and move it to correct position:
      tmppos <- grep(pattern=paste0("^Abundance F[0-9]+ ",channels[i],"[A-z0-9,. -]+$"), names(data), value=FALSE)
      collength <- length(names(data))
      data <- data[ ,c(1:(nread+2+i),tmppos,(nread+3+i):collength)]
    }
  }

  # Unique Peptides & PSMs
  tmppos <- grep(pattern="^# Unique Peptides$", names(data), value=FALSE)
  collength <- length(names(data))
  if (PD21 & abdread) {
    data <- data[,c(1:(2*nread+3),tmppos,(2*nread+4):collength)]
  } else {
    data <- data[,c(1:(nread+3),tmppos,(nread+4):collength)]
  }

  tmppos <- grep(pattern="^# PSMs$", names(data), value=FALSE)
  collength <- length(names(data))
  if (PD21 & abdread) {
    data <- data[ ,c(1:(2*nread+4),tmppos,(2*nread+5):collength)]
  } else {
    data <- data[ ,c(1:(nread+4),tmppos,(nread+5):collength)]
  }

  # Abundance counts
  tmppos <- grep(pattern=paste0("^Abundances Count [A-z0-9,. -]+",refchannel,"[A-z0-9,. -]+$"), names(data), value=FALSE)
  if (length(tmppos)) {
    collength <- length(names(data))
    if (PD21 & abdread) {
      data <- data[ ,c(1:(2*nread+5),tmppos)]
    } else {
      data <- data[ ,c(1:(nread+5),tmppos)]
    }
  } else {
    collength <- length(names(data))
    if (PD21 & abdread) {
      data <- data[ ,c(1:(2*nread+5))]
    } else {
      data <- data[ ,c(1:(nread+5))]
    }
  }

  # Rename the channel to correct treatment
  for (i in 4:(nread+3)) {
    names(data)[i] <- treatment[i-3]
  }

  if (PD21 & abdread) {
    for (i in (nread+4):(2*nread+3)) {
      names(data)[i] <- paste0("Abundance_",treatment[i-(nread+3)])
    }
  }

  # Sum unique peptides, psms & Accession name correction
  names(data) <- gsub(pattern="^# Unique Peptides$","sumUniPeps", names(data))
  names(data) <- gsub(pattern="^Description$","description", names(data))
  names(data) <- gsub(pattern="^# PSMs$","sumPSMs", names(data))
  names(data) <- gsub(pattern=paste0("^Abundances Count [A-z0-9,. -]+",refchannel,"[A-z0-9,. -]+$"),"countNum", names(data))
  names(data) <- gsub(pattern="^Accession$","id", names(data))

  # data$id <- as.character(data$id)
  # data$description <- as.character(data$description)

  return(data)
}
