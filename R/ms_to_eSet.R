#' ms_to_eSet
#'
#' Function to reformat dataframe or matrix data into an expressionset (eSet) format
#'
#' @param data dataset to transform into an expressionset (eSet) format
#' @param matrixonly whether the input dataset is a matrix format, default set to FALSE
#' @param proteininfo a dataframe containing the correponding protein id-
#' protein description information, which is only required when the input data is a matrix
#' @param nread number of columns with protein abundance information
#' @param pdata a dataframe containing the experimental phenotype data information
#' @param refeset a reference eset
#' @param logtransform whether the numeric data input should be log2 transformed
#'
#' @keywords internal
#'
#' @import dplyr Biobase
#' @return a dataframe
#' @examples \dontrun{
#' }
#'
#'

ms_to_eSet <- function(data, matrixonly=FALSE, proteininfo=NULL, nread=NULL,
                       pdata=NULL, refeset=NULL, logtransform=FALSE) {

  # The typical input data should contains id, description and quanti channels
  if (length(nread)==0) { nread <- nrow(pdata) }
  if (length(refeset)==0 & length(pdata)==0) {
    stop("Meta data missing!")
  } else if (length(refeset)>0) {
    pdata <- Biobase::pData(refeset)
    proteininfo <- Biobase::fData(refeset)
    proteininfo <- tibble::rownames_to_column(proteininfo, "id")[ ,c("id","description")]
    #print(head(proteininfo))
    #print(length(unique(proteininfo$id)))
  }
  stopifnot(nread == nrow(pdata)) # make sure all the numeric column has proper phenotype labels

  if (matrixonly) {
    data <- data.frame(data)
    data <- tibble::rownames_to_column(data, var="id")
    #print(head(data))
    #print(length(unique(data$id)))
    #return(list(data=data, proteininfo=proteininfo))

    if (length(proteininfo)==0) {
      stop ("missing a support dataframe with matched id and protein information.")
    }
    if (length(setdiff(data$id, proteininfo$id))) {
      stop ("make sure the provided dataframe is complete.")
    }

    data <- merge(proteininfo, data, by="id")
    #print(head(data))
  }
  #return(data)

  matrix.please <- function(x) {
    m <- as.matrix(x[ ,-c(1:2)])
    rownames(m) <- x$id
    return(m)
  }

  exprs <- matrix.please(data[ ,c(1:(nread+2))])
  if (logtransform) {exprs = log2(exprs)}
  #print(class(exprs))
  #print(paste0("The dimension of input dataset is ", dim(exprs)[1], " X ", dim(exprs)[2]))
  #print(colnames(exprs))
  #print(head(exprs))
  #print(dim(pdata))
  #print(rownames(pdata))
  if (sum(grepl("^X[0-9]+C",colnames(exprs)))==nread) {
    colnames(exprs) <- gsub("^X","",colnames(exprs))
  }
  #print(summary(pdata))
  #print(sapply(pdata, class))
  stopifnot(all(rownames(pdata)==colnames(exprs)))
  metadata <- data.frame(labelDescription=names(pdata),
                         row.names=names(pdata))
  phenodata <- new("AnnotatedDataFrame", data=pdata, varMetadata=metadata)
  #print(phenodata)

  if (matrixonly==FALSE) {
    fdata <- data[ ,-c(1,3:(nread+2))]
    row.names(fdata) <- data$id
  } else {
    fdata <- fData(refeset)
  }

  metadata <- data.frame(labelDescription=names(fdata),
                         row.names=names(fdata))
  featuredata <- new("AnnotatedDataFrame", data=fdata, varMetadata=metadata)
  #print(featuredata)

  experimentData <- new("MIAME",
                        name="mineCETSA",
                        lab="PN Lab",
                        contact="PNlab@ntu.edu.sg",
                        title="CETSA Experiment",
                        abstract="ExpressionSet format",
                        other=list(notes=""))

  eSet <- ExpressionSet(assayData=exprs,
                        phenoData=phenodata,
                        featureData=featuredata,
                        experimentData=experimentData)

  return(eSet)
}

#' ms_eSet_to_df
#'
#' Function to reformat an expressionset (eSet) to a dataframe format
#'
#' @param expressionset expression dataset to transform into a dataframe
#' @param columnorder a vector for specifying the order of columns
#'
#' @keywords internal
#'
#' @import Biobase
#' @return a dataframe
#' @examples \dontrun{
#' }
#'
#'
ms_eSet_to_df <- function(expressionset, columnorder=NULL) {

  df <- cbind(Biobase::fData(expressionset), Biobase::exprs(expressionset))
  df$id <- row.names(expressionset)
  df$description <- as.character(df$description)
  row.names(df) <- NULL
  numberofcolumn <- ncol(df)
  if (length(columnorder)) {
    df <- df[ ,columnorder]
  } else {
    df <- df[ ,c(numberofcolumn,1,2:(numberofcolumn-1))]
  }
  return(df)
}
