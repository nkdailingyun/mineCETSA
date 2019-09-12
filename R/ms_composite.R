#' ms_composite_ID_Gene
#'
#' Functions to make a composite Uniprot ID/Gene name from the Uniprot description
#'
#' @param data a dataframe with id and description column
#' @param pfdatabase whether the input is from P.falciparum database
#'
#' @export
#'
#' @return a dataframe with a compositie id column
#' @examples \dontrun{
#' }

ms_composite_ID_Gene <- function(data) {
  data <- data %>% rowwise() %>%
    mutate(description = getGeneName(description)) %>%
    mutate(id = paste(id, description, sep="\n"))
  data$description<-NULL
  return(data)
}

#' ms_composite_ID_Protein
#'
#' Functions to make a composite Uniprot ID/Protein name from the Uniprot description
#'
#' @param data a dataframe with id and description column
#' @param pfdatabase whether the input is from P.falciparum database
#'
#' @export
#'
#' @return a dataframe with a compositie id column
#' @examples \dontrun{
#' }

ms_composite_ID_Protein <- function(data, pfdatabase=FALSE) {
  data <- data %>% rowwise() %>%
    mutate(description = getProteinName(description, pfdatabase)) %>%
    mutate(id = paste(id, description, sep="\n"))
  data$description<-NULL
  return(data)
}

#' ms_composite_ID_Gene_Protein
#'
#' Functions to make a composite Uniprot ID/Gene name/Protein name from the Uniprot description
#'
#' @param data a dataframe with id and description column
#' @param pfdatabase whether the input is from P.falciparum database
#'
#' @export
#'
#' @return a dataframe with a compositie id column
#' @examples \dontrun{
#' }
#'

ms_composite_ID_Gene_Protein <- function(data, pfdatabase=FALSE) {
  data <- data %>% rowwise() %>% mutate(description1 = getProteinName(description, pfdatabase)) %>%
    mutate(description2 = getGeneName(description)) %>%
    mutate(id = paste(id, description1, description2, sep="\n"))
  data$description1<-NULL
  data$description2<-NULL
  data$description<-NULL
  return(data)
}

#' ms_getGene
#'
#' Functions to add one extra column with gene information from the Uniprot description
#'
#' @param data a dataframe with id and description column
#'
#' @export
#'
#' @return a dataframe with an extra gene column
#' @examples \dontrun{
#' }
#'
ms_getGene <- function(data) {
  data <- data %>% rowwise() %>%
    mutate(gene = getGeneName(description))
  return(data)
}
