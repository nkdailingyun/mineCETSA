#' getName
#'
#' Functions to extract Gene/Protein/Both names from the Uniprot description
#'
#' @param x Uniprot description input
#' @param pfdatabase whether the input is from P.falciparum database
#'
#' @keywords internal
#'
#' @return character
#' @examples \dontrun{
#' }


getGeneName <- function(x) {
  gene=strsplit(strsplit(x, "GN=")[[1]][2], " ")[[1]][1]
  if (length(gene)==0) {return(" ")} else {return(gene)}
}

getProteinName <- function(x, pfdatabase=FALSE) {
  if (pfdatabase) {
    # protein=gsub("product=", "", strsplit(x, "\\|")[[1]][2])
    protein=gsub("gene_product=", "", strsplit(x, "\\|")[[1]][4])
  } else {
    protein=strsplit(x, " OS=")[[1]][1]
  }
  if (length(protein)==0) {return(" ")} else {return(protein)}
}
