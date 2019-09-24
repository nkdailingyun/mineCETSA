#' Scaled dataset of NADPH treated K562 cell lysate in ITDR format.
#'
#' A dataset containing two technical replicates of the isothermal dose responses
#' of K562 lysate proteins when treated with NADPH compound.
#'
#' @format A data frame with 12,248 rows and 16 variables:
#' \describe{
#'   \item{id}{Protein Uniprot ID}
#'   \item{description}{Protein description}
#'   \item{condition}{Sample condition}
#'   \item{sumUniPeps}{Number of associated Unique Peptides}
#'   \item{sumPSMs}{Number of total Peptide Spectrum Matches}
#'   \item{countNum}{Number of Quantifying PSMs}
#'   ...
#' }
#' @source \url{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0208273}
"NADPH_scaled"
