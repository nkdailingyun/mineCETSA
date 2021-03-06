% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ms_subsetting.R
\name{ms_subsetting}
\alias{ms_subsetting}
\title{ms_subsetting}
\usage{
ms_subsetting(
  data,
  hitidlist,
  isfile = TRUE,
  pfdatabase = FALSE,
  allisoform = TRUE,
  revsel = FALSE
)
}
\arguments{
\item{data}{dataset to be subsetted}

\item{hitidlist}{a list of hit UniprotID to parse against the dataset, which
could be a vector of IDs or a txt file, see below}

\item{isfile}{whether the provided hitidlist is in a txt file under
column name "id", default set to TRUE}

\item{pfdatabase}{whether the data is a malaria dataset, default set to FALSE}

\item{allisoform}{whether to retrieve other isoforms of the
same parental Uniprot ID, default set to TRUE}

\item{revsel}{short for reverse selection, when set to TRUE, the output
is the dataset after removing the targets, default set to FALSE}
}
\value{
a dataframe
}
\description{
Function to parse data against a list to subset the data
}
\examples{
\dontrun{
IITDRdata_subset <- ms_subsetting(ITDRdata_f[[1]], hitidlist="hit_list.txt", isfile=TRUE)
ITDRdata_subset <- ms_subsetting(ITDRdata_f[[1]], hitidlist=c("P00000", "P12345-6"), isfile=FALSE)
}


}
