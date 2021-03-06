% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ms_ggplotting.R
\name{ms_ggplotting}
\alias{ms_ggplotting}
\title{ms_ggplotting}
\usage{
ms_ggplotting(
  data,
  legenddata = NULL,
  nread = 10,
  remsinglecondprot = TRUE,
  fitremout = FALSE,
  ctrlcond = NULL,
  bottomcutoff = 0.4,
  topcutoff = 0.8,
  orderAUCdiff = TRUE,
  simpleAUC = TRUE,
  nreplicate = 1,
  topasone = TRUE,
  normTop = TRUE,
  dotconnect = FALSE,
  pfdatabase = FALSE,
  printBothName = TRUE,
  printGeneName = FALSE,
  printcount = TRUE,
  annotypos = 0.5,
  annotyinterval = 0.08,
  presetcolor = TRUE,
  colorpanel = NULL,
  extraidtocomplete = NULL,
  plotfitremout = TRUE,
  withset = FALSE,
  commonlegend = TRUE,
  layout = c(5, 5),
  external = TRUE,
  toplabel = "CETSA data plotting_curve fitting",
  leftlabel = "Non-denatured protein fraction",
  bottomlabel = "Temperature",
  returnplots = FALSE,
  pdfname = "ggplotting.pdf",
  pdfheight = 12,
  pdfwidth = 12
)
}
\arguments{
\item{data}{isothermal dataset to plot}

\item{legenddata}{dataset used for condition extraction, at least one protein
inside this dataset should contains the full set of experiment conditions}

\item{nread}{number of reading channels or sample treatements}

\item{remsinglecondprot}{whether orphan proteins to be plotted,
default value is TRUE, soto exclude them from plotting}

\item{fitremout}{whether to segregate the proteins with messy melt curves}

\item{ctrlcond}{if necessary, could used to specify what conditions to be
referred as control conditions, Ctrl/Control/DMSO is included as default keyword}

\item{bottomcutoff}{the average of the last three points should be lower than
specified bottom cutoff value, which is 0.4 by default}

\item{topcutoff}{the average of the first three points should be higher than
specified bottom cutoff value, which is 0.8 by default}

\item{orderAUCdiff}{whether to order plots by AUC difference among different
treatment for same protein, default set to TRUE}

\item{simpleAUC}{whether to perform a simple calculation of AUC, default set to TRUE}

\item{nreplicate}{number of replicates, default value is 1}

\item{topasone}{whether the top plateau has to be fixed, i.e., 1.0}

\item{normTop}{whether to normalize the AUC based on Top three readings}

\item{dotconnect}{whether to simply dot connect the readings for each curve}

\item{printcount}{whether to annotate the plots with PSM and uniPeptide number}

\item{annotypos}{the starting y-axis position of textual annotation, default value 0.5}

\item{annotyinterval}{the interval on y-axis for textual annotation, default valule 0.08}

\item{presetcolor}{whether to use the pre-defined color scheme}

\item{colorpanel}{a vector of customizable color scheme provided by the user}

\item{plotfitremout}{whether to plot out messy melt curves}

\item{commonlegend}{whether to use one common legend for whole page of plots}

\item{layout}{a vector indicating the panel layout for multi-panel plots per page}

\item{toplabel}{textual label at the top part of the page}

\item{leftlabel}{textual label at the left side of the page}

\item{bottomlabel}{textual label at the bottom part of the page}

\item{pdfname}{name for the pdf plots file}

\item{pdfheight}{a number indicate the height of pdf file, default value 12}

\item{pdfwidth}{a number indicate the width of pdf file, default value 12}
}
\value{
a list of ggplot2 object
}
\description{
Function to generate pdf files with multipanel ggplots for melt curve data
}
\examples{
\dontrun{
ms_ggplotting(LY_scaled)
}


}
