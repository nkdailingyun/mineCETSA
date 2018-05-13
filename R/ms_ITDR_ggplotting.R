#' ms_ITDR_ggplotting
#'
#' Function to generate pdf files with multipanel ggplots for ITDR data
#'
#' @param data ITDR dataset to plot
#' @param legenddata dataset used for condition extraction, at least one protein
#' inside this dataset should contains the full set of experiment conditions,
#' most of the time there is no need to supply, just use the data instead
#' @param levelvector a vector of experimental conditions, not complusory for isothermal functions
#' @param nread number of reading channels or sample treatements, default value 10
#' @param remsinglecondprot whether orphan proteins to be plotted,
#' default to exclude them from plotting
#' @param robustfitting whether to apply robust fitting, which is new since v.0.3.7
#' @param PSMcutoff whether to apply PSM threshold cutoff on hit selection,
#' default set to FALSE
#' @param PSMthreshold the threshold of averaged PSM numbers to consider protein
#' quantification as reliable, default value is 3
#' @param nreplicate number of replicates, default value is 2
#' @param minireplicate number of replicates to keep in final data, default set
#' to NULL, this is an added-on feature to subset enough replicated samples
#' especially for the format with error bar option as discussed below
#' @param barplotformat whether to plot in a bar graph format, default to FALSE
#' @param witherrorbar whether to plot in a mean +/- error bar(se) graph format, default to FALSE
#' @param usegradient whether to plot the bar plot in a color-gradient format,
#' which is only applicable to 1-sample format, default set to TRUE
#' @param orderEC whether to order the plots based on EC (Effective concentration)
#' @param orderAUC whether to order the plots based on AUC (Area Under the Curve)
#' @param orderRep whether to order the plots based on the measurement replicate numbers
#' @param plotseq a vector of plots arragement sequence (composite ID)
#' @param loess whether to perform curve fitting using loess model
#' @param dotconnect whether to simply dot connect the readings for each curve
#' @param printcount whether to annotate the plots with PSM and uniPeptide number
#' @param pfdatabase whether it is Plasmodium falciparum dataset, default to FALSE
#' @param printBothName whether to print both Protein and Gene names, default to TRUE
#' @param printGeneName whether to print only Gene names, default to FALSE, when
#' both printBothName and printBothName are FALSE, the protein name will be print out
#' @param unit textual annotation for the dose unit, default is "mM"
#' @param xlinear whether the x-axis should be in linear scale
#' @param xlog10 whether the x-axis should be in log10 scale
#' @param xsqrt whether the x-axis should be in square-root transformed scale
#' @param xcubert whether the x-axis should be in cube-root transformed scale
#' @param xinterval a number indicating the numerical interval for linear x-axis
#' @param fixedy whether the y-axis should use a fixed range scale
#' @param presetcolor whether to use the pre-defined color scheme
#' @param colorpanel a vector of customizable color scheme provided by the user
#' @param layout a vector indicating the panel layout for multi-panel plots
#' per page, default value is c(5,5)
#' @param toplabel textual label at the top of the page
#' @param bottomlabel textual label at the bottom of the page
#' @param leftlabel textual label at the left side of the page
#' @param shadearea a four element vector indicating the xmin, xmax, ymin, ymax for x and y axis for shading
#' @param shadeoutlinecolor color used to shade the area, default value is black
#' @param shadefillcolor color used to shade the area, default value is gray90
#' @param pdfheight a number indicate the height of pdf file, default value 12
#' @param pdfwidth a number indicate the width of pdf file, default value 12
#'
#'
#' @import tidyr
#' @import dplyr
#' @import RColorBrewer
#' @import ggplot2
#'
#' @export
#'
#' @return a list of ggplot2 object
#' @examples \dontrun{
#' ms_ITDR_ggplotting(ITDRdata_filtered[[1]], orderAUC=TRUE)
#' }
#'
#'
ms_ITDR_ggplotting <- function(data, legenddata=NULL, levelvector=NULL, nread=10,
                               remsinglecondprot=TRUE, robustfitting=FALSE,
                               PSMcutoff=FALSE, PSMthreshold=3,
                               nreplicate=2, minireplicate=NULL,
                               orderEC=FALSE, orderAUC=FALSE, orderRep=FALSE, plotseq=NULL,
                               barplotformat=FALSE, witherrorbar=FALSE, usegradient=TRUE,
                               loess=FALSE, dotconnect=FALSE, pfdatabase=FALSE,
                               printBothName=TRUE, printGeneName=FALSE, printcount=FALSE, unit="mM",
                               xlinear=FALSE, xlog10=TRUE, xsqrt=FALSE, xcubert=FALSE,
                               xinterval=NULL, fixedy=FALSE, layout=c(5,5),
                               presetcolor=TRUE, colorpanel=NULL,
                               toplabel="ITDR CETSA data plotting", bottomlabel=NULL,
                               leftlabel="Non-denatured protein fraction",
                               shadearea=NULL,shadeoutlinecolor="black",shadefillcolor="gray90",
                               pdfname="ITDR_ggplotting.pdf", external=TRUE,
                               pdfheight=12, pdfwidth=12) {

  # legenddata is any dataset containing the full levels of conditions, same as data
  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)

  # To select the proteins with more than 3 PSM (average)
  if (PSMcutoff) {
    PSMcol <- grep("PSM", names(data), value=FALSE)
    PSMname <- names(data)[PSMcol]
    names(data)[PSMcol] <- "PSM"
    PSMkeep <- data %>% group_by(id) %>%
      summarize(PSMmean=mean(PSM)) %>%
      filter(PSMmean>PSMthreshold)
    fkeep <- which(data$id %in% PSMkeep$id)
    names(data)[PSMcol] <- PSMname
    data_PSMsmall <- data[-fkeep, ]
    data <- data[fkeep, ]
  }
  #return(data)
  if (any(duplicated(data[ ,c(1,3)]))) {
    print("Warning for duplicated protein name entries")
    print("Double check the following proteins for duplicated entries:")
    print(paste0(data[duplicated(data[, c(1,3)]), ]$id," in ",
                 data[duplicated(data[, c(1,3)]), ]$condition))
    stop("Remove the duplicated entires from original dataset then start again!")
  }

  # remove single condition proteins if desired
  if (remsinglecondprot) {
    counttable <- data %>% group_by(id) %>% summarize(count=n()) %>% filter(count > 1)
    fkeep <- which(data$id %in% counttable$id)
    data <- data[fkeep, ]
  }

  # to concatenate id and description
  nrowdata <- nrow(data)
  if ( nrowdata==0 ) {
    print("Make sure there are more than one experimental condition in dataset.")
    stop("Otherwise specify remsinglecondprot==FALSE !")
  }

  if (printBothName) {
    data <- data %>% rowwise() %>% mutate(description1 = getProteinName(description, pfdatabase)) %>%
      mutate(description2 = getGeneName(description)) %>%
      mutate(id = paste(id, description1, description2, sep="\n"))
    data$description1<-NULL
    data$description2<-NULL
  } else if (printGeneName) {
    data <- data %>% rowwise() %>% mutate(description = getGeneName(description)) %>%
      mutate(id = paste(id, description, sep="\n"))
  } else {
    data <- data %>% rowwise() %>% mutate(description = getProteinName(description, pfdatabase)) %>%
      mutate(id = paste(id, description, sep="\n"))
  }
  data$description<-NULL

  if (any(duplicated(data[, c(1,2)]))) {
    print("Warning for duplicated protein name entries")
    print("Double check the following proteins for duplicated entries:")
    print(paste0(data[duplicated(data[, c(1,2)]), ]$id," in ",
                 data[duplicated(data[, c(1,2)]), ]$condition))
    stop("Remove the duplicated entires from original dataset then start again!")
  }

  if (printcount) {
    PSMcol <- grep("PSM", names(data), value=FALSE)
    PSM_annod <- data[ ,c(1,2,PSMcol)]
    PSM_annod <- spread(PSM_annod, condition, sumPSMs, drop=FALSE)
    #PSM_annod <- dcast(PSM_annod, id~condition, value.var="sumPSMs", drop=FALSE)
    Pepcol <- grep("Pep", names(data), value=FALSE)
    Pep_annod <- data[ ,c(1,2,Pepcol)]
    Pep_annod <- spread(Pep_annod, condition, sumUniPeps, drop=FALSE)
    #Pep_annod <- dcast(Pep_annod, id~condition, value.var="sumUniPeps", drop=FALSE)
    #return(list(PSM=PSM_annod, Pep=Pep_annod))
  } else {
    PSM_annod <- NULL
    Pep_annod <- NULL
  }

  if (orderEC) {
    plotseq <- unique(data[order(data$EC, decreasing=FALSE, na.last=TRUE), ]$id)
  }

  if (orderAUC) {
    data$AUC <- rowSums(data[ ,c(3:(nread+2))], na.rm=TRUE)
    plotseq <- names(sort(with(data, tapply(AUC, id, mean)), decreasing=TRUE))
    data$AUC <- NULL
    #plotseq <- names(sort(tapply(data_l$reading, data_l$id, mean), decreasing=TRUE))
  }

  if (orderRep) {
    plotseq <- dplyr::count(data, id, sort=T)$id
  }

  if (orderEC | orderAUC | orderRep | length(plotseq)) {
    data$id <- factor(data$id, levels=plotseq)
  } else {
    data$id <- factor(data$id)
  }

  if (length(legenddata)==0) { legenddata <- data }
  if (length(bottomlabel)==0) {
    bottomlabel <- paste0("Compound concentration(", unit, ")")
  }
  #plotlegend <- ms_legend(data, nread, colorpanel)

  if (external) { external_graphs(T) }

  if (barplotformat) {
    pl <- ms_isothermal_bar_innerplot(data, legenddata, levelvector, nread,
                                      minireplicate, witherrorbar, usegradient, fixedy,
                                      presetcolor, colorpanel, layout,
                                      toplabel, leftlabel, bottomlabel)
  } else {
    pl <- ms_isothermal_line_innerplot(data, legenddata, levelvector, nread, robustfitting,
                                  nreplicate, minireplicate, witherrorbar, loess,
                                  dotconnect, printcount, PSM_annod, Pep_annod,
                                  xlinear, xlog10, xsqrt, xcubert, xinterval,
                                  fixedy, presetcolor, colorpanel, layout,
                                  toplabel, leftlabel, bottomlabel,
                                  shadearea, shadeoutlinecolor, shadefillcolor)
  }

  if (length(outdir)) {
    ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"), pdfname),
           pl, height=pdfheight, width=pdfwidth)
  } else {
    ggsave(file=paste0(format(Sys.time(), "%y%m%d_%H%M"), pdfname), pl,
           height=pdfheight, width=pdfwidth)
  }

  if (external) { external_graphs(F) } # switch off the external graphs
  print("ITDR plot file generated successfully.")
}
