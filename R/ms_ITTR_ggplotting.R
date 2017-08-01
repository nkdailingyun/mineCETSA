#' ms_ITTR_ggplotting
#'
#' Function to generate pdf files with multipanel ggplots for ITTR data
#'
#' @param data ITTR dataset to plot
#' @param legenddata dataset used for condition extraction, at least one protein
#' inside this dataset should contains the full set of experiment conditions,
#' most of the time there is no need to supply, just use the data instead
#' @param nread number of reading channels or sample treatements, default value 10
#' @param remsinglecondprot whether orphan proteins to be plotted,
#' default to exclude them from plotting
#' @param PSMcutoff whether to apply PSM threshold cutoff on hit selection,
#' default set to FALSE
#' @param PSMthreshold the threshold of averaged PSM numbers to consider protein
#' quantification as reliable, default value is 3
#' @param nreplicate number of replicates, default value is 2
#' @param orderET whether to order the plots based on ET (Effective concentration)
#' @param orderAUC whether to order the plots based on AUC (Area Under the Curve)
#' @param plotseq a vector of plots arragement sequence (composite ID)
#' @param loess whether to perform curve fitting using loess model
#' @param dotconnect whether to simply dot connect the readings for each curve
#' @param PSManno whether to annotate the plots with PSM and uniPeptide number
#' @param unit textual annotation for the dose unit, default is "min"
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
#' ms_ITTR_ggplotting(ITTRdata_filtered[[1]], orderAUC=TRUE)
#' }
#'
#'
ms_ITTR_ggplotting <- function(data, legenddata=NULL, nread=10, remsinglecondprot=TRUE,
                               PSMcutoff=FALSE, PSMthreshold=3, nreplicate=2,
                               remfragment=FALSE, remribosomal=FALSE,
                               orderET=FALSE, orderAUC=FALSE, plotseq=NULL,
                               loess=FALSE, dotconnect=FALSE, pfdatabase=FALSE,
                               printGeneName=FALSE, PSManno=FALSE, unit="min",
                               xlinear=FALSE, xlog10=FALSE, xsqrt=TRUE, xcubert=FALSE,
                               xinterval=NULL, fixedy=FALSE, layout=c(5,5),
                               presetcolor=TRUE, colorpanel=NULL,
                               top_label="ITTR CETSA data plotting",
                               pdfname="ITTR_ggplotting.pdf", external=TRUE){

  # legenddata is any dataset containing the full levels of conditions, same as data
  dataname <- deparse(substitute(data))
  outdir <- data$outdir[1]
  data$outdir <- NULL

  # to prevent the change of sub-directory folder
  if (!length(outdir)) {
    outdir <- paste0(dataname,"_",format(Sys.time(), "%y%m%d_%H%M"))
    dir.create(outdir)
  } else if (dir.exists(outdir)==FALSE) {
    dir.create(outdir)
  }

  if (remfragment) {
    if(length(grep("Fragment", data$description))) {
      data <- data[-grep("Fragment", data$description), ]
    }
  }

  if (remribosomal) {
    if (length(grep("ribosomal", data$description))) {
      data_ribo <- data[grep("ribosomal", data$description), ]
      data <- data[-grep("ribosomal", data$description), ]
    } else{
      data_ribo <- data[0, ]
    }
  }

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
  if (any(duplicated(data[, c(1,3)]))) {
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
  getGeneName <- function(x) {return (strsplit(strsplit(x, "GN=")[[1]][2], " ")[[1]][1])}
  getProteinName <- function(x) {return (strsplit(x, " OS=")[[1]][1])}
  if (pfdatabase) {
    getProteinName <- function(x) {return (gsub("product=", "", strsplit(x, "\\|")[[1]][2]))}
  }
  if (printGeneName) {
    data <- data %>% rowwise() %>% mutate(description = getGeneName(description)) %>%
      mutate(id = paste(id, description, sep="\n"))
  } else {
    data <- data %>% rowwise() %>% mutate(description = getProteinName(description)) %>%
      mutate(id = paste(id, description, sep="\n"))
  }

  # for(i in 1:nrowdata){
  #   if(geneName){
  #     data[i, 'description']<-strsplit(strsplit(data[i, 'description'], "GN=")[[1]][2], " ")[[1]][1]
  #     data[i,'id']<-paste(data[i, 'id'], data[i, 'description'], sep="\n")
  #   }else{
  #     cid<-paste(data[i, 'id'], data[i, 'description'], sep="\n")
  #     data[i,'id']<-strsplit(cid, " OS")[[1]][1]
  #   }
  # }
  data$description<-NULL

  if (any(duplicated(data[, c(1,2)]))) {
    print("Warning for duplicated protein name entries")
    print("Double check the following proteins for duplicated entries:")
    print(paste0(data[duplicated(data[, c(1,2)]), ]$id," in ",
                 data[duplicated(data[, c(1,2)]), ]$condition))
    stop("Remove the duplicated entires from original dataset then start again!")
  }

  if (PSManno) {
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

  if (orderET) {
    plotseq <- unique(data[order(data$ET, decreasing=FALSE, na.last=TRUE), ]$id)
  }

  if (orderAUC) {
    data$AUC <- rowSums(data[ ,c(3:(nread+2))])
    plotseq <- names(sort(with(data, tapply(AUC, id, mean)), decreasing=TRUE))
    data$AUC <- NULL
    #plotseq <- names(sort(tapply(data_l$reading, data_l$id, mean), decreasing=TRUE))
  }

  if (orderET | orderAUC | length(plotseq)) {
    data$id <- factor(data$id, levels=plotseq)
  } else {
    data$id <- factor(data$id)
  }

  if (!length(legenddata)) { legenddata <- data }
  bottom_label <- paste0("Treatment time(", unit, ")")
  #plotlegend <- ms_legend(data, nread, colorpanel)

  if (external) { external_graphs(T) }

  pl <- ms_isothermal_innerplot(data, legenddata, nread, nreplicate, loess,
                                dotconnect, PSManno, PSM_annod, Pep_annod,
                                xlinear, xlog10, xsqrt, xcubert, xinterval,
                                fixedy, presetcolor, colorpanel, layout,
                                top_label, bottom_label)

  if (length(outdir)) {
    ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"), pdfname),
           pl, height=12, width=12)
  } else {
    ggsave(file=paste0(format(Sys.time(), "%y%m%d_%H%M"), pdfname),
           pl, height=12, width=12)
  }

  if (external) { external_graphs(F) } # switch off the external graphs
  print("ITTR plot file generated successfully.")
}
