#' ms_ggplotting
#'
#' Function to generate pdf files with multipanel ggplots for melt curve data
#'
#' @param data isothermal dataset to plot
#' @param legenddata dataset used for condition extraction, at least one protein
#' inside this dataset should contains the full set of experiment conditions
#' @param nread number of reading channels or sample treatements
#' @param remsinglecondprot whether orphan proteins to be plotted,
#' default value is TRUE, soto exclude them from plotting
#' @param orderAUCdiff whether to order plots by AUC difference among different
#' treatment for same protein, default set to TRUE
#' @param nreplicate number of replicates, default value is 1
#' @param topasone whether the top plateau has to be fixed, i.e., 1.0
#' @param normTop whether to normalize the AUC based on Top three readings
#' @param dotconnect whether to simply dot connect the readings for each curve
#' @param PSManno whether to annotate the plots with PSM and uniPeptide number
#' @param presetcolor whether to use the pre-defined color scheme
#' @param colorpanel a vector of customizable color scheme provided by the user
#' @param commonlegend whether to use one common legend for whole page of plots
#' @param layout a vector indicating the panel layout for multi-panel plots per page
#' @param pdfname name for the pdf plots file
#'
#'
#' @import tidyr
#' @import RColorBrewer
#' @import scales
#' @importFrom plyr . dlply
#' @importFrom grid unit.c
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @import ggplot2
#'
#' @export
#'
#' @return a list of ggplot2 object
#' @examples \dontrun{
#' ms_ggplotting(LY_scaled)
#' }
#'
#'

ms_ggplotting <- function(data, legenddata=NULL, nread=10, remsinglecondprot=TRUE,
                          orderAUCdiff=TRUE, nreplicate=1, topasone=TRUE, normTop=TRUE,
                          dotconnect=FALSE, pfdatabase=FALSE, printGeneName=FALSE,
                          PSManno=TRUE, presetcolor=TRUE, colorpanel=NULL, commonlegend=TRUE,
                          layout=c(5,5), pdfname="ggplotting.pdf") {

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

  if (any(duplicated(data[, c(1,3)]))) {
    print("Warning for duplicated protein name entries")
    print("Double check the following proteins for duplicated entries:")
    print(paste0(data[duplicated(data[, c(1,3)]), ]$id," in ",
                 data[duplicated(data[, c(1,3)]), ]$condition))
    stop("1.Remove the duplicated entires from original dataset then start again!")
  }

  # remove single condition proteins if desired
  if (remsinglecondprot) {
    counttable <- data %>% group_by(id) %>% summarize(count=n()) %>% filter(count > 1)
    fkeep <- which(data$id %in% counttable$id)
    data <- data[fkeep, ]
  }
  if ( nrow(data)==0 ) {
    print("Make sure there are more than one experimental condition in dataset.")
    stop("Otherwise specify remsinglecondprot==FALSE !")
  }

  nametempvector <- names(data[4:(nread+3)])
  ncond <- length(unique(data$condition))
  # if (ncond>8 & ncond<=12) { colorpanel <- brewer.pal(ncond, "Paired") }
  # if (ncond>12) { stop("12 colors are the current limit for plotting in single plot") }
  if (presetcolor & length(colorpanel)==0) {
    if (nreplicate==2 & ncond<=8) {colorpanel <- c("blue1", "blue4", "red1", "red4", "green1", "green4", "darkgrey", "black")}
    else if (nreplicate==2 & ncond<=12) {colorpanel <- brewer.pal(ncond, "Paired")}
    else if (nreplicate==1 & ncond<=8) {colorpanel <- c("black", "red", "blue", "orange", "green", "purple", "yellow", "cyan")}
    else if (nreplicate==1 & ncond<=12) {colorpanel <- brewer.pal(ncond, "Set3")}
    else if (nreplicate==3 & ncond<=12) {colorpanel <- c("blue1", "blue4", "blueviolet", "red1", "red3", "brown",
                                                   "green1", "green3", "forestgreen", "dimgray", "gray30", "black")}
    else if (nreplicate==4 & ncond<=8) {colorpanel <- c("#E0E0E0", "#BABABA", "#878787", "#4D4D4D", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")}
    else {stop("The number of conditions in dataset exceeds the preset number of colors, pls provide a vector of colors in colorpanel")}
  } else if (length(colorpanel) < ncond){
    stop("The number of conditions in dataset exceeds the provided number of colors, pls provide enough colors in colorpanel")
  }
  #colorpanel <- c("blue1", "blue4", "red1", "red4", "green1", "green4", "yellow", "yellow4", "gray30", "black")
  #colorpanel <- brewer_pal("qual")(length(levels(data_l$condition)))

  # to concatenate id and description
  nrowdata <- nrow(data)
  getGeneName <- function(x) {return (strsplit(strsplit(x, "GN=")[[1]][2], " ")[[1]][1])}
  getProteinName <- function(x) {return (strsplit(x, " OS=")[[1]][1])}
  if (pfdatabase) {
    getProteinName <- function(x) {return (gsub("product=", "", strsplit(x, "\\|")[[1]][2]))}
  }
  if (printGeneName) {
    data <- data %>% rowwise() %>%
      mutate(description = getGeneName(description)) %>%
      mutate(id = paste(id, description, sep="\n"))
  } else {
    data <- data %>% rowwise() %>%
      mutate(description = getProteinName(description)) %>%
      mutate(id = paste(id, description, sep="\n"))
  }
  data$description<-NULL

  if(any(duplicated(data[, c(1,2)]))){
    print("Warning for duplicated protein name entries")
    print("Double check the following proteins for duplicated entries:")
    print(paste0(data[duplicated(data[, c(1,2)]), ]$id," in ",
                 data[duplicated(data[, c(1,2)]), ]$condition))
    stop("2.Remove the duplicated entires from original dataset then start again!")
  }

  if (PSManno) {
    PSMcol <- grep("PSM", names(data), value=FALSE)
    PSM_annod <- data[ ,c(1,2,PSMcol)]
    PSM_annod <- tidyr::spread(PSM_annod, condition, sumPSMs, drop=FALSE)
    #PSM_annod <- dcast(PSM_annod, id~condition, value.var="sumPSMs", drop=FALSE)
    Pepcol <- grep("Pep", names(data), value=FALSE)
    Pep_annod <- data[ ,c(1,2,Pepcol)]
    Pep_annod <- tidyr::spread(Pep_annod, condition, sumUniPeps, drop=FALSE)
    #Pep_annod <- dcast(Pep_annod, id~condition, value.var="sumUniPeps", drop=FALSE)
    #return(list(PSM=PSM_annod, Pep=Pep_annod))
  } else {
    PSM_annod <- NULL
    Pep_annod <- NULL
  }
  if (orderAUCdiff) {
    # Ranking based by Area under the curve(AUC) difference
    if (!("^AUC" %in% colnames(data))) {
      data$AUC <- rowSums(data[ ,c(3:(nread+2))])
    }

    if (normTop) {
      data$Top <- rowMeans(data[ ,c(3:5)])
    } else {
      data$Top <- rep(1.0, nrowdata)
    }

    data <- mutate(data, para = AUC/Top)
    delta <- function(dat) {
      max(dat$para) - min(dat$para) # simple ranking based on AUC of raw data
    }
    rank = sapply(split(data[, c('id', 'para')], factor(data$id), drop=T), delta)
    ord <- order(rank, decreasing=T)
    plotseq <- as.factor(names(rank[ord]))
    data$id <- factor(data$id, levels=plotseq)
  } else {
    data$id <- factor(data$id)
  }

  if (!length(legenddata)) { legenddata <- data }
  plotlegend <- ms_melt_legend(legenddata, nread, colorpanel)

  pl <- ms_melt_innerplot(data, nread, topasone, dotconnect,
                          PSManno, PSM_annod, Pep_annod,
                          colorpanel, plotlegend, commonlegend, layout)

  if (dotconnect) { pdfname=paste("simple",pdfname,sep="_") }
  else { pdfname=paste("fitted",pdfname,sep="_") }

  if (length(outdir)){
    ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"), pdfname),
           pl, height=12, width=12)
  } else {
    ggsave(file=paste0(format(Sys.time(), "%y%m%d_%H%M_"), pdfname),
           pl, height=12, width=12)
  }

  print("ggplotting plot file generated successfully.")
}
