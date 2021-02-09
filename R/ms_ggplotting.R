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
#' @param fitremout whether to segregate the proteins with messy melt curves
#' @param ctrlcond if necessary, could used to specify what conditions to be
#' referred as control conditions, Ctrl/Control/DMSO is included as default keyword
#' @param bottomcutoff the average of the last three points should be lower than
#' specified bottom cutoff value, which is 0.4 by default
#' @param topcutoff the average of the first three points should be higher than
#' specified bottom cutoff value, which is 0.8 by default
#' @param orderAUCdiff whether to order plots by AUC difference among different
#' treatment for same protein, default set to TRUE
#' @param simpleAUC whether to perform a simple calculation of AUC, default set to TRUE
#' @param nreplicate number of replicates, default value is 1
#' @param topasone whether the top plateau has to be fixed, i.e., 1.0
#' @param normTop whether to normalize the AUC based on Top three readings
#' @param dotconnect whether to simply dot connect the readings for each curve
#' @param printcount whether to annotate the plots with PSM and uniPeptide number
#' @param annotypos the starting y-axis position of textual annotation, default value 0.5
#' @param annotyinterval the interval on y-axis for textual annotation, default valule 0.08
#' @param presetcolor whether to use the pre-defined color scheme
#' @param plotfitremout whether to plot out messy melt curves
#' @param colorpanel a vector of customizable color scheme provided by the user
#' @param commonlegend whether to use one common legend for whole page of plots
#' @param layout a vector indicating the panel layout for multi-panel plots per page
#' @param toplabel textual label at the top part of the page
#' @param leftlabel textual label at the left side of the page
#' @param bottomlabel textual label at the bottom part of the page
#' @param pdfname name for the pdf plots file
#' @param pdfheight a number indicate the height of pdf file, default value 12
#' @param pdfwidth a number indicate the width of pdf file, default value 12
#'
#'
#' @importFrom MESS auc
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
                          fitremout=FALSE, ctrlcond=NULL, bottomcutoff=0.4, topcutoff=0.8,
                          orderAUCdiff=TRUE, simpleAUC=TRUE, nreplicate=1, topasone=TRUE, normTop=TRUE,
                          dotconnect=FALSE, pfdatabase=FALSE, printBothName=TRUE, printGeneName=FALSE,
                          printcount=TRUE, annotypos=0.5, annotyinterval=0.08,
                          presetcolor=TRUE, colorpanel=NULL,
                          extraidtocomplete=NULL, plotfitremout=TRUE, withset=FALSE,
                          commonlegend=TRUE, layout=c(5,5), external=TRUE,
                          toplabel="CETSA data plotting_curve fitting",
                          leftlabel="Non-denatured protein fraction",
                          bottomlabel="Temperature", returnplots=FALSE,
                          pdfname="ggplotting.pdf", pdfheight=12, pdfwidth=12) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (any(duplicated(data[ ,c(1,3)]))) {
    cat("Warning for duplicated protein name entries\n")
    cat("Double check the following proteins for duplicated entries:\n")
    cat(data[duplicated(data[ ,c(1,3)]), ]$id,"in",data[duplicated(data[ ,c(1,3)]), ]$condition,"\n")
    stop("1.Remove the duplicated entires from original dataset then start again!")
  }

  # remove single condition proteins if desired
  if (remsinglecondprot) {
    counttable <- data %>% group_by(id) %>% summarize(count=n()) %>% filter(count > 1)
    fkeep <- which(data$id %in% counttable$id)
    data <- data[fkeep, ]
  }
  if (nrow(data)==0) {
    cat("Make sure there are more than one experimental condition in dataset.\n")
    stop("Otherwise specify remsinglecondprot==FALSE !")
  }

  nametempvector <- names(data)[4:(nread+3)]
  numtempvector <- as.numeric(nametempvector)
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
  names(colorpanel) <- sort(unique(data$condition))

  # to concatenate id and description
  nrowdata <- nrow(data)
  if (length(extraidtocomplete)) {
    fkeep <- NULL
    for (i in seq_along(extraidtocomplete)) {
      hits <- grep(paste0("^", extraidtocomplete[i], "$"), data$id, value=FALSE)
      fkeep <- c(fkeep, hits)
    }
    data_extra <- data[fkeep, ]
  }

  if (printBothName & !pfdatabase) {
    data <- data %>% rowwise() %>% mutate(description1 = getProteinName(description, pfdatabase)) %>%
      mutate(description2 = getGeneName(description)) %>%
      mutate(id = paste(id, description1, description2, sep="\n"))
    data$description1<-NULL
    data$description2<-NULL
  } else if (printGeneName & !pfdatabase) {
    data <- data %>% rowwise() %>%
      mutate(description = getGeneName(description)) %>%
      mutate(id = paste(id, description, sep="\n"))
  } else {
    data <- data %>% rowwise() %>%
      mutate(description = getProteinName(description, pfdatabase)) %>%
      mutate(id = paste(id, description, sep="\n"))
  }
  data$description<-NULL

  if (length(extraidtocomplete)) {
    if (printGeneName) {
      data_extra <- data_extra %>% rowwise() %>%
        mutate(description = getGeneName(description)) %>%
        mutate(id = paste(id, description, sep="\n"))
    } else {
      data_extra <- data_extra %>% rowwise() %>%
        mutate(description = getProteinName(description)) %>%
        mutate(id = paste(id, description, sep="\n"))
    }
    data_extra$description<-NULL
  }

  if(any(duplicated(data[ ,c(1,2)]))){
    cat("Warning for duplicated protein name entries\n")
    cat("Double check the following proteins for duplicated entries:\n")
    cat(data[duplicated(data[ ,c(1,2)]), ]$id,"in",data[duplicated(data[ ,c(1,2)]), ]$condition,"\n")
    stop("2.Remove the duplicated entires from original dataset then start again!")
  }

  if (printcount) {
    PSMcol <- grep("PSM", names(data), value=FALSE)
    PSM_annod <- data[ ,c(1,2,PSMcol)]
    PSM_annod <- tidyr::spread(PSM_annod, condition, sumPSMs, drop=FALSE)
    Pepcol <- grep("Pep", names(data), value=FALSE)
    Pep_annod <- data[ ,c(1,2,Pepcol)]
    Pep_annod <- tidyr::spread(Pep_annod, condition, sumUniPeps, drop=FALSE)
    #return(list(PSM=PSM_annod, Pep=Pep_annod))
  } else {
    PSM_annod <- NULL
    Pep_annod <- NULL
  }

  # look for outlier proteins based on melting behavior in controls
  outliers <- NULL
  if (fitremout) {
    cat("Make sure you provide fitted data with Tm and R2 values for this option!\n")
    ctrllist1 <- unique(grep("[Cc][Tt][Rr][Ll]", data$condition, value=TRUE))
    ctrllist2 <- unique(grep("[Cc][Oo][Nn][Tt][Rr][Oo][Ll]", data$condition, value=TRUE))
    ctrllist3 <- unique(grep("[Dd][Mm][Ss][Oo]", data$condition, value=TRUE))
    ctrllist <- c(ctrllist1, ctrllist2, ctrllist3, ctrlcond)
    if (length(ctrllist)==0) {
      stop("Name your control conditions with prefix [Ctrl] or [Control] or [DMSO], or specify in ctrlcond argument")
    }
    data$Top <- rowMeans(data[, c(3:5)])
    data$Bottom <- rowMeans(data[ ,c(nread:(nread+2))])

    tmp1 <- which(is.na(data$Tm))
    tmp2 <- which(data$Slope<=0)
    #tmp3 <- which(data$R2<=0.8)
    R2table <- plyr::ddply(data, .(id), summarize, meanR2=median(R2, na.rm=TRUE))
    R2table <- subset(R2table, meanR2<0.8)
    tmp3 <- which(data$id %in% R2table$id)
    #tmp4 <- which(data$condition %in% ctrllist & data$Tm>=100);
    tmp5 <- which(data$condition %in% ctrllist & data$Top<data$Bottom)
    tmp6 <- which(data$condition %in% ctrllist & data$Top<topcutoff)
    tmp7 <- which(data$condition %in% ctrllist & data$Bottom>bottomcutoff)
    tmp <- c(tmp1, tmp2, tmp3, tmp5, tmp6, tmp7)#, tmp4)
    tmp <- unique(tmp)
    data$Top <- NULL
    data$Bottom <- NULL
    if (length(tmp)>0) {
      cat(length(tmp), "measurements were messy in melting behavior and removed.\n")
      outlierid <- data[tmp, ]$id
      outlierid1 <- which(data$id %in% outlierid)
      outliers <- data[outlierid1, ]
      ms_filewrite(outliers, "Messy proteins.txt", outdir=outdir)
      #print(paste0("The outlier protein list is ", outlierids))
      data <- data[-outlierid1, ]
    }
  }

  if (orderAUCdiff) {
    # Ranking based by Area under the curve(AUC) difference
    if (!("^AUC" %in% colnames(data))) {
      if (simpleAUC==TRUE) {
        data$AUC <- rowMeans(data[ ,c(3:(nread+2))],na.rm=T)
        # data <- dplyr::mutate(data, AUC=AUC/nread)
      } else {
        data$AUC <- apply(data, 1, function(x) auc(numtempvector, x[c(3:(nread+2))], type="linear"))
      }
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
    plotseq <- names(rank[ord])
    data$AUC <- NULL
    data$Top <- NULL
    data$para <- NULL
    data$id <- factor(data$id, levels=plotseq)
    if (length(extraidtocomplete)) {
      data <- rbind(data_extra, data)
      data <- data[!duplicated(data), ]
      data$id <- factor(data$id, levels=unique(c(unique(data_extra$id), plotseq)))
    }
  } else {
    if (length(extraidtocomplete)) {
      data <- rbind(data_extra, data)
      data <- data[!duplicated(data), ]
      data$id <- factor(data$id, levels=unique(c(unique(data_extra$id), unique(data$id))))
    } else {
      data$id <- factor(data$id)
    }
  }

  if (!length(legenddata)) { legenddata <- data }
  plotlegend <- ms_melt_legend(legenddata, nread, colorpanel)

  if (external & !returnplots) { external_graphs(T) }

  pl <- ms_melt_innerplot(data, nread, topasone, dotconnect,
                          printcount, PSM_annod, Pep_annod,
                          annotypos, annotyinterval,
                          colorpanel, plotlegend, commonlegend,
                          toplabel, leftlabel, bottomlabel,
                          withset, layout, returnplots, outdir)
  if (returnplots) {
    if (external) { external_graphs(F) } # switch off the external graphs
    return(pl)
  }

  if (dotconnect) { pdfname=paste("simple",pdfname,sep="_") }
  else { pdfname=paste("fitted",pdfname,sep="_") }

  if (length(outdir)){
    ggsave(filename=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"), pdfname),
           pl, height=pdfheight, width=pdfwidth)
  } else {
    ggsave(filename=paste0(format(Sys.time(), "%y%m%d_%H%M_"), pdfname),
           pl, height=pdfheight, width=pdfwidth)
  }

  if (external) { external_graphs(F) } # switch off the external graphs
  cat("ggplotting plot file generated successfully.\n")
}
