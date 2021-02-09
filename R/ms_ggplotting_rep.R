#' ms_ggplotting_rep
#'
#' Function to generate pdf files with multipanel ggplots for melt curve data,
#' particularly for dataset with replicated runs
#'
#' @param data CETSA melt curve dataset to plot
#' @param legenddata dataset used for condition extraction, at least one protein
#' inside this dataset should contains the full set of experiment conditions
#' @param levelvector a vector of conditions, preferably starting with Ctrl conditions
#' @param nread number of reading channels or sample treatements
#' @param remsinglecondprot whether orphan proteins to be plotted,
#' default value is TRUE, so to exclude them from plotting
#' @param withset whether there is set column to perform facet_grid, default set to FALSE
#' @param PSMcutoff whether to apply PSM threshold cutoff on hit selection,
#' default set to FALSE
#' @param PSMthreshold the threshold of averaged PSM numbers to consider protein
#' quantification as reliable, default value is 3
#' @param fitremout whether to segregate the proteins with messy melt curves, default set to FALSE
#' @param ctrlcond if necessary, could used to specify what conditions to be
#' referred as control conditions, for example "Vech"; by default "DMSO","Ctrl" and "Control"
#' have been automatically included as the keyword for control condition
#' @param bottomcutoff the average of the last three points should be lower than
#' specified bottom cutoff value, which is 0.4 by default
#' @param topcutoff the average of the first three points should be higher than
#' specified bottom cutoff value, which is 0.8 by default
#' @param variancecutoff whether to segregate the proteins with large inter-replicate variance
#' @param nMAD_var the number of MADs to set the significance cutoff about variance distribution,
#' default value is 2.5
#' @param simpleAUC whether to perform a simple calculation of AUC, default set to TRUE
#' @param topasone whether the top plateau has to be fixed, i.e., 1.0
#' @param normTop whether to normalize the AUC based on Top three readings
#' @param dotconnect whether to simply dot connect the readings for each curve
#' @param printcount whether to annotate the plots with PSM and uniPeptide number
#' @param annotypos the starting y-axis position of textual annotation, default value 0.5
#' @param annotyinterval the interval on y-axis for textual annotation, default valule 0.08
#' @param presetcolor whether to use the pre-defined color scheme, default set to TRUE
#' @param plotfitremout whether to plot out messy melt curves, default set to TURE
#' @param plotvarremout whether to plot out large inter-replicate variance melt curves, default set to TRUE
#' @param colorpanel a vector of customizable color scheme provided by the user
#' @param commonlegend whether to use one common legend for whole page of plots, default set to TRUE
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
#' @importFrom plyr . dlply ddply
#' @importFrom grid unit.c
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @import ggplot2
#'
#' @export
#'
#' @return a list of ggplot2 object
#' @examples \dontrun{
#' ms_ggplotting_rep(LY_scaled,
#'   levelvector=c("Ctrl.1", "Ctrl.2", "Treatment.1","Treatment.2"))
#' }
#'
#'

ms_ggplotting_rep <- function(data, legenddata=NULL, levelvector=NULL, nread=10, withset=FALSE,
                              remsinglecondprot=TRUE, PSMcutoff=FALSE, PSMthreshold=3,
                              fitremout=FALSE, ctrlcond=NULL, bottomcutoff=0.4, topcutoff=0.8,
                              variancecutoff=FALSE, nMAD_var=2.5, simpleAUC=TRUE,
                              topasone=TRUE, normTop=TRUE, dotconnect=FALSE,
                              pfdatabase=FALSE, printBothName=TRUE, printGeneName=FALSE,
                              printcount=TRUE, annotypos=0.5, annotyinterval=0.08,
                              presetcolor=TRUE, extraidtocomplete=NULL,
                              plotfitremout=TRUE, plotvarremout=TRUE,
                              colorpanel=NULL, commonlegend=TRUE, layout=c(5,5), external=TRUE,
                              toplabel="CETSA data plotting_curve fitting",
                              leftlabel="Non-denatured protein fraction",
                              bottomlabel="Temperature",
                              pdfname="ggplotting.pdf", pdfheight=12, pdfwidth=12, returnplots=FALSE) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  png(filename=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                      dataname, "_curve replicate distribution.png"))
  barplot(table(table(data$id)),
          main="Melting curve replicate distribution",
          xlab="Number of curves per unique protein in replicate exp")
  dev.off()

  checkpos <- c(1,3)
  if (withset) {
    setpos <- grep("^set", names(data))
    checkpos <- c(checkpos, setpos)
  }
  if (any(duplicated(data[ ,checkpos]))) {
    cat("Warning for duplicated protein name entries\n")
    cat("Double check the following proteins for duplicated entries:\n")
    cat(data[duplicated(data[ ,checkpos]), ]$id,"in",data[duplicated(data[ ,checkpos]), ]$condition,"\n")
    stop("1.Remove the duplicated entires from original dataset then start again!")
  }

  ncond <- length(unique(data$condition)) # unique condition names including replicate info
  if(length(levelvector)!=ncond | !setequal(levelvector, unique(data$condition))) {
    stop("Make sure you provide the correct number of conditions,
         first Controls then Treatments")
  }

  # to concatenate id and description
  nrowdata <- nrow(data)
  if (length(extraidtocomplete)) {
    fkeep <- NULL
    for (i in 1:length(extraidtocomplete)){
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
    if (printBothName & !pfdatabase) {
      data_extra <- data_extra %>% rowwise() %>% mutate(description1 = getProteinName(description, pfdatabase)) %>%
        mutate(description2 = getGeneName(description)) %>%
        mutate(id = paste(id, description1, description2, sep="\n"))
      data_extra$description1<-NULL
      data_extra$description2<-NULL
    } else if (printGeneName & !pfdatabase) {
      data_extra <- data_extra %>% rowwise() %>%
        mutate(description = getGeneName(description)) %>%
        mutate(id = paste(id, description, sep="\n"))
    } else {
      data_extra <- data_extra %>% rowwise() %>%
        mutate(description = getProteinName(description, pfdatabase)) %>%
        mutate(id = paste(id, description, sep="\n"))
    }
    data_extra$description<-NULL
  }

  checkpos <- c(1,2)
  if (withset) {
    setpos <- grep("^set", names(data))
    checkpos <- c(checkpos, setpos)
  }
  if (any(duplicated(data[ ,checkpos]))) {
    cat("Warning for duplicated protein name entries\n")
    cat("Double check the following proteins for duplicated entries:\n")
    cat(data[duplicated(data[ ,checkpos]), ]$id,"in",data[duplicated(data[ ,checkpos]), ]$condition,"\n")
    stop("Remove the duplicated entires from original dataset then start again!")
  }

  if (printcount & !withset) {
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
    printcount <- FALSE
  }

  # look for outlier proteins based on melting behavior in controls
  outliers <- NULL
  if (fitremout) {
    cat("Make sure you provide fitted data with Tm and R2 values for this option!\n")
    ctrllist1 <- unique(grep("[Cc][Tt][Rr][Ll]", data$condition, value=TRUE))
    ctrllist2 <- unique(grep("[Cc][Oo][Nn][Tt][Rr][Oo][Ll]", data$condition, value=TRUE))
    ctrllist3 <- unique(grep("[Dd][Mm][Ss][Oo]", data$condition, value=TRUE))
    if (length(ctrlcond)==1) { ctrlcond <- unique(grep(ctrlcond, data$condition, value=TRUE)) }
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
      data <- data[-outlierid1, ]
    }
  }

  # remove single condition proteins if desired
  if (remsinglecondprot) {
    counttable <- data %>% group_by(id) %>% summarize(count=n()) %>% filter(count>1)
    fkeep <- which(data$id %in% counttable$id)
    data <- data[fkeep, ]
  }
  if ( nrow(data)==0 ) {
    message("Make sure there are more than one experimental condition in dataset.")
    stop("Otherwise specify remsinglecondprot==FALSE !")
  }

  # To select the proteins with at least 3 PSM (average)
  if (PSMcutoff) {
    PSMcol <- grep("PSM", names(data), value=FALSE)
    PSMname <- names(data)[PSMcol]
    names(data)[PSMcol] <- "PSM"
    PSMkeep <- data %>% group_by(id) %>%
      summarize(PSMmean=mean(PSM,na.rm=T)) %>%
      filter(PSMmean>=PSMthreshold)
    fkeep <- which(data$id %in% PSMkeep$id)
    names(data)[PSMcol] <- PSMname
    data_PSMsmall <- data[-fkeep, ]
    data <- data[fkeep, ]
  }

  data$condition <- factor(data$condition, levels=levelvector)
  data1 <- tidyr::separate(data, condition, into=c("condition", "replicate"), sep="\\.")
  ncondition <- length(unique(data1$condition))
  nreplicate <- length(unique(data1$replicate))
  if (withset) {
    uniquecond <- unique(data1[ ,c("set", "condition", "replicate")])
  } else {
    uniquecond <- unique(data1[ ,c("condition", "replicate")])
  }
  row.names(uniquecond) <- NULL
  uniquecondlength <- nrow(uniquecond)
  cat("Replicates information were extracted as follows:\n")
  print(as.data.frame(uniquecond))
  nametempvector <- names(data1)[4:(nread+3)]
  numtempvector <- as.numeric(nametempvector)
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
  names(colorpanel) <- levelvector

  if (withset) {
    if (external & !returnplots) { external_graphs(T) }
    # subsetting complete replicates for legend generation:
    keep1 <- which(table(data$id) == uniquecondlength)
    if (length(keep1)==0) {
      #stop("No proteins contain complete replicate")
      keep1 <- which(table(data$id) > ncond)
      nkeep1 <- names(keep1)[1]
      fkeep1 <- which(data$id %in% nkeep1)
      data_complete <- data[fkeep1, ]
      data_complete <- data_complete[duplicated(data_complete$condition), ]
    } else {
      nkeep1 <- names(keep1)
      fkeep1 <- which(data$id %in% nkeep1)
      data_complete <- data[fkeep1, ]
    }

    plotlegend <- ms_melt_legend(data_complete, nread, colorpanel)

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
    ggsave(filename=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                       length(unique(data$id)),
                       "_whole_set_", pdfname), pl, height=pdfheight, width=pdfwidth)
    if (external) { external_graphs(F) } # switch off the external graphs
    message("Whole set plot file generated successfully.")
  } else {
    returnplots <- FALSE # don't allow plot to be returned when splitted...
    # subsetting complete replicates:
    keep1 <- which(table(data$id) == ncond)
    if(length(keep1)==0){ stop("No proteins contain complete replicate") }
    nkeep1 <- names(keep1)
    fkeep1 <- which(data$id %in% nkeep1)
    data_complete <- data[fkeep1, ]
    cat("The number of proteins with complete replicates is:",
                 length(unique(data_complete$id)), "\n")
    cat("The percentage of proteins with complete replicates is:",
                 round(length(unique(data_complete$id))/length(unique(data$id)),3)*100,"%\n")

    # Calculate euclidean distance metrics
    dism <- plyr::ddply(data_complete, "id", calEDscore, nread, nreplicate)
    png(filename=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                        dataname, "Euclidean distance score.png"))
    nf <- layout(mat = matrix(c(1,2),2,1, byrow=TRUE), height = c(3,1))
    par(mar=c(3.1, 3.1, 1.1, 2.1))
    hist(dism$EDscore, n=100, xlim=c(0,2), main="Euclidean distance score distribution")
    boxplot(dism$EDscore, horizontal=T, frame=F, width=0.75, ylim=c(0,2), col="black")
    dev.off()

    dism <- dism[order(dism$EDscore, decreasing=T), ]
    cat("The range of ED score is between",
        round(range(dism$EDscore)[1],3), "and", round(range(dism$EDscore)[2],3),"\n")
    ms_filewrite(dism, "Euclidean distance score table.txt", outdir=outdir)

    # Turkey boxplot [Q1-c*IQD, Q3+c*IQD]
    # pos10 <- which(dism$EDscore >= quantile(dism$EDscore, probs=0.9))
    # pos10id <- dism[pos10, ]
    # ms_filewrite(pos10id,"ED score Top10%_id.txt", outdir=outdir)

    plotlegend <- ms_melt_legend(data_complete, nread, colorpanel)

    data_largevar <- NULL
    if (variancecutoff) {
      cutoff <- NULL
      # MAD guided significance
      message("use MAD scheme to assign variance cutoff...\n")
      cutoff <- round(median(dism$variance)+nMAD_var*mad(dism$variance), 3)
      cat("The variance cutoff limit (", nMAD_var, "* mad ) sets at", cutoff, "\n")
      dism <- subset(dism, variance<cutoff)
      fkeep1 <- which(data_complete$id %in% dism$id)
      data_largevar <- data_complete[-fkeep1, ]
      data_complete <- data_complete[fkeep1, ]
      cat("The number of reproducible proteins with complete replicates is:",
          length(unique(data_complete$id)), "\n")
    }

    if (external & !returnplots) { external_graphs(T) }

    # print out the PSM small data file
    if (PSMcutoff) {
      pl <- ms_melt_innerplot(data_PSMsmall, nread, topasone, dotconnect,
                              printcount, PSM_annod, Pep_annod,
                              annotypos, annotyinterval,
                              colorpanel, plotlegend, commonlegend,
                              toplabel, leftlabel, bottomlabel,
                              withset, layout, returnplots, outdir)
      ggsave(filename=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                         length(unique(data_PSMsmall$id)),
                         "_PSMsmall proteins_", pdfname), pl, height=pdfheight, width=pdfwidth)
      message("PSMsmall plot file generated successfully.")
    }

    # print out the outliers data file
    if (fitremout & plotfitremout) {
      pl <- ms_melt_innerplot(outliers, nread, topasone, dotconnect,
                              printcount, PSM_annod, Pep_annod,
                              annotypos, annotyinterval,
                              colorpanel, plotlegend, commonlegend,
                              toplabel, leftlabel, bottomlabel,
                              withset, layout, returnplots, outdir)
      ggsave(filename=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                         length(unique(outliers$id)),
                         "_Messy_proteins_", pdfname), pl, height=pdfheight, width=pdfwidth)
      message("Messy plot file generated successfully.")
    }

    # print out the large variance data file
    if (variancecutoff & plotvarremout) {
      pl <- ms_melt_innerplot(data_largevar, nread, topasone, dotconnect,
                              printcount, PSM_annod, Pep_annod,
                              annotypos, annotyinterval,
                              colorpanel, plotlegend, commonlegend,
                              toplabel, leftlabel, bottomlabel,
                              withset, layout, returnplots, outdir)
      ggsave(filename=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                         length(unique(data_largevar$id)),
                         "_Non_reproducible_proteins_", pdfname), pl, height=pdfheight, width=pdfwidth)
      message("Non reproducible plot file generated successfully.")
    }

    message("Generating first complete plot file, pls wait...")
    #seq1 <- as.factor(dism$id)
    if (length(extraidtocomplete)) {
      data_complete <- rbind(data_extra, data_complete)
      data_complete <- data_complete[!duplicated(data_complete), ]
      data_complete$id <- factor(data_complete$id, levels=c(unique(data_extra$id), dism$id))
    } else {
      data_complete$id <- factor(data_complete$id, levels=dism$id)
    }

    pl <- ms_melt_innerplot(data_complete, nread, topasone, dotconnect,
                            printcount, PSM_annod, Pep_annod,
                            annotypos, annotyinterval,
                            colorpanel, plotlegend, commonlegend,
                            toplabel, leftlabel, bottomlabel,
                            withset, layout, returnplots, outdir)
    ggsave(filename=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                       length(unique(data_complete$id)),
                       "_Complete replicates_", pdfname), pl, height=pdfheight, width=pdfwidth)

    message("first complete plot file generated successfully.")

    # subsetting incomplete replicates:
    keep2 <- which(table(data$id) < ncond)
    if(length(keep2)==0){ stop("No proteins contain incomplete replicate") }
    nkeep2 <- names(keep2)
    fkeep2 <- which(data$id %in% nkeep2)
    data_incomp <- data[fkeep2, ]
    cat("The number of proteins without complete replicates is:",
        length(unique(data_incomp$id)), "\n")

    message("Generating second incomplete plot file, pls wait...")

    if (!("AUC" %in% colnames(data_incomp))) {
      if (simpleAUC==TRUE) {
        data_incomp$AUC <- rowMeans(data_incomp[ ,c(3:(nread+2))], na.rm=T)
      } else {
        interval <- max(numtempvector,na.rm=T)-min(numtempvector,na.rm=T)
        data_incomp$AUC <- apply(data_incomp, 1, function(x)
          MESS::auc(numtempvector, x[c(3:(nread+2))], type="linear")/interval)
      }
    }
    if (normTop) {
      data_incomp$Top <- rowMeans(data_incomp[ ,c(3:5)])
    } else {
      data_incomp$Top <- rep(1.0, nrow(data_incomp))
    }
    # This section is for un-fitted data
    data_incomp$Top <- rowMeans(data_incomp[ ,c(3:5)])
    data_incomp <- mutate(data_incomp, para=AUC/Top)
    delta <- function(dat) {
      max(dat$para)-min(dat$para) # simple ranking based on AUC of raw data
    }
    rank <- sapply(split(data_incomp[, c('id', 'para')], factor(data_incomp$id), drop=T), delta)
    ord <- order(rank, decreasing=T)
    plotseq <- as.factor(names(rank[ord]))
    data_incomp$id <- factor(data_incomp$id, levels=plotseq)

    pl <- ms_melt_innerplot(data_incomp, nread, topasone, dotconnect,
                            printcount, PSM_annod, Pep_annod,
                            annotypos, annotyinterval,
                            colorpanel, plotlegend, commonlegend,
                            toplabel, leftlabel, bottomlabel,
                            withset, layout, returnplots, outdir)
    ggsave(filename=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                       length(unique(data_incomp$id)),
                       "_Non_replicates_", pdfname), pl, height=pdfheight, width=pdfwidth)

    if (external) { external_graphs(F) } # switch off the external graphs
    message("second incomplete plot file generated successfully.")
  }
}




#' calEDscore
#'
#' Internal function to calculate euclidean distance metrics
#' for CETSA melt curve data
#'
#' @param data melt curve dataset (in complete replicates)
#' @param nread number of reading channels or sample treatments
#' @param nreplicate number of replicates
#'
#' @keywords internal
#'
#' @return a dataframe
#' @examples \dontrun{
#' }
#'
calEDscore <- function(data, nread, nreplicate) {
  # Calculate euclidean distance metrics
  data <- data[order(data$condition), ]
  dm <- as.matrix(stats::dist(data[,c(3:(nread+2))]))
  if (nreplicate==2) {
    Intra1<-dm[1,3]
    Intra2<-dm[2,4]
    InterCtrl<-dm[1,2]
    InterTreatment<-dm[3,4]
    variance<-InterCtrl+InterTreatment
    distance<-Intra1+Intra2-variance
    EDscore<-(Intra1+Intra2)/(10^(variance))
    eucl<-c(Intra1, Intra2, InterCtrl, InterTreatment, variance, distance, EDscore)
    dism <- data.frame(Intra1=eucl[1], Intra2=eucl[2], InterCtrl=eucl[3],
                       InterTreatment=eucl[4], variance=eucl[5],
                       distance=eucl[6], EDscore=eucl[7])
  } else if (nreplicate==3) {
      Intra1<-dm[1,4]
      Intra2<-dm[2,5]
      Intra3<-dm[3,6]
      InterCtrl1<-dm[1,2]
      InterCtrl2<-dm[1,3]
      InterCtrl3<-dm[2,3]
      InterTreatment1<-dm[4,5]
      InterTreatment2<-dm[5,6]
      InterTreatment3<-dm[4,6]
      variance<-(InterCtrl1+InterCtrl2+InterCtrl3+
                   InterTreatment1+InterTreatment2+InterTreatment3)/2
      distance<-Intra1+Intra2+Intra3-variance
      EDscore<-(Intra1+Intra2+Intra3)/(10^(variance))
      eucl<-c(Intra1, Intra2, Intra3, InterCtrl1, InterCtrl2, InterCtrl3,
              InterTreatment1, InterTreatment2, InterTreatment3, variance, distance, EDscore)
      dism <- data.frame(Intra1=eucl[1], Intra2=eucl[2], Intra3=eucl[3],
                 InterCtrl1=eucl[4], InterCtrl2=eucl[5], InterCtrl3=eucl[6],
                 InterTreatment1=eucl[7], InterTreatment2=eucl[8], InterTreatment3=eucl[9],
                 variance=eucl[10], distance=eucl[11], EDscore=eucl[12])
  } else if (nreplicate==4) {
      Intra1<-dm[1,5]
      Intra2<-dm[2,6]
      Intra3<-dm[3,7]
      Intra4<-dm[4,8]
      InterCtrl1<-dm[1,2]
      InterCtrl2<-dm[1,3]
      InterCtrl3<-dm[1,4]
      InterCtrl4<-dm[2,3]
      InterCtrl5<-dm[2,4]
      InterCtrl6<-dm[3,4]
      InterTreatment1<-dm[5,6]
      InterTreatment2<-dm[5,7]
      InterTreatment3<-dm[5,8]
      InterTreatment4<-dm[6,7]
      InterTreatment5<-dm[6,8]
      InterTreatment6<-dm[7,8]
      variance<-(InterCtrl1+InterCtrl2+InterCtrl3+
                   InterCtrl4+InterCtrl5+InterCtrl6+
                   InterTreatment1+InterTreatment2+InterTreatment3+
                   InterTreatment4+InterTreatment5+InterTreatment6)/3
      distance<-Intra1+Intra2+Intra3+Intra4-variance
      EDscore<-(Intra1+Intra2+Intra3+Intra4)/(10^(variance))
      eucl<-c(Intra1, Intra2, Intra3, Intra4,
              InterCtrl1, InterCtrl2, InterCtrl3,
              InterCtrl4, InterCtrl5, InterCtrl6,
              InterTreatment1, InterTreatment2, InterTreatment3,
              InterTreatment4, InterTreatment5, InterTreatment6, variance, distance, EDscore)
      dism <- data.frame(Intra1=eucl[1], Intra2=eucl[2], Intra3=eucl[3], Intra4=eucl[4],
                 InterCtrl1=eucl[5], InterCtrl2=eucl[6], InterCtrl3=eucl[7],
                 InterCtrl4=eucl[8], InterCtrl5=eucl[9], InterCtrl6=eucl[10],
                 InterTreatment1=eucl[11], InterTreatment2=eucl[12], InterTreatment3=eucl[13],
                 InterTreatment4=eucl[14], InterTreatment5=eucl[15], InterTreatment6=eucl[16],
                 variance=eucl[17], distance=eucl[18], EDscore=eucl[19])
  }
  return(dism)
}
