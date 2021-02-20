#' ms_IDIT_R2AUC
#'
#' Function to generate a R2 vs AUC plot for the ITDR or ITTR dataset
#' The R2 (x-axis value), AUC (y-axis) and MDT/MTT (size of dot) are from the input data
#'
#' @param data ITDR or ITTR dataset to generate R2 vs AUC plot, the provided
#' data should be associated with fitting parameters, i.e., after ms_ITDR_fitting
#' or ms_ITTR_fitting
#' @param nread number of reading channels or sample treatements, default value 10
#' @param printBothName whether to print both gene name and protein name, default set to FALSE
#' @param printGeneName whether to print only gene name, default set to TRUE, when
#' both printBothName and printGeneName are FALSE, only the protein name is shown
#' @param pfdatabase whether the data is a malaria dataset, default set to FALSE
#' @param rep the rep indicator used in the naming of experimental conditions,
#' such as "r", or "rep" or ""
#' @param onlyshowstabilized whether to show only stablized hits, ie,
#' direct targets, default set to FALSE
#' @param preaveraged whether the data has already been averaged from replicates,
#' default set to FALSE
#' @param simpleAUC whether to perform a simple calculation of AUC, default set to TRUE;
#' note it is a bit "tricky" to interpret the AUC value in logarithmic space for
#' most of ITDR data, suggest to just use simpleAUC by default
#' @param normalizedAUC whether to use the AUC values against a reference control
#' such as 37C readings, default set to FALSE
#' @param refkeyword a keyword used for the indication of reference control, as
#' most of the time the reference is 37C expression level, default value is 37C
#' @param log2scale whether to transform the readings into log2 scale, default set to FALSE
#' @param nMAD level of significance of AUC values, default set at 2.5
#' @param baselineMAD MAD of baseline variation, default value is 0
#' @param fcthreshold short for fold change threshold, similar parameter as in
#' ms_ITD/TR_filter() function, however this is an added-on selection feature
#' in this R2AUC plot, and applied to the averaged (and/or normalized) data,
#' so default value is 0
#' @param checkpointposition referring to the positions of dose points to check
#' whether their readings surpass fcthreshold, default value is NULL, which
#' would automatically check the highest 3 dose points
#' @param keepreplicate whether to only keep the curves that are measured with
#' all replicates under at least one experimental conditions, default set to
#' FALSE, when there is only one experimental condition in the dataset, this
#' means only keep the full replicates, ie, same as keepfullrep=TRUE
#' @param keepfullrep whether to only keep the curves that are measured with
#' all replicates under all experimental conditions, default set to FALSE, this
#' is useful for the dataset containing more than one experimental condition
#' @param stableref whether to only keep the curves with stable/flat reference
#' (most of the time, 37C) curves, default set to FALSE
#' @param stableref_nMAD level of significance for control of the reference
#' (most of the time, 37C) sample stability, default set at 3.5
#' @param colorbackground whether to color the background with color and shape
#' according to the condition, default set to TRUE
#' @param graybackground whether to use gray color for the background but
#' different shapes according to the condition, default set to FALSE
#' @param nodesizebyMDT whether to set node size according to MDT/MTT value,
#' ie the EC or ET column, default set to TRUE
#' @param PSMcutoff whether to apply PSM threshold cutoff on hit selection,
#' default set to TRUE
#' @param PSMthreshold the minimal threshold of averaged PSM numbers to consider
#' protein quantification as reliable, default value is 3
#' @param PSMcutoffbycondition whether to apply PSM threshold cutoff on each
#' experimental condition, otherwise on all the measured conditions, default set to FALSE
#' @param colornodes whether to fill the nodes in different colors if possible,
#' default set to TRUE
#' @param labelnodes whether to text-label the nodes, default set to TRUE
#' @param yscale a two-element vector to indicate the y-axis scale for plotting
#' if provided
#' @param idtoshow a vector containing the uniprot IDs to show on plot
#' regardless they are hits or not
#' @param idtoexclude a vector containing the uniprot IDs to exclude from plot
#' @param plottitle a character used to label the title of the plot
#' @param returnplot whether to retrieve the plot object, default set to FALSE
#'
#'
#' @importFrom MESS auc
#' @import tidyr
#' @importFrom plyr . dlply ddply
#' @import dplyr
#' @import RColorBrewer
#' @import ggplot2
#' @import ggrepel
#'
#' @export
#'
#' @return NULL
#' @examples \dontrun{
#' ms_IDIT_R2AUC(ITDRdata_fitted)
#' ms_IDIT_R2AUC(ITDRdata_fitted, normalizedAUC=TRUE, stableref=TRUE, yscale=c(0.5,2))
#' }
#'
#'


ms_IDIT_R2AUC <-  function(data, nread=10, printBothName=FALSE, printGeneName=TRUE,
                           pfdatabase=FALSE, rep="r", onlyshowstabilized=FALSE,
                           simpleAUC=TRUE, preaveraged=FALSE,
                           normalizedAUC=FALSE, refkeyword="37C",
                           AUC_uplimit=NULL, AUC_lowlimit=NULL, log2scale=FALSE,
                           nMAD=2.5, baselineMAD=0, fcthreshold=0, checkpointposition=NULL,
                           keepreplicate=FALSE, keepfullrep=FALSE,
                           stableref=FALSE, stableref_nMAD=3.5,
                           colorbackground=TRUE, graybackground=FALSE, nodesizebyMDT=TRUE,
                           colornodes=TRUE, labelnodes=TRUE, yscale=NULL,
                           PSMcutoff=TRUE, PSMthreshold=3, PSMcutoffbycondition=FALSE,
                           idtoshow=NULL, idtoexclude=NULL, plottitle=NULL, returnplot=FALSE) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
  }

  # #assign a very small R2 to negative R2 values
  # if (length(refkeyword)) {
  #   refcond <- unique(data$condition)[grep(refkeyword, unique(data$condition))]
  #   if (nrow(subset(data, condition %in% refcond & R2<=0))>0) {
  #     subset(data, condition %in% refcond & R2<=0)[ ,"R2"] <- 0.0099
  #   }
  #   data <- subset(data, R2>0)
  # }

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

  numdosevector <- as.numeric(names(data)[c(3:(nread+2))])

  if (sum(is.na(data$id))>0) { cat("There were", sum(is.na(data$id)),
                                   "entries not properly parsed...\n") }
  else { cat("Protein labels were successfully parsed.\n") }

  #to keep the data proper replicates as specified
  #to separate condition into condition and replicates
  nset <- length(unique(data$condition))
  data1 <- tidyr::separate(data, condition, into=c("condition", "replicate"), sep="\\.")
  ncondition <- length(unique(data1$condition))
  nreplicate <- length(unique(data1$replicate))
  uniquecond <- unique(data1[ ,c("condition", "replicate")])
  row.names(uniquecond) <- NULL
  cat("Replicates information were extracted as follows:\n")
  print(as.data.frame(uniquecond))
  conditionrep <- dplyr::count(uniquecond, condition)
  data_freq <- dplyr::count(data1, id, condition)
  if (keepfullrep) {
    id_keep1 <- unique(data_freq$id)
    for (i in 1:nrow(conditionrep)) {
      id_keep1 <- intersect(id_keep1, subset(data_freq, condition==conditionrep$condition[i] & n==conditionrep$n[i])$id)
    }
    data1 <- subset(data1, id %in% id_keep1)
    cat(nrow(data1), "measurements were measured with fully complete replicates in all conditions.\n")
  } else if (keepreplicate) {
    id_keep <- data_freq[0, ]
    for (i in 1: nrow(conditionrep)) {
      id_keep <- rbind(id_keep, subset(data_freq, condition==conditionrep$condition[i] & n==conditionrep$n[i]))
    }
    data1 <- merge(data1, id_keep, by=c("id","condition"), all=FALSE)
    data1$n <- NULL
    cat(nrow(data1), "measurements were measured with complete replicates in at least one condition.\n")
  }
  data <- data1

  # To select the proteins with more than 3 PSM (average)
  if (PSMcutoff) {
    PSMcol <- grep("PSM", names(data), value=F)
    PSMname <- names(data)[PSMcol]
    names(data)[PSMcol] <- "PSM"
    if (PSMcutoffbycondition) {
      PSMtable <- plyr::ddply(data, plyr::.(id,condition), summarize, PSMmean=mean(PSM))
    } else {
      PSMtable <- plyr::ddply(data, plyr::.(id), summarize, PSMmean=mean(PSM))
    }
    PSMtable <- subset(PSMtable, PSMmean>=PSMthreshold)
    names(data)[PSMcol] <- PSMname
    data <- merge(data, PSMtable)
    data$PSMmean <- NULL
    cat("PSM cutoff was applied, at least",as.character(PSMthreshold),"PSMs were required.\n")
    cat(length(unique(data$id)), "proteins remain after checking on PSM.\n")
  }

  if (length(idtoexclude)) {
    fkeep <- NULL
    for (i in 1:length(idtoexclude)) {
      hits <- grep(paste0("^", idtoexclude[i]), data$id, value=FALSE)
      fkeep <- c(fkeep, hits)
    }
    if (length(fkeep)>0) { data <- data[-fkeep, ] }
  }

  if (log2scale) {
    # print(head(data))
    data.m <- data.frame(log2(as.matrix(data[ ,c(4:(nread+3))])))
    names(data.m) <- names(data)[c(4:(nread+3))]
    data <- cbind(data[ ,c(1:3)], data.m, data[ ,c((nread+4):ncol(data))])
  }
  data$STD <- apply(data[ ,c(4:(nread+3))], 1, sd, na.rm=T)
  data$Mean <- apply(data[ ,c(4:(nread+3))], 1, mean, na.rm=T)
  data <- dplyr::mutate(data, CV=STD/Mean)
  if (simpleAUC==TRUE) {
    data$AUC <- rowMeans(data[ ,c(4:(nread+3))],na.rm=T)
    # data <- dplyr::mutate(data, AUC=rowSums(data[ ,c(4:(nread+3))],na.rm=T))
    # data <- dplyr::mutate(data, AUC=AUC/nread)
  } else {
    data$AUC <- apply(data, 1, function(x) auc(numdosevector, x[c(4:(nread+3))], type="linear"))
  }
  data$STD <- NULL
  data$Mean <- NULL
  # if (length(attr(data,"outdir"))==0 & length(outdir)>0) {
  #   attr(data,"outdir") <- outdir
  # }
  # ms_filewrite(data, paste0(dataname, "_data_selected_plus_AUC.txt"), outdir=outdir)
  colorder <- names(data)
  colorder1 <- setdiff(colorder, "replicate")
  #return(data)

  if (!preaveraged) {
    # to unite condition and replicates into condition
    data <- tidyr::unite(data, condition, condition, replicate, sep=".")
    # to average the data
    d1 <- tidyr::gather(data, parameter, reading, -id, -condition) #[ ,c(1:ncol(data))]
    d1$condition <- gsub(paste0("\\.", rep, "[0-9]+"), "", d1$condition)
    d1$condition <- gsub("\\.[0-9]", "", d1$condition)

    datac <- plyr::ddply(d1, plyr::.(id,condition,parameter), .drop=TRUE,
                         .fun = function(xx, col) {
                           c(mean = mean   (xx[[col]], na.rm=TRUE)
                             #sd   = sd     (xx[[col]], na.rm=na.rm)
                           )
                         },
                         "reading"
    )
    averageddata <- tidyr::spread(datac, parameter, mean)
    averageddata <- averageddata[ ,colorder1]
  } else {
    data <- tidyr::separate(data, condition, into=c("condition", "replicate"), sep="\\.")
    data$replicate <- NULL
    data["replicate"] <- NULL
    averageddata <- data
  }
  names(averageddata) <- gsub("E[CT]$", "MDT", names(averageddata))
  colorder2 <- gsub("E[CT]$", "MDT", colorder1)
  # return(averageddata)

  if (normalizedAUC) {
    # averageddata1 <- averageddata[ ,c("id", "condition", "AUC", "CV")]
    averageddata1 <- averageddata[ ,!(names(averageddata) %in% c("sumUniPeps","sumPSMs","countNum","Slope","MDT","R2","EC50","ET50"))]
    denominator <- unique(averageddata1$condition)[grep(refkeyword, unique(averageddata1$condition))][1]
    cat("The reference condition is", denominator, ",for AUC normalization.\n")
    averageddata1 <- averageddata1 %>% rowwise() %>% mutate(role=ifelse(condition==denominator,"denominator","numerator"))
    if (stableref) {
      averageddata_withbase <- filter(averageddata1, role=="denominator")
      CVcutoff <- median(averageddata_withbase$CV)+stableref_nMAD*mad(averageddata_withbase$CV)
      #print(CVcutoff)
      averageddata_withbase <- filter(averageddata_withbase, CV < CVcutoff)
      averageddata1 <- subset(averageddata1, id %in% averageddata_withbase$id)
    } else {
      averageddata_withbase <- filter(averageddata1, role=="denominator")
      averageddata1 <- subset(averageddata1, id %in% averageddata_withbase$id)
    }
    #print(head(averageddata1))
    averageddata1$CV <- NULL
    averageddata1$role <- NULL
    averageddata2 <- tidyr::gather(averageddata1, parameter, reading, -id, -condition)
    averageddata2 <- plyr::ddply(averageddata2, plyr::.(id,parameter), function(data) {
      base=data[(data$condition==denominator), ]$reading
      if (log2scale) {
        data <- plyr::mutate(data, reading=reading-base)
      } else {
        data <- plyr::mutate(data, reading=reading/base)
      }
    })
    averageddata1 <- averageddata[ ,-c(3:(nread+2))]
    averageddata1$AUC <- NULL
    averageddata2 <- subset(averageddata2, condition!=denominator) # to remove the denominator
    averageddata2 <- tidyr::spread(averageddata2, parameter, reading)
    averageddata <- merge(averageddata2, averageddata1)
    averageddata <- averageddata[ ,colorder2]
  }
  #print(head(averageddata))
  # return(averageddata)

  if (length(AUC_uplimit)==0) {
    AUC_uplimit <- round(median(averageddata$AUC,na.rm=T)+nMAD*mad(averageddata$AUC,na.rm=T),3)
  }
  if (length(AUC_lowlimit)==0) {
    AUC_lowlimit <- round(median(averageddata$AUC,na.rm=T)-nMAD*mad(averageddata$AUC,na.rm=T),3)
  }
  cat("The upper AUC cutoff set at", AUC_uplimit, "\n")
  cat("The lower AUC cutoff set at", AUC_lowlimit, "\n")
  nrowdata <- nrow(averageddata)
  cat("The number of total proteins to consider is:", nrowdata, "\n")
  # return(averageddata)

  if (onlyshowstabilized) {
    data_changed <- subset(averageddata, R2>0.8 & (AUC>AUC_uplimit))
  } else {
    data_changed <- subset(averageddata, R2>0.8 & (AUC>AUC_uplimit | AUC<AUC_lowlimit))
  }
  data_changed <- tidyr::drop_na(data_changed, MDT) #data.table::na.omit(data_changed, cols="MDT")
  maxdose <- max(as.numeric(names(data_changed)[c(3:(nread+2))]))
  data_changed <- subset(data_changed, MDT>0 & MDT<maxdose)
  nrowdata <- nrow(data_changed)
  cat("The number of proteins passing the AUC cutoff is:", nrowdata, "\n")

  if (fcthreshold>0) { # similar cutoff as the ms_ITDR_filter
    cutoff_high <- round((1+nMAD*baselineMAD)*(1+fcthreshold), 3)
    cutoff_low <- round((1-nMAD*baselineMAD)/(1+fcthreshold), 3)
    if (log2scale) {
      cutoff_high <- round(log2(cutoff_high),3)
      cutoff_low <- round(log2(cutoff_low),3)
    }
    cat("The upper fold change cutoff threshold for shift set at", cutoff_high, "\n")
    cat("The lower fold change cutoff threshold for shift set at", cutoff_low, "\n")
  } else if (fcthreshold==0) {
    cutoff_high <- 1.05
    cutoff_low <- 0.95
    if (log2scale) {
      cutoff_high <- round(log2(cutoff_high),3)
      cutoff_low <- round(log2(cutoff_low),3)
    }
  }

  if (length(checkpointposition)==0) { # default to check the last three points
    checkpointposition <- c(8,9,10)-(10-nread)+2
  } else {
    checkpointposition <- checkpointposition+2
  }
  fkeep <- NULL
  for (i in 1:nrowdata) {
    highest <- max(as.numeric(data_changed[i,checkpointposition]), na.rm=T)
    lowest <- min(as.numeric(data_changed[i,checkpointposition]), na.rm=T)
    if (highest >= cutoff_high & data_changed[i,"Slope"] >0) {
      fkeep <- c(fkeep, i)
    } else if (lowest <= cutoff_low & data_changed[i,"Slope"] <0) {
      fkeep <- c(fkeep, i)
    }
  }
  if (length(fkeep)) {
    data_changed <- data_changed[fkeep, ]
  } else {
    stop("Opps, no proteins passed the set fold change threshold!")
  }
  cat("The number of potential hits is:", nrow(data_changed), "\n")
  cat("They are:\n")
  print(data_changed$id)

  data_extra <- averageddata[0, ]
  if (length(idtoshow)) {
    fkeep <- NULL
    for (i in 1:length(idtoshow)) {
      hits <- grep(paste0("^", idtoshow[i]), averageddata$id, value=FALSE)
      fkeep <- c(fkeep, hits)
    }
    if (length(fkeep)>0) {
      data_extra <- averageddata[fkeep, ]
      data_extra <- dplyr::setdiff(data_extra, data_changed)
    }
    if (nrow(data_extra)) {
      cat("These proteins would be shown in plot, although not as hits:\n")
      print(data_extra$id)
    }
  }

  if (onlyshowstabilized) {
    data_remaining <- subset(averageddata, R2<=0.8 | (AUC<=AUC_uplimit))
  } else {
    data_remaining <- subset(averageddata, R2<=0.8 | (AUC<=AUC_uplimit & AUC>=AUC_lowlimit))
  }
  cat("The number of background proteins is:", nrow(data_remaining), "\n")

  averageddata <- tidyr::separate(averageddata, id, into=c("id", "name"), sep="\n")
  averageddata <- merge(proteininfo, averageddata[ ,-2])
  if (length(attr(averageddata,"outdir"))==0 & length(outdir)>0) {
    attr(averageddata,"outdir") <- outdir
  }
  ms_filewrite(averageddata, paste0(dataname, "_data_mapped_in_R2AUC_plot.txt"),
               outdir=outdir)
  #return(data_changed)

  data_changed1 <- data_changed # save a copy
  if (printGeneName & !pfdatabase) {
    data_remaining <- data_remaining %>% rowwise() %>%
      mutate(id = strsplit(id,"\n")[[1]][2])
    data_changed <- data_changed %>% rowwise() %>%
      mutate(id = strsplit(id,"\n")[[1]][2])
    if (nrow(data_extra)) {
      data_extra <- data_extra %>% rowwise() %>%
        mutate(id = strsplit(id,"\n")[[1]][2])
    }
  }

  colorpanel <- unique(c(brewer.pal(8, "Set1")[c(1:5,7:8)], brewer.pal(7, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(11,"Spectral")))
  # colorpanel <- c("#FF7F00", "#FDBF6F", "#FB9A99", "#F781BF", "#E7298A", "#CAB2D6", "#E6AB02",
  #                 "#E41A1C", "#E31A1C", "#D95F02", "#D53E4F", "#F0027F", "#BF5B17", "#BEAED4",
  #                 "#B2DF8A", "#A6CEE3", "#B15928", "#A6761D", "#A65628", "#9E0142", "#377EB8",
  #                 "#7FC97F", "#7570B3", "#6A3D9A", "#5E4FA2", "#666666", "#66A61E", "#4DAF4A",
  #                 "#386CB0", "#984EA3", "#3288BD", "#1F78B4", "#33A02C", "#1B9E77")
  colorpanel <- c(colorpanel, colorpanel, colorpanel, colorpanel)

  if (colorbackground) {
    n1 <- ggplot() + geom_point(data=data_remaining, aes(x=R2, y=AUC, shape=condition, colour=condition),size=2,alpha=0.5)
  } else if (graybackground) {
    n1 <- ggplot() + geom_point(data=data_remaining, aes(x=R2, y=AUC, shape=condition),colour="gray",size=2,alpha=0.6)
  } else {
    n1 <- ggplot() + geom_point(data=data_remaining, aes(x=R2, y=AUC), colour="gray",size=2,alpha=0.6)
  }
  if (nrow(data_extra)) {
    n1 <- n1 + geom_point(data=data_extra, aes(x=R2, y=AUC, shape=condition), colour="blue", size=2)
    if (labelnodes) {
      n1 <- n1 + geom_text_repel(data=data_extra, aes(x=R2, y=AUC, label=id), colour="blue", size=3, max.overlaps=50)
    }
  }
  if (colornodes & nrow(data_changed)>0 & nrow(data_changed)<=144) {
    if (nodesizebyMDT) {
      n1 <- n1 + geom_point(data=data_changed, aes(x=R2, y=AUC, shape=condition, size=-log10(MDT), colour=id))
    } else {
      n1 <- n1 + geom_point(data=data_changed, aes(x=R2, y=AUC, shape=condition, colour=id), size=4)
    }
    if (labelnodes) {
      n1 <- n1 + geom_text_repel(data=data_changed, aes(x=R2, y=AUC, label=id, colour=id), size=3, max.overlaps=50)
    }
    n1 <- n1 + guides(colour=FALSE)
  } else if (nrow(data_changed)) {
    if (nodesizebyMDT) {
      n1 <- n1 + geom_point(data=data_changed, aes(x=R2, y=AUC, shape=condition, size=-log10(MDT)), colour="orange")
    } else {
      n1 <- n1 + geom_point(data=data_changed, aes(x=R2, y=AUC, shape=condition), size=4, colour="orange")
    }
    if (labelnodes) {
      n1 <- n1 + geom_text_repel(data=data_changed, aes(x=R2, y=AUC, label=id), size=3, colour="orange", max.overlaps=50)
    }
  }
  n1 <- n1 + scale_colour_manual(drop=FALSE, values=colorpanel) +
    geom_vline(xintercept=0.8, colour="gray50", linetype="longdash") +
    geom_hline(yintercept=AUC_uplimit, colour="gray50", linetype="longdash") +
    geom_hline(yintercept=AUC_lowlimit, colour="gray50", linetype="longdash") +
    scale_x_continuous(breaks=seq(0,1,0.2),limits=c(0,1)) + scale_size(range=c(2,5))
  if (length(plottitle)) { n1 <- n1 + ggtitle(plottitle) + theme(plot.title=element_text(hjust=0.5, size=rel(2), face="bold"))}
  if (length(yscale)) { n1 <- n1 + coord_cartesian(ylim = yscale) }
  n2 <- n1 + labs(x="Dose-Response trend (Reading R-squared)", y="Relative Shift (Reading Area under the curve)") +
    theme_bw() + theme(text=element_text(size=12),
                       axis.text=element_text(size=12),
                       legend.position="bottom",
                       legend.text=element_text(size=12), aspect.ratio=1)
  ggsave(filename=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),dataname,"_R2_vs_AUC.pdf"), n2, height=11.69, width=8.27*2)

  if (returnplot) {
    n2 <- n1 + labs(x="Dose-Response trend (Reading R-squared)",y="Relative Shift (Reading coefficient of variation)") +
      theme_bw() + theme(text = element_text(size=12), axis.text=element_text(size=12), legend.position="none", aspect.ratio=1)
    if (length(attr(n2,"outdir"))==0 & length(outdir)>0) {
      attr(n2,"outdir") <- outdir
    }
    return(n2)
  }

  if (printBothName) {
    data_changed1 <- tidyr::separate(data_changed1, id, into=c("id","name1","name2"), sep="\n")
    data_changed1 <- merge(proteininfo, data_changed1[ ,-c(2:3)])
    if (length(attr(data_changed1,"outdir"))==0 & length(outdir)>0) {
      attr(data_changed1,"outdir") <- outdir
    }
  } else {
    data_changed1 <- tidyr::separate(data_changed1, id, into=c("id","name"), sep="\n")
    data_changed1 <- merge(proteininfo, data_changed1[ ,-2])
    if (length(attr(data_changed1,"outdir"))==0 & length(outdir)>0) {
      attr(data_changed1,"outdir") <- outdir
    }
  }
  return(list(hits=data_changed1, total=averageddata))
}
