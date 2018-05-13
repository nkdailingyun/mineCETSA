#' ms_IDIT_R2AUC
#'
#' Function to generate a R2 vs AUC plot for the ITDR or ITTR dataset
#' The R2 (x-axis value), AUC (y-axis) and MDT/MTT (size of dot) are from the input data
#'
#' @param data ITDR or ITTR dataset to generate R2 vs AUC plot, the provided
#' data should be associated with fitting parameters, i.e., after ms_ITDR_fitting
#' or ms_ITTR_fitting
#' @param nread number of reading channels or sample treatements, default value 10
#' @param printBothName whether to print both gene name and protein name, default set to TRUE
#' @param printGeneName whether to print only gene name, default set to FALSE, when
#' both printBothName and printGeneName are FALSE, only the protein name is shown
#' @param pfdatabase whether the data is a malaria dataset, default set to FALSE
#' @param rep the rep indicator used in the naming of experimental conditions,
#' such as "r", or "rep" or ""
#' @param onlyshowstabilized whether to show only stablized hits, ie,
#' direct targets, default set to FALSE
#' @param standardizedAUC whether to use the standardized AUC values, ie, 1.0,
#' default set to TRUE
#' @param preaveraged whether the data has already been averaged from replicates,
#' default set to FALSE
#' @param normalizedAUC whether to use the AUC values against a reference control
#' such as 37C readings, default set to FALSE
#' @param refkeyword a keyword used for the indication of reference control, as
#' most of the time the reference is 37C expression level, default value is 37C
#' @param nMAD level of significance, default set at 2.5
#' @param baselineMAD MAD of baseline variation, default value is 0; if not
#' provided, it will be calculated based on the readings from the lowest two
#' dose groups
#' @param fcthreshold short for fold change threshold, similar parameter as in
#' ms_ITD/TR_filter() function, however this is an added-on selection feature
#' in this R2AUC plot, and applied to the averaged (and/or normalized) data,
#' so default value is 0
#' @param checkpointposition refering to the positions of dose points to check
#' whether their readings surpass fcthreshold, default value is NULL, which
#' would automatically check the highest 3 dose points
#' @param keepreplicate whether to only keep the curves that are measured with
#' all replicates under at least one experimental conditions, default set to
#' FALSE, when there is only one experimental condition in the dataset, this
#' means only keep the full replicates
#' @param stableref whether to only keep the curves with stable/flat reference
#' (most of the time, 37C) curves, default set to FALSE
#' @param stableref_nMAD level of significance for control of the reference
#' (most of the time, 37C) sample stability, default set at 3.5
#' @param colorbackground whether to color the background with color and shape
#' according to the condition, default set to TRUE
#' @param graybackground whether to use gray color for the background but
#' different shapes according to the condition, default set to FALSE
#' @param nodesizebyMDT whether to set node size according to MDT value,
#' default set to TRUE
#' @param PSMcutoff whether to apply PSM threshold cutoff on hit selection,
#' default set to TRUE
#' @param PSMthreshold the threshold of averaged PSM numbers to consider protein
#' quantification as reliable, default value is 3
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
#'
#'
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


ms_IDIT_R2AUC <-  function(data, nread=10, printBothName=TRUE, printGeneName=FALSE,
                           pfdatabase=FALSE, rep="r", onlyshowstabilized=FALSE,
                           byAUC=TRUE, standardizedAUC=TRUE, preaveraged=FALSE,
                           normalizedAUC=FALSE, refkeyword="37C",
                           AUC_uplimit=NULL, AUC_lowlimit=NULL,
                           nMAD=2.5, baselineMAD=0, fcthreshold=0, checkpointposition=NULL,
                           keepreplicate=FALSE, stableref=FALSE, stableref_nMAD=3.5,
                           colorbackground=TRUE, graybackground=FALSE, nodesizebyMDT=TRUE,
                           colornodes=TRUE, labelnodes=TRUE, yscale=NULL,
                           PSMcutoff=TRUE, PSMthreshold=3, PSMcutoffbycondition=FALSE,
                           idtoshow=NULL, idtoexclude=NULL, plottitle=NULL) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)

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

  if (printBothName) {
    data <- data %>% rowwise() %>% mutate(description1 = getProteinName(description,pfdatabase)) %>%
      mutate(description2 = getGeneName(description)) %>%
      mutate(id = paste(id, description1, description2, sep="\n"))
    data$description1<-NULL
    data$description2<-NULL
  } else if (printGeneName) {
    data <- data %>% rowwise() %>%
      mutate(description = getGeneName(description)) %>%
      mutate(id = paste(id, description, sep="\n"))
  } else {
    data <- data %>% rowwise() %>%
      mutate(description = getProteinName(description,pfdatabase)) %>%
      mutate(id = paste(id, description, sep="\n"))
  }
  data$description<-NULL

  if (sum(is.na(data$id))>0) { print(paste0("There were ", sum(is.na(data$id)),
                                            " entries not properly parsed..."))}
  else { print("Protein labels were successfully parsed.") }

  #to separate condition into condition and replicates
  nset <- length(unique(data$condition))
  data <- tidyr::separate(data, condition, into=c("condition", "replicate"), sep="\\.")
  ncondition <- length(unique(data$condition))
  nreplicate <- length(unique(data$replicate))
  uniquecond <- unique(data[ ,c("condition", "replicate")])
  row.names(uniquecond) <- NULL
  print("Replicates information were extracted as follows:")
  print(as.data.frame(uniquecond))

  # to keep the data in at least one condition with full replicates
  if (keepreplicate) {
    conditionrep <- dplyr::count(uniquecond, condition)
    data_freq <- dplyr::count(data, id, condition)
    id_keep <- data_freq[0, ]
    for (i in 1: nrow(conditionrep)) {
      id_keep <- rbind(id_keep, subset(data_freq, condition==conditionrep$condition[i] & n==conditionrep$n[i]))
    }
    data <- merge(data, id_keep)
    data$n <- NULL
    print(paste(nrow(data), "measurements were measured with complete replicates in at least one condition.", sep=" "))
  }
  #return(data)
  # To select the proteins with more than 3 PSM (average)
  if (PSMcutoff) {
    PSMcol <- grep("PSM", names(data), value=F)
    PSMname <- names(data)[PSMcol]
    names(data)[PSMcol] <- "PSM"
    if (PSMcutoffbycondition) {
      PSMtable <- plyr::ddply(data, plyr::.(id,condition), summarize, PSMmean=mean(PSM))
      PSMtable <- PSMtable %>% filter(PSMmean>=PSMthreshold) %>% count(id) %>% filter(n==ncondition)
    } else {
      PSMtable <- plyr::ddply(data, plyr::.(id), summarize, PSMmean=mean(PSM))
      PSMtable <- subset(PSMtable, PSMmean>=PSMthreshold)
    }
    nkeep <- PSMtable$id
    fkeep <- which(data$id %in% nkeep)
    names(data)[PSMcol] <- PSMname
    data_PSMsmall <- data[-fkeep, ]
    data <- data[fkeep, ]
    print("PSM cutoff was applied.")
  }

  if (length(idtoexclude)) {
    fkeep <- NULL
    for (i in 1:length(idtoexclude)) {
      hits <- grep(paste0("^", idtoexclude[i]), data$id, value=FALSE)
      fkeep <- c(fkeep, hits)
    }
    if (length(fkeep)>0) { data <- data[-fkeep, ] }
  }

  data$STD <- apply(data[,c(4:(nread+3))], 1, sd, na.rm=T)
  data$Mean <- apply(data[,c(4:(nread+3))], 1, mean, na.rm=T)
  data <- dplyr::mutate(data, CV=STD/Mean, AUC=rowSums(data[, c(4:(nread+3))],na.rm=T))
  data$STD <- NULL
  data$Mean <- NULL
  if (standardizedAUC) {
    data <- dplyr::mutate(data, AUC=AUC/nread)
  }
  #return(data)

  if (!preaveraged) {
    # to unite condition and replicates into condition
    data <- tidyr::unite(data, condition, condition, replicate, sep=".")
    # to average the data
    d1 <- tidyr::gather(data[ ,c(1:ncol(data))], parameter, reading, -id, -condition)
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
    a <- which(!is.na(as.numeric(names(averageddata)))==TRUE)
    a <- c(names(averageddata)[c(1:(a[1]-1))],
           gtools::mixedsort(names(averageddata)[c(a[1]:a[length(a)])]),
           names(averageddata)[(c(a[length(a)]+1):ncol(averageddata))])
    averageddata <- averageddata[ ,a]
  } else {
    data <- tidyr::separate(data, condition, into=c("condition", "replicate"), sep="\\.")
    data$replicate <- NULL
    data["replicate"] <- NULL
    averageddata <- data
  }
  names(averageddata) <- gsub("E[CT]", "MDT", names(averageddata))
  #return(averageddata)

  if (normalizedAUC) {
    averageddata1 <- averageddata[ ,c("id", "condition", "AUC", "CV")]
    denominator <- unique(averageddata1$condition)[grep(refkeyword, unique(averageddata1$condition))][1]
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
    averageddata1 <- plyr::ddply(averageddata1, plyr::.(id), function(data) {
      base_AUC=data[(data$role=="denominator"), ]$AUC
      #print(base_AUC)
      data <- plyr::mutate(data, AUC=AUC/base_AUC)
    })
    averageddata$AUC <- NULL
    averageddata <- merge(averageddata, averageddata1[ ,c("id","condition","AUC")])
    averageddata <- subset(averageddata, condition!=denominator)
    a <- which(!is.na(as.numeric(names(averageddata)))==TRUE)
    a <- c(names(averageddata)[c(1:(a[1]-1))],
           gtools::mixedsort(names(averageddata)[c(a[1]:a[length(a)])]),
           names(averageddata)[(c(a[length(a)]+1):ncol(averageddata))])
    averageddata <- averageddata[ ,a]
  }
  #print(head(averageddata))
  #return(averageddata)

  if (length(AUC_uplimit)==0) {
    AUC_uplimit <- round(median(averageddata$AUC,na.rm=T)+nMAD*mad(averageddata$AUC,na.rm=T),3)
  }
  if (length(AUC_lowlimit)==0) {
    AUC_lowlimit <- round(median(averageddata$AUC,na.rm=T)-nMAD*mad(averageddata$AUC,na.rm=T),3)
  }
  print(paste0("The upper AUC cutoff set at ", AUC_uplimit, "."))
  print(paste0("The lower AUC cutoff set at ", AUC_lowlimit, "."))

  if (byAUC) {
    nrowdata <- nrow(averageddata)
    print(nrowdata)
    #return(averageddata)
    if (onlyshowstabilized) {
      data_changed <- subset(averageddata, R2>0.8 & (AUC>AUC_uplimit))
    } else {
      data_changed <- subset(averageddata, R2>0.8 & (AUC>AUC_uplimit | AUC<AUC_lowlimit))
    }
    data_changed <- tidyr::drop_na(data_changed, MDT) #data.table::na.omit(data_changed, cols="MDT")
    maxdose <- max(as.numeric(names(data_changed)[c(3:(nread+2))]))
    data_changed <- subset(data_changed, MDT>0 & MDT<maxdose)
    if (fcthreshold>0) { # similar cutoff as the ms_ITDR_filter
      cutoff_high <- round((1+nMAD*baselineMAD)*(1+fcthreshold), 3)
      cutoff_low <- round((1-nMAD*baselineMAD)/(1+fcthreshold), 3)
      print(paste0("The upper cutoff threshold for shift set at ", cutoff_high, "."))
      print(paste0("The lower cutoff threshold for shift set at ", cutoff_low, "."))
      nrowdata <- nrow(data_changed)
      print("The number of proteins passing the AUC cutoff is:")
      print(nrowdata)
      if (length(checkpointposition)==0) { # default to check the last three points
        checkpointposition <- c(8,9,10)-(10-nread)+2
      } else {
        checkpointposition <- checkpointposition+2
      }
      fkeep <- NULL
      for (i in 1:nrowdata) {
        # if (excludelastpoint) {
        #   highest <- max(as.numeric(data_changed[i,c((nread+2-ncheckpoint):
        #                                                 (nread+1))]), na.rm=T)
        #   lowest <- min(as.numeric(data_changed[i,c((nread+2-ncheckpoint):
        #                                                (nread+1))]), na.rm=T)
        # } else {
        #   highest <- max(as.numeric(data_changed[i,c((nread+3-ncheckpoint):
        #                                                (nread+2))]), na.rm=T)
        #   lowest <- min(as.numeric(data_changed[i,c((nread+3-ncheckpoint):
        #                                               (nread+2))]), na.rm=T)
        # }
        highest <- max(as.numeric(data_changed[i,checkpointposition]), na.rm=T)
        lowest <- min(as.numeric(data_changed[i,checkpointposition]), na.rm=T)
        if (highest >= cutoff_high & data_changed[i,"Slope"] <0) {
          fkeep <- c(fkeep, i)
        } else if (lowest <= cutoff_low & data_changed[i,"Slope"] >0) {
          fkeep <- c(fkeep, i)
        }
      }
      if (length(fkeep)) {
        data_changed <- data_changed[fkeep, ]
      } else {
        stop("Opps, no proteins passed the set fold change threshold!")
      }
    }
    print("The number of potential hits is:")
    print(nrow(data_changed))
    print("They are:")
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
        print("These proteins would be shown in plot, although not as hits: ")
        print(data_extra$id)
      }
    }

    if (onlyshowstabilized) {
      data_remaining <- subset(averageddata, R2<=0.8 | (AUC<=AUC_uplimit))
    } else {
      data_remaining <- subset(averageddata, R2<=0.8 | (AUC<=AUC_uplimit & AUC>=AUC_lowlimit))
    }
    print("The number of background proteins is:")
    print(nrow(data_remaining))
  }

  averageddata <- tidyr::separate(averageddata, id, into=c("id", "name"), sep="\n")
  averageddata <- merge(proteininfo, averageddata[ ,-2])
  if (length(attr(averageddata,"outdir"))==0 & length(outdir)>0) {
    attr(averageddata,"outdir") <- outdir
  }
  ms_filewrite(averageddata, paste0(dataname, "_data_averaged_R2AUC.txt"),
               outdir=outdir)
  #return(data_changed)

  data_changed1 <- data_changed # save a copy
  if (printGeneName) {
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
  if (byAUC) {
    if (colorbackground) {
      n1 <- ggplot() + geom_point(data=data_remaining, aes(x=R2, y=AUC, shape=condition, colour=condition),size=2,alpha=0.05)
    } else if (graybackground) {
      n1 <- ggplot() + geom_point(data=data_remaining, aes(x=R2, y=AUC, shape=condition),colour="gray",size=2,alpha=0.6)
    } else {
      n1 <- ggplot() + geom_point(data=data_remaining, aes(x=R2, y=AUC), colour="gray",size=2,alpha=0.6)
    }
    if (nrow(data_extra)) {
      n1 <- n1 + geom_point(data=data_extra, aes(x=R2, y=AUC, shape=condition), colour="blue", size=2)
      if (labelnodes) {
        n1 <- n1 + geom_text_repel(data=data_extra, aes(x=R2, y=AUC, label=id), colour="blue", size=3)
      }
    }
    if (colornodes & nrow(data_changed)>0 & nrow(data_changed)<=144) {
      if (nodesizebyMDT) {
        n1 <- n1 + geom_point(data=data_changed, aes(x=R2, y=AUC, shape=condition, size=-log10(MDT), colour=id))
      } else {
        n1 <- n1 + geom_point(data=data_changed, aes(x=R2, y=AUC, shape=condition, colour=id), size=4)
      }
      if (labelnodes) {
        n1 <- n1 + geom_text_repel(data=data_changed, aes(x=R2, y=AUC, label=id, colour=id), size=3)
      }
    } else if (nrow(data_changed)) {
      if (nodesizebyMDT) {
        n1 <- n1 + geom_point(data=data_changed, aes(x=R2, y=AUC, shape=condition, size=-log10(MDT)), colour="orange")
      } else {
        n1 <- n1 + geom_point(data=data_changed, aes(x=R2, y=AUC, shape=condition), size=4, colour="orange")
      }
      if (labelnodes) {
        n1 <- n1 + geom_text_repel(data=data_changed, aes(x=R2, y=AUC, label=id), size=3, colour="orange")
      }
    }
    n1 <- n1 + scale_colour_manual(drop=FALSE, values=colorpanel) +
      geom_vline(xintercept=0.8, colour="gray50", linetype = "longdash") +
      geom_hline(yintercept=AUC_uplimit, colour="gray50", linetype = "longdash") +
      geom_hline(yintercept=AUC_lowlimit, colour="gray50", linetype = "longdash") +
      scale_x_continuous(breaks=seq(0,1,0.2)) + scale_size(range=c(3,5)) + coord_cartesian(xlim = c(0.0,1.0))
    if (length(plottitle)) { n1 <- n1 + ggtitle(plottitle) + theme(plot.title=element_text(hjust=0.5, size=rel(2), face="bold"))}
    if (length(yscale)) { n1 <- n1 + coord_cartesian(ylim = yscale) }
    n2 <- n1 + labs(x="Dose-Response trend (Reading R-squared)",y="Relative Shift (Reading Area under the curve)") +
      theme_bw() + theme(text=element_text(size=12), axis.text=element_text(size=8), legend.position="bottom", legend.text=element_text(size=4), aspect.ratio=1)
  }

  ggsave(filename = paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M"),"_R2_vs_AUC.pdf"), n2, height=11.69, width=8.27*2)

  n2 <- n1 + labs(x="Dose-Response trend (Reading R-squared)",y="Relative Shift (Reading coefficient of variation)") +
    theme_bw() + theme(text = element_text(size=12), axis.text=element_text(size=8), legend.position="none", aspect.ratio=1)
  if (length(attr(n2,"outdir"))==0 & length(outdir)>0) {
    attr(n2,"outdir") <- outdir
  }
  if (printBothName) {
    data_changed1 <- tidyr::separate(data_changed1, id, into=c("id", "name1", "name2"), sep="\n")
    data_changed1 <- merge(proteininfo, data_changed1[ ,-c(2:3)])
    if (length(attr(data_changed1,"outdir"))==0 & length(outdir)>0) {
      attr(data_changed1,"outdir") <- outdir
    }
  } else {
    data_changed1 <- tidyr::separate(data_changed1, id, into=c("id", "name"), sep="\n")
    data_changed1 <- merge(proteininfo, data_changed1[ ,-2])
    if (length(attr(data_changed1,"outdir"))==0 & length(outdir)>0) {
      attr(data_changed1,"outdir") <- outdir
    }
  }
  return(list(hits=data_changed1, total=averageddata, plot=n2))
}
