#' ms_ggplotting_rep
#'
#' Function to generate pdf files with multipanel ggplots for melt curve data,
#' particularly for dataset with replicate runs
#'
#' @param data isothermal dataset to plot
#' @param legenddata dataset used for condition extraction, at least one protein
#' inside this dataset should contains the full set of experiment conditions
#' @param levelvector a vector of conditions, preferably starting with Ctrl sets
#' @param nread number of reading channels or sample treatements
#' @param remsinglecondprot whether orphan proteins to be plotted,
#' default value is TRUE, so to exclude them from plotting
#' @param PSMcutoff whether to apply PSM threshold cutoff on hit selection,
#' default set to FALSE
#' @param PSMthreshold the threshold of averaged PSM numbers to consider protein
#' quantification as reliable, default value is 3
#' @param nreplicate number of replicates, default value is 2, change the value
#' accordingly, this is relavent to automatic coloring scheme,
#' up to 4 is possible for now
#' @param fitremout whether to segregate the proteins with messy melt curves
#' @param ctrlcond if necessary, could used to specify what conditions to be
#' referred as control conditions
#' @param bottomcutoff the average of the last three points should be lower than
#' specified bottom cutoff value, which is 0.4 by default
#' @param topcutoff the average of the first three points should be higher than
#' specified bottom cutoff value, which is 0.8 by default
#' @param variancecutoff whether to segregate the proteins with large inter-replicate variance
#' @param nMAD_var the number of MADs to set the significance cutoff about variance distribution,
#' default value is 2.5
#' @param topasone whether the top plateau has to be fixed, i.e., 1.0
#' @param normTop whether to normalize the AUC based on Top three readings
#' @param dotconnect whether to simply dot connect the readings for each curve
#' @param PSManno whether to annotate the plots with PSM and uniPeptide number
#' @param PSMannoypos the starting y postion of PSM/uniPeptide number annotation from top,
#' default value is 0.5
#' @param PSMannoyinterval the interval of PSM/uniPeptide number annotation per line,
#' default value is 0.08
#' @param presetcolor whether to use the pre-defined color scheme
#' @param plotfitremout whether to plot out messy melt curves
#' @param plotvarremout whether to plot out large inter-replicate variance melt curves
#' @param colorpanel a vector of customizable color scheme provided by the user
#' @param withset whether there is set column to perform facet_grid
#' @param commonlegend whether to use one common legend for whole page of plots
#' @param layout a vector indicating the panel layout for multi-panel plots per page
#' @param pdfname name for the pdf plots file
#'
#'
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

ms_ggplotting_rep <- function(data, legenddata=NULL, levelvector=NULL, nread=10,
                              remsinglecondprot=TRUE, PSMcutoff=FALSE, PSMthreshold=3,
                              remfragment=FALSE, remribosomal=FALSE,
                              fitremout=FALSE, ctrlcond=NULL, bottomcutoff=0.4, topcutoff=0.8,
                              variancecutoff=FALSE, nMAD_var=2.5,
                              nreplicate=2, topasone=TRUE, normTop=TRUE, dotconnect=FALSE,
                              pfdatabase=FALSE, printBothName=TRUE, printGeneName=FALSE,
                              PSManno=TRUE, PSMannoypos=0.5, PSMannoyinterval=0.08,
                              presetcolor=TRUE, extraidtocomplete=NULL,
                              plotfitremout=TRUE, plotvarremout=TRUE,
                              colorpanel=NULL, withset=FALSE, commonlegend=TRUE,
                              layout=c(5,5), pdfname="ggplotting.pdf", external=TRUE) {

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

  png(filename=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                      dataname, "_curve replicate distribution.png"))
  barplot(table(table(data$id)),
          main="Melting curve replicate distribution",
          xlab="Number of curves per unique protein in replicate exp")
  dev.off()

  if (nreplicate>4) {stop("Only up to four replicates are supported")}
  checkpos <- c(1,3)
  if (withset) {
    setpos <- grep("^set", names(data))
    checkpos <- c(checkpos, setpos)
  }
  if (any(duplicated(data[, checkpos]))) {
    print("Warning for duplicated protein name entries")
    print("Double check the following proteins for duplicated entries:")
    print(paste0(data[duplicated(data[, checkpos]), ]$id," in ",
                 data[duplicated(data[, checkpos]), ]$condition))
    stop("1.Remove the duplicated entires from original dataset then start again!")
  }

  ncond <- length(unique(data$condition))
  if(length(levelvector)!=ncond | !setequal(levelvector, unique(data$condition))) {
    stop("Make sure you provide the correct number of conditions,
         first Controls then Treatments")
    #levelvector=sort(unique(data$condition));
  }

  # to concatenate id and description
  nrowdata <- nrow(data)
  getGeneName <- function(x) {return (strsplit(strsplit(x, "GN=")[[1]][2], " ")[[1]][1])}
  getProteinName <- function(x) {return (strsplit(x, " OS=")[[1]][1])}
  if (pfdatabase) {
    getProteinName <- function(x) {return (gsub("product=", "", strsplit(x, "\\|")[[1]][2]))}
  }

  if (length(extraidtocomplete)) {
    fkeep <- NULL
    for (i in 1:length(extraidtocomplete)){
      hits <- grep(paste0("^", extraidtocomplete[i], "$"), data$id, value=FALSE)
      fkeep <- c(fkeep, hits)
    }
    data_extra <- data[fkeep, ]
  }

  if (printBothName) {
    data <- data %>% rowwise() %>% mutate(description1 = getProteinName(description)) %>%
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
      mutate(description = getProteinName(description)) %>%
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
  # for(i in 1:nrowdata){
  #   if(geneName){
  #     data[i, 'description']<-strsplit(strsplit(data[i, 'description'], "GN=")[[1]][2], " ")[[1]][1]
  #     data[i,'id']<-paste(data[i, 'id'], data[i, 'description'], sep="\n")
  #   }else{
  #     cid<-paste(data[i, 'id'], data[i, 'description'], sep="\n")
  #     data[i,'id']<-strsplit(cid, " OS")[[1]][1]
  #   }
  # }

  checkpos <- c(1,2)
  if (withset) {
    setpos <- grep("^set", names(data))
    checkpos <- c(checkpos, setpos)
  }
  if (any(duplicated(data[, checkpos]))) {
    print("Warning for duplicated protein name entries")
    print("Double check the following proteins for duplicated entries:")
    print(paste0(data[duplicated(data[, checkpos]), ]$id," in ",
                 data[duplicated(data[, checkpos]), ]$condition))
    stop("Remove the duplicated entires from original dataset then start again!")
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

  # look for outlier proteins based on melting behavior in controls
  outliers <- NULL
  if (fitremout) {
    print("Make sure you provide fitted data with Tm and R2 values for this option!")
    ctrllist1 <- unique(grep("[Cc][Tt][Rr][Ll]", data$condition, value=TRUE))
    ctrllist2 <- unique(grep("[Cc][Oo][Nn][Tt][Rr][Oo][Ll]", data$condition, value=TRUE))
    ctrllist3 <- unique(grep("[Dd][Mm][Ss][Oo]", data$condition, value=TRUE))
    ctrllist <- c(ctrllist1, ctrllist2, ctrllist3, ctrlcond)
    if(length(ctrllist)==0) {
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
      print(paste0(length(tmp), " measurements were messy in melting behavior and removed."))
      outlierid <- data[tmp, ]$id
      outlierid1 <- which(data$id %in% outlierid)
      outliers <- data[outlierid1, ]
      ms_filewrite(outliers, "Messy proteins.txt", outdir=outdir)
      #print(paste0("The outlier protein list is ", outlierids))
      data <- data[-outlierid1, ]
    }
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

  # To select the proteins with more than 3 PSM (average)
  if (PSMcutoff & !withset) {
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

  data$condition <- factor(data$condition, levels=levelvector)

  nametempvector <- names(data[4:(nread+3)])
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

  if (withset) {
    if (external) { external_graphs(T) }
    # subsetting complete replicates:
    keep1 <- which(table(data$id) == ncond)
    if(length(keep1)==0){ stop("No proteins contain complete replicate") }
    nkeep1 <- names(keep1)
    fkeep1 <- which(data$id %in% nkeep1)
    data_complete <- data[fkeep1, ]

    plotlegend <- ms_melt_legend(data_complete, nread, colorpanel)

    pl <- ms_melt_innerplot(data, nread, topasone, dotconnect,
                            PSManno, PSM_annod, Pep_annod, PSMannoypos, PSMannoyinterval,
                            colorpanel, plotlegend, commonlegend, withset, layout)
    ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                       length(unique(data$id)),
                       "_whole_set_", pdfname), pl, height=12, width=12)
    if (external) { external_graphs(F) } # switch off the external graphs
    print("whole set plot file generated successfully.")
    stop("Done")
  }


  # subsetting complete replicates:
  keep1 <- which(table(data$id) == ncond)
  if(length(keep1)==0){ stop("No proteins contain complete replicate") }
  nkeep1 <- names(keep1)
  fkeep1 <- which(data$id %in% nkeep1)
  data_complete <- data[fkeep1, ]
  print(paste0("The number of proteins with complete replicates is: ",
               length(unique(data_complete$id))))
  print(paste0("The percentage of proteins with complete replicates is: ",
               round(length(unique(data_complete$id))/length(unique(data$id)),3)*100, "%"))

  # Calculate euclidean distance metrics
  if (nreplicate==2) {
    dism <- plyr::ddply(data_complete, "id", function(data) {
      data<-data[order(data$condition), ]
      dm<-as.matrix(dist(data[,c(3:(nread+2))]))
      Intra1<-dm[1,3]
      Intra2<-dm[2,4]
      InterCtrl<-dm[1,2]
      InterTreatment<-dm[3,4]
      variance<-InterCtrl+InterTreatment
      distance<-Intra1+Intra2-variance
      EDscore<-(Intra1+Intra2)/(10^(variance))
      eucl<-c(Intra1, Intra2, InterCtrl, InterTreatment, variance, distance, EDscore)
      data.frame(Intra1=eucl[1], Intra2=eucl[2], InterCtrl=eucl[3],
                 InterTreatment=eucl[4], variance=eucl[5],
                 distance=eucl[6], EDscore=eucl[7])
    })
  } else if (nreplicate==3) {
    dism <- plyr::ddply(data_complete, "id", function(data) {
      data<-data[order(data$condition), ]
      dm<-as.matrix(dist(data[,c(3:(nread+2))]))
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
      data.frame(Intra1=eucl[1], Intra2=eucl[2], Intra3=eucl[3],
                 InterCtrl1=eucl[4], InterCtrl2=eucl[5], InterCtrl3=eucl[6],
                 InterTreatment1=eucl[7], InterTreatment2=eucl[8], InterTreatment3=eucl[9],
                 variance=eucl[10], distance=eucl[11], EDscore=eucl[12])
    })
  } else if (nreplicate==4) {
    dism <- plyr::ddply(data_complete, "id", function(data){
      data<-data[order(data$condition), ]
      dm<-as.matrix(dist(data[,c(3:(nread+2))]))
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
      data.frame(Intra1=eucl[1], Intra2=eucl[2], Intra3=eucl[3], Intra4=eucl[4],
                 InterCtrl1=eucl[5], InterCtrl2=eucl[6], InterCtrl3=eucl[7],
                 InterCtrl4=eucl[8], InterCtrl5=eucl[9], InterCtrl6=eucl[10],
                 InterTreatment1=eucl[11], InterTreatment2=eucl[12], InterTreatment3=eucl[13],
                 InterTreatment4=eucl[14], InterTreatment5=eucl[15], InterTreatment6=eucl[16],
                 variance=eucl[17], distance=eucl[18], EDscore=eucl[19])
    })
  }

  png(filename=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                      dataname, "Euclidean distance score.png"))
  nf <- layout(mat = matrix(c(1,2),2,1, byrow=TRUE), height = c(3,1))
  par(mar=c(3.1, 3.1, 1.1, 2.1))
  hist(dism$EDscore, n=100, xlim=c(0,2), main="Euclidean distance score distribution")
  boxplot(dism$EDscore, horizontal=T, frame=F, width=0.75, ylim=c(0,2), col="black")
  dev.off()

  dism <- dism[order(dism$EDscore, decreasing=T), ]
  print(paste0("The range of ED score is between ",
               round(range(dism$EDscore)[1], 3), " and ",
               round(range(dism$EDscore)[2], 3)))
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
    print("use MAD scheme to assign variance cutoff...")
    cutoff <- round(median(dism$variance)+nMAD_var*mad(dism$variance), 3)
    print(paste0("The variance cutoff limit (", nMAD_var, "*mad) sets at ", cutoff))
    dism <- subset(dism, variance<cutoff)
    fkeep1 <- which(data_complete$id %in% dism$id)
    data_largevar <- data_complete[-fkeep1, ]
    data_complete <- data_complete[fkeep1, ]
    print(paste0("The number of reproducible proteins with complete replicates is: ",
                 length(unique(data_complete$id))))
  }

  if (external) { external_graphs(T) }

  # print out the PSM small data file
  if (PSMcutoff) {
    pl <- ms_melt_innerplot(data_PSMsmall, nread, topasone, dotconnect,
                            PSManno, PSM_annod, Pep_annod, PSMannoypos, PSMannoyinterval,
                            colorpanel, plotlegend, commonlegend, withset, layout)
    ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                       length(unique(data_PSMsmall$id)),
                       "_PSMsmall proteins", pdfname), pl, height=12, width=12)
    print("PSMsmall plot file generated successfully.")
  }

  # print out the outliers data file
  if (fitremout & plotfitremout) {
    pl <- ms_melt_innerplot(outliers, nread, topasone, dotconnect,
                            PSManno, PSM_annod, Pep_annod, PSMannoypos, PSMannoyinterval,
                            colorpanel, plotlegend, commonlegend, withset, layout)
    ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                       length(unique(outliers$id)),
                       "_Messy_proteins_", pdfname), pl, height=12, width=12)
    print("Messy plot file generated successfully.")
  }

  # print out the large variance data file
  if (variancecutoff & plotvarremout) {
    pl <- ms_melt_innerplot(data_largevar, nread, topasone, dotconnect,
                            PSManno, PSM_annod, Pep_annod, PSMannoypos, PSMannoyinterval,
                            colorpanel, plotlegend, commonlegend, withset, layout)
    ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                       length(unique(data_largevar$id)),
                       "_Non_reproducible_proteins_", pdfname), pl, height=12, width=12)
    print("Non reproducible plot file generated successfully.")
  }

  print("Generating first complete plot file, pls wait.")
  #seq1 <- as.factor(dism$id)
  if (length(extraidtocomplete)) {
    data_complete <- rbind(data_extra, data_complete)
    data_complete <- data_complete[!duplicated(data_complete), ]
    data_complete$id <- factor(data_complete$id, levels=c(unique(data_extra$id), dism$id))
  } else {
    data_complete$id <- factor(data_complete$id, levels=dism$id)
  }

  pl <- ms_melt_innerplot(data_complete, nread, topasone, dotconnect,
                          PSManno, PSM_annod, Pep_annod, PSMannoypos, PSMannoyinterval,
                          colorpanel, plotlegend, commonlegend, withset, layout)
  ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                     length(unique(data_complete$id)),
                     "_Complete replicates_", pdfname), pl, height=12, width=12)

  print("first complete plot file generated successfully.")

  # subsetting incomplete replicates:
  keep2 <- which(table(data$id) < ncond)
  nkeep2 <- names(keep2)
  fkeep2 <- which(data$id %in% nkeep2)
  data_incomp <- data[fkeep2, ]
  print(paste0("The number of proteins without complete replicates is: ",
               length(unique(data_incomp$id))))

  print("Generating second incomplete plot file, pls wait.")

  if (!("AUC" %in% colnames(data_incomp))){
    data_incomp$AUC <-rowSums(data_incomp[ ,c(3:(nread+2))])
  }
  if (normTop) {
    data_incomp$Top <- rowMeans(data_incomp[ ,c(3:5)])
  }else {
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
                          PSManno, PSM_annod, Pep_annod, PSMannoypos, PSMannoyinterval,
                          colorpanel, plotlegend, commonlegend, withset, layout)
  ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                     length(unique(data_incomp$id)),
                     "_Non_replicates_", pdfname), pl, height=12, width=12)
  if (external) { external_graphs(F) } # switch off the external graphs
  print("second incomplete plot file generated successfully.")
}
