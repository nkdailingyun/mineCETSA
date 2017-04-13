#' ms_IDIT_QC
#'
#' Function to generate a QC report for the ITDR or ITTR dataset
#'
#' @param data ITDR or ITTR dataset to perform QC analysis
#' @param foldername name for the QC result folder
#' @param reportname name for the QC report
#' @param nread number of reading channels or sample treatements, default value 10
#' @param isdatafitted whether the provided data is associated with
#' fitting parameters, i.e., after ms_ITDR_fitting or ms_ITTR_fitting
#' @param PSMcutoff whether to apply PSM threshold cutoff on hit selection,
#' default set to TRUE
#' @param PSMthreshold the threshold of averaged PSM numbers to consider protein
#' quantification as reliable, default value is 3
#'
#' @import Nozzle.R1
#' @import VennDiagram
#' @import tidyr
#' @importFrom plyr dlply ddply
#' @import dplyr
#' @import RColorBrewer
#' @import ggplot2
#'
#' @export
#'
#' @return NULL
#' @examples \dontrun{
#' ms_IDIT_QC(ITDRdata_fitted)
#' ms_IDIT_QC(ITDRdata_scaled, isdatafitted=FALSE, foldername="QC1")
#' }
#'
#'


ms_IDIT_QC <- function(data, foldername=NULL, reportname=NULL, nread=10,
                       PSMcutoff=TRUE, PSMthreshold=3, isdatafitted=TRUE) {
  # provide the fitted data
  print("Make sure you provide the scaled data with fitting parameters as input
        for maximal utilization of this QC function")

  dataname <- deparse(substitute(data))
  outdir <- data$outdir[1]
  data$outdir <- NULL

  if (!length(foldername)) {
    foldername <- paste0(dataname, "_QCreports")
  }
  if (!length(reportname)) {
    reportname <- foldername
  }
  dir.create( paste0(outdir,"/",foldername), showWarnings=FALSE )

  #to separate condition into condition and replicates
  data <- tidyr::separate(data, condition, into=c("condition", "replicate"), sep="\\.")
  # list <- strsplit(as.character(data$condition), "\\.")
  # df <- ldply(list)
  # colnames(df) <- c("condition", "replicate")
  # data <- cbind(data[,c(1,2)], df, data[,-c(1:3)])
  ncondition <- length(unique(data$condition))
  nreplicate <- length(unique(data$replicate))
  nset <- ncondition * nreplicate
  if (nset > 5) {print("To plot venn diagram, maximum five unique data sets.")}

  data_N1 <- plyr::ddply(data, .(condition, replicate), summarise, Number_run=length(unique(id)))
  data_N2 <- plyr::ddply(data, .(condition), summarise, Number_sample=length(unique(id)))
  data_N <- merge(data_N1, data_N2)

  drawvennplot <- function (data, numberofgroup, foldername, plotname) {
    vennplot <- venn.diagram(
      x = data,
      euler.d = TRUE,
      scaled = TRUE,
      height = 8,
      width = 8,
      unit = "in",
      resolution = 500,
      filename = NULL,
      lwd = 1,
      cex = 0.8,
      cat.cex = 1,
      cat.col = brewer.pal(8,"Dark2")[c(1:numberofgroup)],
      fill = brewer.pal(8,"Set2")[c(1:numberofgroup)],
      margin = 0.1
    )
    png(paste0(outdir,"/",foldername,"/",plotname), width=8, height=8, units="in", res=200)
    grid.draw(vennplot)
    dev.off()
  }

  if (nset>1 & nset<=5) {
    data_id1 <- plyr::dlply(data, .(condition, replicate), function(x) n=unique(x$id))
    drawvennplot(data_id1, nset, foldername, "Quantified_IDs_run.png")
  }
  if (ncondition > 1 & ncondition<=5) {
    data_id2 <- plyr::dlply(data, .(condition), function(x) condition=unique(x$id))
    drawvennplot(data_id2, ncondition, foldername, "Quantified_IDs_sample.png")
  }

  d1 <- tidyr::gather(data[,c(1,3,5:(4+nread))], treatment, reading, -id, -condition)
  d1$treatment <- factor(d1$treatment, levels=sort(as.numeric(unique(d1$treatment)), decreasing=FALSE))
  #d1 <- melt(data[,c(1,3,5:(4+nread))], id.vars=c("id", "condition"), variable.name="treatment", value.name="reading");
  q <- ggplot(d1, aes(x = treatment, y = reading)) +
    coord_cartesian(ylim = c(0,2)) +
    scale_y_continuous(breaks=c(0,0.8,1,1.25,2))
  q <- q + labs(x="Treatment", y="Normalized Ratio") +
    geom_boxplot(aes(fill=condition)) +
    geom_hline(yintercept=c(0.8,1.25), linetype = "longdash") +
    theme_bw() + theme(text = element_text(size=12),
                       axis.text.x = element_text(angle = 45,hjust = 1),
                       aspect.ratio=1)
  ggsave(filename = paste0(outdir,"/",foldername,"/","wholeset_trend.png"), q, height=8, width=8)


  data$STD <- apply(data[ ,c(5:(4+nread))], 1, sd)
  data <- mutate(data, CV=STD/rowMeans(data[,c(5:(4+nread))]))
  limit <- quantile(data$CV, 0.999)
  m <- ggplot(data, aes(x=CV, fill=condition)) + coord_cartesian(xlim = c(0,limit))
  m <- m + geom_histogram(binwidth = limit/100, alpha=0.5, position="identity") +
    labs(x="Reading coefficient of variation (CV)") + theme_bw() +
    theme(text = element_text(size=12), aspect.ratio=1)
  ggsave(filename = paste0(outdir,"/",foldername,"/","CV_distribution.png"), m, height=8, width=8)

  data$AUC <- rowSums(data[ ,c(5:(4+nread))])
  limit <- quantile(data$AUC, 0.999)
  m <- ggplot(data, aes(x=AUC, fill=condition)) + coord_cartesian(xlim = c(0,limit))
  m <- m + geom_histogram(binwidth = limit/100, alpha=0.5, position="identity") +
    labs(x="Area under the curve (AUC)") + theme_bw() +
    theme(text = element_text(size=12), aspect.ratio=1)
  ggsave(filename = paste0(outdir,"/",foldername,"/","AUC_distribution.png"), m, height=8, width=8)

  limity <- quantile(data$CV, 0.999)
  limitx <- quantile(data$AUC, 0.999)
  n1 <- ggplot(data, aes(x=AUC, y=CV)) + geom_point(aes(colour=condition),size=0.5,alpha=0.3) +
    coord_cartesian(xlim = c(0,limitx), ylim = c(0,limity))
  n1 <- n1 + facet_grid(.~condition+replicate) +
    labs(x="Reading Area Under the Curve",y="Reading coefficient of variation") +
    theme_bw() + theme(
      text = element_text(size=12),
      axis.text=element_text(size=6),
      legend.position="bottom", aspect.ratio=1)
  ggsave(filename = paste0(outdir,"/",foldername,"/","AUC_CV_distribution1.png"), n1, height=4, width=8)

  if (isdatafitted) {
    n2 <- ggplot(data, aes(x=R2, y=CV)) + geom_point(aes(colour=condition),size=0.5,alpha=0.3) +
      geom_vline(xintercept=0.8, colour="blue", linetype = "longdash") +
      coord_cartesian(ylim = c(0,limity))
    n2 <- n2 + labs(x="Reading R-squared",y="Reading coefficient of variation") +
      facet_grid( .~condition+replicate) + theme_bw() +
      theme(text = element_text(size=12), axis.text=element_text(size=6),
            legend.position="bottom", aspect.ratio=1)
    ggsave(filename = paste0(outdir,"/",foldername,"/","R2_CV_distribution1.png"), n2, height=4, width=8)

    a<-data.frame(t(as.matrix(summary(data$CV))), stringsAsFactors=F, check.names=F)
    a$parameter <- "CV"
    b<-data.frame(t(as.matrix(summary(data$AUC))), stringsAsFactors=F, check.names=F)
    b$parameter <- "AUC"
    c<-data.frame(t(as.matrix(summary(data$R2))), stringsAsFactors=F, check.names=F)
    c$parameter <- "R2"
    R2CV<-rbind(a,b)
    R2CV<-rbind(R2CV,c)
    R2CV<-R2CV[ ,c(7,1:6)]
  } else {
    a<-data.frame(t(as.matrix(summary(data$CV))), stringsAsFactors=F, check.names=F)
    a$parameter <- "CV"
    b<-data.frame(t(as.matrix(summary(data$AUC))), stringsAsFactors=F, check.names=F)
    b$parameter <- "AUC"
    R2CV<-rbind(a,b)
    R2CV<-R2CV[ ,c(7,1:6)]
  }
  #print(R2CV)

  if (PSMcutoff) {# To select the proteins with more than 3 PSM (average)
    PSMcol <- grep("PSM", names(data), value=FALSE)
    PSMname <- names(data)[PSMcol]
    names(data)[PSMcol] <- "PSM"
    PSMtable <- plyr::ddply(data, .(id), summarize, PSMmean=mean(PSM))

    PSMplot <- ggplot(PSMtable, aes(x = PSMmean)) + geom_histogram(fill="darkgray") + geom_vline(xintercept=3, colour="blue", linetype = "longdash")
    PSMplot <- PSMplot + labs(x="Average PSM number for each protein") + scale_x_log10() +
      theme_bw() + theme(text = element_text(size=12), axis.text.x = element_text(angle = 45,hjust = 1), aspect.ratio=1);
    ggsave(filename = paste0(outdir,"/",foldername,"/","PSM_distribution.png"), PSMplot, height=8, width=8);

    PSMtable <- subset(PSMtable, PSMmean>3)
    nkeep <- PSMtable$id
    fkeep <- which(data$id %in% nkeep)
    names(data)[PSMcol] <- PSMname
    data_PSMsmall <- data[-fkeep, ]
    data_PSMbig <- data[fkeep, ]
  }

  data_N2 <- plyr::ddply(data_PSMbig, .(condition), summarise, Number_sample=length(unique(id)))
  data_N_PSMbig <- merge(data_N1, data_N2)
  limity <- quantile(data_PSMbig$CV, 0.999)
  limitx <- quantile(data_PSMbig$AUC, 0.999)
  n1 <- ggplot(data_PSMbig, aes(x=AUC, y=CV)) +
    geom_point(aes(colour=condition),size=0.5,alpha=0.3) +
    coord_cartesian(xlim = c(0,limitx), ylim = c(0,limity))
  n1 <- n1 + labs(x="Reading Area Under the Curve",y="Reading coefficient of variation") +
    facet_grid( .~condition+replicate) + theme_bw() +
    theme(text = element_text(size=12), axis.text=element_text(size=6),
          legend.position="bottom", aspect.ratio=1)
  ggsave(filename = paste0(outdir,"/",foldername,"/","AUC_CV_distribution2.png"), n1, height=4, width=8)

  if (isdatafitted) {
    n2 <- ggplot(data_PSMbig, aes(x=R2, y=CV)) +
      geom_point(aes(colour=condition),size=0.5,alpha=0.3) +
      geom_vline(xintercept=0.8, colour="blue", linetype = "longdash") +
      coord_cartesian(ylim = c(0,limity))
    n2 <- n2 + labs(x="Reading R-squared",y="Reading coefficient of variation") +
      facet_grid(.~condition+replicate) + theme_bw() +
      theme(text = element_text(size=12), axis.text=element_text(size=6),
            legend.position="bottom", aspect.ratio=1)
    ggsave(filename = paste0(outdir,"/",foldername,"/","R2_CV_distribution2.png"), n2, height=4, width=8)

    a<-data.frame(t(as.matrix(summary(data_PSMbig$CV))), stringsAsFactors=F, check.names=F)
    a$parameter <- "CV"
    b<-data.frame(t(as.matrix(summary(data_PSMbig$AUC))), stringsAsFactors=F, check.names=F)
    b$parameter <- "AUC"
    c<-data.frame(t(as.matrix(summary(data_PSMbig$R2))), stringsAsFactors=F, check.names=F)
    c$parameter <- "R2"
    R2CV_PSMbig<-rbind(a,b)
    R2CV_PSMbig<-rbind(R2CV_PSMbig,c)
    R2CV_PSMbig<-R2CV_PSMbig[ ,c(7,1:6)]
  }else{
    a<-data.frame(t(as.matrix(summary(data_PSMbig$CV))), stringsAsFactors=F, check.names=F)
    a$parameter <- "CV"
    b<-data.frame(t(as.matrix(summary(data_PSMbig$AUC))), stringsAsFactors=F, check.names=F)
    b$parameter <- "AUC"
    R2CV_PSMbig<-rbind(a,b)
    R2CV_PSMbig<-R2CV_PSMbig[ ,c(7,1:6)]
  }
  #print(R2CV_PSMbig)

  r <- newCustomReport( reportname )
  ss1 <- newSection( "Quantified Protein Numbers" )
  ss2 <- newSection( "Quantified Protein Overlapping" )
  ss3 <- newSection( "The distributions of readings (grouped by 10plex-curve)" )
  ss4 <- newSection( "The distributions of readings (grouped by 10plex-curve) in PSM>3 group" )

  if (nset>1 & nset<=5) {
    fig1 <- newFigure( "Quantified_IDs_run.png", "The overlapping of proteins from each run.")
  } else {
    fig1 <- NULL
  }
  if(ncondition>1 & ncondition<=5) {
    fig2 <- newFigure( "Quantified_IDs_sample.png", "The overlapping of proteins from each sample.")
  } else {
    fig2 <- NULL
  }
  fig3 <- newFigure( "wholeset_trend.png", "The distribution of readings in each channel.")
  fig4 <- newFigure( "AUC_distribution.png", "The distribution of AUCs of readings per curve in each condition.")
  fig5 <- newFigure( "CV_distribution.png", "The distribution of CVs of readings per curve in each condition.")

  fig6 <- newFigure( "AUC_CV_distribution1.png", "The distribution of AUC vs CV of readings per curve in each condition.")
  fig7 <- newFigure( "AUC_CV_distribution2.png", "The distribution of AUC vs CV of readings per curve in each condition from PSM>3 group.")
  if (isdatafitted) {
    fig8 <- newFigure( "R2_CV_distribution1.png", "The distribution of R2 vs CV of readings per curve in each condition.")
    fig9 <- newFigure( "R2_CV_distribution2.png", "The distribution of R2 vs CV of readings per curve in each condition from PSM>3 group.")
  } else {
    fig8 <- NULL
    fig9 <- NULL
  }
  fig10 <- newFigure( "PSM_distribution.png", "The distribution of PSMs for each quantified protein")

  t1 <- newTable( data_N, paste0("The number of quantified proteins for in total ", length(unique(data$id))), " unique proteins." )
  t2 <- newTable( R2CV, "The distribution of CV and R2 of readings (grouped by 10plex-curve)." )
  t3 <- newTable( data_N_PSMbig, paste0("The number of quantified proteins in PSM>3 group for in total ", length(unique(data_PSMbig$id))), " unique proteins." )
  t4 <- newTable( R2CV_PSMbig, "The distribution of CV and R2 of readings (grouped by 10plex-curve) in PSM>3 group." )
  p1 <- newParagraph( "The coefficient of variation (CV) is defined as the ratio of standard deviation to mean of ten readings in one TMT10 set defined curve.")
  p2 <- newParagraph( "Note that to prevent the extreme values compressing the plot, only up to 99.9% percentile of CV distribution in the data was plotted." )
  p3 <- newParagraph( "PSM>3 group is the subset of quantified proteins with on average more than 3 Peptide Spectrum Match (PSM) support." );

  ss1 <- addTo( ss1, t1)
  if(ncondition>1){
    ss2 <- addTo( ss2, fig1, fig2 ) # parent, child_1, ..., child_n
  }else{
    ss2 <- addTo( ss2, fig1)
  }
  if(isdatafitted){
    ss3 <- addTo( ss3, t2, p1, fig3, fig4, fig5, fig6, fig8, p2 )
    ss4 <- addTo( ss4, p3, fig10, t3, t4, fig7, fig9, p2 )
  }else{
    ss3 <- addTo( ss3, t2, p1, fig3, fig4, fig5, fig6, p2 )
    ss4 <- addTo( ss4, p3, fig10, t3, t4, fig7, p2 )
  }
  r <- addTo( r, ss1, ss2, ss3, ss4 )
  writeReport( r, filename=paste0(outdir,"/",foldername, "/", reportname )) # w/o extension
  rm(p1, p2, p3, r, ss1, ss2, ss3, ss4, t1, t2, t3, t4, fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10)

}
