#' ms_ggplotting_errorbar
#'
#' Function to generate pdf files with multipanel ggplots for replicated
#' (ideally at least three) CETSA melt curve data, in the format of a mean +/- errorbar
#'
#' @param data isothermal dataset to plot
#' @param legenddata dataset used for condition extraction, at least one protein
#' inside this dataset should contains the full set of experiment conditions,
#' most of the time there is no need to supply, just use the data instead
#' @param levelvector a vector of experimental conditions, not complusory for isothermal functions
#' @param nread number of reading channels or sample treatements
#' @param rep the replicate indicator used in the naming of experimental conditions,
#' such as "r", or "rep" or ""
#' @param minireplicate ideally at least three replicates of measurments for the
#' same proteins are required to plot in a box plot format
#' @param fitremout whether to segregate the proteins with messy melt curves
#' @param bottomcutoff the average of the last three points should be lower than
#' specified bottom cutoff value, which is 0.4 by default
#' @param topcutoff the average of the first three points should be higher than
#' specified bottom cutoff value, which is 0.8 by default
#' @param topasone whether the top plateau has to be fixed, i.e., 1.0
#' @param printTm whether to annotate the plots with Tm values, however in this version,
#' when there are one samples, the Tm annotation would be switched off
#' @param annotypos the starting y-axis position of textual annotation, default value 0.5
#' @param annotyinterval the interval on y-axis for textual annotation, default valule 0.08
#' @param layout a vector indicating the panel layout for multi-panel plots per page
#' @param presetcolor whether to use the pre-defined color scheme
#' @param colorpanel a vector of customizable color scheme provided by the user
#' @param toplabel textual label at the top part of the page
#' @param leftlabel textual label at the left side of the page
#' @param bottomlabel textual label at the bottom part of the page
#' @param pdfname name for the pdf plots file
#' @param pdfheight a number indicate the height of pdf file, default value 12
#' @param pdfwidth a number indicate the width of pdf file, default value 12
#'
#'
#' @import tidyr
#' @import RColorBrewer
#' @import scales
#' @import drc
#' @importFrom plyr . dlply ddply
#' @importFrom grid unit.c
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @import ggplot2
#'
#' @export
#'
#' @return a list of ggplot2 object
#' @examples \dontrun{
#' ms_ggplotting_errorbar(LY_fitted)
#' }
#'
#'

ms_ggplotting_errorbar <- function(data, legenddata=NULL, levelvector=NULL,
                              nread=10, rep="r", minireplicate=NULL,
                              fitremout=FALSE, bottomcutoff=0.4, topcutoff=0.8,
                              topasone=TRUE, pfdatabase=FALSE,
                              printBothName=TRUE, printGeneName=FALSE,
                              printTm=TRUE, annotypos=0.5, annotyinterval=0.08,
                              layout=c(5,5), external=TRUE,
                              presetcolor=TRUE, colorpanel=NULL,
                              toplabel="CETSA data plotting_barplot",
                              leftlabel="Non-denatured protein fraction",
                              bottomlabel="Temperature",
                              pdfname="ggplotting.pdf", pdfheight=12, pdfwidth=12) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
  }

  png(filename=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                      dataname, "_curve replicate distribution.png"))
  barplot(table(table(data$id)),
          main="Melting curve replicate distribution",
          xlab="Number of curves per unique protein in replicate exp")
  dev.off()

  checkpos <- c(1,3)
  if (any(duplicated(data[, checkpos]))) {
    print("Warning for duplicated protein name entries")
    print("Double check the following proteins for duplicated entries:")
    print(paste0(data[duplicated(data[ ,checkpos]), ]$id," in ",
                 data[duplicated(data[ ,checkpos]), ]$condition))
    stop("1.Remove the duplicated entires from original dataset then start again!")
  }

  # to concatenate id and description
  nrowdata <- nrow(data)
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

  checkpos <- c(1,2)
  if (any(duplicated(data[, checkpos]))) {
    print("Warning for duplicated protein name entries")
    print("Double check the following proteins for duplicated entries:")
    print(paste0(data[duplicated(data[ ,checkpos]), ]$id," in ",
                 data[duplicated(data[ ,checkpos]), ]$condition))
    stop("Remove the duplicated entires from original dataset then start again!")
  }

  # look for outlier proteins based on melting behavior in controls
  outliers <- NULL
  if (fitremout) {
    print("Make sure you provide fitted data with Tm and R2 values for this option!")
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
      data <- data[-outlierid1, ]
    }
  }

  # remove proteins sampled with too few replicates
  if (length(minireplicate)>0) {
    counttable <- data %>% group_by(id) %>% summarize(count=n()) %>% filter(count>=minireplicate)
    fkeep <- which(data$id %in% counttable$id)
    data <- data[fkeep, ]
  }

  if ( nrow(data)==0 ) {
    stop("Make sure there are enough experimental conditions in dataset.")
  } else {
    print(paste0(length(unique(data$id)), " proteins pass the measurement replicates cutoff."))
  }

  if (length(legenddata)==0) { legenddata <- data }

  data_l <- tidyr::gather(data[ ,c(1:(nread+2))], temperature, reading, -id, -condition)
  data_l <- tidyr::separate(data_l, condition, into=c("sample","rep"), sep="\\.")
  cdata <- plyr::ddply(data_l, c("id", "sample", "temperature"), summarise,
                       N    = length(reading),
                       mean = mean(reading),
                       sd   = sd(reading),
                       se   = sd / sqrt(N)
  )

  nsample <- length(unique(cdata$sample))
  if (presetcolor & length(colorpanel)==0) {
    if (nsample<=9) {colorpanel <- brewer.pal(9, "Set1")}
    else {stop("The number of conditions in dataset exceeds the preset number of colors, pls provide a vector of colors in colorpanel")}
  } else if (length(colorpanel) < nsample){
    stop("The number of conditions in dataset exceeds the provided number of colors, pls provide enough colors in colorpanel")
  }
  # when there is more than one sample, don't print out Tm annotation
  if (nsample > 1) { printTm <- FALSE }

  data1 <- tidyr::separate(data, condition, into=c("sample","rep"), sep="\\.")
  data1$rep <- NULL
  Tmpos <- grep("^Tm", names(data1), value=F)
  R2pos <- grep("^R2", names(data1), value=F)
  Slopepos <- grep("^Slope", names(data1), value=F)
  RSEpos <- grep("^RSE", names(data1), value=F)
  pos <- c(Tmpos, R2pos, Slopepos, RSEpos)
  if (length(pos) > 1) {
    Tmtable <- data1[ ,c(1,2,pos)] %>% group_by(id, sample) %>%
      summarize(repnumber=n(),
                Tm_mean=mean(Tm, na.rm=T),
                Tm_sd=sd(Tm, na.rm=T),
                Tm_se=Tm_sd/sqrt(repnumber),
                R2_median=median(R2, na.rm=T),
                Slope_median=median(Slope, na.rm=T))
  }
  errortable <- cdata %>% group_by(id, sample) %>% summarize(sum_error=sqrt(sum(sd^2)))
  errortable <- merge(Tmtable, errortable)
  errortable <- mutate(errortable, sum_error=sum_error/R2_median)
  errortable <- errortable[order(errortable$sum_error, decreasing=FALSE), ]

  errortable1 <- tidyr::separate(errortable, id, into=c("id", "name"), sep="\n")
  errortable1 <- merge(proteininfo, errortable1[ ,-2])
  ms_filewrite(errortable1, paste0(dataname, "_error_metrics.txt"), outdir=outdir)

  seq <- as.factor(errortable$id)
  cdata <- tidyr::complete(cdata, id, sample, temperature)
  cdata$id <- factor(cdata$id, levels=seq)

  if (external) { external_graphs(T) }

  y <- c(1.0, 0.98, 0.95, 0.8, 0.7, 0.5, 0.5, 0.2, 0.15, 0.1)
  x <- c(37.0, 40.0, 43.0, 46.0, 49.0, 52.0, 55.0, 58.0, 61.0, 64.0)
  df <- data.frame(x, y)
  t1 <- drc::drm(y~x, data=df, fct=LL.4())
  #print(t1)

  if (topasone==TRUE) { top <- 1.0 } else { top <- NA }

  plotting <- function(d1, printTm=TRUE) {

    d1$temperature <- as.numeric(as.character(d1$temperature))
    xmin <- min(d1$temperature)
    xmax <- max(d1$temperature)
    xpos <- (xmax-xmin)/4 + xmin

    q <- ggplot(d1, aes(x=temperature, y=mean, group=sample, colour=sample)) +
      geom_point(aes(colour=sample), shape=20, size=2) +
      geom_smooth(method=drc::"drm", formula=t1, method.args=list(fct=drc::LL.4(fixed=c(NA,NA,top,NA))), se=FALSE, size=0.5) +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(0.05)) +
      scale_x_continuous(breaks=sort(unique(d1$temperature), decreasing=FALSE)) +
      coord_cartesian(xlim=c(xmin,xmax), ylim = c(0,1.2)) +
      scale_y_continuous(breaks=seq(0,1,0.2)) + labs(x=" ", y=" ") +
      ggtitle(as.character(unique(d1$id)))

    if (printTm) {
      Tm_label <- errortable[errortable$id==unique(d1$id), ]
      q <- q + annotate("text", x=xpos, y=annotypos, label=paste("No. of replicate:", unique(d1$N), sep=" "), size=2)
      q <- q + annotate("text", x=xpos, y=annotypos-annotyinterval, label=paste0("Tm: ", round(Tm_label$Tm_mean,1), "+/-", round(Tm_label$Tm_se,1)), size=2)
      q <- q + annotate("text", x=xpos, y=annotypos-2*annotyinterval, label=paste0("R2 median: ", round(Tm_label$R2_median,2)), size=2)
    }

    q <- q + theme_classic() +
      theme(
        text = element_text(size=10),
        strip.text.x = element_text(size = 5),
        plot.title = element_text(hjust = 0.5, size = rel(0.7)),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position="none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        aspect.ratio = 1
      )
    return(q)
  }
  plots <- plyr::dlply(cdata, plyr::.(id), .fun=plotting, printTm=printTm)
  #plots[[1]]
  plotlegend <- ms_isothermal_legend2(legenddata, levelvector,
                                      nread, nsample, presetcolor, colorpanel)
  legend <- plotlegend$legend
  lheight <- plotlegend$lheight

  params <- list(nrow=layout[1], ncol=layout[2])
  n <- with(params, nrow*ncol)
  ## add one page if division is not complete
  pages <- length(plots) %/% n + as.logical(length(plots) %% n)
  groups <- split(seq_along(plots), gl(pages, n, length(plots)))

  pl <- lapply(names(groups), function(i){
    gridExtra::grid.arrange(
      do.call(arrangeGrob,
              c(plots[groups[[i]]], params, top=toplabel,
                left=leftlabel, bottom=bottomlabel)),
      legend,
      ncol = 1,
      heights = grid::unit.c(unit(1, "npc") - lheight, lheight))
  })

  class(pl) <- c("arrangelist", "ggplot", class(pl))
  ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                     length(unique(cdata$id)), pdfname),
                     pl, height=pdfheight, width=pdfwidth)

  if (external) { external_graphs(F) } # switch off the external graphs
  print("CETSA melt curve errorbar plot file generated successfully.")
}
