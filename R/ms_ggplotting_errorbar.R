#' ms_ggplotting_errorbar
#'
#' Function to generate pdf files with multipanel ggplots for replicated
#' (ideally at least three) CETSA melt curve data, in the format of a mean +/- errorbar
#'
#' @param data dataset to plot, ideally with the fitting information but not complusory
#' @param extradata dataset with extrapolated readings, default to NULL
#' @param legenddata dataset used for condition extraction, at least one protein
#' inside this dataset should contains the full set of experiment conditions,
#' most of the time there is no need to supply, just use the data instead
#' @param levelvector a vector of experimental conditions
#' @param nread number of reading channels or sample treatements
#' @param minireplicate ideally at least three replicates of measurments for the
#' same proteins are required to plot in a box plot format
#' @param witherrorbar whether to plot in a mean +/- error bar(se) graph format, default to TRUE
#' @param fitremout whether to segregate the proteins with messy melt curves, default set to FALSE
#' @param bottomcutoff the average of the last three points should be lower than
#' specified bottom cutoff value, which is 0.4 by default
#' @param topcutoff the average of the first three points should be higher than
#' specified bottom cutoff value, which is 0.8 by default
#' @param topasone whether the top plateau has to be fixed, i.e., 1.0
#' @param keepcommonprotein whether to only keep the proteins also present in the
#' provided extradata set, default set to TRUE
#' @param printTm whether to annotate the plots with Tm values, however in this version,
#' when there is more than one sample or no replicates, the Tm annotation would be switched off
#' @param annotypos the starting y-axis position of textual annotation, default value 0.5
#' @param annotyinterval the interval on y-axis for textual annotation, default valule 0.08
#' @param layout a vector indicating the panel layout for multi-panel plots per page
#' @param withpoints whether to fit curve with the reading points, default set to TRUE
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
#' @keywords internal
#'
#' @return a list of ggplot2 object
#' @examples \dontrun{
#' ms_ggplotting_errorbar(LY_fitted)
#' }
#'
#'

ms_ggplotting_errorbar <- function(data, extradata=NULL, legenddata=NULL, levelvector=NULL,
                              nread=10, minireplicate=NULL, witherrorbar=TRUE,
                              fitremout=FALSE, bottomcutoff=0.4, topcutoff=0.8,
                              topasone=TRUE, pfdatabase=FALSE, keepcommonprotein=TRUE,
                              printBothName=TRUE, printGeneName=FALSE, printcount=TRUE,
                              printTm=TRUE, annotypos=0.5, annotyinterval=0.08,
                              layout=c(5,5), external=TRUE, withpoints=TRUE,
                              presetcolor=TRUE, colorpanel=NULL, nshape=5, nsize=2,
                              toplabel="CETSA data plotting_error bar",
                              leftlabel="Non-denatured protein fraction",
                              bottomlabel="Temperature", returnplots=FALSE,
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
  if (any(duplicated(data[ ,checkpos]))) {
    cat("Warning for duplicated protein name entries\n")
    cat("Double check the following proteins for duplicated entries:\n")
    cat(data[duplicated(data[ ,checkpos]), ]$id,"in",data[duplicated(data[ ,checkpos]), ]$condition,"\n")
    stop("1.Remove the duplicated entires from original dataset then start again!")
  }

  # to concatenate id and description
  nrowdata <- nrow(data)
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

  if (class(extradata)!="NULL") {
    if (printBothName & !pfdatabase) {
      extradata <- extradata %>% rowwise() %>% mutate(description1 = getProteinName(description, pfdatabase)) %>%
        mutate(description2 = getGeneName(description)) %>%
        mutate(id = paste(id, description1, description2, sep="\n"))
      extradata$description1<-NULL
      extradata$description2<-NULL
    } else if (printGeneName & !pfdatabase) {
      extradata <- extradata %>% rowwise() %>%
        mutate(description = getGeneName(description)) %>%
        mutate(id = paste(id, description, sep="\n"))
    } else {
      extradata <- extradata %>% rowwise() %>%
        mutate(description = getProteinName(description, pfdatabase)) %>%
        mutate(id = paste(id, description, sep="\n"))
    }
    extradata$description<-NULL
    if (length(grep("pvalue", names(extradata)))) {
      # should not use the dataframe with readings difference and deltaAUC
      extradata <- extradata[ ,!(names(extradata) %in% c("deltaAUC","deltaAUC.z","deltaAUC.p",
                                                         "t","pvalue","fdr"))]
    }

    # only keep the commonly presented proteins
    if (keepcommonprotein) {
      commonid <- intersect(extradata$id, data$id)
      extradata <- subset(extradata, id %in% commonid)
      data <- subset(data, id %in% commonid)
    }
  }

  checkpos <- c(1,2)
  if (any(duplicated(data[ ,checkpos]))) {
    cat("Warning for duplicated protein name entries\n")
    cat("Double check the following proteins for duplicated entries:\n")
    cat(data[duplicated(data[ ,checkpos]), ]$id,"in",data[duplicated(data[ ,checkpos]), ]$condition,"\n")
    stop("Remove the duplicated entires from original dataset then start again!")
  }

  # look for outlier proteins based on melting behavior in controls
  outliers <- NULL
  if (fitremout) {
    cat("Make sure you provide fitted data with Tm and R2 values for this option!\n")
    data$Top <- rowMeans(data[ ,c(3:5)])
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

  # remove proteins sampled with too few replicates if specified
  if (length(minireplicate)>0) {
    counttable <- data %>% group_by(id) %>% summarize(count=n()) %>% filter(count>=minireplicate)
    fkeep <- which(data$id %in% counttable$id)
    data <- data[fkeep, ]
    if ( nrow(data)==0 ) {
      stop("Make sure there are enough experimental conditions in dataset.")
    } else {
      cat(length(unique(data$id)), "proteins pass the measurement replicates cutoff.\n")
    }
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
  if (length(levelvector)) {
    cdata$sample <- factor(cdata$sample, levels=levelvector)
  } else {
    cdata$sample <- factor(cdata$sample)
  }

  if (class(extradata)!="NULL") {
    extradata_l <- tidyr::gather(extradata, temperature, reading, -id, -condition)
    extradata_l$temperature <- gsub("C", "", extradata_l$temperature)
    extradata_l$temperature <- as.numeric(extradata_l$temperature)
  }

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
  multirep <- length(unique(data1$rep))>1
  data1$rep <- NULL
  pos <- which(names(data1) %in% c("Tm","R2","Slope","RSE"))
  errortable <- NA
  if (length(pos) > 1 & multirep) {
    Tmtable <- data1[ ,c(1,2,pos)] %>% group_by(id, sample) %>%
      summarize(repnumber=n(),
                Tm_mean=mean(Tm, na.rm=T),
                Tm_sd=sd(Tm, na.rm=T),
                Tm_se=Tm_sd/sqrt(repnumber),
                R2_median=median(R2, na.rm=T),
                Slope_median=median(Slope, na.rm=T))

    errortable <- cdata %>% group_by(id, sample) %>% summarize(sum_error=sqrt(sum(sd^2)))
    errortable <- merge(Tmtable, errortable)
    errortable <- mutate(errortable, sum_error=sum_error*exp(R2_median))
    errortable <- mutate(errortable, sum_error=sum_error*Tm_sd)
    errortable <- errortable[order(errortable$sum_error, decreasing=FALSE, na.last=TRUE), ]

    errortable1 <- tidyr::separate(errortable, id, into=c("id", "name"), sep="\n")
    errortable1 <- merge(proteininfo, errortable1[ ,-2])
    ms_filewrite(errortable1, paste0(dataname, "_error_metrics.txt"), outdir=outdir)
  }

  if (printcount) {
    pos <- which(names(data1)=="countNum")
    counttable <- data1[ ,c(1,2,pos)] %>% group_by(id, sample) %>%
      summarize(countnumber=round(median(countNum, na.rm=T),1))
    counttable <- spread(counttable, sample, countnumber)
  }

  if (is.na(errortable)) {
    cdata <- tidyr::complete(cdata, id, sample, temperature)
    cdata$id <- factor(cdata$id)
  } else { # rank by errortable metrics
    if (length(unique(data1$sample))==2) {
      Tm_diff <- tidyr::spread(errortable[,c(1,2,4)],sample,Tm_mean)
      Tm_diff <- mutate(Tm_diff, difference=eval(parse(text=names(Tm_diff)[2]))-eval(parse(text=names(Tm_diff)[3])))
      Tm_diff <- Tm_diff[order(Tm_diff$difference, decreasing=FALSE, na.last=TRUE), ]
      seq <- unique(as.factor(Tm_diff$id))
    } else {
      seq <- unique(as.factor(errortable$id))
    }
    cdata <- tidyr::complete(cdata, id, sample, temperature)
    cdata$id <- factor(cdata$id, levels=seq)
  }
  # print(head(cdata))
  if (external & !returnplots) { external_graphs(T) }

  y <- c(1.0, 0.98, 0.95, 0.8, 0.7, 0.5, 0.5, 0.2, 0.15, 0.1)
  x <- c(37.0, 40.0, 43.0, 46.0, 49.0, 52.0, 55.0, 58.0, 61.0, 64.0)
  df <- data.frame(x, y)
  t1 <- drc::drm(y~x, data=df, fct=drc::LL.4())
  #print(t1)

  if (topasone==TRUE) { top <- 1.0 } else { top <- NA }

  plotting <- function(d1) {

    if (printcount) {
      countlabel <- counttable[counttable$id==unique(d1$id), ][ ,c(2:ncol(counttable))]
    }

    d1$temperature <- as.numeric(as.character(d1$temperature))
    xmin <- min(d1$temperature)
    xmax <- max(d1$temperature)
    xpos <- (xmax-xmin)/4 + xmin

    q <- ggplot(d1, aes(x=temperature, y=mean, group=sample, colour=sample))
    if (withpoints) {
      q <- q + geom_point(aes(colour=sample), shape=20, size=2)
    }
    q <- q + geom_smooth(method=drc::"drm", formula=t1, method.args=list(fct=drc::LL.4(fixed=c(NA,NA,top,NA))), se=FALSE, size=0.5) +
      scale_x_continuous(breaks=sort(unique(d1$temperature), decreasing=FALSE)) +
      coord_cartesian(xlim=c(xmin,xmax), ylim = c(0,1.2)) +
      scale_y_continuous(breaks=seq(0,1,0.2)) + labs(x=" ", y=" ") +
      scale_colour_manual(drop=FALSE, values=colorpanel) +
      ggtitle(as.character(unique(d1$id)))

    if (witherrorbar) { q <- q + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=2,
                                               position=position_dodge(0.05)) }

    if (class(extradata)!="NULL") {
      data1 <- subset(extradata_l, id==unique(d1$id))
      # print(head(data1))
      if (nrow(data1)) {
        q <- q + geom_point(data=data1, aes(x=temperature, y=reading, group=condition, colour=condition), shape=nshape, size=nsize)
      }
    }

    if (printTm & !is.na(errortable)) {
      Tm_label <- errortable[errortable$id==unique(d1$id), ]
      q <- q + annotate("text", x=xpos, y=annotypos, label=paste("No. of replicate:", unique(d1$N), sep=" "), size=2)
      q <- q + annotate("text", x=xpos, y=annotypos-annotyinterval, label=paste0("Tm: ", round(Tm_label$Tm_mean,1),
                                                                                 "+/-", round(Tm_label$Tm_se,1)), size=2)
      q <- q + annotate("text", x=xpos, y=annotypos-2*annotyinterval, label=paste0("R2 median: ", round(Tm_label$R2_median,2)), size=2)
    }

    if (printcount) {
      q <- q + annotate("text", x=xpos, y=annotypos-0.1, label="#qPSM", size=2)
      for (i in 1:length(countlabel)) {
        q <- q + annotate("text", x=xpos, y=annotypos-0.1-annotyinterval*i,
                          label=paste0(names(countlabel)[i], ":  ", as.character(countlabel)[i], " "),
                          size=2, hjust=1)
        i <- i + 1
      }
    }

    q <- q + theme_classic() +
      theme(
        text = element_text(size=10),
        strip.text.x = element_text(size=5),
        plot.title = element_text(hjust=0.5, size=rel(0.7)),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        aspect.ratio = 1
      )
    return(q)
  }
  plots <- plyr::dlply(cdata, plyr::.(id), .fun=plotting)
  #plots[[1]]
  if (returnplots) {
    if (external) { external_graphs(F) } # switch off the external graphs
    return(plots)
  }
  plotlegend <- ms_isothermal_legend2(legenddata, levelvector, nread,
                                      nsample, presetcolor, colorpanel)
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
  ggsave(filename=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),
                     length(unique(cdata$id)), "_", pdfname),
                     pl, height=pdfheight, width=pdfwidth)

  if (external) { external_graphs(F) } # switch off the external graphs
  message("CETSA melt curve errorbar plot file generated successfully.")
}
