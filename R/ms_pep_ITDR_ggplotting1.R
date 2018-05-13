#' ms_pep_ITDR_ggplotting1
#'
#' Function to generate pdf files of peptide level ggplots for isothermal CETSA data


#' @param data peptide level ITDR dataset to plot, for easy readout, it is suggested
#' to merge with protein level (Uniprot) id-description pairs beforehand
#' @param nread number of reading channels or sample treatements, default value 10
#' @param setminimalrep whether to set a criteria of minimal replicates number, defalut set to FALSE
#' @param minimalrep the number of a mininal replicates for the annotated sequence to be plot
#' @param unit textual annotation for the dose unit, default is "mM"
#' @param colorpanel a vector of customizable color scheme provided by the user
#' @param layout a vector indicating the panel layout for multi-panel plots
#' per page, default value is c(5,5)
#'
#' @import scales
#' @importFrom plyr . dlply
#' @importFrom grid unit.c
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @import tidyr
#' @import dplyr
#' @import RColorBrewer
#' @import ggplot2
#'
#' @export
#'
#' @return a list of ggplot2 object
#' @examples \dontrun{
#' ms_pep_ITDR_ggplotting1(ITDR_pep_data_with_desciption)
#' }
#'
#'

ms_pep_ITDR_ggplotting1 <- function(data, legenddata=NULL, nread=10, pfdatabase=FALSE, printGeneName=FALSE,
                                    setminimalrep=FALSE, minimalrep=1,nreplicate=2,
                                    orderEC=FALSE, orderAUC=FALSE, plotseq=NULL,
                                    loess=FALSE, dotconnect=FALSE,
                                    xlinear=FALSE, xlog10=TRUE, xsqrt=FALSE, xcubert=FALSE,
                                    xinterval=NULL, fixedy=FALSE, unit="mM", external=TRUE,
                                    presetcolor=TRUE, colorpanel=NULL, layout=c(5,5),
                                    top_label="CETSA peptide level plotting", pdfname="_ggplotting.pdf") {


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

  getGeneName <- function(x) {return (strsplit(strsplit(x, "GN=")[[1]][2], " ")[[1]][1])}
  getProteinName <- function(x) {return (strsplit(x, " OS=")[[1]][1])}
  if (pfdatabase) {
    getProteinName <- function(x) {return (gsub("product=", "", strsplit(x, "\\|")[[1]][2]))}
  }
  if (printGeneName) {
    data <- data %>% rowwise() %>% mutate(description = getGeneName(description)) %>%
      mutate(MasterProteinAccessions = paste(MasterProteinAccessions, description, sep="\n")) %>%
      mutate(AnnotatedSequence = paste(AnnotatedSequence, PhosphorylatedAA, sep="-")) %>%
      mutate(AnnotatedSequence = paste(AnnotatedSequence, MasterProteinAccessions, sep="\n"))
  } else {
    data <- data %>% rowwise() %>% mutate(description = getProteinName(description)) %>%
      mutate(MasterProteinAccessions = paste(MasterProteinAccessions, description, sep="\n")) %>%
      mutate(AnnotatedSequence = paste(AnnotatedSequence, PhosphorylatedAA, sep="-")) %>%
      mutate(AnnotatedSequence = paste(AnnotatedSequence, MasterProteinAccessions, sep="\n"))
  }

  # remove proteins without enough replicates if desired
  if (setminimalrep) {
    counttable <- data %>% group_by(AnnotatedSequence) %>% summarize(count=n()) %>% filter(count >= minimalrep)
    fkeep <- which(data$AnnotatedSequence %in% counttable$AnnotatedSequence)
    data <- data[fkeep, ]
  }

  data$PhosphorylatedAA <- NULL
  data$MasterProteinAccessions <- NULL
  names(data)[c(1:2)] <- c("id","condition")
  data <- data[ ,c(1:(nread+2))]

  if (xlinear | xsqrt | xcubert) {
    xmin <- as.numeric(names(data)[3])
    xmax <- as.numeric(names(data)[(nread+2)])
    print("The x axis scale range is: ")
    print(c(xmin, xmax)) # scale used for plotting
  } else if (xlog10) {
    first <- log10(as.numeric(names(data)[3]))
    last <- log10(as.numeric(names(data)[(nread+2)]))
    xmin <- (floor(first) + ceiling(first))/2
    if ( xmin > first ) { xmin <- floor(first) }
    if ( xmin == -Inf ) {
      # to hardcode zero dose to two log units smaller than the minimal dose
      xmin <- floor(log10(as.numeric(names(data)[4])/10))
      sub <- as.numeric(names(data)[4])/100
      print(paste0("In log scale, the zero was substituted as ", sub, " for regression."))
      names(data) <- gsub("^0$", sub, names(data))
    }
    xmax <- (floor(last) + ceiling(last))/2
    if (xmax < last) { xmax <- ceiling(last) }

    print("The log10 x axis scale range is: ")
    print(c(xmin, xmax)) # scale used for plotting
  }

  nametreatmentvector <- names(data)[3:(nread+2)]
  print("The points are as follows: ")
  print(as.numeric(nametreatmentvector))

  if (xlinear | xsqrt | xcubert) {
    treat_min <- xmin
    treat_max <- xmax
  } else if (xlog10) {
    treat_min <- log(min(as.numeric(nametreatmentvector)),10)-0.05 # extra 10% room for buffering
    treat_max <- log(max(as.numeric(nametreatmentvector)),10)+0.05
  }

  bottom_label <- paste0("Compound concentration(", unit, ")")

  if (orderEC) {
    plotseq <- unique(data[order(data$EC, decreasing=FALSE, na.last=TRUE), ]$id)
  }

  if (orderAUC) {
    data$AUC <- rowSums(data[ ,c(3:(nread+2))], na.rm=T)
    plotseq <- names(sort(with(data, tapply(AUC, id, mean)), decreasing=TRUE))
    data$AUC <- NULL
    #plotseq <- names(sort(tapply(data_l$reading, data_l$id, mean), decreasing=TRUE))
  }

  if (orderEC | orderAUC | length(plotseq)) {
    data$id <- factor(data$id, levels=plotseq)
  } else {
    data$id <- factor(data$id)
  }

  data_l <- tidyr::gather(data, treatment, reading, -id, -condition)

  data_l$condition <- factor(data_l$condition)
  data_l$treatment <- as.numeric(as.character(data_l$treatment))

  ncond <- length(unique(data_l$condition))
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
  names(colorpanel) <- levels(data_l$condition)
  #print(colorpanel)

  #plotlegend <- ms_melt_legend(data, nread, colorpanel)

  print("Generating fitted plot file, pls wait.")

  if (external) { external_graphs(T) }

  pep_isothermal_innerplotting <- function(d1) {

    q <- ggplot()

    if (loess) {
      q <- ggplot(d1, aes(x = treatment, y = reading, group = condition, colour = condition)) +
        geom_point(aes(colour=condition), shape=20, size=2) +
        geom_smooth(method = "loess", aes(colour=condition), se=FALSE, size=0.5) +
        scale_colour_manual(drop=FALSE, values=colorpanel) + labs(x=" ", y=" ")
    } else if (dotconnect) {
      q <- ggplot(d1, aes(x = treatment, y = reading, group = condition, colour = condition)) +
        geom_point(aes(colour=condition), shape=20, size=2) +
        geom_line(aes(colour=condition), size=0.5) +
        scale_colour_manual(drop=FALSE, values=colorpanel) + labs(x=" ", y=" ")
    } else {
      d1 <- droplevels(d1)
      d1_list <- split(d1, d1$condition)
      for (j in names(d1_list)) { # use levels of conditions, but not length of condition for individual proteins (length(d1_list))
        #print(colorpanel[j])
        q <- q + geom_point(data = d1_list[[j]], aes(x = treatment, y = reading), colour=colorpanel[j], shape=20, size=2)
        d1_list[[j]] <- na.omit(d1_list[[j]])
        model.drm <- try(drm (reading ~ treatment, data = d1_list[[j]], fct = LL.4()))
        if(class(model.drm) != "try-error"){
          if (xlinear | xsqrt | xcubert) {
            mda <- data.frame(treatment=seq(treat_min,treat_max,length.out=100))
            mda$reading <- predict(model.drm, mda)
          } else if (xlog10) {
            mda <- data.frame(treatment=10^(seq(treat_min,treat_max,length.out=100)))
            mda$reading <- predict(model.drm,mda)
          }
          q <- q + geom_line(data = mda, aes(x = treatment, y = reading), colour=colorpanel[j], size=0.5)
        } else {
          q <- q + geom_smooth(data = d1_list[[j]], method = "lm", formula = y~poly(x,1), aes(x = treatment, y = reading), colour=colorpanel[j], se=FALSE, size=0.5)
        }
      }
    }

    if (xlinear) {
      if (!length(xinterval)) {
        q <- q + coord_cartesian(xlim=c(xmin,xmax)) + labs(x=" ", y=" ")
      } else {
        q <- q + coord_cartesian(xlim=c(xmin,xmax)) + scale_x_continuous(breaks=seq(xmin, xmax, xinterval)) + labs(x=" ", y=" ")
      }
    } else if (xlog10) {
      q <- q + coord_cartesian(xlim=c(10^xmin,10^xmax)) + scale_x_log10(limits=c(10^xmin,10^xmax), breaks = trans_breaks("log10", function(x) 10^x),
                                                                        labels = trans_format("log10", math_format(10^.x))) + labs(x=" ", y=" ")
    } else if (xsqrt) {
      q <- q + coord_cartesian(xlim=c(0,xmax)) + scale_x_continuous(trans=sqrt_trans(), breaks=trans_breaks("sqrt", function(x) x^2)) + labs(x=" ", y=" ")
    } else if (xcubert) {
      cuberoot_trans = function() trans_new("cuberoot", function(x) x^(1/3), function(x) x^3)
      q <- q + coord_cartesian(xlim=c(0,xmax)) + scale_x_continuous(trans=cuberoot_trans(), breaks=seq(0,xmax,xinterval)) + labs(x=" ", y=" ")
    }

    if (min(abs(d1[,4]), na.rm=T) <= 0.5) {
      mi = 0
    } else {
      mi = 0.5
    }

    if (fixedy | max(abs(d1[,4]), na.rm=T) < 2.0 ) {
      q <- q + coord_cartesian(ylim = c(mi-0.05,2)) + scale_y_continuous(breaks=seq(mi, 2, 0.5))
    } else {
      ma <- ceiling(2*max(abs(d1[,4]), na.rm=T))/2
      if (ma<=3.0) {
        q <- q + coord_cartesian(ylim = c(mi-0.05,ma+0.1)) + scale_y_continuous(breaks=seq(mi, ma, 0.5))
      } else if (ma<=5.0) {
        q <- q + coord_cartesian(ylim = c(mi-0.05,ma+0.1)) + scale_y_continuous(breaks=seq(mi, ma, 1))
      } else {
        q <- q + coord_cartesian(ylim = c(mi-0.05,ma+0.1)) + scale_y_continuous(breaks=seq(mi, ma, 2))
      }
    }

    q <- q + ggtitle(as.character(unique(d1$id)))

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

  # conditions proteins
  plots <- plyr::dlply(data_l, plyr::.(id), pep_isothermal_innerplotting)

  if (!length(legenddata)) { legenddata <- data }
  plotlegend <- ms_isothermal_legend(legenddata, nread, nreplicate, presetcolor, colorpanel)
  legend <- plotlegend$legend
  lheight <- plotlegend$lheight

  params <- list(nrow=layout[1], ncol=layout[2])
  n <- with(params, nrow*ncol)
  ## add one page if division is not complete
  pages <- length(plots) %/% n + as.logical(length(plots) %% n)
  groups <- split(seq_along(plots), gl(pages, n, length(plots)))

  # pl <- lapply(names(groups), function(g){
  #   do.call(gridExtra::arrangeGrob, c(plots[groups[[g]]], params, top="CETSA peptide level plotting",
  #                          left="Non-denatured protein fraction", bottom=bottom_label))
  # })

  pl <- lapply(names(groups), function(g) {
    gridExtra::grid.arrange(
      do.call(gridExtra::arrangeGrob,
              c(plots[groups[[g]]], params,
                top=top_label,
                left="Non-denatured protein fraction",
                bottom=bottom_label)),
      #bottom=expression(Temperature~(R^{o}~C))
      legend,
      ncol = 1,
      heights = grid::unit.c(unit(1, "npc") - lheight, lheight))
  })

  class(pl) <- c("arrangelist", "ggplot", class(pl))
  # print(outdir)
  ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M"),"_peptide",pdfname), pl, height=12, width=12)
  if (external) { external_graphs(F) } # switch off the external graphs
  print("peptide level plots generated.")
}
