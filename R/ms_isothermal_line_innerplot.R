#' ms_isothermal_line_innerplot
#'
#' Internal function to generate ggplot objects for isothermal data in a curve plot format
#'
#' @param data isothermal dataset to plot
#' @param legenddata dataset used for condition extraction, at least one protein
#' inside this dataset should contains the full set of experiment conditions
#' @param levelvector a vector of experimental conditions, not complusory for isothermal functions
#' @param nread number of reading channels or sample treatements
#' @param robustfitting whether to apply robust fitting, which is new since v.0.3.7
#' @param nreplicate number of replicates
#' @param minireplicate number of replicates to keep in final data
#' @param loess whether to perform curve fitting using loess model
#' @param dotconnect whether to simply dot connect the readings for each curve
#' @param printcount whether to annotate the plots with PSM and uniPeptide number
#' @param PSM_annod dataframe containing the PSM number information
#' @param Pep_annod dataframe containing the uniPeptide number information
#' @param xliner whether the x-axis should be in linear scale
#' @param xlog10 whether the x-axis should be in log10 scale
#' @param xsqrt whether the x-axis should be in square-root transformed scale
#' @param xcubert whether the x-axis should be in cube-root transformed scale
#' @param xinterval a number indicating the numerical interval for linear x-axis
#' @param fixedy whether the y-axis should use a fixed range scale
#' @param presetcolor whether to use the pre-defined color scheme
#' @param colorpanel a vector of customizable color scheme provided by the user
#' @param layout a vector indicating the panel layout for multi-panel plots per page
#' @param toplabel textual label at the top of the page
#' @param leftlabel textual label at the left side of the page
#' @param bottomlabel textual label at the bottom of the page
#' @param shadearea a four element vector indicating the xmin, xmax, ymin, ymax for x and y axis for shading
#' @param shadeoutlinecolor color used to shade the area, default value is black
#' @param shadefillcolor color used to shade the area, default value is gray90
#'
#'
#' @import tidyr
#' @import RColorBrewer
#' @import scales
#' @importFrom plyr . dlply
#' @importFrom grid unit.c
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @import ggplot2
#' @keywords internal
#'
#' @return a list of ggplot2 object
#' @examples \dontrun{
#' }
#'
#'

ms_isothermal_line_innerplot <- function(data, legenddata, levelvector, nread,
                                    robustfitting, nreplicate, minireplicate,
                                    witherrorbar, loess, dotconnect,
                                    printcount, PSM_annod, Pep_annod,
                                    xlinear, xlog10, xsqrt, xcubert, xinterval,
                                    fixedy, presetcolor, colorpanel, layout,
                                    toplabel, leftlabel, bottomlabel,
                                    shadearea, shadeoutlinecolor, shadefillcolor) {

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

  print("Generating fitted plot file, pls wait.")

  if (witherrorbar) {
    data_l <- tidyr::gather(data[ ,c(1:(nread+2))], treatment, reading, -id, -condition)
    #data_l <- melt(data[ ,c(1:(nread+2))], id.vars=c("id", "condition"), variable.name="treatment", value.name="reading")
    data_l <- tidyr::separate(data_l, condition, into=c("sample","rep"), sep="\\.")
    cdata <- plyr::ddply(data_l, c("id", "sample", "treatment"), summarise,
                         N    = length(reading),
                         mean = mean(reading),
                         sd   = sd(reading),
                         se   = sd / sqrt(N)
    )
    if (length(minireplicate)>0) {
      cdata <- subset(cdata, N>=minireplicate)
    }
    #cdata  <- na.omit(cdata) # to remove the measurements without enough replicates
    #cdata <- tidyr::complete(cdata, id, sample, treatment)
    nsample <- length(unique(data_l$sample))

    if (presetcolor & length(colorpanel)==0) {
      if (nsample<=9) {colorpanel <- brewer.pal(9, "Set1")}
      else {stop("The number of conditions in dataset exceeds the preset number of colors, pls provide a vector of colors in colorpanel")}
    } else if (length(colorpanel) < nsample){
      stop("The number of conditions in dataset exceeds the provided number of colors, pls provide enough colors in colorpanel")
    }

    if (length(levelvector)) {
      cdata$sample<-factor(as.character(cdata$sample), levels=levelvector)
    } else {
      cdata$sample<-factor(cdata$sample, levels=unique(cdata$sample))
    }
    cdata$treatment <- as.numeric(as.character(cdata$treatment))
    names(colorpanel) <- levels(cdata$sample)

  } else {
    data_l <- tidyr::gather(data[ ,c(1:(nread+2))], treatment, reading, -id, -condition)
    #data_l <- melt(data[ ,c(1:(nread+2))], id.vars=c("id", "condition"), variable.name="treatment", value.name="reading")

    if (length(levelvector)) {
      data_l$condition<-factor(as.character(data_l$condition), levels=levelvector)
    } else {
      data_l$condition<-factor(data_l$condition, levels=sort(unique(data_l$condition)))
    }
    data_l$treatment <- as.numeric(as.character(data_l$treatment))
    #data_l<-data.frame(data_l)
    #return(data_l)

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
  }
  #print(colorpanel)
  plotting1 <- function(d1) {

    if (printcount) {
      PSM_label <- PSM_annod[PSM_annod$id==unique(d1$id), ][ ,c(2:ncol(PSM_annod))]
      Pep_label <- Pep_annod[Pep_annod$id==unique(d1$id), ][ ,c(2:ncol(Pep_annod))]
    }

    q <- ggplot()

    if (loess) {
      q <- ggplot(d1, aes(x=treatment, y=reading, group=condition, colour=condition)) +
        geom_point(aes(colour=condition), shape=20, size=2) +
        geom_smooth(method="loess", aes(colour=condition), se=FALSE, size=0.5) +
        scale_colour_manual(drop=FALSE, values=colorpanel) + labs(x=" ", y=" ")
    } else if (dotconnect) {
      q <- ggplot(d1, aes(x=treatment, y=reading, group=condition, colour=condition)) +
        geom_point(aes(colour=condition), shape=20, size=2) +
        geom_line(aes(colour=condition), size=0.5) +
        scale_colour_manual(drop=FALSE, values=colorpanel) + labs(x=" ", y=" ")
    } else {
      d1 <- droplevels(d1)
      d1_list <- split(d1, d1$condition)
      for (j in names(d1_list)) { # use levels of conditions, but not length of condition for individual proteins (length(d1_list))
        #print(colorpanel[j])
        q <- q + geom_point(data=d1_list[[j]], aes(x=treatment, y=reading), colour=colorpanel[j], shape=20, size=2)
        d1_list[[j]] <- na.omit(d1_list[[j]])
        if (robustfitting) {
          model.drm <- try(drm(reading~treatment, data=d1_list[[j]], fct=LL.4(), robust="median"))
        } else {
          model.drm <- try(drm(reading~treatment, data=d1_list[[j]], fct=LL.4()))
        }
        if(class(model.drm) != "try-error"){
          if (xlinear | xsqrt | xcubert) {
            mda <- data.frame(treatment=seq(treat_min,treat_max,length.out=100))
            mda$reading <- predict(model.drm, mda)
          } else if (xlog10) {
            mda <- data.frame(treatment=10^(seq(treat_min,treat_max,length.out=100)))
            mda$reading <- predict(model.drm,mda)
          }
          q <- q + geom_line(data=mda, aes(x=treatment, y=reading), colour=colorpanel[j], size=0.5)
        } else {
          q <- q + geom_smooth(data=d1_list[[j]], method="lm", formula=y~poly(x,1), aes(x=treatment, y=reading), colour=colorpanel[j], se=FALSE, size=0.5)
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
      q <- q + coord_cartesian(xlim=c(10^xmin,10^xmax*1.5)) + scale_x_log10(limits=c(10^xmin,10^xmax*1.5), breaks = trans_breaks("log10", function(x) 10^x),
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
      q <- q + coord_cartesian(ylim=c(mi-0.05,2)) + scale_y_continuous(breaks=seq(mi, 2, 0.5))
    } else {
      ma <- ceiling(2*max(abs(d1[,4]), na.rm=T))/2
      if (ma<=3.0) {
        q <- q + coord_cartesian(ylim=c(mi-0.05,ma+0.1)) + scale_y_continuous(breaks=seq(mi,ma,0.5))
      } else if (ma<=5.0) {
        q <- q + coord_cartesian(ylim=c(mi-0.05,ma+0.1)) + scale_y_continuous(breaks=seq(mi,ma,1))
      } else {
        q <- q + coord_cartesian(ylim=c(mi-0.05,ma+0.1)) + scale_y_continuous(breaks=seq(mi,ma,2))
      }
    }

    if (printcount & fixedy) {
      q <- q + annotate("text", x=15, y=0.5, label="  #PSM  #Peptides", size=1.5)
      for (i in 1:length(PSM_label)) {
        q <- q + annotate("text", x=15, y=0.5-0.1*i, label=paste0(names(PSM_label)[i], ":  ", as.character(PSM_label)[i], " "), size=1.5, hjust=1)
        i=i+1
      }
      for (i in 1:length(Pep_label)) {
        q <- q + annotate("text", x=45, y=0.5-0.1*i, label=as.character(Pep_label)[i], size=1.5, hjust=1)
        i=i+1
      }
    }
    if (length(shadearea)) {
      q <- q + geom_rect(mapping=aes(xmin=shadearea[1], xmax=shadearea[2],
                                     ymin=shadearea[3], ymax=shadearea[4]),
                         color=shadeoutlinecolor, fill=shadefillcolor, alpha=0.1)
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

  plotting2 <- function(d1) {

    q <- ggplot()

    if (loess) {
      q <- ggplot(d1, aes(x=treatment, y=mean, group=sample, colour=sample)) +
        geom_point(aes(colour=sample), shape=20, size=2) +
        geom_smooth(method="loess", aes(colour=sample), se=FALSE, size=0.5) +
        scale_colour_manual(drop=FALSE, values=colorpanel) + labs(x=" ", y=" ")
      q <- q + geom_errorbar(d1, aes(ymin=mean-se, ymax=mean+se, colour=sample), width=.2,
                             position=position_dodge(.9))
    } else if (dotconnect) {
      q <- ggplot(d1, aes(x=treatment, y=mean, group=sample, colour=sample)) +
        geom_point(aes(colour=sample), shape=20, size=2) +
        geom_line(aes(colour=sample), size=0.5) +
        scale_colour_manual(drop=FALSE, values=colorpanel) + labs(x=" ", y=" ")
      q <- q + geom_errorbar(d1, aes(ymin=mean-se, ymax=mean+se, colour=sample), width=.2,
                             position=position_dodge(.9))
    } else {
      d1 <- droplevels(d1)
      d1_list <- split(d1, d1$sample)
      for (j in names(d1_list)) {
        # use levels of conditions, but not length of condition for individual proteins (length(d1_list))
        q <- q + geom_point(data=d1_list[[j]], aes(x=treatment, y=mean), colour=colorpanel[j], shape=20, size=2)
        d1_list[[j]] <- na.omit(d1_list[[j]])
        if (robustfitting) {
          model.drm <- try(drm(mean~treatment, data=d1_list[[j]], fct=LL.4(), robust="median"))
        } else {
          model.drm <- try(drm (mean~treatment, data=d1_list[[j]], fct=LL.4()))
        }
        if (class(model.drm) != "try-error") {
          if (xlinear | xsqrt | xcubert) {
            mda <- data.frame(treatment=seq(treat_min,treat_max,length.out=100))
            mda$reading <- predict(model.drm, mda)
          } else if (xlog10) {
            mda <- data.frame(treatment=10^(seq(treat_min,treat_max,length.out=100)))
            mda$reading <- predict(model.drm,mda)
          }
          q <- q + geom_line(data=mda, aes(x=treatment, y=reading), colour=colorpanel[j], size=0.5)
        } else {
          q <- q + geom_smooth(data=d1_list[[j]], method="lm", formula=y~poly(x,1),
                               aes(x=treatment, y=mean), colour=colorpanel[j], se=FALSE, size=0.5)
        }

        q <- q + geom_errorbar(data=d1_list[[j]], aes(x=treatment, ymin=mean-se, ymax=mean+se),
                               colour=colorpanel[j], width=.2, # Width of the error bars
                               position=position_dodge(.9))
      }
    }

    if (xlinear) {
      if (!length(xinterval)) {
        q <- q + coord_cartesian(xlim=c(xmin,xmax)) + labs(x=" ", y=" ")
      } else {
        q <- q + coord_cartesian(xlim=c(xmin,xmax)) + scale_x_continuous(breaks=seq(xmin, xmax, xinterval)) + labs(x=" ", y=" ")
      }
    } else if (xlog10) {
      q <- q + coord_cartesian(xlim=c(10^xmin,10^xmax*1.5)) + scale_x_log10(limits=c(10^xmin,10^xmax*1.5), breaks = trans_breaks("log10", function(x) 10^x),
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

    if (fixedy | max(abs(d1[,5]), na.rm=T) < 2.0) {
      q <- q + coord_cartesian(ylim=c(mi-0.05,2)) + scale_y_continuous(breaks=seq(mi,2,0.5))
    } else {
      ma <- ceiling(2*max(abs(d1[,5]), na.rm=T))/2
      if (ma<=3.0) {
        q <- q + coord_cartesian(ylim=c(mi-0.05,ma+0.1)) + scale_y_continuous(breaks=seq(mi,ma,0.5))
      } else if (ma<=5.0) {
        q <- q + coord_cartesian(ylim=c(mi-0.05,ma+0.1)) + scale_y_continuous(breaks=seq(mi,ma,1))
      } else {
        q <- q + coord_cartesian(ylim=c(mi-0.05,ma+0.1)) + scale_y_continuous(breaks=seq(mi,ma,2))
      }
    }

    if (length(shadearea)) {
      q <- q + geom_rect(mapping=aes(xmin=shadearea[1], xmax=shadearea[2],
                                     ymin=shadearea[3], ymax=shadearea[4]),
                         color=shadeoutlinecolor, fill=shadefillcolor, alpha=0.1)
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

  if (witherrorbar) {
    plots <- plyr::dlply(cdata, .(id), .fun = plotting2)
    plotlegend <- ms_isothermal_legend2(legenddata, levelvector,
                                        nread, nsample, presetcolor, colorpanel)
  } else {
    plots <- plyr::dlply(data_l, .(id), .fun = plotting1)
    plotlegend <- ms_isothermal_legend(legenddata, levelvector,
                                       nread, nreplicate, presetcolor, colorpanel)
  }
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
  return(pl)
}