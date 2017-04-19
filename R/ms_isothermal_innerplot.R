#' ms_isothermal_innerplot
#'
#' Internal function to generate ggplot objects for isothermal data
#'
#' @param data isothermal dataset to plot
#' @param legenddata dataset used for condition extraction, at least one protein
#' inside this dataset should contains the full set of experiment conditions
#' @param nread number of reading channels or sample treatements
#' @param nreplicate number of replicates
#' @param loess whether to perform curve fitting using loess model
#' @param dotconnect whether to simply dot connect the readings for each curve
#' @param PSManno whether to annotate the plots with PSM and uniPeptide number
#' @param PSM_annod dataframe containing the PSM number information
#' @param Pep_annod dataframe containing the uniPeptide number information
#' @param xliner whether the x-axis should be in linear scale
#' @param xlog10 whether the x-axis should be in log10 scale
#' @param xsqrt whether the x-axis should be in square-root transformed scale
#' @param xinterval a number indicating the numerical interval for linear x-axis
#' @param fixedy whether the y-axis should use a fixed range scale
#' @param presetcolor whether to use the pre-defined color scheme
#' @param colorpanel a vector of customizable color scheme provided by the user
#' @param layout a vector indicating the panel layout for multi-panel plots per page
#' @param toplabel textual label at the top of the page
#' @param bottomlabel textual label at the bottom of the page
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

ms_isothermal_innerplot <- function(data, legenddata, nread, nreplicate, loess,
                                    dotconnect, PSManno, PSM_annod, Pep_annod,
                                    xlinear, xlog10, xsqrt, xinterval, fixedy,
                                    presetcolor, colorpanel, layout,
                                    top_label, bottom_label) {

  if (xlinear | xsqrt) {
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

  if (xlinear | xsqrt) {
    treat_min <- xmin
    treat_max <- xmax
  } else if (xlog10) {
    treat_min <- log(min(as.numeric(nametreatmentvector)),10)-0.05 # extra 10% room for buffering
    treat_max <- log(max(as.numeric(nametreatmentvector)),10)+0.05
  }

  print("Generating fitted plot file, pls wait.")

  data_l <- tidyr::gather(data[ ,c(1:(nread+2))], treatment, reading, -id, -condition)
  #data_l <- melt(data[ ,c(1:(nread+2))], id.vars=c("id", "condition"), variable.name="treatment", value.name="reading")

  data_l$condition <- factor(data_l$condition)
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

  plotting <- function(d1) {

    if (PSManno) {
      PSM_label <- PSM_annod[PSM_annod$id==unique(d1$id), ][ ,c(2:ncol(PSM_annod))]
      Pep_label <- Pep_annod[Pep_annod$id==unique(d1$id), ][ ,c(2:ncol(Pep_annod))]
    }

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
          if (xlinear | xsqrt) {
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
    }

    if (fixedy | max(abs(d1[,4]), na.rm=T)<=2.0) {
      q <- q + coord_cartesian(ylim = c(0,2)) + scale_y_continuous(breaks=seq(0, 2, 0.5))
    } else {
      ma <- (floor(max(abs(d1[,4]), na.rm=T)) + ceiling(max(abs(d1[,4]), na.rm=T)))/2
      if (ma<=3.0) {
        q <- q + coord_cartesian(ylim = c(0,ma)) + scale_y_continuous(breaks=seq(0, ma, 0.5))
      } else {
        q <- q + coord_cartesian(ylim = c(0,ma)) + scale_y_continuous(breaks=seq(0, ma, 1))
      }
    }

    if (PSManno & fixedy) {
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

    q <- q + ggtitle(as.character(unique(d1$id)))

    q <- q + theme_classic() +
      theme(
        text = element_text(size=10),
        strip.text.x = element_text(size = 5),
        plot.title = element_text(hjust = 0.5, size = rel(0.7)),
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

  plots <- plyr::dlply(data_l, .(id), .fun = plotting)
  plotlegend <- ms_isothermal_legend(legenddata, nread, nreplicate, presetcolor, colorpanel)
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
              c(plots[groups[[i]]], params, top=top_label,
                left="Non-denatured protein fraction",bottom=bottom_label)),
      legend,
      ncol = 1,
      heights = grid::unit.c(unit(1, "npc") - lheight, lheight))

  })

  class(pl) <- c("arrangelist", "ggplot", class(pl))
  return(pl)
}
