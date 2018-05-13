#' ms_isothermal_bar_innerplot
#'
#' Internal function to generate ggplot objects for isothermal data in a bar plot format
#'
#' @param data isothermal dataset to plot
#' @param legenddata dataset used for condition extraction, at least one protein
#' inside this dataset should contains the full set of experiment conditions
#' @param levelvector a vector of experimental conditions, not complusory for isothermal functions
#' @param nread number of reading channels or sample treatements
#' @param minireplicate number of replicates to keep in final data, default to NULL
#' @param fixedy whether the y-axis should use a fixed range scale
#' @param presetcolor whether to use the pre-defined color scheme
#' @param colorpanel a vector of customizable color scheme provided by the user
#' @param layout a vector indicating the panel layout for multi-panel plots per page
#' @param toplabel textual label at the top of the page
#' @param leftlabel textual label at the left side of the page
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

ms_isothermal_bar_innerplot <- function(data, legenddata, levelvector, nread,
                                        minireplicate, witherrorbar, usegradient,
                                        fixedy, presetcolor, colorpanel, layout,
                                        toplabel, leftlabel, bottomlabel) {

  nametreatmentvector <- names(data)[3:(nread+2)]
  print("The points are as follows: ")
  print(nametreatmentvector)
  #print(as.numeric(nametreatmentvector))

  print("Generating fitted plot file, pls wait.")

  data_l <- tidyr::gather(data[ ,c(1:(nread+2))], treatment, reading, -id, -condition)
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
  cdata <- tidyr::complete(cdata, id, sample, treatment)
  nsample <- length(unique(cdata$sample))

  if (presetcolor & length(colorpanel)==0) {
    if (nsample<=9) {colorpanel <- brewer.pal(9, "Set1")}
    else {stop("The number of conditions in dataset exceeds the preset number of colors, pls provide a vector of colors in colorpanel")}
  } else if (length(colorpanel) < nsample){
    stop("The number of conditions in dataset exceeds the provided number of colors, pls provide enough colors in colorpanel")
  }
  #colorpanel <- c("blue1", "blue4", "red1", "red4", "green1", "green4", "yellow", "yellow4", "gray30", "black")
  #colorpanel <- brewer_pal("qual")(length(levels(data_l$condition)))
  if (length(levelvector)) {
    cdata$sample<-factor(as.character(cdata$sample), levels=levelvector)
  } else {
    cdata$sample<-factor(cdata$sample, levels=unique(cdata$sample))
  }
  cdata$treatment<-factor(as.character(cdata$treatment), levels=nametreatmentvector)

  plotting <- function(d1, minreading=0.5, maxreading=2, witherrorbar=TRUE, usegradient=TRUE) {

    legendscale = c(min(round(min(d1$mean, na.rm=T)-0.1,2), minreading), max(round(max(d1$mean, na.rm=T)+0.1,2), maxreading))

    if (nsample==1) {
      if (usegradient) {
        #print("use dichromatic color scheme")
        q <- ggplot(d1, aes(x=treatment, y=mean, fill=mean)) +
          geom_bar(stat="identity") + coord_cartesian(ylim=legendscale) +
          scale_fill_gradient2(limits=legendscale, low="#4575B4", mid="ivory", high="#D73027",
                                              midpoint=1, na.value="gray90", guide=guide_colorbar(""))
      } else {
        #print("use monochromatic color scheme")
        q <- ggplot(d1, aes(x=treatment, y=mean, fill=treatment)) +
          geom_bar(stat="identity", position=position_dodge()) + coord_cartesian(ylim=legendscale) +
          scale_fill_manual(drop=FALSE, values=colorpanel)
      }
    } else {
      q <- ggplot(d1, aes(x=treatment, y=mean, fill=sample, color=sample)) +
        geom_bar(stat="identity", position=position_dodge()) + coord_cartesian(ylim=legendscale) +
        scale_fill_manual(drop=FALSE, values=colorpanel)
    }
    if (witherrorbar) { q <- q + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, # Width of the error bars
                                                 position=position_dodge(.9)) }


    q <- q + labs(x=" ", y=" ")

    if (fixedy | max(abs(d1$mean), na.rm=T) < 2.0 ) {
      q <- q + coord_cartesian(ylim=c(legendscale[1]-0.05,2)) + scale_y_continuous(breaks=seq(legendscale[1], 2, 0.5))
    } else {
      ma <- ceiling(2*max(abs(d1$mean), na.rm=T))/2
      if (ma<=3.0) {
        q <- q + coord_cartesian(ylim=c(legendscale[1]-0.05,legendscale[2]+0.1)) + scale_y_continuous(breaks=seq(legendscale[1], legendscale[2], 0.5))
      } else if (ma<=5.0) {
        q <- q + coord_cartesian(ylim=c(legendscale[1]-0.05,legendscale[2]+0.1)) + scale_y_continuous(breaks=seq(legendscale[1], legendscale[2], 1))
      } else {
        q <- q + coord_cartesian(ylim=c(legendscale[1]-0.05,legendscale[2]+0.1)) + scale_y_continuous(breaks=seq(legendscale[1], legendscale[2], 2))
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

  plots <- plyr::dlply(cdata, .(id), .fun = plotting, witherrorbar=witherrorbar, usegradient=usegradient)
  if(nsample>1) { # legend is only for multiple sample groups
    plotlegend <- ms_isothermal_legend2(legenddata, levelvector, nread, nsample, presetcolor, colorpanel)
    legend <- plotlegend$legend
    lheight <- plotlegend$lheight
  }

  params <- list(nrow=layout[1], ncol=layout[2])
  n <- with(params, nrow*ncol)
  ## add one page if division is not complete
  pages <- length(plots) %/% n + as.logical(length(plots) %% n)
  groups <- split(seq_along(plots), gl(pages, n, length(plots)))

  if(nsample>1) {
    pl <- lapply(names(groups), function(i){
      gridExtra::grid.arrange(
        do.call(arrangeGrob,
                c(plots[groups[[i]]], params, top=toplabel,
                  left=leftlabel,bottom=bottomlabel)),
      legend,
      ncol = 1,
      heights = grid::unit.c(unit(1, "npc") - lheight, lheight))
    })
  } else {
    pl <- lapply(names(groups), function(i){
      gridExtra::grid.arrange(
        do.call(arrangeGrob,
                c(plots[groups[[i]]], params, top=toplabel,
                  left=leftlabel,bottom=bottomlabel)))
    })
  }

  class(pl) <- c("arrangelist", "ggplot", class(pl))
  return(pl)
}
