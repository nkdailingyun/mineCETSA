#' ms_legend
#'
#' Internal function to generate legend for ggplotting
#'
#' @param legenddata dataset used for condition extraction, at least one protein
#' inside this dataset should contains the full set of experiment conditions
#' @param nread number of reading channels or sample treatements
#' @param nreplicate number of replicates
#' @param presetcolor whether to use the pre-defined color scheme
#' @param colorpanel a vector of customizable color scheme provided by the user
#'
#' @import tidyr
#' @import RColorBrewer
#' @import ggplot2
#' @keywords internal
#'
#' @return a legend object
#' @examples \dontrun{
#' }
#'
#'

ms_legend <- function(legenddata, nread, nreplicate, presetcolor, colorpanel) {

  #The input data must contain at least one full condition
  protein_with_fulllevels = names(sort(with(legenddata,tapply(condition,id,function(u)length(unique(u)))),decreasing=T)[1])
  legenddata <- legenddata[legenddata$id == protein_with_fulllevels, ]
  if (!length(grep("description", names(legenddata)))) {
    d1 <- tidyr::gather(legenddata[ ,c(1:(nread+2))], treatment, reading, -id, -condition)
  } else {
    d1 <- tidyr::gather(legenddata[ ,c(1,3:(nread+3))], treatment, reading, -id, -condition)
  }
  #d1 <- melt(legenddata[, c(1:(nread+2))], id.vars=c("id", "condition"), variable.name="treatment", value.name="reading");
  d1$condition <- factor(d1$condition)
  ncond <- length(unique(d1$condition))
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
  names(colorpanel) <- levels(d1$condition)

  demo <- ggplot(data = d1, aes(x=treatment,y=reading,group=condition,colour=condition))
  demo <- demo + geom_point() + scale_colour_manual(name='', values=colorpanel)
  g <- ggplotGrob(demo + theme(legend.position="bottom", legend.background=element_rect(fill=NULL)))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  return(list(legend=legend, lheight=lheight))
}
