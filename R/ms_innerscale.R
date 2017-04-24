#' ms_innerscale
#'
#' Internal function to generate median or sum of reading data, according to condition group
#'
#' @param data dataset to be calculated
#' @param nread number of reading channels or sample treatements, default value 10
#' @param sumf whether to calculate sum or not, default to FALSE, i.e., to calculate median
#'
#' @importFrom gtools mixedsort
#' @keywords internal
#'
#' @return a dataframe
#' @examples \dontrun{
#' }
#'
#'

ms_innerscale <- function(data, nread, sumf=FALSE) {

  # create vector of unique ids and prepare to iterate
  uniqueCondition <- unique(data$condition)
  niter <- length(uniqueCondition)
  uniques <- which(!duplicated(data$condition))
  # create outdataset
  outdata <- data[uniques, c(1:(nread+3))]
  # calculate median
  for (i in 1:niter) {
    pattern <- grep(paste0("^",uniqueCondition[i],"$"), data$condition, value=FALSE)
    tmpdata <- data[pattern, ]
    for (j in 4:(nread+3)) {
      if (sumf) {
        outdata[i,j] <- sum(tmpdata[[j]], na.rm=TRUE)
      } else {
        outdata[i,j] <- median(tmpdata[[j]], na.rm=TRUE)
      }
    }
  }
  return(outdata)
}


#' ms_innerplotbox
#'
#' Internal function to generate box plot on reading data, according to condition group
#'
#' @param data dataset to be calculated
#' @param filename filename to be specified
#' @param xlabel label on x axis to be specified
#' @param ylabel label on y axis to be specified
#' @param isratio whether the reading data are ratio or abundance
#' @param isothermalstyle whether the reading data are from Isothermal or CETSA melt curve data
#' @param outdir the subfolder name to save the plot
#'
#' @keywords internal
#'
#' @import tidyr
#' @import dplyr
#' @import ggplot2
#'
#' @return NULL
#' @examples \dontrun{
#' }
#'
#'

ms_innerplotbox <- function(data, filename, xlabel, ylabel, isratio, isothermalstyle, outdir=NULL) {

  data <- tidyr::gather(data, treatment, reading, -id, -condition)
  #d1 <- melt(data, id.vars=c("id", "condition"), variable.name="treatment", value.name="reading")
  if (isratio) {
    data$treatment <- factor(data$treatment, levels=gtools::mixedsort(unique(data$treatment)))
    q <- ggplot(data, aes(x = treatment, y = reading))
    if (isothermalstyle) {
      q <- q + coord_cartesian(ylim=c(0,2)) + scale_y_continuous(breaks=c(0,0.9,1,1.1,2)) +
        geom_hline(yintercept=c(0.9,1.1), color="blue", size=0.5, linetype=c(2,2))
    } else {
      q <- q + coord_cartesian(ylim=c(0,1.2)) + scale_y_continuous(breaks=c(0,0.5,1))
    }
    q <- q + labs(x=xlabel, y=ylabel) + geom_boxplot(aes(fill=condition)) +
      theme(text=element_text(size=20), axis.text.x = element_text(angle=45,hjust=1), aspect.ratio=1)
  } else {
    data$treatment <- factor(data$treatment, levels=gtools::mixedsort(unique(data$treatment)))
    q <- ggplot(data, aes(x=treatment, y = log10(reading)))
    q <- q + labs(x=xlabel, y=ylabel) + geom_boxplot(aes(fill=condition)) +
      theme(text=element_text(size=20), axis.text.x = element_text(angle=45,hjust=1), aspect.ratio=1)
  }
  if (length(outdir)) {
    ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),filename), q, height=12, width=12)
  } else {
    ggsave(file=paste0(format(Sys.time(), "%y%m%d_%H%M_"),filename), q, height=12, width=12)
  }
}

#' ms_innerplotmedian
#'
#' Internal function to generate median line plot on reading data,
#' comparing pre-scaling and post-scaling median change, according to condition group
#'
#' @param mediandata pre-scaling and post-scaling median data, generated from
#' \code{ms_innerscale} and with a set category factor
#' @param filename filename to be specified
#' @param xlabel label on x axis to be specified
#' @param ylabel label on y axis to be specified
#' @param isothermalstyle whether the reading data are from Isothermal or CETSA melt curve data
#' @param outdir the subfolder name to save the plot
#'
#' @keywords internal
#'
#' @import dplyr
#' @import ggplot2
#'
#' @return NULL
#' @examples \dontrun{
#' }
#'
#'

ms_innerplotmedian <- function(mediandata, filename, xlabel, ylabel, isothermalstyle, outdir=NULL) {

  # mdata$set <- "Pre-normalization"
  # scalemdata$set <- "Post-normalization"
  # mdata <- tidyr::gather(mdata, treatment, reading, -id, -set, -condition)
  # scalemdata <- tidyr::gather(scalemdata, treatment, reading, -id, -set, -condition)
  #mdata <- melt(mdata, id.vars=c("id", "set", "condition"), variable.name="treatment", value.name="reading")
  #scalemdata <- melt(scalemdata, id.vars=c("id", "set", "condition"), variable.name="treatment", value.name="reading")
  # data <- rbind(mdata, scalemdata)
  # data$set <- factor(data$set, levels=c("Pre-normalization", "Post-normalization"))
  # mediandata$treatment <- factor(mediandata$treatment,
  #                                levels=sort(as.numeric(unique(mediandata$treatment)), decreasing=FALSE))
  q <- ggplot(mediandata, aes(x = treatment, y = reading))
  if (isothermalstyle) {
    q <- q + coord_cartesian(ylim=c(0,2)) + scale_y_continuous(breaks=seq(0,2,0.5))
  } else {
    q <- q + coord_cartesian(ylim=c(0,1.2)) + scale_y_continuous(breaks=seq(0,1,0.5))
  }
  q <- q + labs(x=xlabel, y=ylabel) + geom_point(aes(group = condition, colour = condition)) +
    geom_line(aes(group = condition, colour = condition)) + facet_grid( .~set) +
    theme(text=element_text(size=20), axis.text.x = element_text(angle=45,hjust=1),
          aspect.ratio=1, legend.position="bottom")
  if (length(outdir)) {
    ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"),filename), q, height=6, width=12)
  } else {
    ggsave(file=paste0(format(Sys.time(), "%y%m%d_%H%M_"),filename), q, height=6, width=12)
  }
}
