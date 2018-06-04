#' ms_2D_barplotting
#'
#' Function to generate pdf files with multipanel bar plots for 2D-CETSA data
#'
#' @param data dataset after ms_2D_caldiff to plot
#' @param treatmentlevel a vector of treatment labels, such as c("DMSO","TNFa","AT26533")
#' the order determines the arrangement, so in this case DMSO group would be the first group
#' @param setlevel a vector of set information if any, such as c("M13","M16")
#' @param colorpanel a vector of customizable color scheme provided by the user
#' @param log2scale whether the yscales should be in log2 scale, default set to TRUE
#' @param layout a vector indicating the panel layout for multi-panel plots
#' per page, default value is c(2,3) for set containing data, otherwise c(4,3)
#' @param toplabel textual label at the top part of the page
#' @param leftlabel textual label at the left side of the page
#' @param bottomlabel textual label at the bottom part of the page
#' @param pdfheight a number indicate the height of pdf file, default value 12
#' @param pdfwidth a number indicate the width of pdf file, default value 12
#'
#'
#' @import tidyr
#' @import dplyr
#' @import RColorBrewer
#' @import ggplot2
#'
#' @export
#'
#' @return a list of ggplot2 object
#' @examples \dontrun{
#'   ms_2D_barplotting(MOLM, treatmentlevel=c("DMSO","TNFa","AT26533"), setlevel=c("M13","M16"))
#' }
#'
#'

ms_2D_barplotting <- function(data, treatmentlevel=NULL, setlevel=NULL,
                              printBothName=TRUE, printGeneName=FALSE, pfdatabase=FALSE,
                              witherrorbar=TRUE, colorpanel=c("gray","blue","orange"),
                              layout=NULL, external=TRUE, log2scale=TRUE,
                              toplabel="2D-CETSA bar plotting", leftlabel="", bottomlabel="",
                              pdfname="bar_ggplotting.pdf", pdfheight=12, pdfwidth=12) {

  # legenddata is any dataset containing the full levels of conditions, same as data
  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  # to concatenate id and description
  nrowdata <- nrow(data)
  if ( nrowdata==0 ) {
    print("Make sure there are more than one experimental condition in dataset.")
    stop("Otherwise specify remsinglecondprot==FALSE !")
  }

  if (printBothName) {
    data <- data %>% rowwise() %>% mutate(description1=getProteinName(description)) %>%
      mutate(description2=getGeneName(description)) %>%
      mutate(id=paste(id, description1, description2, sep="\n"))
    data$description1<-NULL
    data$description2<-NULL
  } else if (printGeneName) {
    data <- data %>% rowwise() %>% mutate(description=getGeneName(description)) %>%
      mutate(id=paste(id, description, sep="\n"))
  } else {
    data <- data %>% rowwise() %>% mutate(description=getProteinName(description)) %>%
      mutate(id=paste(id, description, sep="\n"))
  }
  data$description<-NULL

  data1 <- tidyr::gather(data[ ,!(names(data) %in% c("sumUniPeps","sumPSMs","countNum"))], condition, reading, -id)
  if (!log2scale) {
    data1 <- mutate(data1, reading=2^reading)
  }
  a <- data1$condition[1]
  if (length(unlist(strsplit(a, "_")))==4) {
    withset <- TRUE
    data1 <- tidyr::separate(data1, condition, into=c("set","temperature","replicate","treatment"), sep="_")
    temperature <- sort(unique(data1$temperature))
    cdata <- plyr::ddply(data1, c("id", "set", "temperature", "treatment"), summarise,
                         N    = length(na.omit(reading)),
                         mean = mean(reading, na.rm=T),
                         sd   = sd(reading, na.rm=T),
                         se   = sd / sqrt(N)
    )
    if (length(layout)==0) {layout <- c(2,3)}
  } else if (length(unlist(strsplit(a, "_")))==3) {
    withset <- FALSE
    data1 <- tidyr::separate(data1, condition, into=c("temperature","replicate","treatment"), sep="_")
    temperature <- sort(unique(data1$temperature))
    cdata <- plyr::ddply(data1, c("id", "temperature", "treatment"), summarise,
                         N    = length(na.omit(reading)),
                         mean = mean(reading, na.rm=T),
                         sd   = sd(reading, na.rm=T),
                         se   = sd / sqrt(N)
    )
    if (length(layout)==0) {layout <- c(4,3)}
  } else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }

  cdata <- cdata %>% rowwise() %>% mutate(condition=paste(temperature, treatment, sep="_"))
  if (withset) { cdata$set <- factor(as.character(cdata$set), levels=setlevel) }
  cdata$treatment <- factor(as.character(cdata$treatment), levels=treatmentlevel)
  cdata$condition <- factor(as.character(cdata$condition),
                            levels = apply(expand.grid(temperature,treatmentlevel), 1, paste, collapse="_"))

  print("Generating fitted plot file, pls wait.")

  barplotting <- function(d1, withset=FALSE) {

    # print(max(d1$reading))
    # print(min(d1$reading))
    if (withset) {
      d1 <- droplevels(d1)
      d1_list <- split(d1, d1$set)
      q_list <- list()
      for (j in names(d1_list)) {
        #print(d1_list[[j]])
        if (nrow(d1_list[[j]])>0) {
          d2 <- d1_list[[j]]
          if (!log2scale) {
            minreading=0.5
            maxreading=2
            legendscale = c(min(max(min(d2$mean, na.rm=T)-0.5, 0), minreading), max(max(d2$mean, na.rm=T)+0.5, maxreading))
          } else {
            minreading=-0.5
            maxreading=0.5
            legendscale = c(min(min(d2$mean, na.rm=T)-0.1, minreading), max(max(d2$mean, na.rm=T)+0.1, maxreading))
          }
          q <- ggplot(d2, aes(x = condition, y = mean, fill=treatment)) +
            geom_bar(stat = "identity") + coord_cartesian(ylim = legendscale) +
            scale_fill_manual(drop=FALSE, values=colorpanel)
          if (witherrorbar) { q <- q + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, # Width of the error bars
                                                       position=position_dodge(.9)) }
          if (log2scale) {
            q <- q + ylab("fold change(log2)") + ggtitle(paste(j, as.character(unique(d2$id)), sep="\n"))
          } else {
            q <- q + ylab("fold change") + ggtitle(paste(j, as.character(unique(d2$id)), sep="\n"))
          }

          q <- q + theme_classic() +
            theme(
              text = element_text(size=10),
              strip.text.x = element_text(size = 5),
              plot.title = element_text(hjust=0.5, size = rel(0.8)),
              legend.background=element_rect(fill=NULL),
              legend.key.height=unit(0.5, "cm"),
              legend.key.width=unit(0.15, "cm"),
              legend.title = element_text(face = "bold"),
              legend.text = element_text(size = rel(0.7)),
              legend.justification="center",
              #legend.position="none", #c(0.2,0.8),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              #strip.text = element_blank(),
              axis.line.x = element_line(),
              axis.line.y = element_line(),
              axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.7)),
              aspect.ratio = 0.6
            )
          q_list[[j]] <- q
        } else {
          q <- ggplot()
          q_list[[j]] <- q
        }
      }
      q_list <- gridExtra::grid.arrange(grobs=q_list, ncol=1)
      return(q_list)
    } else {
      if (!log2scale) {
        minreading=0.5
        maxreading=2
        legendscale = c(min(max(min(d1$mean, na.rm=T)-0.5, 0), minreading), max(max(d1$mean, na.rm=T)+0.5, maxreading))
      } else {
        minreading=-0.5
        maxreading=0.5
        legendscale = c(min(min(d1$mean, na.rm=T)-0.1, minreading), max(max(d1$mean, na.rm=T)+0.1, maxreading))
      }
      q <- ggplot(d1, aes(x = condition, y = mean, fill=treatment)) +
        geom_bar(stat = "identity") + coord_cartesian(ylim = legendscale) +
        scale_fill_manual(drop=FALSE, values=colorpanel)
      if (witherrorbar) { q <- q + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, # Width of the error bars
                                                 position=position_dodge(.9)) }
      if (log2scale) {
        q <- q + ylab("fold change(log2)") + ggtitle(as.character(unique(d1$id)))
      } else {
        q <- q + ylab("fold change") + ggtitle(as.character(unique(d1$id)))
      }

      q <- q + theme_classic() +
        theme(
          text = element_text(size=10),
          strip.text.x = element_text(size = 5),
          plot.title = element_text(hjust=0.5, size = rel(0.8)),
          legend.background=element_rect(fill=NULL),
          legend.key.height=unit(0.5, "cm"),
          legend.key.width=unit(0.15, "cm"),
          legend.title = element_text(face = "bold"),
          legend.text = element_text(size = rel(0.7)),
          legend.justification="center",
          #legend.position="none", #c(0.2,0.8),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          strip.background = element_blank(),
          #strip.text = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.7)),
          aspect.ratio = 0.6
        )
      return(q)
    }
  }

  if (external) { external_graphs(T) }

  plots <- plyr::dlply(cdata, plyr::.(id), .fun=barplotting, withset=withset)

  params <- list(nrow=layout[1], ncol=layout[2])
  n <- with(params, nrow*ncol)
  ## add one page if division is not complete
  pages <- length(plots) %/% n + as.logical(length(plots) %% n)
  groups <- split(seq_along(plots), gl(pages, n, length(plots)))

  pl <- lapply(names(groups), function(i){

    gridExtra::grid.arrange(
      do.call(gridExtra::arrangeGrob,
              c(plots[groups[[i]]], params, top=toplabel,
                left=leftlabel,bottom=bottomlabel)))
  })

  class(pl) <- c("arrangelist", "ggplot", class(pl))
  pdfname <- gsub("/", " ", pdfname)
  if (length(outdir)) {
    ggsave(file=paste0(outdir,"/",format(Sys.time(), "%y%m%d_%H%M_"), pdfname),
           pl, height=pdfheight, width=pdfwidth)
  } else {
    ggsave(file=paste0(format(Sys.time(), "%y%m%d_%H%M"), pdfname), pl,
           height=pdfheight, width=pdfwidth)
  }

  if (external) { external_graphs(F) } # switch off the external graphs
  print("2D-CETSA bar plot file generated successfully.")
}
