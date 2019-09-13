#' ms_melt_innerplot
#'
#' Internal function to generate ggplot objects for CETSA melt curve data
#'
#' @param data melt curve dataset to plot
#' @param nread number of reading channels or sample treatments
#' @param topasone whether the top plateau has to be fixed, i.e., 1.0
#' @param dotconnect whether to simply dot connect the readings for each curve
#' @param printcount whether to annotate the plots with PSM and uniPeptide number
#' @param PSM_annod dataframe containing the PSM number information
#' @param Pep_annod dataframe containing the uniPeptide number information
#' @param annotypos the starting y-axis position of textual annotation
#' @param annotyinterval the interval on y-axis for textual annotation
#' @param colorpanel a vector of customizable color scheme provided by the user
#' @param plotlegend a legend object from ms_melt_legend
#' @param commonlegend whether to use one common legend for whole page of plots
#' @param toplabel textual label at the top of the page
#' @param leftlabel textual label at the left side of the page
#' @param bottomlabel textual label at the bottom of the page
#' @param withset whether there is set column to perform facet_grid
#' @param layout a vector indicating the panel layout for multi-panel plots per page
#'
#'
#' @import tidyr
#' @import scales
#' @import drc
#' @importFrom plyr . dlply
#' @importFrom grid unit.c
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @import ggplot2
#'
#' @keywords internal
#'
#' @return a list of ggplot2 object
#' @examples \dontrun{
#' }
#'
#'

ms_melt_innerplot <- function(data, nread, topasone, dotconnect, printcount,
                              PSM_annod, Pep_annod, annotypos, annotyinterval,
                              colorpanel, plotlegend, commonlegend,
                              toplabel, leftlabel, bottomlabel, withset,
                              layout, returnplots, outdir) {

  if (withset) {
    setpos <- grep("^set", names(data))
    if ( length(setpos)==0 ) {stop("There is no set column to perform facet_grid")}
    d1 <- tidyr::gather(data[ ,c(1:(nread+2),setpos)], temperature, reading, -id, -set, -condition)
    d1 <- d1 %>% complete(id, set, condition, temperature)
  } else {
    d1 <- tidyr::gather(data[ ,c(1:(nread+2))], temperature, reading, -id, -condition)
  }

  d1$condition <- as.factor(d1$condition)
  d1$temperature <- as.numeric(as.character(d1$temperature))
  xmin <- min(d1$temperature)
  xmax <- max(d1$temperature)
  xpos <- (xmax-xmin)/4 + xmin

  y <- c(1.0, 0.98, 0.95, 0.8, 0.7, 0.5, 0.5, 0.2, 0.15, 0.1)
  x <- c(37.0, 40.0, 43.0, 46.0, 49.0, 52.0, 55.0, 58.0, 61.0, 64.0)
  df <- data.frame(x, y)
  t1 <- drc::drm(y~x, data=df, fct=LL.4())
  #print(t1)

  if (topasone==TRUE) { top <- 1.0 } else { top <- NA }

  zz <- file(paste0(outdir,"/","curve_fitting_record.txt"), open="wt")
  sink(zz, type = "message")
  on.exit(sink(type="message"))

  plotting <- function(d1) {

    q <- ggplot()

    if (printcount) {
      PSM_label <- PSM_annod[PSM_annod$id==unique(d1$id), ][ ,c(2:ncol(PSM_annod))]
      Pep_label <- Pep_annod[Pep_annod$id==unique(d1$id), ][ ,c(2:ncol(Pep_annod))]
    }

    if (dotconnect==FALSE) {
      # q <- q + geom_smooth(method = "drm", formula = t1,
      #                      fct = LL.4(fixed=c(NA,NA,top,NA)),
      #                      aes(colour=condition), se=FALSE, size=0.5)
      # q <- q + geom_smooth(method=drc::"drm", formula=t1,
      #                      method.args=list(fct=LL.4(fixed=c(NA,NA,top,NA))),
      #                      aes(colour=condition), se=FALSE, size=0.5)
      d1 <- droplevels(d1)
      d1_list <- split(d1, d1$condition)
      for (j in names(d1_list)) {
        q <- q + geom_point(data=d1_list[[j]], aes(x=temperature, y=reading), colour=colorpanel[j], shape=20, size=2)
        d1_list[[j]] <- na.omit(d1_list[[j]])
        model.drm <- try(drc::drm(reading~temperature, data=d1_list[[j]], fct=LL.4(fixed=c(NA,NA,top,NA)),
                          na.action=na.omit, control=drmc(noMessage=TRUE)), silent=TRUE)#, outFile=zz)
        if (class(model.drm) != "try-error") {
          mda <- data.frame(temperature1=seq(xmin,xmax,length.out=100))
          mda$reading <- predict(model.drm, mda)
          q <- q + geom_line(data=mda, aes(x=temperature1, y=reading), colour=colorpanel[j], size=0.5)
        } else {
          q <- q + geom_smooth(data=d1_list[[j]], method="lm", formula=y~poly(x,1), aes(x=temperature, y=reading),
                               colour=colorpanel[j], se=FALSE, size=0.5)
        }
      }
    } else {
      q <- ggplot(d1, aes(x=temperature,y=reading,group=condition,colour=condition)) +
        geom_point(aes(colour=condition), shape=20, size=2) + geom_line(aes(colour=condition), size=0.5) +
        scale_colour_manual(drop=FALSE, values=colorpanel)
    }

    q <- q + scale_x_continuous(breaks=sort(unique(d1$temperature), decreasing=FALSE)) +
      coord_cartesian(xlim=c(xmin,xmax), ylim = c(0,1.2)) + theme_classic() +
      scale_y_continuous(breaks=seq(0,1,0.2)) + labs(x=" ", y=" ") +
      ggtitle(as.character(unique(d1$id)))

    if (withset) { q <- q + facet_grid(set~., drop=FALSE) }

    if (printcount) {
      q <- q + annotate("text", x=xpos, y=annotypos, label="  #PSM  #Peptides", size=1.5)
      for (i in 1:length(PSM_label)) {
        q <- q + annotate("text", x=xpos, y=annotypos-annotyinterval*i,
                          label=paste0(names(PSM_label)[i], ":  ", as.character(PSM_label)[i], " "),
                          size=1.5, hjust=1)
        i <- i + 1
      }
      for (i in 1:length(Pep_label)) {
        q <- q + annotate("text", x=xpos+2, y=annotypos-annotyinterval*i,
                          label=as.character(Pep_label)[i], size=1.5)
        i <- i + 1
      }
    }

    if (commonlegend==FALSE) {
      q <- q + theme(
        text = element_text(size=10),
        plot.title = element_text(hjust=0.5, size=rel(0.7)),
        legend.background = element_rect(fill=NULL),
        legend.key.height = unit(0.05, "cm"),
        legend.key.width = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.4)),
        legend.justification = "center",
        legend.position = c(0.8,0.9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        aspect.ratio = 1
      )
    } else {
      q <- q + theme(
        text = element_text(size=10),
        plot.title = element_text(hjust=0.5, size=rel(0.7)),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        aspect.ratio = 1
      )
    }
    return(q)
  }

  plots <- plyr::dlply(d1, .(id), .fun=plotting)
  #plots[[1]]
  zz <- paste0(outdir,"/","curve_fitting_record.txt")
  if (file.exists(zz)) { file.remove(zz) }
  if (returnplots) { return(plots) }
  legend <- plotlegend$legend
  lheight <- plotlegend$lheight
  params <- list(nrow=layout[1], ncol=layout[2])
  n <- with(params, nrow*ncol)
  ## add one page if division is not complete
  pages <- length(plots) %/% n + as.logical(length(plots) %% n)
  groups <- split(seq_along(plots), gl(pages, n, length(plots)))

  if (commonlegend==FALSE) {
    pl <- lapply(names(groups), function(g) {
      do.call(gridExtra::arrangeGrob, c(plots[groups[[g]]], params,
                                        top=toplabel,
                                        left=leftlabel,
                                        bottom=bottomlabel))
    })
  } else {
    pl <- lapply(names(groups), function(g) {
      gridExtra::grid.arrange(
        do.call(gridExtra::arrangeGrob,
                c(plots[groups[[g]]], params,
                  top=toplabel,
                  left=leftlabel,
                  bottom=bottomlabel)),
                  #bottom=expression(Temperature~(R^{o}~C))
        legend,
        ncol = 1,
        heights = grid::unit.c(unit(1, "npc") - lheight, lheight))
    })
  }
  class(pl) <- c("arrangelist", "ggplot", class(pl))
  return(pl)
}


#' ms_complex_melt_innerplot
#'
#' Internal function to generate ggplot objects for CETSA melt curve data grouped by complex
#'
#' @param data melt curve dataset with complex information to plot
#' @param nread number of reading channels or sample treatements
#' @param topasone whether the top plateau has to be fixed, i.e., 1.0
#' @param colorpanel a vector of customizable color scheme provided by the user
#' @param bycondition whether to color by condition
#' @param toplabel textual label at the top of the page
#' @param leftlabel textual label at the left side of the page
#' @param bottomlabel textual label at the bottom of the page
#' @param layout a vector indicating the panel layout for multi-panel plots per page
#'
#'
#' @import tidyr
#' @import scales
#' @import drc
#' @import RColorBrewer
#' @importFrom plyr . dlply
#' @importFrom grid unit.c
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @import ggplot2
#'
#' @keywords internal
#'
#' @return a list of ggplot2 object
#' @examples \dontrun{
#' }
#'
#'

ms_complex_melt_innerplot <- function(data, nread, topasone, printpvalue, colorpanel, bycondition,
                                      toplabel, leftlabel, bottomlabel, layout) {

  if (length(unique(data$condition))>1) {
    data <- data %>% rowwise() %>% mutate(name=paste(paste0("Complex ID: ",ComplexID),
                                                     ComplexName,paste0("Total subunit number: ", subunitsNum),sep="\n"))
  } else {
    data <- data %>% rowwise() %>% mutate(name=paste(paste0("Complex ID: ",ComplexID),
                                                     ComplexName,paste0("Total subunit number: ", subunitsNum),
                                                     paste0("Identified subunit number: ", subunitsIdentifiedNum),sep="\n"))
  }

  colinclusion <- c("name","ComplexID","ComplexName","subunitsUniProt_IDs","condition",
                    "subunitsNum","subunitsIdentifiedNum","subunitsIdentifiedPerc",
                    "Tm.mean","Tm.sd","Tm.cv","distance.mean","distance.sd","distance.cv","p.value")
  colinclusion <- intersect(colnames(data), colinclusion)
  complexinfo <- unique(data[,colinclusion])
  complexinfo <- complexinfo[order(complexinfo$p.value,decreasing=F), ]
  proteininfo <- unique(data[,c("id","description")])

  colinclusion <- c("ComplexID","ComplexName","subunitsUniProt_IDs","subunitsNum",
                    "id","description","sumUniPeps","sumPSMs","countNum",
                    "subunitsIdentifiedNum","subunitsIdentifiedPerc",
                    "Tm.mean","Tm.sd","Tm.cv","distance.mean","distance.sd",
                    "distance.cv","p.value","Tm","R2","Slope","RSE")
  data1 <- data[ ,!(names(data) %in% colinclusion)]
  # print(head(data1))
  d1 <- tidyr::gather(data1, temperature, reading, -name, -condition, -gene)
  if (length(unique(data$condition))>1) {
    d1 <- d1 %>% complete(name, condition, temperature)
  }
  # print(head(d1))
  d1$condition <- as.factor(d1$condition)
  d1$temperature <- as.numeric(as.character(d1$temperature))
  d1$name <- factor(d1$name, levels=unique(complexinfo$name))

  xmin <- min(d1$temperature)
  xmax <- max(d1$temperature)
  xpos <- (xmax-xmin)/4 + xmin

  annotypos <- 0.3
  if (length(colorpanel)==0) {
    colorpanel <- unique(c(brewer.pal(8, "Set1")[c(1:5,7:8)], brewer.pal(7, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(11,"Spectral")))
    colorpanel <- c(colorpanel, colorpanel, colorpanel, colorpanel)
  }

  y <- c(1.0, 0.98, 0.95, 0.8, 0.7, 0.5, 0.5, 0.2, 0.15, 0.1)
  x <- c(37.0, 40.0, 43.0, 46.0, 49.0, 52.0, 55.0, 58.0, 61.0, 64.0)
  df <- data.frame(x, y)
  t1 <- drc::drm(y~x, data=df, fct=drc::LL.4())
  #print(t1)

  if (topasone==TRUE) { top <- 1.0 } else { top <- NA }
  plotting <- function(d1) {

    if (printpvalue & !bycondition) {
      pvalue_label <- complexinfo[complexinfo$name==unique(d1$name), ]
      if (length(unique(pvalue_label$condition))>1) {
        pvalue_label <- pvalue_label[order(pvalue_label$condition), ]
      }
    }

    if (bycondition) {
      q <- ggplot(d1, aes(x=temperature,y=reading,group=interaction(gene,condition),colour=condition)) +
        geom_point(aes(colour=condition), shape=20, size=2)
    } else {
      q <- ggplot(d1, aes(x=temperature,y=reading,group=gene,colour=gene)) +
        geom_point(aes(colour=gene), shape=20, size=2)
    }
    q <- q + scale_x_continuous(breaks=sort(unique(d1$temperature), decreasing=FALSE)) +
      coord_cartesian(xlim=c(xmin,xmax), ylim = c(0,1.2)) + theme_classic() +
      scale_y_continuous(breaks=seq(0,1,0.2)) + labs(x=" ", y=" ") +
      scale_colour_manual(drop=FALSE, values=colorpanel) +
      ggtitle(as.character(unique(d1$name)))

    if (length(unique(d1$condition))>1 & !bycondition) { q <- q + facet_grid(condition~., drop=FALSE) }

    if (bycondition) {
      q <- q + geom_smooth(method=drc::"drm", formula=t1,
                           method.args=list(fct=drc::LL.4(fixed=c(NA,NA,top,NA))),
                           aes(colour=condition), se=FALSE, size=0.5)
    } else{
      q <- q + geom_smooth(method=drc::"drm", formula=t1,
                           method.args=list(fct=drc::LL.4(fixed=c(NA,NA,top,NA))),
                           aes(colour=gene), se=FALSE, size=0.5)
    }

    if (printpvalue & !bycondition) {
      q <- q + annotate("text",x=xpos,y=annotypos,label=paste0("p value: ",
                       formatC(pvalue_label$p.value,digits=2,format="e")),size=2)
    }

    q <- q + theme(
      text = element_text(size=8),
      plot.title = element_text(hjust=0.5, size=rel(1)),
      legend.background=element_rect(fill=NULL),
      legend.key.height=unit(0.05, "cm"),
      legend.key.width=unit(0.1, "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = rel(1)),
      legend.justification="center",
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      axis.line.x = element_line(),
      axis.line.y = element_line(),
      aspect.ratio = 1
    )
    return(q)
  }

  plots <- plyr::dlply(d1, plyr::.(name), .fun=plotting)
  #plots[[1]]
  params <- list(nrow=layout[1], ncol=layout[2])
  n <- with(params, nrow*ncol)
  ## add one page if division is not complete
  pages <- length(plots) %/% n + as.logical(length(plots) %% n)
  groups <- split(seq_along(plots), gl(pages, n, length(plots)))

  pl <- lapply(names(groups), function(g) {
    do.call(gridExtra::arrangeGrob, c(plots[groups[[g]]], params,
                                      top=toplabel,
                                      left=leftlabel,
                                      bottom=bottomlabel))
  })
  class(pl) <- c("arrangelist", "ggplot", class(pl))
  return(pl)
}
