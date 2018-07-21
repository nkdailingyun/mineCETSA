#' ms_2D_score
#'
#' Function to calculate the relative abundance score and thermal stability score
#'
#' @param data dataset after calculating the relative protein abundance differences
#' @param basetemp the character indicating the baseline temperature for
#' expression levels, default value is 37C
#' @param allzscore when calculating the z score, whether to use all the readings from
#' all the different treatment groups, default set to TRUE
#' @param fdrthreshold the significance level of global fdr, default value is 0.01
#' @param labelnodes whether to text-label the selected nodes, default set to FALSE
#' @param labelcategory the categories of nodes to label, default value is c("CC","NC","CN")
#'
#' @import dplyr fdrtool
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- ms_2D_score(MOLM)
#' }
#'
#'


ms_2D_score <- function(data, basetemp="37C", allzscore=TRUE, fdrthreshold=0.01,
                        labelnodes=FALSE, labelcategory=c("CC","NC","CN"),
                        xrange=NULL, yrange=NULL) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
    data$description <- NULL
  }
  if (length(grep("countNum", names(data)))) {
    countinfo <- unique(data[ ,c("id","sumUniPeps","sumPSMs","countNum")])
    data <- data[ ,!(names(data) %in% c("sumUniPeps","sumPSMs","countNum"))]
  }

  refcol <- which(apply(data[,-1],2,sum,na.rm=T)==0)
  data <- gather(data[ ,-(refcol+1)], condition, reading, -id)
  if (length(unlist(strsplit(data$condition[1], "_")))==4) {
    data <- tidyr::separate(data, condition, into=c("set","temperature","replicate","treatment"), sep="_")
    data_abd <- data %>% group_by(id, set, treatment) %>%
      filter(temperature==basetemp) %>%
      summarise(abundance.score=mean(reading,na.rm=T))
    data_thermal <- data %>% left_join(data_abd) %>% rowwise() %>%
      mutate(stability=reading-abundance.score) %>% filter(temperature!=basetemp)
    data_thermal1 <- data_thermal %>% group_by(id,set,treatment,temperature) %>%
      summarise(stability.score=mean(stability,na.rm=T))
    data_thermal1t <- na.omit(data_thermal) %>% group_by(id,set,treatment,temperature) %>%
      summarise(stability.score.t=t.test(stability)$statistic)
    data_thermal1 <- merge(data_thermal1, data_thermal1t,all.x=T)
    data_thermal2 <- data_thermal1 %>% group_by(id,set,treatment) %>%
      summarise(stability.score.mean=mean(stability.score,na.rm=T),stability.score.sd=sd(stability.score,na.rm=T))
    data_thermal1.score <- mutate(data_thermal1, temperature=paste("stability.score",temperature,sep="."))
    data_thermal1.score <- tidyr::spread(data_thermal1.score[,-5], temperature, stability.score)
    data_thermal1.t <- mutate(data_thermal1, temperature=paste0("stability.score.",temperature,".t"))
    data_thermal1.t <- tidyr::spread(data_thermal1.t[,-4], temperature, stability.score.t)
    data_thermal3 <- merge(data_thermal1.score, data_thermal1.t)
    data_thermal3 <- merge(data_thermal3, data_thermal2)
    data_score <- merge(data_abd, data_thermal3)
    #data1 <- tidyr::unite(data1, condition, set, temperature, replicate, treatment, sep="_")
  } else if (length(unlist(strsplit(data$condition[1], "_")))==3) {
    data <- tidyr::separate(data, condition, into=c("temperature","replicate","treatment"), sep="_")
    data_abd <- data %>% group_by(id, treatment) %>%
      filter(temperature==basetemp) %>%
      summarise(abundance.score=mean(reading,na.rm=T))
    data_thermal <- data %>% left_join(data_abd) %>% rowwise() %>%
      mutate(stability=reading-abundance.score) %>% filter(temperature!=basetemp)
    data_thermal1 <- data_thermal %>% group_by(id,treatment,temperature) %>%
      summarise(stability.score=mean(stability,na.rm=T))
    data_thermal1t <- na.omit(data_thermal) %>% group_by(id,treatment,temperature) %>%
      summarise(stability.score.t=t.test(stability)$statistic)
    data_thermal1 <- merge(data_thermal1, data_thermal1t,all.x=T)
    data_thermal2 <- data_thermal1 %>% group_by(id,treatment) %>%
      summarise(stability.score.mean=mean(stability.score,na.rm=T),stability.score.sd=sd(stability.score,na.rm=T))
    data_thermal1.score <- mutate(data_thermal1, temperature=paste("stability.score",temperature,sep="."))
    data_thermal1.score <- tidyr::spread(data_thermal1.score[,-5], temperature, stability.score)
    data_thermal1.t <- mutate(data_thermal1, temperature=paste0("stability.score.",temperature,".t"))
    data_thermal1.t <- tidyr::spread(data_thermal1.t[,-4], temperature, stability.score.t)
    data_thermal3 <- merge(data_thermal1.score, data_thermal1.t)
    data_thermal3 <- merge(data_thermal3, data_thermal2)
    data_score <- merge(data_abd, data_thermal3)
    #data_score <- tidyr::unite(data_score, condition, temperature, replicate, treatment, sep="_")
  } else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }

  # z-transformation of each score
  data_score_z <- NULL

  if (allzscore) {
    data_score_nocontrol <- data_score[0,]
    for (treat in unique(data_score$treatment)) {
      data_score1 <- subset(data_score, treatment==treat)
      #print(treat)
      if (sum(abs(data_score1$abundance.score),na.rm=T)==0) {next}
      else { data_score_nocontrol <- rbind(data_score_nocontrol, data_score1) }
    }
    col.names <- setdiff(names(data_score_nocontrol),c("id","treatment"))
    col.names <- col.names[-grep("\\.sd", col.names)]
    for (cn in col.names) {
      if (grepl("\\.t", cn)) {
        next
      } else {
        ratios <- data_score_nocontrol[ ,cn]
        quants <- quantile(ratios,probs=c(0.1587,0.5,0.8413),na.rm=T)
        ratios[is.na(ratios)] <- 0
        pratios <- ratios>=0
        nratios <- ratios<0
        z <- ratios
        z[nratios] <- ratios[nratios]/as.numeric(diff(quants)[1])
        z[pratios] <- ratios[pratios]/as.numeric(diff(quants)[2])
        data_score_nocontrol[ ,paste0(cn,".z")] <- z
        data_score_nocontrol[ ,paste0(cn,".fdr")] <- fdrtool(z,plot=F,verbose=F)$qval
      }
    }
      data_score_z <- data_score_nocontrol
  } else {
    for (treat in unique(data_score$treatment)) {
      data_score1 <- subset(data_score, treatment==treat)
      #print(treat)
      if (sum(abs(data_score1$abundance.score),na.rm=T)==0) {next}
      col.names <- setdiff(names(data_score1),c("id","treatment"))
      col.names <- col.names[-grep("\\.sd", col.names)]
      for (cn in col.names) {
        if (grepl("\\.t", cn)) {
          next
        } else {
          ratios <- data_score1[ ,cn]
          quants <- quantile(ratios,probs=c(0.1587,0.5,0.8413),na.rm=T)
          ratios[is.na(ratios)] <- 0
          pratios <- ratios>=0
          nratios <- ratios<0
          z <- ratios
          z[nratios] <- ratios[nratios]/as.numeric(diff(quants)[1])
          z[pratios] <- ratios[pratios]/as.numeric(diff(quants)[2])
          data_score1[ ,paste0(cn,".z")] <- z
          data_score1[ ,paste0(cn,".fdr")] <- fdrtool(z,plot=F,verbose=F)$qval
        }
      }
      data_score_z <- rbind(data_score_z, data_score1)
    }
  }
  data_score_z <- data_score_z %>% rowwise() %>% mutate(abundance.hit=ifelse(abundance.score.fdr<fdrthreshold,T,F))
  data_score_z <- data_score_z %>% rowwise() %>% mutate(stability.hit=ifelse(stability.score.mean.fdr<fdrthreshold,T,F))
  data_score_z <- data_score_z %>% mutate(category=paste0(as.character(abundance.hit),as.character(stability.hit)))
  data_score_z$category <- gsub("FALSE","N",data_score_z$category)
  data_score_z$category <- gsub("TRUE","C",data_score_z$category)
  data_score <- merge(proteininfo, data_score_z)
  print("The number of proteins in each category is as follows:")
  print(table(data_score$category, data_score$treatment))

  ms_filewrite(data_score, paste0(dataname, "_Abundance_Stability_score.txt"), outdir=outdir)

  data_score <- data_score %>% rowwise() %>% mutate(gene=getGeneName(description))
  q <- ggpubr::ggscatter(data_score, x = "abundance.score", y = "stability.score.mean",
                         color = "category", shape=20, alpha=0.9, facet.by="treatment",
                         palette = c("#FC4E07", "#00AFBB", "#E7B800", "gray"),
                         xlab = "Protein Abundance score",
                         ylab = "Protein Thermal stability score")
  if (labelnodes) {
    q <- q + ggrepel::geom_text_repel(data=subset(data_score, category %in% labelcategory),
                                      aes(label=gene))
  }
  q <- q + theme(text = element_text(size=12), plot.title = element_text(hjust=0.5, size=rel(1.2)), aspect.ratio=1)
  # q <- q + coord_cartesian(xlim=xrange, ylim=yrange)
  ggsave(file=paste0(format(Sys.time(), "%y%m%d_%H%M"), "_Abundance_Stability_score", ".pdf"), q, width=11.69, height=8.27)

  q <- ggpubr::ggscatter(data_score, x = "abundance.score.z", y = "stability.score.mean.z",
                         color = "category", shape=20, alpha=0.9, facet.by="treatment",
                         palette = c("#FC4E07", "#00AFBB", "#E7B800", "gray"),
                         xlab = "Protein Abundance z-score",
                         ylab = "Protein Thermal stability z-score")
  if (labelnodes) {
    q <- q + ggrepel::geom_text_repel(data=subset(data_score, category %in% labelcategory),
                                      aes(label=gene))
  }
  q <- q + theme(text = element_text(size=12), plot.title = element_text(hjust=0.5, size=rel(1.2)), aspect.ratio=1)
  # q <- q + coord_cartesian(xlim=xrange, ylim=yrange)
  ggsave(file=paste0(format(Sys.time(), "%y%m%d_%H%M"), "_Abundance_Stability_score_z", ".pdf"), q, width=11.69, height=8.27)

  return(data_score)
}
