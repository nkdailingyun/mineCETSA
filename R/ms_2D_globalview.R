#' ms_2D_globalview
#'
#' Function to have a global overview of the 2D CETSA dataset, regarding to
#' expression level changes and thermal shifts
#'
#' @param data dataset after ms_2D_caldiff() function
#' @param set a single character to indicate the sample name
#' @param treatment a single character to indicate the sample name
#' @param cvthreshold the CV threshold value for subsetting reproducible
#' measurements, default value is 0.1
#' @param corrthreshold the Correlation threshold value for subsetting
#' reproducible measurements, default value is 0.5
#' @param basetemp the character indicating the baseline temperature for
#' expression levels, default value is 37C
#' @param expressionchange_nMAD the number of MADs to set the significance
#' cutoff on protein expression level change, default value is 2.5
#' @param thermalchange_nMAD the number of MADs to set the significance cutoff
#' on protein thermal shift, default value is 2.5
#' @param labelnodes whether to text-label the selected nodes, default set to FALSE
#' @param labelcategory the categories of nodes to label, default value is c("CC","NC","CN")

#'
#' @import dplyr Biobase
#' @import ggpubr
#' @import ggrepel
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- ms_2D_globalview(MOLM, set="M13", treatment="AT26533")
#' }
#'
ms_2D_globalview <- function(data, set=NULL, treatment=NULL, cvthreshold=0.1, corrthreshold=0.5,
                             basetemp="37C", expressionchange_nMAD=2.5, thermalchange_nMAD=2.5,
                             labelnodes=FALSE, labelcategory=c("CC","NC","CN")) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)

  if (length(set)>0) { data <- data[ ,c(1,2,grep(set,names(data)))] }
  if (length(treatment)!=1) {stop("Provide only one treatment keyword for globalview")}
  else { data <- data[ ,c(1,2,grep(treatment,names(data)))] }

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
  }

  datal <- tidyr::gather(data[ ,-2], treatment, reading, -id)
  if (length(set)>0) {
    datal1 <- tidyr::separate(datal, treatment, into=c("set","temperature","replicate","treatment"), sep="_")
    datal <- datal1 %>% group_by(id, set, treatment, temperature) %>%
      summarize(mreading=mean(reading,na.rm=T), sdreading=sd(reading, na.rm=T))
    reproducible1 <- datal %>% group_by(id) %>% summarise(cvreading=sqrt(exp(mean(sdreading,na.rm=T)^2)-1)) %>%
      filter(cvreading<=cvthreshold)
    reproducible2 <- tidyr::spread(datal1, replicate, reading)
    reproducible2 <- plyr::ddply(reproducible2, "id", function(data) {
      a <- cor(data[ ,-c(1:4)], use="complete.obs")
      data.frame(corr=mean(a[lower.tri(a)]))
    })
    reproducible2 <- subset(reproducible2, corr>=corrthreshold)

    datal <- subset(datal, id %in% unique(c(reproducible1$id,reproducible2$id)))
    datal_Cchange <- datal %>% group_by(id, treatment) %>% summarize(change=max(mreading,na.rm=T)-min(mreading,na.rm=T))
    datal_Echange <- datal[grep(basetemp, datal$temperature), c(1,2,3,5)]
    datal_change <- na.omit(merge(datal_Echange, datal_Cchange))
    names(datal_change)[c(4,5)] <- c("expressionchange", "thermalchange")
  } else {
    datal1 <- tidyr::separate(datal, treatment, into=c("temperature","replicate","treatment"), sep="_")
    datal <- datal1 %>% group_by(id, treatment, temperature) %>%
      summarize(mreading=mean(reading,na.rm=T), sdreading=sd(reading, na.rm=T))
    reproducible1 <- datal %>% group_by(id) %>% summarise(cvreading=sqrt(exp(mean(sdreading,na.rm=T)^2)-1)) %>%
      filter(cvreading<=cvthreshold)
    reproducible2 <- tidyr::spread(datal1, replicate, reading)
    reproducible2 <- plyr::ddply(reproducible2, "id", function(data) {
      a <- cor(data[ ,-c(1:3)], use="complete.obs")
      data.frame(corr=mean(a[lower.tri(a)]))
    })
    reproducible2 <- subset(reproducible2, corr>=corrthreshold)

    datal <- subset(datal, id %in% unique(c(reproducible1$id,reproducible2$id)))
    datal_Cchange <- datal %>% group_by(id, treatment) %>% summarize(change=max(mreading,na.rm=T)-min(mreading,na.rm=T))
    datal_Echange <- datal[grep(basetemp, datal$temperature), c(1,2,4)]
    datal_change <- na.omit(merge(datal_Echange, datal_Cchange))
    names(datal_change)[c(3,4)] <- c("expressionchange", "thermalchange")
  }
    expressioncutoff <- median(datal_change$expressionchange) + expressionchange_nMAD*mad(datal_change$expressionchange)
    thermalcutoff <- median(datal_change$thermalchange) + expressionchange_nMAD*mad(datal_change$thermalchange)

    datal_change <- datal_change %>% rowwise() %>%
      mutate(expression = ifelse(abs(expressionchange)<expressioncutoff, "N", "C")) %>%
      #mutate(basechangedir = ifelse(mreading < 0, "-", "+")) %>%
      mutate(thermal = ifelse(thermalchange<thermalcutoff, "N", "C")) %>%
      mutate(category=paste0(expression, thermal))

    datal_change <- proteininfo %>% rowwise() %>% mutate(gene=getGeneName(description)) %>%
      inner_join(datal_change) %>% arrange(category)

    print(paste0("The category of expression level change and thermal shift are as follows: "))
    print(table(datal_change$category))

  q <- ggpubr::ggscatter(datal_change, x = "expressionchange", y = "thermalchange",
                         color = "category", shape=20, alphla=0.1,
                         palette = c("#FC4E07", "#00AFBB", "#E7B800", "gray"),
                         title = paste0("Changes in ",treatment),
                         xlab = "37C expression level change",
                         ylab = "Delta fold change across temperature")
  if (labelnodes) {
    q <- q + ggrepel::geom_text_repel(data=subset(datal_change, category %in% labelcategory),
                                      aes(label=gene))
  }
  q <- q + geom_hline(yintercept=thermalcutoff, linetype="dashed", color="black") +
    geom_vline(xintercept=-expressioncutoff, linetype="dashed", color="black") +
    geom_vline(xintercept=expressioncutoff, linetype="dashed", color="black") +
    theme(text = element_text(size=8), plot.title = element_text(hjust=0.5, size=rel(1.2)), aspect.ratio=0.5)

  ggsave(file=paste0(format(Sys.time(), "%y%m%d_%H%M"), "_Changes_in_", treatment, ".pdf"), q, width=11.69, height=8.27)
  ms_filewrite(datal_change, paste0("Changes_in_", treatment, ".txt"), outdir=outdir)
  return(datal_change)
}
