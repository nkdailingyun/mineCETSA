#' ms_2D_reproducible
#'
#' Function to retrieve a reproducible subset of the 2D CETSA dataset, regarding to
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
#'
#' @import dplyr Biobase
#' @import ggpubr
#' @import ggrepel
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- ms_2D_reproducible(MOLM)
#' }
#'

ms_2D_reproducible <- function(data, set=NULL, treatment=NULL,
                             cvthreshold=0.1, corrthreshold=0.5,
                             basetemp="37C") {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (length(grep("countNum", names(data)))) {
    countinfo <- unique(data[ ,c("id","sumUniPeps","sumPSMs","countNum")])
    data <- data[ ,!(names(data) %in% c("sumUniPeps","sumPSMs","countNum"))]
  }

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
  }

  if (length(set)>0) { data <- data[ ,c(1,2,grep(set,names(data)))] }

  datal <- tidyr::gather(data[ ,-2], treatment, reading, -id)
  if (length(set)>0) {
    datal1 <- tidyr::separate(datal, treatment, into=c("set","temperature","replicate","treatment"), sep="_")
    datal2 <- list()
    for (i in unique(datal1$treatment)) {
      datal_temp <- subset(datal1, treatment==i)
      if (sum(abs(datal_temp$reading),na.rm=T)==0) {next}
      print(paste0("in condition: ", i))
      datal_score <- datal_temp %>% group_by(id, set, treatment, temperature) %>%
        summarize(mreading=mean(reading,na.rm=T), sdreading=sd(reading, na.rm=T))
      reproducible1 <- datal_score %>% group_by(id) %>% summarise(cvreading=sqrt(exp(mean(sdreading,na.rm=T)^2)-1)) %>%
        filter(cvreading<=cvthreshold)
      print(paste0(nrow(reproducible1), " proteins passed the cv cutoff..."))
      reproducible2 <- tidyr::spread(datal_temp, replicate, reading)
      reproducible2 <- plyr::ddply(reproducible2, "id", function(data) {
        a <- cor(data[ ,-c(1:4)], use="complete.obs")
        data.frame(corr=mean(a[lower.tri(a)]))
      })
      reproducible2 <- subset(reproducible2, corr>=corrthreshold)
      print(paste0(nrow(reproducible2), " proteins passed the correlation cutoff..."))
      datal_temp <- subset(datal_temp, id %in% unique(c(reproducible1$id,reproducible2$id)))
      datal_temp <- tidyr::unite(datal_temp, condition, temperature, replicate, treatment, sep="_")
      datal_temp <- tidyr::spread(datal_temp, condition, reading)
      datal_temp <- merge(datal_temp, countinfo)
      datal_temp <- merge(proteininfo, datal_temp)
      if (length(attr(datal_temp,"outdir"))==0 & length(outdir)>0) {
        attr(datal_temp,"outdir") <- outdir
      }
      datal2[[i]] <- datal_temp
    }
  } else {
    datal1 <- tidyr::separate(datal, treatment, into=c("temperature","replicate","treatment"), sep="_")
    datal2 <- list()
    for (i in unique(datal1$treatment)) {
      datal_temp <- subset(datal1, treatment==i)
      if (sum(abs(datal_temp$reading),na.rm=T)==0) {next}
      print(paste0("in condition: ", i))
      datal_score <- datal_temp %>% group_by(id, treatment, temperature) %>%
        summarize(mreading=mean(reading,na.rm=T), sdreading=sd(reading, na.rm=T))
      reproducible1 <- datal_score %>% group_by(id) %>% summarise(cvreading=sqrt(exp(mean(sdreading,na.rm=T)^2)-1)) %>%
        filter(cvreading<=cvthreshold)
      print(paste0(nrow(reproducible1), " proteins passed the cv cutoff..."))
      reproducible2 <- tidyr::spread(datal_temp, replicate, reading)
      reproducible2 <- plyr::ddply(reproducible2, "id", function(data) {
        a <- cor(data[ ,-c(1:3)], use="complete.obs")
        data.frame(corr=mean(a[lower.tri(a)]))
      })
      reproducible2 <- subset(reproducible2, corr>=corrthreshold)
      print(paste0(nrow(reproducible2), " proteins passed the correlation cutoff..."))
      datal_temp <- subset(datal_temp, id %in% unique(c(reproducible1$id,reproducible2$id)))
      datal_temp <- tidyr::unite(datal_temp, condition, temperature, replicate, treatment, sep="_")
      datal_temp <- tidyr::spread(datal_temp, condition, reading)
      datal_temp <- merge(datal_temp, countinfo)
      datal_temp <- merge(proteininfo, datal_temp)
      if (length(attr(datal_temp,"outdir"))==0 & length(outdir)>0) {
        attr(datal_temp,"outdir") <- outdir
      }
      datal2[[i]] <- datal_temp
    }
  }
  return(datal2)
}
