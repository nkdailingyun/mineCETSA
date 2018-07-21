#' ms_2D_rearrange
#'
#' Function to rearrange 2D-CETSA data into a wide one-unique-protein-per-row format
#'
#' @param data dataset to transform into a wide format
#' @param nread number of reading channels, default value 9 in 2D-CETSA scheme
#' @param repthreshold the minimal percentage threshold of protein being sampled
#' from multiple runs, default value is 0.75
#' @param with37Creading whether the kept proteins should have readings at 37C
#' @param averagecount whether to average the supporting PSM/peptide/count numbers,
#' default set to TRUE
#' @param countthreshold the minimal threshold number of associated abundance
#' count of proteins, default value is 2
#'
#' @import dplyr tidyr
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'  MOLM <- ms_rearrange(MOLM, nread=9, countthreshold = 3)
#' }
#'
#'


ms_2D_rearrange <- function(data, nread=9, repthreshold=0.75, with37Creading=TRUE,
                            averagecount=TRUE, countthreshold=2) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (repthreshold>0) {
    counttable <- count(data, id) %>% mutate(freq=n/max(n)) %>% filter(freq>=repthreshold)
    data <- subset(data, id %in% counttable$id)
    print(paste0(nrow(counttable), " proteins pass the measurement replicates cutoff ", repthreshold*100, "%."))
  }

  if (with37Creading) {
    cond37C <- unique(data$condition)[grep("37C", unique(data$condition))]
    data37C <- subset(data, condition==cond37C)
    print(paste0(nrow(data37C), " proteins are measured at 37C and they are kept."))
    data <- subset(data, id %in% data37C$id)
  }

  d1 <- tidyr::gather(data[ ,c(1,3,4:(3+nread))], treatment, reading, -id, -condition)
  d1 <- tidyr::unite(d1, combinedcol, condition, treatment, sep="_")
  d1 <- tidyr::spread(d1, combinedcol, reading)
  data1 <- unique(data[ ,c(1,2)]) %>% inner_join(d1) #%>% rowwise() %>% mutate(gene=getGeneName(description))

  if (averagecount) {
    counttable <- NULL
    peppos <- grep("^sumUniPeps", names(data), value=F)
    psmpos <- grep("^sumPSMs", names(data), value=F)
    countpos <- grep("^countNum", names(data), value=F)
    pos <- c(peppos, psmpos, countpos)
    if (length(pos) > 1) {
      counttable <- group_by(data[ ,c(1,pos)], id) %>%
        summarize(sumUniPeps=median(sumUniPeps), sumPSMs=median(sumPSMs), countNum = median(countNum))
      if (countthreshold > 0) {
        counttable <- subset(counttable, countNum>=countthreshold)
        print(paste0(nrow(counttable), " proteins pass the count number cutoff ", countthreshold, "."))
      }
    }
    data1 <- merge(data1, counttable)
  }

  if (length(attr(data1,"outdir"))==0 & length(outdir)>0) {
    attr(data1,"outdir") <- outdir
  }
  mineCETSA::ms_filewrite(data1, paste0(dataname, "_data_pre_normalization.txt"),
                          outdir=outdir)
  return(data1)
}
