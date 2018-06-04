#' ms_2D_DE
#'
#' Function to have an overview of the 2D CETSA dataset, regarding to
#' expression level changes only, to select out the differentially expressed ones
#'
#' @param data dataset after ms_2D_caldiff() function
#' @param set a single character to indicate the sample name
#' @param treatment a single character to indicate the sample name
#' @param cvthreshold the CV threshold value for subsetting reproducible
#' measurements, default value is 0.1
#' @param basetemp the character indicating the baseline temperature for
#' expression levels, default value is 37C
#' @param expressionchange_nMAD the number of MADs to set the significance
#' cutoff on protein expression level change, default value is 2.5
#' @param expressionchange_cutoff a two numeric element vector to indicate
#'the significance cutoff on protein expression level change
#' @param densityformat whether plot the histogram in density plot format
#' @param binwidth the width of bin used in histogram plot
#' @param xlimit a two numeric element vector to indicate the limit of x-axis

#'
#' @import dplyr Biobase
#' @import ggplot2
#' @import ggrepel
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- ms_2D_DE(MOLM, set="M13", treatment="AT26533")
#' }
#'
ms_2D_DE <- function(data, set=NULL, treatment=NULL, cvthreshold=0.1, basetemp="37C",
                     expressionchange_nMAD=2.5, expressionchange_cutoff=NULL,
                     densityformat=FALSE, binwidth=NULL, xlimit=NULL) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (length(set)>0) { data <- data[ ,c(1,2,grep(set,names(data)))] }
  if (length(treatment)!=1) {stop("Provide only one treatment keyword for globalview")}
  else { data <- data[ ,c(1,2,grep(treatment,names(data)))] }

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
  }

  datal <- tidyr::gather(data[ ,-2], treatment, reading, -id)
  if (length(set)>0) {
    datal <- na.omit(datal[grep(basetemp, datal$treatment), ])
    datal1 <- tidyr::separate(datal, treatment, into=c("set","temperature","replicate","treatment"), sep="_")
    datal <- datal1 %>% group_by(id, set, treatment, temperature) %>%
      summarize(mreading=mean(reading,na.rm=T), sdreading=sd(reading, na.rm=T))
    reproducible1 <- datal %>% group_by(id) %>% summarise(cvreading=sqrt(exp(mean(sdreading,na.rm=T)^2)-1)) %>%
      filter(cvreading<=cvthreshold)
    nrem <- length(unique(datal$id))-length(unique(reproducible1$id))
    print(paste0(nrem, " Proteins did not pass the CV threshold and were removed from downstream analysis."))
    datal_change <- merge(datal, reproducible1)[ ,c(1,2,3,4,5,7)]
    names(datal_change)[5] <- "expressionchange"
    names(datal_change)[6] <- "expressionchange_cv"
  } else {
    datal <- na.omit(datal[grep(basetemp, datal$treatment), ])
    datal1 <- tidyr::separate(datal, treatment, into=c("temperature","replicate","treatment"), sep="_")
    datal <- datal1 %>% group_by(id, treatment, temperature) %>%
      summarize(mreading=mean(reading,na.rm=T), sdreading=sd(reading, na.rm=T))
    reproducible1 <- datal %>% group_by(id) %>% summarise(cvreading=sqrt(exp(mean(sdreading,na.rm=T)^2)-1)) %>%
      filter(cvreading<=cvthreshold)
    nrem <- length(unique(datal$id))-length(unique(reproducible1$id))
    print(paste0(nrem, " Proteins did not pass the CV threshold and were removed from downstream analysis."))
    datal_change <- merge(datal, reproducible1)[ ,c(1,2,3,4,6)]
    names(datal_change)[4] <- "expressionchange"
    names(datal_change)[5] <- "expressionchange_cv"
  }
  if (length(expressionchange_cutoff)==0) {
    expressioncutoff <- median(datal_change$expressionchange,na.rm=T) + expressionchange_nMAD*mad(datal_change$expressionchange,na.rm=T)
    print(paste0("The cutoffs at ", expressionchange_nMAD, "*MAD are ", -round(expressioncutoff,2), " and ", round(expressioncutoff,2)))
    datal_change <- datal_change %>% rowwise() %>%
      mutate(expression = ifelse(abs(expressionchange)<expressioncutoff, "N", "C"))
  } else if (length(expressionchange_cutoff)==2) {
    expressioncutoff <- expressionchange_cutoff
    datal_change <- datal_change %>% rowwise() %>%
      mutate(expression = ifelse(expressionchange>expressioncutoff[1] & expressionchange<expressioncutoff[2], "N", "C"))
  } else {
    stop ("pls provide a two numeric element vector of expression level cutoff.")
  }

  datal_change <- proteininfo %>% rowwise() %>% mutate(gene=getGeneName(description)) %>%
    inner_join(datal_change) %>% arrange(expression)

  print(paste0("The category of expression level change are as follows: "))
  print(table(datal_change$expression))

  if (densityformat) {
    q <- ggplot2::ggplot(datal_change, aes(x=expressionchange)) +
      geom_histogram(aes(y=..density..), color="gray", fill="white", binwidth=binwidth) +
      geom_density(alpha=0.2, fill='#FF6666')
  } else {
    q <- ggplot2::ggplot(datal_change, aes(x=expressionchange)) +
      geom_histogram(binwidth=binwidth)
  }

  if (length(expressionchange_cutoff)==0) {
    q <- q + geom_vline(xintercept=-expressioncutoff, linetype="dashed", color="black") +
      geom_vline(xintercept=expressioncutoff, linetype="dashed", color="black") +
      theme(text = element_text(size=8), plot.title = element_text(hjust=0.5, size=rel(1.2)), aspect.ratio=1)
  } else if (length(expressionchange_cutoff)==2) {
    q <- q + geom_vline(xintercept=expressioncutoff[1], linetype="dashed", color="black") +
      geom_vline(xintercept=expressioncutoff[2], linetype="dashed", color="black") +
      theme(text = element_text(size=8), plot.title = element_text(hjust=0.5, size=rel(1.5)), aspect.ratio=1)
  }
  q <- q + ggtitle(paste0("Changes in ",treatment)) + xlab("37C expression level change")
  if (length(xlimit)==2) { q <- q + coord_cartesian(xlim=xlimit) }

  ggsave(file=paste0(format(Sys.time(), "%y%m%d_%H%M"), "_Differential_Expression_in_", treatment, ".pdf"), q, width=8.27, height=8.27)
  ms_filewrite(datal_change, paste0("Changes_in_", treatment, ".txt"), outdir=outdir)
  return(datal_change)
}
