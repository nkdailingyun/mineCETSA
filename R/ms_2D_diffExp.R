#' ms_2D_diffExp
#'
#' Function to have an overview of the 2D CETSA dataset, regarding to
#' expression level changes only, to select out the differentially expressed proteins
#'
#' @param data dataset after ms_2D_normalization function, readings in log2 format
#' @param set a single character to indicate the sample name
#' @param basetemp the character indicating the baseline temperature for
#' expression levels, default value is 37C
#' @param contrast a character to indicate the contrasting treatment conditions
#' @param logFC_threshold the threshold value for log fold changes, default set at 0.2
#' @param adjp_threshold the threshold value for adjusted p values, default set at 0.01
#' @param labelnodes whether to label the proteins with significant differential expression
#' @param xlimit a two numeric element vector to indicate the limit of x-axis
#' @param ylimit a two numeric element vector to indicate the limit of y-axis

#'
#' @import dplyr Biobase
#' @import ggplot2
#' @import ggrepel
#' @export
#' @return a dataframe
#' @examples \dontrun{
#'   MOLM <- ms_2D_diffExp(MOLM, set="M13", contrast="TNFa-DMSO")
#' }
#'
ms_2D_diffExp <- function(data, set=NULL, basetemp="37C", contrast=NULL,
                          logFC_threshold=0.2, adjp_threshold=0.01,
                          labelnodes=TRUE, xlimit=NULL, ylimit=NULL, returneset=FALSE) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)$outdir
  data <- ms_directory(data, dataname)$data

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
  }
  if (length(grep("countNum", names(data)))) {
    countinfo <- unique(data[ ,c("id","sumUniPeps","sumPSMs","countNum")])
  }

  if (length(set)>0) { data <- data[ ,c(1,2,grep(set,names(data)))] }
  # if (length(treatment)!=1) {stop("Provide only one treatment keyword for globalview")}
  # else { data <- data[ ,c(1,2,grep(treatment,names(data)))] }
  if (sum(grepl(basetemp,names(data)))) { data <- data[ ,c(1,2,grep(basetemp,names(data)))] }
  else {stop("Make sure the basetemp info is embeded in the column names.")}

  cname <- setdiff(names(data), c("id","description","sumUniPeps","sumPSMs","countNum"))
  if (length(unlist(strsplit(cname[1], "_")))==3) {
    temperature <- unlist(lapply(strsplit(cname, "_"),`[`,1))
    replicate <- unlist(lapply(strsplit(cname, "_"),`[`,2))
    treatment <- unlist(lapply(strsplit(cname, "_"),`[`,3))
    pdata1 <- data.frame(temperature=temperature, replicate=replicate, treatment=treatment)
    row.names(pdata1) <- cname
    nread <- nrow(pdata1)
    #print(pdata1)
    data <- merge(data, countinfo)
    data_eset <- ms_to_eSet(data=data, nread=nread, pdata=pdata1)
  } else if (length(unlist(strsplit(cname[1], "_")))==4) {
    set <- unlist(lapply(strsplit(cname, "_"),`[`,1))
    temperature <- unlist(lapply(strsplit(cname, "_"),`[`,2))
    replicate <- unlist(lapply(strsplit(cname, "_"),`[`,3))
    treatment <- unlist(lapply(strsplit(cname, "_"),`[`,4))
    pdata1 <- data.frame(set=set, temperature=temperature, replicate=replicate, treatment=treatment)
    row.names(pdata1) <- cname
    nread <- nrow(pdata1)
    #print(pdata1)
    data <- merge(data, countinfo)
    data_eset <- ms_to_eSet(data=data, nread=nread, pdata=pdata1)
  } else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }

  if (returneset) {return(data_eset)}

  ct <- pdata1$treatment
  design <- model.matrix(~0+ct)
  colnames(design) <- levels(ct)
  #print(design)

  if (length(contrast)) {
    contrast.matrix <- makeContrasts(contrasts=contrast, levels=design)
  } else {
    stop("pls specify a contrast expression, such as 'TNFa-DMSO'.")
  }
  fit <- lmFit(data_eset, design)
  fit1 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit1)#, proportion=0.2, trend=T)
  fittoptable <- NULL
  for (i in 1:length(contrast)) {
    contrast.name <- colnames(fit2$contrasts)[i]
    top <- topTable(fit2, coef=i, number=Inf, adjust="BH", sort.by="p")
    top <- tibble::rownames_to_column(top, "id")
    top$contrast <- contrast.name
    top$category <- ifelse(abs(top$logFC)>=logFC_threshold & top$adj.P.Val<=adjp_threshold, "C", "N")
    print(paste0("The category of expression level change in ", contrast.name, " are as follows: "))
    print(table(top$category))
    write.csv(top, paste0(format(Sys.time(), "%y%m%d_%H%M_"), "_", dataname, "_", contrast.name, "_eBays.csv"), row.names=F)
    fittoptable <- rbind(fittoptable, top)
  }

  fittoptable <- fittoptable %>% rowwise() %>%
    mutate(gene = getGeneName(description), log10p=-log10(adj.P.Val))

  q <- ggpubr::ggscatter(fittoptable, x = "logFC", y = "log10p",
                         color = "category", shape=20, alphla=0.1,
                         palette = c("#FC4E07", "gray"),
                         title = paste0("Expression Changes of ",contrast),
                         xlab = "fold change of 37C expression level [log2]",
                         ylab = "adjusted p values [-log10]")
  if (length(xlimit)) { q <- q + coord_cartesian(xlim=xlimit)}
  if (length(ylimit)) { q <- q + coord_cartesian(ylim=ylimit)}
  if (labelnodes) {
    q <- q + ggrepel::geom_text_repel(data=subset(fittoptable, category=="C"),
                                      aes(label=gene))
  }
  q <- q + theme(text = element_text(size=12), plot.title = element_text(hjust=0.5, size=rel(1.2)), aspect.ratio=1)

  ggsave(file=paste0(format(Sys.time(), "%y%m%d_%H%M"), "_Changes_in_", contrast, ".pdf"), q, width=8.27, height=11.69)

  return(fittoptable)
}
