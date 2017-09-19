#' ms_ITTR_filter
#'
#' Function to perform ITTR data segregation and positive hits selection based
#' on desired fc threshold level and the associated fitting parameters
#'
#' @param data dataset to be filtered
#' @param nread number of reading channels or sample treatements, default value
#' is 10
#' @param fcthreshold short for fold change threshold, indicate the threshold
#' fold change compared to baseline (i.e.,1.0) above which the readings can be
#' considered significantly stabilized or destabilized, default value 0.3
#' @param R2threshold minimal R2 level to consider time-reponse fitting as reliable
#' @param nbaseline the number of points (count from lowest concentration) used
#' for calculate the baseline MAD, default is 2
#' @param baselineMAD MAD of baseline variation, default value is 0; if not
#' provided, it will be calculated based on the readings from the lowest two dose groups
#' @param nMAD the significance level of MAD cutoff, default value is 2.5
#' @param checkreplicate whether to check replicated measurements, when set to
#' TRUE, make sure the proteins with all the replicated measurements from at least
#' one condition pass fold change threshold will be segregated, default set to FALSE
#' @param treatcondition specify the condition names to perform replicate check,
#' no need to specify the replicate part of the full condition names, default value is NULL
#' @param allpresence whether require the proteins be present (i.e., successfully
#' quantified) in all the different conditions, default set to FALSE
#' @param ncheckpoint the number of dose groups (count down from highest
#' concentration) to check whether their readings surpass fold change threshold
#' @param nreplicate number of replicates
#' @param ncondition number of treatment conditions
#' @param PSMcutoff whether to apply PSM threshold cutoff on hit selection,
#' default set to FALSE
#' @param PSMthreshold the threshold of averaged PSM numbers to consider protein
#' quantification as reliable, default value is 3
#' @param remfragment whether to remove fragment proteins, default set to FALSE
#' @param remribosomal whether to remove ribosomal proteins, default set to FALSE
#'
#' @import dplyr
#' @export
#' @return a list of segregated dataframes
#' @examples \dontrun{
#' ITTRdata_filtered <- ms_ITTR_filter(ITTRdata_fitted, keeptwice=TRUE,
#' PSMcutoff=TRUE)
#' }
#'
#'



ms_ITTR_filter <- function(data, nread=10, fcthreshold=0.3, R2threshold=0.8,
                           nbaseline=2, baselineMAD=0, nMAD=2.5,
                           checkreplicate=FALSE, treatcondition=NULL,
                           allpresence=FALSE, ncheckpoint=3, nreplicate=2, ncondition=2,
                           PSMcutoff=FALSE, PSMthreshold=3,
                           remfragment=FALSE, remribosomal=FALSE) {

  dataname <- deparse(substitute(data))
  # outdir <- data$outdir[1]
  # data$outdir <- NULL
  if (remfragment) {
    if(length(grep("Fragment", data$description))) {
      data <- data[-grep("Fragment", data$description), ]
    }
  }
  if (remribosomal) {
    if (length(grep("ribosomal", data$description))) {
      data_ribo <- data[grep("ribosomal", data$description), ]
      data <- data[-grep("ribosomal", data$description), ]
    } else{
      data_ribo <- data[0, ]
    }
  }

  if (PSMcutoff) {# To select the proteins with more than 3 PSM (average)
    PSMcol <- grep("PSM", names(data), value=FALSE)
    PSMname <- names(data)[PSMcol]
    names(data)[PSMcol] <- "PSM"
    PSMkeep <- data %>% group_by(id) %>%
      summarize(PSMmean=mean(PSM)) %>%
      filter(PSMmean>PSMthreshold)
    fkeep <- which(data$id %in% PSMkeep$id)
    names(data)[PSMcol] <- PSMname
    data_PSMsmall <- data[-fkeep, ]
    data <- data[fkeep, ]
  }
  #return(list(PSMsmall_protein=data_PSMsmall, PSMlarge_protein=data))

  pattern <- which( !is.na(data$R2) & data$R2>R2threshold )
  if (length(pattern)) {
    data_good_tem <- data[pattern, ]
    data_rem <- data[-pattern, ]
  } else {
    data_good_tem <- data
    data_rem <- data[0, ]
  }

  if (!length(baselineMAD)) {
    baselineMAD <- round(mad(unlist(data[ ,c(4:(3+nbaseline))]), na.rm=T), 4)
    print(paste0("The baseline variance (based on the first ", nbaseline,
                 " points) is ", baselineMAD, "."))
  }
  cutoff_high <- round((1+nMAD*baselineMAD)*(1+fcthreshold), 3)
  cutoff_low <- round((1-nMAD*baselineMAD)/(1+fcthreshold), 3)
  print(paste0("The upper cutoff threshold for shift set at ", cutoff_high, "."))
  print(paste0("The lower cutoff threshold for shift set at ", cutoff_low, "."))

  fkeep <- NULL
  nrowdata <- nrow(data_good_tem)
  for (i in 1:nrowdata) {
    highest <- max(as.numeric(data_good_tem[i,c((nread+4-ncheckpoint):
                                                  (nread+3))]), na.rm=T)
    lowest <- min(as.numeric(data_good_tem[i,c((nread+4-ncheckpoint):
                                                 (nread+3))]), na.rm=T)
    if (highest >= cutoff_high & data_good_tem[i,"Slope"] <0) {
      fkeep <- c(fkeep, i)
    }else if (lowest <= cutoff_low & data_good_tem[i,"Slope"] >0) {
      fkeep <- c(fkeep, i)
    }
  }

  if (length(fkeep)) {
    data_rem <- rbind(data_rem, data_good_tem[-fkeep, ])
    data_good_tem <- data_good_tem[fkeep, ]
  } else {
    stop("Opps, no proteins passed the selection threshold!")
  }

  if (checkreplicate) {
    # Focus only on the treatment or heat challenged conditions
    if (length(treatcondition)==1) {
      pattern <- grep(pattern=paste0("^", treatcondition, "\\."), data_good_tem$condition)
      data_good_tem <- data_good_tem[pattern, ]
    } else if(length(treatcondition)>1) {
      allpattern <- NULL
      for (i in 1:length(treatcondition)) {
        pattern <- grep(pattern=paste0("^", treatcondition[i], "\\."), data_good_tem$condition)
        allpattern <- c(allpattern, pattern)
      }
      data_good_tem <- data_good_tem[allpattern, ]
    }
    if (nrow(data_good_tem)==0) {
      stop("Opps, no proteins remained in your treament groups! Pls double check!")
    }

    data_good_copy <- data_good_tem
    data_good_copy$condition <- gsub("\\.[Rr][Ee][Pp][1-9]+", "", data_good_copy$condition)
    data_good_copy$condition <- gsub("\\.[Rr][1-9]+", "", data_good_copy$condition)
    data_good_copy$condition <- gsub("\\.[1-9]+", "", data_good_copy$condition)
    protein_rep <- which(with(data_good_copy,
                              tapply(condition,id,function(u) length(unique(u)))
                              < tapply(condition,id,function(u) length(u))))
    name_protein_rep <- names(protein_rep)
    name_protein_solo <- setdiff(data_good_copy$id, name_protein_rep)

    fkeep1 <- NULL
    fkeep2 <- NULL
    if (length(name_protein_rep)!=0) {
      fkeep1 <- which(data$id %in% name_protein_rep)
      data_good_p <- data[fkeep1, ]
      if (allpresence) { # whether the final hits need to be found in complete set
        protein_all <- which(with(data_good_p,
                                  tapply(condition,id,function(u)length(unique(u)))
                                  ==nreplicate*ncondition))
        nall <- names(protein_all)
        fkeep_all <- which(data_good_p$id %in% nall)
        data_good_p1 <- data_good_p[-fkeep_all, ]
        data_good_p <- data_good_p[fkeep_all, ]
      }
    } else {
      stop("Opps, no proteins with good replicative shifts were retrieved!")
    }

    if (length(name_protein_solo)!=0) {
      fkeep2 <- which(data$id %in% name_protein_solo)
      data_good_single_p <- data[fkeep2, ]
      if (allpresence) { # whether the final hits need to be found in complete set
        data_good_single_p <- rbind(data_good_single_p, data_good_p1)
      }
    } else {
      data_good_single_p <- data[0, ]
    }
    data_rem_p <- data[-c(fkeep1,fkeep2), ]
    #print(intersect(data_good_p$id, data_good_single_p$id))
    print(paste0("The number of proteins with replicative significant shifts in ", dataname, " :"))
    print(length(unique(data_good_p$id)))

    print(paste0("The number of proteins with solo significant shifts in ", dataname, " :"))
    print(length(unique(data_good_single_p$id)))
  } else {
    nkeep1 <- data_good_tem$id
    fkeep1 <- which(data$id %in% nkeep1)
    data_good_p <- data[fkeep1, ]
    data_rem_p <- data[-fkeep1, ]
  }

  if (checkreplicate & PSMcutoff) {
    return(list(good_protein=data_good_p, single_protein=data_good_single_p,
                PSMsmall_protein=data_PSMsmall, removed_protein=data_rem_p))
  } else if (checkreplicate) {
    return(list(good_protein=data_good_p, single_protein=data_good_single_p,
                removed_protein=data_rem_p))
  } else if (PSMcutoff) {
    return(list(good_protein=data_good_p, PSMsmall_protein=data_PSMsmall,
                removed_protein=data_rem_p))
  } else {
    return(list(good_protein=data_good_p, removed_protein=data_rem_p))
  }
}
