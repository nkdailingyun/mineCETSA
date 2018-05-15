#' ms_ITDR_filter
#'
#' Function to perform ITDR data segregation and positive hits selection based
#' on desired fc threshold level and the associated fitting parameters
#'
#' @param data dataset to be filtered
#' @param nread number of reading channels or sample treatements, default value
#' is 10
#' @param checkpointposition refering to the positions of dose points to check
#' whether their readings surpass fcthreshold, default value is NULL, which
#' would automatically check the highest 3 dose points
#' @param fcthreshold short for fold change threshold, indicate the threshold
#' fold change compared to baseline (i.e.,1.0) above which the readings can be
#' considered significantly stabilized or destabilized, default value is 0.3
#' @param R2threshold minimal R2 to consider dose-reponse fitting as reliable
#' @param baselineMAD MAD of baseline variation, default value is 0; if not
#' provided, it will be calculated based on the readings from the lowest few
#' dose points, specified by nbaseline value
#' @param nbaseline the number of points (count from lowest concentration) used
#' for calculate the baseline MAD, default is 3
#' @param nMAD the significance level of MAD cutoff, default value is 2.5
#' @param checkreplicate whether to check replicated measurements, when set to
#' TRUE, further segregate the proteins that all the replicated measurements from
#' at least one experimental condition surpass fold change threshold, default
#' set to TRUE
#' @param treatcondition specify the condition names to perform replicate check,
#' no need to specify the replicate part of the full condition names, useful when
#' multiple experimental conditions are contained in the dataset
#' @param PSMcutoff whether to apply PSM threshold cutoff on hit selection,
#' default set to TRUE
#' @param PSMthreshold the threshold of averaged PSM numbers to consider protein
#' quantification as reliable, default value is 3
#'
#' @import dplyr
#' @export
#' @return a list of segregated dataframes
#' @examples \dontrun{
#' ITDRdata_filtered <- ms_ITDR_filter(ITDRdata_fitted, fcthreshold=0.3,
#'                      checkreplicate=TRUE, PSMcutoff=TRUE)
#' }
#'
#'

ms_ITDR_filter <- function(data, nread=10, checkpointposition=NULL,
                           fcthreshold=0.3, R2threshold=0.8,
                           baselineMAD=0, nbaseline=3, nMAD=2.5,
                           checkreplicate=TRUE, treatcondition=NULL,
                           PSMcutoff=TRUE, PSMthreshold=3) {

  dataname <- deparse(substitute(data))
  outdir <- ms_directory(data, dataname)

  if (PSMcutoff) { # To select the proteins with more than 3 PSM (average)
    PSMcol <- grep("PSM", names(data), value=FALSE)
    PSMname <- names(data)[PSMcol]
    names(data)[PSMcol] <- "PSM"
    PSMkeep <- data %>% group_by(id) %>%
      summarize(PSMmean=mean(PSM)) %>%
      filter(PSMmean>=PSMthreshold)
    fkeep <- which(data$id %in% PSMkeep$id)
    names(data)[PSMcol] <- PSMname
    data_PSMsmall <- data[-fkeep, ]
    data <- data[fkeep, ]
    if (length(attr(data_PSMsmall,"outdir"))==0 & length(outdir)>0) {
      attr(data_PSMsmall,"outdir") <- outdir
    }
  }
  #return(list(PSMsmall_protein=data_PSMsmall, PSMlarge_protein=data))

  #to keep the data proper replicates as specified
  #to separate condition into condition and replicates
  nset <- length(unique(data$condition))
  data1 <- tidyr::separate(data, condition, into=c("condition", "replicate"), sep="\\.")
  ncondition <- length(unique(data1$condition))
  nreplicate <- length(unique(data1$replicate))
  uniquecond <- unique(data1[ ,c("condition", "replicate")])
  row.names(uniquecond) <- NULL
  # print("Replicates information were extracted as follows:")
  # print(as.data.frame(uniquecond))

  pattern <- which( !is.na(data$R2) & data$R2>=R2threshold )
  if (length(pattern)) {
    data_good_tem <- data[pattern, ]
    data_rem <- data[-pattern, ]
  } else {
    data_good_tem <- data
    data_rem <- data[0, ]
  }

  if (length(baselineMAD)==0) {
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

  if (length(checkpointposition)==0) { # default to check the last three points
    checkpointposition <- c(8,9,10)-(10-nread)+3
  } else {
    checkpointposition <- checkpointposition+3
  }
  for (i in 1:nrowdata) {
    # if (excludeNA) {
    #' whether to exclude the NA points when counting points to check,
    #' default set to FALSE, when set to TRUE, the exact number of points with readings
    #' would be used as specified by ncheckpoint parameter
    #   data_to_check <- na.omit(as.numeric(data_good_tem[i,c(4:(nread+3))]))
    #   npoint <- length(data_to_check)
    #   highest <- max(data_to_check[c((npoint-ncheckpoint+1):npoint)], na.rm=T)
    #   lowest <- min(data_to_check[c((npoint-ncheckpoint+1):npoint)], na.rm=T)
    # } else if (excludelastpoint) {
    #   highest <- max(as.numeric(data_good_tem[i,c((nread+3-ncheckpoint):
    #                                                 (nread+2))]), na.rm=T)
    #   lowest <- min(as.numeric(data_good_tem[i,c((nread+3-ncheckpoint):
    #                                                (nread+2))]), na.rm=T)
    # } else {
    #   highest <- max(as.numeric(data_good_tem[i,c((nread+4-ncheckpoint):
    #                                                 (nread+3))]), na.rm=T)
    #   lowest <- min(as.numeric(data_good_tem[i,c((nread+4-ncheckpoint):
    #                                                (nread+3))]), na.rm=T)
    # }
    highest <- max(as.numeric(data_good_tem[i,checkpointposition]), na.rm=T)
    lowest <- min(as.numeric(data_good_tem[i,checkpointposition]), na.rm=T)

    if (highest >= cutoff_high & data_good_tem[i,"Slope"] <0) {
      fkeep <- c(fkeep, i)
    } else if (lowest <= cutoff_low & data_good_tem[i,"Slope"] >0) {
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
    if (length(treatcondition)>=1) {
      uniquecond1 <- subset(uniquecond, condition %in% treatcondition)
    } else {
      uniquecond1 <- uniquecond
    }
    conditionrep <- dplyr::count(uniquecond1, condition)
    data_good_tem <- tidyr::separate(data_good_tem, condition, into=c("condition", "replicate"), sep="\\.")
    data_good_tem$replicate <- NULL
    data_good_tem_freq <- dplyr::count(data_good_tem, id, condition)
    data_good_keep <- data_good_tem_freq[0, ]
    for (i in 1:nrow(conditionrep)) {
      data_good_keep <- rbind(data_good_keep, subset(data_good_tem_freq, condition==conditionrep$condition[i] & n==conditionrep$n[i]))
    }
    if (nrow(data_good_keep)==0) {
      stop("Opps, no proteins remained in your treament groups! Pls double check!")
    }
    name_protein_rep <- data_good_keep$id
    name_protein_solo <- setdiff(data_good_tem$id, name_protein_rep)

    fkeep1 <- NULL
    fkeep2 <- NULL
    if (length(name_protein_rep)!=0) {
      fkeep1 <- which(data$id %in% name_protein_rep)
      data_good_p <- data[fkeep1, ]
    } else {
      stop("Opps, no proteins with good replicative shifts were retrieved!")
    }

    if (length(name_protein_solo)!=0) {
      fkeep2 <- which(data$id %in% name_protein_solo)
      data_good_single_p <- data[fkeep2, ]
    } else {
      data_good_single_p <- data[0, ]
    }
    data_rem_p <- data[-c(fkeep1,fkeep2), ]
    if (length(attr(data_good_p,"outdir"))==0 & length(outdir)>0) {
      attr(data_good_p,"outdir") <- outdir
    }
    if (length(attr(data_good_single_p,"outdir"))==0 & length(outdir)>0) {
      attr(data_good_single_p,"outdir") <- outdir
    }
    if (length(attr(data_rem_p,"outdir"))==0 & length(outdir)>0) {
      attr(data_rem_p,"outdir") <- outdir
    }
    print(paste0("The number of proteins with replicative significant shifts in ", dataname, " :"))
    print(length(unique(data_good_p$id)))

    print(paste0("The number of proteins with solo significant shifts in ", dataname, " :"))
    print(length(unique(data_good_single_p$id)))
  } else {
    nkeep1 <- data_good_tem$id
    fkeep1 <- which(data$id %in% nkeep1)
    data_good_p <- data[fkeep1, ]
    data_rem_p <- data[-fkeep1, ]
    if (length(attr(data_good_p,"outdir"))==0 & length(outdir)>0) {
      attr(data_good_p,"outdir") <- outdir
    }
    if (length(attr(data_rem_p,"outdir"))==0 & length(outdir)>0) {
      attr(data_rem_p,"outdir") <- outdir
    }
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
