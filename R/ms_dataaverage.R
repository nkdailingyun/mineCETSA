#' ms_find_replicate
#'
#' Function to perform complete replicated dataset subsetting on dataframe
#'
#' @param data dataset to look for complete set
#' @param nread number of reading channels or sample treatements, default value
#' is 10
#' @param nreplicate number of replicates, default value is 2, change the value
#' accordingly
#' @param subsetonly whether to just perform a simple subsetting of dataset
#' with specified replicate numbers
#' @param variancecutoff whether to segregate the proteins with large
#' inter-replicate variance
#' @param nMAD_var the significance level of MAD cutoff, default value is 2.5
#'
#' @importFrom plyr daply
#' @import tidyr
#' @export
#' @return dataframe containing the full replicated complete dataset
#' @examples \dontrun{
#' data_complete_set <- ms_finc_replicate(data_scaled, nreplicate=2)
#' }
#'
#'

ms_find_replicate <- function(data, nreplicate=2, nread=10, subsetonly=FALSE,
                              variancecutoff=TRUE, nMAD_var=2.5) {

  dataname <- deparse(substitute(data))
  outdir <- data$outdir[1]
  data$outdir <- NULL

  # to prevent the change of sub-directory folder
  if (!length(outdir)) {
    outdir <- paste0(dataname,"_",format(Sys.time(), "%y%m%d_%H%M"))
    dir.create(outdir)
  } else if (dir.exists(outdir)==FALSE) {
    dir.create(outdir)
  }

  #if (nreplicate>4) {stop("Only up to four replicates are supported for now")}

  # if (length(grep("description", names(data)))) {
  #   proteininfo <- unique(data[ ,c("id","description")])
  # }
  # data$description <- NULL

  nrowdata <- nrow(data)
  #distance_matrix <- matrix(nrow = nrowdata, ncol = nread)

  # subsetting complete replicates:
  keep1 <- which(table(data$id) == nreplicate)
  if(length(keep1)==0){ stop("No proteins contain complete replicate") }
  nkeep1 <- names(keep1)
  fkeep1 <- which(data$id %in% nkeep1)
  data_complete <- data[fkeep1, ]
  print(paste0("The number of proteins with complete replicates is: ",
               length(unique(data_complete$id))))
  print(paste0("The percentage of proteins with complete replicates is: ",
               round(length(unique(data_complete$id))/length(unique(data$id)),3)*100, "%"))

  if (subsetonly) {
    return(data_complete)
  }

  if (nreplicate==2) {
    dism <- plyr::daply(data_complete, "id", function(data) {
      data<-data[order(data$condition), ]
      dm<-as.vector(dist(data[ ,c(4:(nread+3))]))
      eucl<-dm
      eucl
    })

    dism <- data.frame(id=names(dism), distance=dism, row.names=NULL)
    #return(dism)
    print("The distances of inter-replicates readings distribution is as follows:")
    print(summary(dism$distance))
  } else {
    dism <- plyr::daply(data_complete, "id", function(data) {
      data<-data[order(data$condition), ]
      dm<-as.vector(dist(data[ ,c(4:(nread+3))]))
      nsd<-sd(dm)
      nmean<-mean(dm,na.rm=T)
      ncv<-nsd/nmean
      eucl<-c(dm, nmean, nsd, ncv)
      eucl
    })

    dism <- data.frame(id=rownames(dism), dism, row.names=NULL)
    nlength <- nreplicate*(nreplicate-1)/2
    colnames(dism) <- c("id", paste0("distance_", 1:nlength), "mean", "sd", "cv")
    #return(dism)
    print("The CVs of inter-replicates readings distribution is as follows:")
    print(summary(dism$cv))
  }

  data_largevar <- NULL
  if (variancecutoff) {
    cutoff <- NULL
    # MAD guided significance
    print("use MAD scheme to assign variance cutoff...")
    if (nreplicate==2) {
      cutoff <- round(median(dism$distance,na.rm=T)+nMAD_var*mad(dism$distance,na.rm=T), 3)
      print(paste0("The distance cutoff limit (", nMAD_var, "*mad) sets at ", cutoff))
      dism1 <- subset(dism, distance<cutoff)
    } else {
      cutoff <- round(median(dism$cv,na.rm=T)+nMAD_var*mad(dism$cv,na.rm=T), 3)
      print(paste0("The variance (CV) cutoff limit (", nMAD_var, "*mad) sets at ", cutoff))
      dism1 <- subset(dism, cv<cutoff)
    }
    fkeep1 <- which(data_complete$id %in% dism1$id)
    data_largevar <- data_complete[-fkeep1, ]
    data_complete <- data_complete[fkeep1, ]
    print(paste0("The number of reproducible proteins with complete replicates is: ",
                 length(unique(data_complete$id))))
    data_complete$outdir <- outdir
    data_largevar$outdir <- outdir
    return(list(data_complete=data_complete, data_largevar=data_largevar, dism=dism))
  } else {
    data_complete$outdir <- outdir
    return(list(data_complete=data_complete, dism=dism))
  }

}


#' ms_data_average
#'
#' Function to perform reading averaging on dataframe
#'
#' @param data dataset to be filtered
#' @param nread number of reading channels or sample treatements, default value
#' is 10
#' @param rep the rep indicator on the naming of experimental conditions
#' @param calPSM whether to calculate arithmatic mean of PSMs, Unique peptides
#' and count numbers
#' @param calsd whether to calculate standard deviation on readings, default FALSE
#'
#' @importFrom plyr . ddply
#' @import dplyr
#' @import tidyr
#' @export
#' @return dataframe containing the averaged readings
#' @examples \dontrun{
#' data_averaged <- ms_data_average(data_complete_set, rep="")
#' }
#'
#'
ms_data_average <- function(data, nread=10, rep="r", weightbycountnum=TRUE, calPSM=FALSE, calsd=FALSE) {
  # The input data is post fitting and filterring

  # calratio=FALSE, numerator="52C", denominator="37C"

  dataname <- deparse(substitute(data))
  outdir <- data$outdir[1]
  data$outdir <- NULL

  # to prevent the change of sub-directory folder
  if (!length(outdir)) {
    outdir <- paste0(dataname,"_",format(Sys.time(), "%y%m%d_%H%M"))
    dir.create(outdir)
  } else if (dir.exists(outdir)==FALSE) {
    dir.create(outdir)
  }

  nrowdata <- nrow(data)

  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
  }
  data$description <- NULL

  d1 <- tidyr::gather(data[ ,c(1:(nread+2))], treatment, reading, -id, -condition)

  # summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
  #datarr <- summarySE(datal, measurevar="reading", groupvars=c("id", "condition", "temperature"))
  if (weightbycountnum) {
    d1 <- merge(d1, data[ ,c("id","condition","sumUniPeps","sumPSMs","countNum")])
    d1$condition <- gsub(paste0("\\.", rep, "[1-9]+"), "", d1$condition)
    d1$condition <- gsub("\\.[1-9]", "", d1$condition)

    datac <- plyr::ddply(d1, plyr::.(id,condition,treatment), .drop=TRUE,
                         .fun = function(xx) {
                           c(reading_mean = weighted.mean(xx[["reading"]], xx[["countNum"]], na.rm=TRUE),
                             sumUniPeps_new = mean(xx[["sumUniPeps"]]),
                             sumPSMs_new = mean(xx[["sumPSMs"]]),
                             countNum_new = mean(xx[["countNum"]])
                           )
                         }
    )
  } else {
    d1$condition <- gsub(paste0("\\.", rep, "[1-9]+"), "", d1$condition)
    d1$condition <- gsub("\\.[1-9]", "", d1$condition)

    datac <- plyr::ddply(d1, plyr::.(id,condition,treatment), .drop=TRUE,
                         .fun = function(xx) {
                           c(reading_mean = mean(xx[["reading"]], na.rm=TRUE)
                             #sd   = sd(xx[[col]], na.rm=na.rm)
                           )
                         }
    )
  }

  #datac <- d1 %>% group_by(id,condition,treatment) %>% summarise(reading_mean=mean(reading, na.rm=TRUE))

  if (calPSM) {
    data_copy <- data
    data_copy$condition <- gsub(paste0("\\.", rep, "[1-9]+"), "", data_copy$condition)
    data_copy$condition <- gsub("\\.[1-9]", "", data_copy$condition)
    if (grep("sumUniPeps", names(data_copy))) {
      pepinfo <- data_copy[ ,c("id","condition","sumUniPeps")] %>% group_by(id,condition) %>% summarise(sumUniPeps_mean=mean(sumUniPeps, na.rm=TRUE))
    } else {
      stop("sumUniPeps column is missing")
    }
    if (grep("sumPSMs", names(data_copy))) {
      PSMinfo <- data_copy[ ,c("id","condition","sumPSMs")] %>% group_by(id,condition) %>% summarise(sumPSMs_mean=mean(sumPSMs, na.rm=TRUE))
    } else {
      stop("sumPSMs column is missing")
    }
    if (grep("countNum", names(data_copy))) {
      countinfo <- data_copy[ ,c("id","condition","countNum")] %>% group_by(id,condition) %>% summarise(countNum_mean=mean(countNum, na.rm=TRUE))
    } else {
      stop("countNums column is missing")
    }
    pepinfo <- merge(pepinfo, PSMinfo, by=c("id","condition"))
    pepinfo <- merge(pepinfo, countinfo, by=c("id","condition"))
  }

  if (calsd) {
    datac1 <- d1 %>% group_by(id,condition,treatment) %>% summarise(reading_sd=sd(reading, na.rm=TRUE))
    datac <- merge(datac, datac1)
    return(datac)
  }

  # if (calratio) {
  #   datac <- mutate(datac, role=ifelse(condition==numerator,"numerator","denominator"))
  #   ratios <- ddply(datac, .(id,concentration), summarize,
  #                   Ratio52to37 = mean[role == "numerator"] / mean[role == "denominator"]);
  #   ratios <- dcast(ratios, id~concentration, value.var="Ratio52to37", drop=FALSE);
  #
  #   list <- strsplit(ratios$id, "\n");
  #   df <- ldply(list);
  #   colnames(df) <- c("id", "Description");
  #   ratios<-cbind(df, ratios[,-1]);
  # }

  if (weightbycountnum) {
    averageddata <- tidyr::spread(datac[ ,c(1:4)], treatment, reading_mean)
    averageddata <- merge(averageddata, unique(datac[ ,c(1,2,5:7)]), all=FALSE)
  } else {
    averageddata <- tidyr::spread(datac, treatment, reading_mean)
  }

  if (calPSM) {
    averageddata <- merge(averageddata, pepinfo, by=c("id","condition"))
  }

  averageddata <- merge(proteininfo, averageddata)
  names(averageddata) <- gsub("_new","", names(averageddata))
  averageddata$outdir <- outdir

  # if(calratio){
  #   return(list(averageddata=averageddata,ratios=ratios))
  # }else{
  #   return(averageddata)
  # }
  return(averageddata)
}
