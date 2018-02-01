#' ms_goodcurve_selection
#'
#' Function to perform good quality melting curve selection on dataframe
#' The main criteria including averaged R2>0.8, with steady plateau
#' at the top and bottom part of the melting curves
#'
#' @param data dataset to look for complete set
#' @param nread number of reading channels or sample treatements, default value
#' is 10
#' @param ctrlcond by default the function only apply on the data containing condition
#' keyword "Ctrl"/"DMSO", the user could add customized condition keyword using this argument
#' @param topCVcutoff the threshold CV value used to define the top plateau cutoff, default value is 0.1
#' @param bottomCVcutoff the threshold CV value used to define the bottom plateau cutoff, default value is 0.1
#' @param bottomcutoff the threshold value to control the maximal bottom plateau, default value is 0.4
#' @param nMAD the significance level of MAD cutoff, default value is 2.5

#' @importFrom plyr ddply
#' @import tidyr
#' @export
#' @return dataframe containing the subset of melting curves with good melting profile
#' @examples \dontrun{
#' data_good <- ms_goodcurve_selection(data_scaled, bottomcutoff=0.3)
#' }
#'
#'
ms_goodcurve_selection <- function(data, nread=10, ctrlcond=NULL,
                                   topCVcutoff=0.1,bottomCVcutoff=0.1,
                                   bottomcutoff=0.4, nMAD=2.5) {

  # look for outlier proteins based on melting behavior in controls
  outliers <- NULL
  print("Make sure you provide fitted data with Tm and R2 values for this option!")
  ctrllist1 <- unique(grep("[Cc][Tt][Rr][Ll]", data$condition, value=TRUE))
  ctrllist2 <- unique(grep("[Cc][Oo][Nn][Tt][Rr][Oo][Ll]", data$condition, value=TRUE))
  ctrllist3 <- unique(grep("[Dd][Mm][Ss][Oo]", data$condition, value=TRUE))
  ctrllist <- c(ctrllist1, ctrllist2, ctrllist3, ctrlcond)
  if(length(ctrllist)==0) {
    stop("Name your control conditions with prefix [Ctrl] or [Control] or [DMSO], or specify in ctrlcond argument")
  }
  data$Top <- rowMeans(data[, c(4:6)],na.rm =T)
  data$Bottom <- rowMeans(data[ ,c((nread+1):(nread+3))],na.rm =T)
  data$TopCV <- apply(data[, c(4:6)],1,sd,na.rm =T)/rowMeans(data[, c(4:6)],na.rm =T)
  data$BottomCV <- apply(data[ ,c((nread+1):(nread+3))],1,sd,na.rm =T)/rowMeans(data[ ,c((nread+1):(nread+3))])
  if (length(topCVcutoff)==0) {
    topCVcutoff <- median(data$TopCV,na.rm =T)+nMAD*mad(data$TopCV,na.rm =T)
    print("The top plateau CV cutoff value based on mad scheme is ", topCVcutoff)
  }
  if (length(bottomCVcutoff)==0) {
    bottomCVcutoff <- median(data$BottomCV,na.rm =T)+nMAD*mad(data$BottomCV,na.rm =T)
    print("The bottom plateau CV cutoff value based on mad scheme is ", bottomCVcutoff)
  }

  tmp <- which(is.na(data$Tm))
  tmp <- c(tmp,which(data$Slope<=0))
  tmp <- c(tmp,which(data$R2<=0))
  R2table <- plyr::ddply(data, .(id), summarize, meanR2=median(R2, na.rm=TRUE))
  R2table <- subset(R2table, meanR2<0.8)
  tmp <- c(tmp,which(data$id %in% R2table$id))
  #tmp4 <- which(data$condition %in% ctrllist & data$Tm>=100);
  tmp <- c(tmp,which(data$condition %in% ctrllist & data$Top<data$Bottom))
  #tmp6 <- which(data$condition %in% ctrllist & data$Top<topcutoff)
  tmp <- c(tmp,which(data$condition %in% ctrllist & data$Bottom>bottomcutoff))
  tmp <- c(tmp,which(data$condition %in% ctrllist & data$TopCV>topCVcutoff))
  tmp <- c(tmp,which(data$condition %in% ctrllist & data$BottomCV>bottomCVcutoff))
  tmp <- unique(tmp)
  data$Top <- NULL
  data$Bottom <- NULL
  data$TopCV <- NULL
  data$BottomCV <- NULL
  if (length(tmp)>0) {
    print(paste0(length(tmp), " measurements were messy in melting behavior and removed."))
    # outlierid <- data[tmp, ]$id
    # outlierid1 <- which(data$id %in% outlierid)
    # outliers <- data[outlierid1, ]
    # ms_filewrite(outliers, "Messy proteins.txt", outdir=outdir)
    #print(paste0("The outlier protein list is ", outlierids))
    data <- data[-tmp, ]
  }
  return(data)
}

#' ms_find_replicate
#'
#' Function to perform complete replicated dataset subsetting on dataframe
#'
#' @param data dataset to look for complete set
#' @param nread number of reading channels or sample treatements, default value
#' is 10
#' @param getcompletesubsetonly whether to just perform a simple subsetting of dataset
#' with specified replicate numbers
#' @param variancecutoff whether to segregate the proteins with large
#' inter-replicate variance, default set to TRUE, only used when there is single condition
#' @param nMAD_var the significance level of MAD cutoff, default value is 2.5
#'
#' @importFrom plyr daply
#' @import tidyr
#' @export
#' @return dataframe containing the full replicated complete dataset
#' @examples \dontrun{
#' data_complete_set <- ms_find_replicate(data_scaled)
#' }
#'
#'

ms_find_replicate <- function(data, nread=10, getcompletesubsetonly=TRUE,
                              variancecutoff=FALSE, nMAD_var=2.5) {

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
  colorders <- names(data)

  # subsetting complete replicates:
  # to keep the data in at least one condition with full replicates
  data1 <- tidyr::separate(data, condition, into=c("condition", "replicate"), sep="\\.")
  ncondition <- length(unique(data1$condition))
  nreplicate <- length(unique(data1$replicate))
  uniquecond <- unique(data1[ ,c("condition", "replicate")])
  #print(as.data.frame(uniquecond))
  conditionrep <- dplyr::count(uniquecond, condition)
  #print(as.data.frame(conditionrep))
  data_freq <- dplyr::count(data1, id, condition)
  #print(head(data_freq, n=20))
  id_keep <- data_freq[0, ]
  for (i in 1: nrow(conditionrep)) {
    # print(conditionrep$condition[i])
    # print(conditionrep$n[i])
    id_keep <- rbind(id_keep, subset(data_freq, condition==conditionrep$condition[i] & n==conditionrep$n[i]))
  }
  data_complete <- merge(data1, id_keep, by=c("id","condition"), all = FALSE)
  data_complete <- tidyr::unite(data_complete, condition, condition, replicate, sep=".")
  data_complete$n <- NULL
  print(paste(nrow(data_complete), "measurements were measured with complete replicates in at least one condition.", sep=" "))

  # keep1 <- which(table(data$id) == nreplicate)
  # if(length(keep1)==0){ stop("No proteins contain complete replicate") }
  # nkeep1 <- names(keep1)
  # fkeep1 <- which(data$id %in% nkeep1)
  # data_complete <- data[fkeep1, ]
  # print(paste0("The number of proteins with complete replicates is: ",
  #              length(unique(data_complete$id))))
  # print(paste0("The percentage of proteins with complete replicates is: ",
  #              round(length(unique(data_complete$id))/length(unique(data$id)),3)*100, "%"))

  if (getcompletesubsetonly) {
    data_complete <- data_complete[ ,colorders]
    data_complete$outdir <- outdir
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
#' @param rep the rep indicator used in the naming of experimental conditions, such as "r", or "rep" or ""
#' @param weightbycountnum whether weight the readings with count information, default to TRUE
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
ms_data_average <- function(data, nread=10, rep="r", weightbycountnum=TRUE) {
  # The input data is post fitting and filterring

  dataname <- deparse(substitute(data))
  outdir <- data$outdir[1]
  data$outdir <- NULL

  # to prevent the change of sub-directory folder
  if (!length(outdir)) {
    print("no output directory information, creating one..")
    outdir <- paste0(dataname,"_",format(Sys.time(), "%y%m%d_%H%M"))
    dir.create(outdir)
  } else if (dir.exists(outdir)==FALSE) {
    dir.create(outdir)
  }

  nrowdata <- nrow(data)
  colorders <- names(data)
  if (length(grep("description", names(data)))) {
    proteininfo <- unique(data[ ,c("id","description")])
  }
  data$description <- NULL

  d1 <- tidyr::gather(data[ ,c(1:(nread+2))], treatment, reading, -id, -condition)
  ##summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
  # datarr <- summarySE(datal, measurevar="reading", groupvars=c("id", "condition", "temperature"))
  if (weightbycountnum) {
    d1 <- merge(d1, data[ ,c("id","condition","countNum")])
    d1$condition <- gsub(paste0("\\.", rep, "[0-9]+"), "", d1$condition)
    d1$condition <- gsub("\\.[0-9]", "", d1$condition)

    datac <- plyr::ddply(d1, plyr::.(id,condition,treatment), .drop=TRUE,
                         .fun = function(xx) {
                           c(reading_mean = weighted.mean(xx[["reading"]], xx[["countNum"]], na.rm=TRUE)
                           )
                         }
    )
  } else {
    d1$condition <- gsub(paste0("\\.", rep, "[0-9]+"), "", d1$condition)
    d1$condition <- gsub("\\.[0-9]", "", d1$condition)

    datac <- plyr::ddply(d1, plyr::.(id,condition,treatment), .drop=TRUE,
                         .fun = function(xx) {
                           c(reading_mean = mean(xx[["reading"]], na.rm=TRUE)
                             #sd   = sd(xx[[col]], na.rm=na.rm)
                           )
                         }
    )
  }
  averageddata <- tidyr::spread(datac, treatment, reading_mean)
  # print(head(averageddata))

  d1 <- tidyr::gather(data[ ,-c(3:(nread+2))], parameter, count, -id, -condition)
  d1$condition <- gsub(paste0("\\.", rep, "[0-9]+"), "", d1$condition)
  d1$condition <- gsub("\\.[0-9]", "", d1$condition)
  # print(head(d1))

  datacount <- plyr::ddply(d1, plyr::.(id,condition,parameter), .drop=TRUE,
                           .fun = function(xx, col) {
                             c(count_mean = mean(xx[[col]], na.rm=TRUE)
                               #sd  = sd(xx[[col]], na.rm=na.rm)
                             )
                           },
                           "count"
  )
  # print(head(datacount))
  datacount <- tidyr::spread(datacount, parameter, count_mean)
  # print(head(datacount))
  averageddata <- merge(averageddata, datacount)
  averageddata <- merge(averageddata, proteininfo)
  averageddata <- averageddata[ ,colorders]

  # if (calsd) {
  #   datac1 <- d1 %>% group_by(id,condition,treatment) %>% summarise(reading_sd=sd(reading, na.rm=TRUE))
  #   datac <- merge(datac, datac1)
  #   return(datac)
  # }

  averageddata$outdir <- outdir
  return(averageddata)
}
