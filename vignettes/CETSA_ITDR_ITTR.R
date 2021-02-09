## ---- message=TRUE------------------------------------------------------------
library("mineCETSA")

## ---- eval=FALSE--------------------------------------------------------------
#  ITDRdata <- ms_ITDR_rawread(c("SN1.txt", "SN2.txt", "ST1.txt", "ST2.txt"), PD21=TRUE, dose=c(20,5,1.25,0.3123,0.078125,0.019531,0.004883,0.001221,0.000305,0))
#  
#  ITDRdata_1file <- ms_ITDR_rawread("SN1.txt", PD21=TRUE, dose=c(20,5,1.25,0.3123,0.078125,0.019531,0.004883,0.001221,0.000305,0))
#  
#  ITTRdata <- ms_ITTR_rawread(c("file1.txt", "file2.txt", "file3.txt", "file4.txt"), PD21=TRUE, time=c(0,5,10,15,20,30,45,60,90,120))

## ---- eval=FALSE--------------------------------------------------------------
#  ITDRdata <- ms_conditionrename(ITDRdata, incondition=c("drugS.3","drugS.4","drugX.1","drugX.2"), outcondition=c("drugS.r1","drugS.r2","drugX.r1","drugX.r2"))

## ---- eval=FALSE--------------------------------------------------------------
#  ITDRdata_c <- ms_clean(ITDRdata)
#  ITTRdata_c <- ms_clean(ITTRdata)

## ---- eval=FALSE--------------------------------------------------------------
#  ITDRdata_c1 <- ms_isoform_resolve(ITDRdata_c)
#  # The next isoform cleaning step is optional. You are encouraged to double check and adjust the match table entries.
#  ITDRdata_c2 <- ms_isoform_consolidate(ITDRdata_c1, matchtable="./subfolder/tobe_consolidated.txt")
#  # When the data contains protein abundance information.
#  ITDRdata_c2 <- ms_isoform_consolidate(ITDRdata_c1, matchtable="./subfolder/tobe_consolidated.txt", withabd=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  ITDRdata_s <- ms_ITDR_scaling(ITDRdata_c2)
#  ITTRdata_s <- ms_ITTR_scaling(ITTRdata_c2)

## ---- eval=FALSE--------------------------------------------------------------
#  ITDRdata_f <- ms_ITDR_fitting(ITDRdata_s)
#  ITTRdata_f <- ms_ITTR_fitting(ITTRdata_s, calMTT=FALSE, fc=0.3)

## ---- eval=FALSE--------------------------------------------------------------
#  ms_IDIT_QC(ITDRdata_f)
#  ms_IDIT_QC(ITDRdata_s, isdatafitted=FALSE, foldername="QC1") #not recommended

## ----eval=FALSE---------------------------------------------------------------
#  ITDRdata_shifted <- ms_IDIT_R2AUC(ITDRdata_f, printBothName=TRUE)
#  ITDRdata_stabilized <- ms_IDIT_R2AUC(ITDRdata_f, printBothName=FALSE, printGeneName=TRUE, keepreplicate=TRUE, onlyshowstabilized=TRUE)
#  ITDRdata_hits <- subset(ITDRdata_f, id %in% ITDRdata_shifted[[1]]$id)
#  ITDRdata_stabilized_hits <- subset(ITDRdata_f, id %in% ITDRdata_stabilized[[1]]$id)

## ---- eval=FALSE--------------------------------------------------------------
#  ITDRdata_filtered <- ms_ITDR_filter(ITDRdata_f, fcthreshold=0.3, checkreplicate=TRUE, PSMcutoff=TRUE)
#  ITTRdata_filtered <- ms_ITTR_filter(ITTRdata_f, fcthreshold=0.3, checkreplicate=TRUE, PSMcutoff=TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  ms_ITDR_ggplotting(ITDRdata_hits, orderAUC=TRUE)
#  ms_ITDR_ggplotting(ITDRdata_stabilized_hits, orderAUC=TRUE)
#  ms_ITTR_ggplotting(ITTRdata_filtered[[1]], orderAUC=TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  ms_filewrite(ITDRdata_hits, "ITDRdata_hits.txt")
#  ms_filewrite(ITDRdata_stabilized_hits, "ITDRdata_stabilized_hits.txt")
#  ms_filewrite(ITDRdata_filtered[[1]], "ITDRdata_good.txt")
#  ms_filewrite(ITDRdata_filtered[[2]], "ITDRdata_single.txt")
#  ms_filewrite(ITDRdata_filtered[[3]], "ITDRdata_PSMsmall.txt")
#  ms_filewrite(ITDRdata_filtered[[4]], "ITDRdata_notselected.txt")

## ---- eval=FALSE--------------------------------------------------------------
#  ITDRdata_selected_good <- ms_fileread("./subfolder_name/ITDRdata_good.txt")
#  # . indicate the current working directory
#  # Here the subfolder_name refer to a first level folder name directly under current working directory
#  # Remember that typically now you are still working in the main working directory
#  # The subfolder should contain the data file to read in

## ---- eval=FALSE--------------------------------------------------------------
#  ITDRdata_subset <- ms_subsetting(ITDRdata_f, hitidlist="hit_list.txt", isfile=TRUE)
#  ITDRdata_subset <- ms_subsetting(ITDRdata_f, hitidlist=c("P00000", "P12345-6"), isfile=FALSE)

