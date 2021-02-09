## ---- message=TRUE------------------------------------------------------------
library("mineCETSA")

## ---- eval=FALSE--------------------------------------------------------------
#  LY <- ms_rawread(c("LY_K562_Ctrl_rep1_PD21_Proteins.txt","LY_K562_Ctrl_rep2_PD21_Proteins.txt","LY_K562_Treatment_rep1_PD21_Proteins.txt","LY_K562_Treatment_rep2_PD21_Proteins.txt"))

## ---- eval=FALSE--------------------------------------------------------------
#  LY1 <- ms_conditionrename(LY, incondition=c("Treatmentrep1.3", "Treatmentrep2.4"), outcondition=c("Treatment.1","Treatment.2"))

## ---- eval=FALSE--------------------------------------------------------------
#  LY_cleaned <- ms_clean(LY1)

## ---- eval=FALSE--------------------------------------------------------------
#  LY_cleaned1 <- ms_isoform_resolve(LY_cleaned)
#  # The next isoform cleaning step is optional. You are encouraged to double check and adjust the match table entries.
#  LY_cleaned2 <- ms_isoform_consolidate(LY_cleaned1, matchtable="./subfolder/tobe_consolidated.txt")

## ---- eval=FALSE--------------------------------------------------------------
#  LY_scaled <- ms_scaling(LY_cleaned2)

## ---- eval=FALSE--------------------------------------------------------------
#  LY_fs_all <- ms_forcescaling(LY_cleaned, refcondition="all")
#  LY_fs <- ms_forcescaling(LY_cleaned, refcondition=c("Ctrl.1", "Ctrl.2"))

## ---- eval=FALSE--------------------------------------------------------------
#  LY_fitted <- ms_fitting(LY_scaled)

## ---- eval=FALSE--------------------------------------------------------------
#  ms_ggplotting(LY_scaled)

## ---- eval=FALSE--------------------------------------------------------------
#  ms_ggplotting_rep(LY_scaled, levelvector=c("Ctrl.1", "Ctrl.2", "Treatment.1","Treatment.2"))
#  ms_ggplotting_rep(LY_fitted, levelvector=c("Ctrl.1", "Ctrl.2", "Treatment.1","Treatment.2"))

## ---- eval=FALSE--------------------------------------------------------------
#  ms_filewrite(LY_scaled, "CETSA_data_scaled.txt")
#  ms_filewrite(LY_fitted, "CETSA_data_with_fitting_parameters.txt")

## ---- eval=FALSE--------------------------------------------------------------
#  MyCETSAdata <- ms_fileread("./subfolder_name/myCETSAdata_to_readin.txt")
#  # . indicate the current working directory
#  # Here the subfolder_name refer to a first level folder name directly under current working directory
#  # Remember that typically now you are still working in the main working directory
#  # The subfolder should contain the data file to read in

## ---- eval=FALSE--------------------------------------------------------------
#  CETSAdata_subset <- ms_subsetting(LY_fitted, hitidlist="hit_list.txt", isfile=TRUE)
#  CETSAdata_subset <- ms_subsetting(LY_fitted, hitidlist=c("P00000", "P12345-6"), isfile=FALSE)

