## ---- message=FALSE------------------------------------------------------
library("mineCETSA")

## ---- eval=FALSE---------------------------------------------------------
#  Chem2D <- ms_2D_rawread(c("20170517_LY_K562_CellCycle_ChemArrest_37C_3bio_rep1_Proteins.txt","20170517_LY_K562_CellCycle_ChemArrest_47C_3bio_rep1_Proteins.txt","20170517_LY_K562_CellCycle_ChemArrest_50C_3bio_rep1_Proteins.txt","20170517_LY_K562_CellCycle_ChemArrest_52C_3bio_rep1_Proteins.txt","20170517_LY_K562_CellCycle_ChemArrest_54C_3bio_rep1_Proteins.txt","20170517_LY_K562_CellCycle_ChemArrest_57C_3bio_rep1_Proteins.txt"), treatment=c("B1_G1S","B1_S","B1_PM","B2_G1S","B2_S","B2_PM","B3_G1S","B3_S","B3_PM","Mix"))

## ---- eval=FALSE---------------------------------------------------------
#  Chem2D_c <- ms_clean(Chem2D)

## ---- eval=FALSE---------------------------------------------------------
#  Chem2D_c <- ms_conditionrename(Chem2D_c, incondition = c("C37C.1","C47C.2","C50C.3","C52C.4","C54C.5","C57C.6"),  outcondition=c("37C","47C","50C","52C","54C","57C"))
#  Chem2D_c <- Chem2D_c[,-13] # to remove the "Mix" channel sample
#  

## ---- eval=FALSE---------------------------------------------------------
#  Chem2D_c1 <- ms_isoform_resolve(Chem2D_c)
#  # The next isoform cleaning is optional. You are encouraged to double check and adjust the match table entries.
#  Chem2D_c2 <- ms_isoform_consolidate(Chem2D_c1, nread=9, matchtable = "./subfolder/tobe_consolidated.txt")
#  

## ----eval=FALSE----------------------------------------------------------
#  Chem2D_c3 <- ms_2D_rearrange(Chem2D_c2, nread=9, repthreshold=0.8, countthreshold=3)

## ----eval=FALSE----------------------------------------------------------
#  Chem2D_s <- ms_2D_normalization(Chem2D_c3)

## ----eval=FALSE----------------------------------------------------------
#  Chem2D_s1 <- ms_2D_caldiff(Chem2D_s, treatmentlevel=c("G1S","S","PM"))

## ---- eval=FALSE---------------------------------------------------------
#  Chem2D_category <- ms_2D_globalview(Chem2D_s1, treatment="PM", labelnodes=TRUE, labelcategory="CC")

## ---- eval=FALSE---------------------------------------------------------
#  ms_2D_barplotting(Chem2D_s1[grep("cyclin",Chem2D_s1$description),], treatmentlevel=c("G1S","S","PM"))
#  
#  Chem2D_CC <- subset(Chem2D_s1, id %in% subset(Chem2D_category, category=="CC")$id)
#  ms_2D_barplotting(Chem2D_CC, treatmentlevel=c("G1S","S","PM"))
#  
#  ms_2D_barplotting(Chem2D_s1, treatmentlevel=c("G1S","S","PM"))

