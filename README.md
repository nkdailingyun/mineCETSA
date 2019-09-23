# mineCETSA
**A package for Processing and Visualization of Proteome-wide MS-CETSA data**


## What is CETSA?  
The Cellular Thermal Shift Assay (CETSA) (orginially described in [Science 341(6141):84-87](http://www.sciencemag.org/lookup/doi/10.1126/science.1233606)) is a biophysical assay based on the principle of ligand-induced thermal stabilization of target proteins, meaning that a protein's melting temperature will change upon ligand interaction.  

By heating samples (lysate, cells or tissue pieces) to different temperatures, and quantifying proteins in the soluble fraction we can detect altered protein interactions after for example drug treatment. This can either be done for selected proteins of interest by using antibody-based assays or on a proteome-wide level by using mass spectrometry.  

CETSA allows direct monitoring of ligand binding to a specific target (target engagement) in lysate, live cells or even tissue pieces. It can also be used to study downstream effects on protein interaction, providing a novel perspective on protein function in situ. For more details please refer to the [CETSA website](https://www.cetsa.org/about).  

## How to install mineCETSA?  
To install mineCETSA from GitHub, first, you need to invoke R, typically by double-clicking the R or RStudio application. Then, type and run the following commands in R console:  
*> install.packages("devtools")*  
*> library("devtools")*  
*> install_github("nkdailingyun/mineCETSA")*  

Depending on the configuration of the R environment in the user computer, a few dependency packages might need to be manually installed or updated for the proper installation of mineCETSA package. More specifically, the following packages are necessary:  

>R (>= 3.3), arrayQualityMetrics, Biobase, dplyr (>= 0.7.0), drc (>= 3.0), fdrtool, ggplot2 (>= 3.0.0), ggpubr, ggrepel (>= 0.7.0), grid, gridExtra, gtools, limma, Nozzle.R1, plyr (>= 1.8.0), RColorBrewer, readr (>=1.2.0), scales, tidyr (>= 0.8.1), VennDiagram

Once it is installed, you can load the mineCETSA package by running the following command:  
*> library("mineCETSA")*  
The successful loading of the mineCETSA package shows a welcome message as below:
*“Welcome to mineCETSA! Now it's the time to start mining your CETSA data!”.*  

The main steps and detailed explanations of using mineCETSA package for analyzing CETSA melt curve or ITDR data are described in the accompanying vignettes.  
To access the vignettes, key in and run: **vignette('CETSA_melt_curve', package='mineCETSA')** or **vignette('CETSA_ITDR_ITTR', package='mineCETSA')**, respectively.   



## sessionInfo()
R version 3.5.3 (2019-03-11)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: OS X El Capitan 10.11.6

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
 [1] zip_2.0.1         Rcpp_1.0.1        cellranger_1.1.0  pillar_1.3.1     
 [5] compiler_3.5.3    forcats_0.4.0     tools_3.5.3       tibble_2.1.1     
 [9] lattice_0.20-38   pkgconfig_2.0.2   rlang_0.4.0       openxlsx_4.1.0   
[13] Matrix_1.2-16     rstudioapi_0.10   yaml_2.2.0        curl_3.3         
[17] mvtnorm_1.0-10    haven_2.1.0       rio_0.5.16        gtools_3.8.1     
[21] drc_3.0-1         hms_0.4.2         grid_3.5.3        data.table_1.12.0
[25] plotrix_3.7-4     readxl_1.3.1      survival_2.43-3   foreign_0.8-71   
[29] multcomp_1.4-10   TH.data_1.0-10    carData_3.0-2     car_3.0-2        
[33] magrittr_1.5      scales_1.0.0      codetools_0.2-16  splines_3.5.3    
[37] MASS_7.3-51.1     abind_1.4-5       colorspace_1.4-1  sandwich_2.5-0   
[41] munsell_0.5.0     crayon_1.3.4      zoo_1.8-4      

## License

This project is covered under the GNU General Public License, version 3.0 (GPL-3.0).
