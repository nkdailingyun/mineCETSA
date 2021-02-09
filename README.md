# mineCETSA
**A package for Processing and Visualization of Proteome-wide MS-CETSA data**


## What is CETSA?  
The Cellular Thermal Shift Assay (CETSA) (orginially described in [Science 341(6141):84-87 (2013)](http://www.sciencemag.org/lookup/doi/10.1126/science.1233606)) is a biophysical assay based on the principle of ligand-induced thermal stabilization of target proteins, meaning that a protein's melting temperature will change upon ligand interaction.  

By heating samples (lysate, cells or tissue pieces) to different temperatures, and quantifying proteins in the soluble fraction we can detect altered protein interactions after for example drug treatment. This can either be done for selected proteins of interest by using antibody-based assays or on a proteome-wide level by using mass spectrometry.  

CETSA allows direct monitoring of ligand binding to a specific target (target engagement) in lysate, live cells or even tissue pieces. It can also be used to study downstream effects on protein interaction, providing a novel perspective on protein function in situ. For more details please refer to the [CETSA website](https://www.cetsa.org/about).  

## How to install mineCETSA?  
To install mineCETSA from GitHub, first, you need to invoke R, typically by double-clicking the R or RStudio application. Then, type and run the following commands in R console:  
*> install.packages("devtools")*  
*> library("devtools")*  
*> install_github("nkdailingyun/mineCETSA")*  

Depending on the configuration of the R environment in the user computer, a few dependency packages might need to be manually installed or updated for the proper installation of mineCETSA package. More specifically, the following packages are necessary:  

>R (>= 4.0), arrayQualityMetrics, Biobase, dplyr (>= 1.0.0), drc (>= 3.0), fdrtool, ggplot2 (>= 3.0.0), ggpubr, ggrepel (>= 0.7.0), grid, gridExtra, gtools, limma, MESS, Nozzle.R1, plyr (>= 1.8.0), RColorBrewer, readr (>=1.2.0), scales, tidyr (>= 1.0.0), VennDiagram

For the ones in bioconductor, the package installation is as follows:
*if (!requireNamespace("BiocManager", quietly = TRUE))*
*install.packages("BiocManager")*

*BiocManager::install("arrayQualityMetrics")*

Once it is installed, you can load the mineCETSA package by running the following command:  
*> library("mineCETSA")*  
The successful loading of the mineCETSA package shows a welcome message as below:
*“Welcome to mineCETSA! Now it's the time to start mining your CETSA data!”.*  

The main steps and detailed explanations of using mineCETSA package for analyzing CETSA melt curve or ITDR data are described in the accompanying vignettes.  
To access the vignettes, key in and run: **vignette('CETSA_melt_curve', package='mineCETSA')** or **vignette('CETSA_ITDR_ITTR', package='mineCETSA')**, respectively.   

For the detailed protocol of installing and using the mineCETSA package, the reader could refer to [Nature Protocols 15(6),1881–1921(2020)](https://www.nature.com/articles/s41596-020-0310-z) 

## sessionInfo()
R version 4.0.3 (2020-10-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

Random number generation:
 RNG:     Mersenne-Twister 
 Normal:  Inversion 
 Sample:  Rounding 
 
locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.6        mvtnorm_1.1-1     lattice_0.20-41   tidyr_1.1.2      
 [5] zoo_1.8-8         gtools_3.8.2      ggforce_0.3.2     R6_2.5.0         
 [9] cellranger_1.1.0  plyr_1.8.6        ggridges_0.5.3    backports_1.2.1  
[13] labelled_2.7.0    ggstance_0.3.5    ggplot2_3.3.3     pillar_1.4.7     
[17] geepack_1.3-2     rlang_0.4.10      curl_4.3          multcomp_1.4-16  
[21] readxl_1.3.1      data.table_1.13.6 car_3.0-10        Matrix_1.3-2     
[25] splines_4.0.3     readr_1.4.0       stringr_1.4.0     foreign_0.8-81   
[29] polyclip_1.10-0   munsell_0.5.0     geeM_0.10.1       broom_0.7.4      
[33] compiler_4.0.3    pkgconfig_2.0.3   tidyselect_1.1.0  mosaicCore_0.9.0 
[37] tibble_3.0.6      gridExtra_2.3     rio_0.5.16        codetools_0.2-18 
[41] crayon_1.4.1      dplyr_1.0.4       MASS_7.3-53       grid_4.0.3       
[45] gtable_0.3.0      lifecycle_0.2.0   ggformula_0.10.1  magrittr_2.0.1   
[49] MESS_0.5.7        scales_1.1.1      zip_2.1.1         stringi_1.5.3    
[53] carData_3.0-4     farver_2.0.3      reshape2_1.4.4    ellipsis_0.3.1   
[57] drc_3.0-1         generics_0.1.0    vctrs_0.3.6       sandwich_3.0-0   
[61] openxlsx_4.2.3    TH.data_1.0-10    tools_4.0.3       forcats_0.5.1    
[65] glue_1.4.2        tweenr_1.0.1      purrr_0.3.4       hms_1.0.0        
[69] abind_1.4-5       plotrix_3.8-1     survival_3.2-7    yaml_2.2.1       
[73] colorspace_2.0-0  haven_2.3.1       

## License

This project is covered under the GNU General Public License, version 3.0 (GPL-3.0).
