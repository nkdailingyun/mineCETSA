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
Depending on the configuration of the R environment in the user computer, a few dependency packages might need to be manually installed or updated for the proper installation of mineCETSA package.  

Once it is installed, you can load the mineCETSA package by running the following command:  
*> library("mineCETSA")*  
The successful loading of the mineCETSA package shows a welcome message as below:
*“Welcome to mineCETSA! Now it's the time to start mining your CETSA data!”.*  

The main steps and detailed explanations of using mineCETSA package for analyzing CETSA melt curve or ITDR data are described in the accompanying vignettes.  
To access the vignettes, key in and run: **vignette('CETSA_melt_curve', package='mineCETSA')** or **vignette('CETSA_ITDR_ITTR', package='mineCETSA')**, respectively.   
