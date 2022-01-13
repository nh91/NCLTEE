# Natural Community Long-Term Evolution Experiment (NCLTEE)

 ### Data and scripts to reproduce results from the NCLTEE experiment 
 
 ## Guide
 
 #### Software
 In 'Scripts/' you will find the script 'Packages.R' that contains all packages used.
 
 #### Data

  All raw sequence reads can be found in BioProject (SRA): PRJNA795303. Additional data (and data generated from the raw reads) can be found in the folder: 'Data/'. Specifically:
 
 - The phyloseq object (PS) 'ps.genus.reduced.RData' is used for figure plots found in 'Scripts'.
 - The PS object 'ps.clean.RData' is used for 'Fig5C_PhaeobacterASVs.R'.
 - The folders 'MetalonDA_*' contains the intermediate figures and results from the MetaLonDA analysis (Figure 4). 
 
 #### Scripts
 All scripts used for analysis can be found in 'Scripts/'. The scripts are named according to their figure number in the manuscript.

The preprocessing and creation of the PS object(s) can be found in 'Scripts' using the following scripts:
 
 - 1.DADA2_pipeline.R
 - 2.CreatePS.R
 - 3.CleaningPS.R

 
  #### SNP analysis (table S3)
  
  
 