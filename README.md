# Natural Community Long-Term Evolution Experiment (NCLTEE)

 ### Data and scripts to reproduce results from the publication 
'The ability of Phaeobacter inhibens to produce tropodithietic acid influences the community dynamics of a microalgal microbiome'
published in ISME Communcations.

 
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
  
  We used the <i>snippy</i> tool (https://github.com/tseemann/snippy) to predicted SNP in the 48 <i>Phaeobacter</i> isolates. 
  
  To reproduce the results follow the installation guide at https://github.com/tseemann/snippy and collect the following input:
  - The raw reads from the 48 isolates from BioProject (SRA): PRJNA795303. 
  - Update the 'input.txt' file (found in 'Data/') with your path to the R1 and R2 reads. 
  - Get the reference genome (genbank file) from <i>Phaobacter inhibens </i> DSM17395 (GenBank assembly accession GCA_000154765.2).
  
  Then run the following code:

snippy-multi input.txt --ref GCA_000154765.2_ASM15476v2_genomic.gbff --cpus 20 > runme.sh
 
