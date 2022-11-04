# Natural Community Long-Term Evolution Experiment (NCLTEE)

 ##### Data and scripts to reproduce results from the publication 
Henriksen, N.N.S.E., Schostag, M.D., Balder, S.R. et al. 'The ability of <i> Phaeobacter inhibens </i> to produce tropodithietic acid influences the community dynamics of a microalgal microbiome.' ISME COMMUN. 2, 109 (2022). 
https://doi.org/10.1038/s43705-022-00193-6
 
 ## Guide
 
 #### Software
 In 'Scripts/' you will find the script 'Packages.R' that contains all packages used.
 
 #### Data

  All raw sequence reads can be found in BioProject (SRA): PRJNA795303. Additional data (and data generated from the raw reads) can be found in the folder: 'Data/'.
 
 #### Scripts
 All scripts used for analysis can be found in 'Scripts/'. The scripts are named according to their figure number in the publication.

The preprocessing and creation of the PS object(s) can be found in 'Scripts' using the following scripts:
 
 - 1.DADA2pipeline.R
 - 2.CreatePS-ASV.R
 - 3.PreprocessPS-ASV.R
 - 4.PhytreeCreation-ASV.R
 
  #### SNP analysis (table S3)
  
  We used the <i>snippy</i> tool (https://github.com/tseemann/snippy) to predicted SNP in the 48 <i>Phaeobacter</i> isolates. 
  
  To reproduce the results follow the installation guide at https://github.com/tseemann/snippy and collect the following input:
  - The raw reads from the 48 isolates from BioProject (SRA): PRJNA795303. 
  - Update the 'input.txt' file (found in 'Data/') with your path to the R1 and R2 reads. 
  - Get the reference genome (genbank file) from <i>Phaobacter inhibens </i> DSM17395 (GenBank assembly accession GCA_000154765.2).
  
  Then run the following code:

snippy-multi input.txt --ref GCA_000154765.2_ASM15476v2_genomic.gbff --cpus 20 > runme.sh
 
