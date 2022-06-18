###########################################################################################################
                                      #### DADA2-denoise #####
                                        #Updated June2022
                                        
###########################################################################################################  
### Reads are demultiplexed with cutadapt - cutadapt removes the barcode (now located in the file name). Primer trimmed with cutadapt and reads reoriented

  #### ---- Packages ---- ####

library("dplyr")
library("tidyr")
library("dada2"); packageVersion("dada2")


  #### ---- Directory ---- ####
# In fastq both forward and reverse are located together

  #### ---- Path and prepare demultiplexed seqs ---- ####
pathF <- "../demultiplexed/FWD2/" # Path to forward reads
pathR <- "../demultiplexed/REV2/" #Path to reverse reads

# Forward and reverse fastq filenames have format: SAMPLENAME_1_trimmed.fastq.gz and SAMPLENAME_2_trimmed.fastq.gz (read oriented)
fnFs <- sort(list.files(pathF, pattern="_1_trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(pathR, pattern="_2_trimmed.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

  #### ---- Quality profiles ---- #### 
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

  #### ---- Filter and trim reads ---- ####

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(pathF, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(pathR, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


# set maxEE to 2 as it is a better measure than average quality score
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2,
                     rm.phix=TRUE, compress=TRUE, 
                     multithread=10, verbose = TRUE)

head(out)

# calculate percentage passing filter and trim (out_perc)
out_perc = data.frame(out)
out_perc$perc_pass = out_perc$reads.out/out_perc$reads.in*100


  ##### ---- Learn error rates ---- ####

#use a parametric error model (err) to learn error rates in dataset.
#Learns by alternating estimation of err rates and inference of sample composition until solution. Starts with initial guess, 
#for which max possible error rate in the data set are used 
#(error rate if only most abundant seq is correct and all the rest are errors)

# Make error models
errF <- learnErrors(filtFs, multithread=10, verbose = TRUE, nbases = 1e9)
errR <- learnErrors(filtRs, multithread=10, verbose = TRUE, nbases = 1e9)


# visualize the estimated error rates
#The error rates for each possible transition (A→C, A→G, …) are shown.
# Points are the observed error rates for each consensus quality score.
#The black line shows the estimated error rates after convergence of the machine-learning algorithm.
#The red line shows the error rates expected under the nominal definition of the Q-score

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
#  If the plotted error model does not look like a good fit, try increasing the nreads parameter to see if the fit improves.


#### --- Dereplication ---- #####

# Combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”
derepF <- derepFastq(filtFs,  verbose=TRUE)
derepR <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepF) <- sample.names
names(derepR) <- sample.names


#### ---- Sample inference ---- ####
dadaFs <- dada(derepF, err=errF, multithread=10)
dadaRs <- dada(derepR, err=errR, multithread=10)


#### ---- Merge ---- ####
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)

#Remove sequences below 380
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 380:434]
dim(seqtab2) #13460 ASVs

#### ---- Get sequencetable with no chim ---- ####
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, "seqtab.nochim.rds")

#### ---- Statistics ---- ####

sum(seqtab.nochim)/sum(seqtab2) #chimeras only account: 0.9651209, around 3%
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
#write.table(track, file = "track_dada2.tsv", sep="\t")


#### ---- Assign taxonomy ---- ####

#Silva database train set
tt = assignTaxonomy(seqtab.nochim, "../silva_assignment_sep2021/silva_nr99_v138.1_train_set.fa.gz", 
                    tryRC = TRUE, verbose = TRUE, multithread = 10)

# Silva database species set
tt.plus = addSpecies(tt, "../silva_assignment_sep2021/silva_species_assignment_v138.1.fa.gz", verbose = TRUE, tryRC = TRUE)

save.image("dada2-17062022.RData")
savehistory(file = "dada2.Rhistory")

rm(list = ls())
