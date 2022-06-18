### Creation of Phyloseq object (phys) with read oriented###


#### ---- Import data ---- ####

#Cleaned metadata
metadata <- read.csv("Data/metadata_clean.csv", row.names = 1)


# load ASV table
seqtab <- readRDS("Data/ReadO/seqtab.nochim.rds") # file from dada2 

#Load taxanomy
#taxa_silva138 <- readRDS("Data/tax_silva138.rds") # file from assign taxonomy in dada2 october2020
taxa_silva138.1 <- readRDS("Data/ReadO/taxtable.rds") # file from assign taxonomy and add species in dada2 with 138.1 version from september 2021

#### ---- Remove short length sequences ---- ####

# Distribution length of asv table
table(nchar(getSequences(seqtab)))
hist(nchar(getSequences(seqtab)))
# abline(v=380) # cut off for short sequences
# 
# # Remove short length sequences from asv table
# seqtab3 <- seqtab2[,nchar(colnames(seqtab2)) %in% 381:441] # remove everything below 381 bp length, they are in columns
# table(nchar(getSequences(seqtab3)))
# hist(nchar(getSequences(seqtab3)))

# Distribution of length of taxonomy and remove short length seqs (remove same length sequences in the taxonomy file) change to tax file you want
# table(nchar(getSequences(taxa_silva138.1)))
# hist(nchar(getSequences(taxa_silva138.1)))
# taxa_silva138.2 <- taxa_silva138.1[nchar(rownames(taxa_silva138.1)) %in% 381:441,] # here they are rows
# table(nchar(getSequences(taxa_silva138.2)))

#### ---- Standard files + add ASV names ---- ####

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")

for (i in 1:dim(seqtab)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# Making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
#write(asv_fasta, "Data/ReadO/ASVs_sequences.fa")

# count table:
asv_tab <- t(seqtab)
row.names(asv_tab) <- sub(">", "", asv_headers)
#write.table(asv_tab, "Data/ReadO/ASVs_tab.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
colnames(taxa_silva138.1) <- ranks
rownames(taxa_silva138.1) <- gsub(pattern=">", replacement="", x=asv_headers)

#write.table(taxa_silva138.2, "Data/ASVs_taxa_silva138.2.tsv", sep = "\t", quote=F, col.names=NA)

#### ---- Phyloseq object ---- ####

ASV = otu_table(asv_tab, taxa_are_rows = TRUE)
TAX = tax_table(taxa_silva138.1)
MET = sample_data(metadata)
REF <- refseq(Biostrings::readDNAStringSet("Data/ReadO/ASVs_sequences.fa", format = "fasta"))

taxa_names(TAX)
taxa_names(ASV)
taxa_names(REF)

sample_names(ASV)
sample_names(MET)


phys <- phyloseq(ASV, TAX, MET, REF)
phys.reado = phys %>%
  prune_taxa(taxa_sums(.) > 0, .) # 1353 ASVs


sample_names(phys.reado)
rank_names(phys.reado)
sample_variables(phys.reado)

#Save phyloseqobject 
saveRDS(phys.reado, "Data/phys-readO.rds")

# Import again
#phys <- readRDS("Data/phys.rds")

