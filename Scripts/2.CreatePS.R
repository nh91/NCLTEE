### Creation of Phyloseq object (phys) ###


#### ---- Import data ---- ####

# Old Metadata
metadata_old <- read.table("Data/metadata.tsv", sep="\t", row.names = 2)
metadatatest = subset(metadata_old, select = c("sample.name", "type", "time", "exp.condition", "bio.rep")) #remove unnesseray variables

#Cleaned metadata - simplified
metadata <- read.csv("Data/metadata_clean.csv", row.names = 1)


# load ASV table
seqtab <- readRDS("Data/seqtab.nochim.rds") # file from dada2 

# remove samples not related to the project
seqtab2 <- seqtab[-c(141:150),]

#Load taxanomy
#taxa_silva138 <- readRDS("Data/tax_silva138.rds") # file from assign taxonomy in dada2 october2020
taxa_silva138.1 <- readRDS("Data/tax_silva138.1_wSpecies.rds") # file from assign taxonomy and add species in dada2 with 138.1 version from september 2021

#### ---- Remove short length sequences ---- ####

# Distribution length of asv table
table(nchar(getSequences(seqtab2)))
hist(nchar(getSequences(seqtab2)))
abline(v=380) # cut off for short sequences

# Remove short length sequences from asv table
seqtab3 <- seqtab2[,nchar(colnames(seqtab2)) %in% 381:441] # remove everything below 381 bp length, they are in columns
table(nchar(getSequences(seqtab3)))
hist(nchar(getSequences(seqtab3)))

# Distribution of length of taxonomy and remove short length seqs (remove same length sequences in the taxonomy file) change to tax file you want
table(nchar(getSequences(taxa_silva138.1)))
hist(nchar(getSequences(taxa_silva138.1)))
taxa_silva138.2 <- taxa_silva138.1[nchar(rownames(taxa_silva138.1)) %in% 381:441,] # here they are rows
table(nchar(getSequences(taxa_silva138.2)))


#### ---- Rep.seqs fasta with all length of seqs --- ####

# get representative sequences with all lengths (before filter out short seqs)
#asv_alllen <- colnames(seqtab2)
#asv_headers_alllen = vector(dim(seqtab2) [2], mode = "character")


#for (i in 1:dim(seqtab2)[2]) {
#asv_headers_alllen[i] <- paste(">ASV", i, sep="_")
#}
#asv_alllen_fasta <- c(rbind(asv_headers_alllen, asv_alllen))
#write(asv_alllen_fasta, "shortlenseqs/ASVs_alllen_sequences.fa")


#### ---- Standard files + add ASV names ---- ####

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab3)
asv_headers <- vector(dim(seqtab3)[2], mode="character")

for (i in 1:dim(seqtab3)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# Making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "Data/ASVs_sequences.fa")

# count table:
asv_tab <- t(seqtab3)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "Data/ASVs_tab.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
colnames(taxa_silva138.2) <- ranks
rownames(taxa_silva138.2) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(taxa_silva138.2, "Data/ASVs_taxa_silva138.2.tsv", sep = "\t", quote=F, col.names=NA)

#### ---- Phyloseq object ---- ####

ASV = otu_table(asv_tab, taxa_are_rows = TRUE)
TAX = tax_table(taxa_silva138.2)
MET = sample_data(metadata)
REF <- refseq(Biostrings::readDNAStringSet("Data/ASVs_sequences.fa", format = "fasta"))

taxa_names(TAX)
taxa_names(ASV)
taxa_names(REF)

sample_names(ASV)
sample_names(MET)


phys <- phyloseq(ASV, TAX, MET, REF)
phys = phys %>%
  prune_taxa(taxa_sums(.) > 0, .) # 2648 ASVs


sample_names(phys)
rank_names(phys)
sample_variables(phys)

#Save phyloseqobject 
saveRDS(phys, "Data/phys.rds")

# Import again
#phys <- readRDS("Data/phys.rds")

