# Cleaning and subset PS

#Raw PS object 
ps <- readRDS("Data/phys.rds")
                                                                                
#### ---- Library sizes - w. negs and true ---- ####
#Library sizes/number of reads in each sample as a function of whether that sample was a sample or neg-control

# put sample_data into a ggplot-friendly data.frame
sample.data = as.data.frame(sample_data(ps))

#Add read counts to sadata as a variable
sample.data$lib.size = sample_sums(ps)

# Order the sadata after the library sizes
sample.data = sample.data[order(sample.data$lib.size),]

#Make index of the order
sample.data$Index = seq(nrow(sample.data))

# Plot library sizes 
ggplot(sample.data, aes(x=Index, y=lib.size, color=type)) +
  geom_point() + 
  theme_bw() +
  scale_y_continuous(breaks = seq(10000,600000,by=100000)) +
  scale_x_continuous(breaks = seq(0,145,by=10))

#The library sizes of the positive samples primarily fall from 10,000 to 510,00 reads, 
#but there are some low-read outliers. The negative control samples have fewer reads as expected. 
#Note: It is important to keep the low-read samples for now, because we want to use those negative controls to help identify contaminants!

#### ---- Contaminants - prevalence method---- ####

# Start by looking at content in negative compared to positive

# phyloseq object without negative samples                   (variable type in META has two levels sample or negative, keep sample)
ps_clean = ps %>%
  subset_samples(type == "sample") %>%
  prune_taxa(taxa_sums(.) > 0, .) # 2504, 144 ASVs unique to the negatives (?)

# ps-object with only negatives
ps_negonly = ps %>%
  subset_samples(type == "negative") %>%
  prune_taxa(taxa_sums(.) > 0, .)
sample_names(ps_negonly) # 219 ASVs found in the negatives


# Compare ps_sub with ps                    (8 negative samples have to be removed)
ps_clean
ps_negonly

# make data.frame with only negatives to allow for manual sorting and observation
ps_negonly_df_genus = ps_negonly %>%
  psmelt()

# Plot only negatives to hav a look
ps_negonly_df_genus2 = ps_negonly %>%
  tax_glom(taxrank = "genus") %>% # gloomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>%  #Transform to relative abundance
  psmelt() %>%                                          # Melt into long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(genus)                                      # Sort data frame alphabetically by phylum

#Subset for only one negative sample (if needed)
only.one.sample = ps_negonly_df_genus2[ which(ps_negonly_df_genus2$Sample=='91'),] # to look only at one negative per time

ggplot(only.one.sample, aes(x = Sample, y = Abundance, fill = genus)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("Relative Abundance > 1% \n") +
  xlab("\n Sample ID")+
  ggtitle("Negative control at genus level")

# Save csv with negatives only (to allow for manual inspection)
#write.csv(ps_negonly_df_genus, file="Data/negatives_samples.csv") #Adjusted the negatives

# Import csv with ASV that are true negatives
true.negatives = read.table("Data/negatives_sorted_true_ASVlist.csv")
true.negatives = true.negatives[ , "V1"]

# Remove True negatives ASVs

alltaxa = taxa_names(ps_clean)
alltaxa.t = alltaxa[!(alltaxa %in% true.negatives)]
ps_clean_nc = prune_taxa(alltaxa.t, ps_clean) %>%
  prune_taxa(taxa_sums(.) > 0, .) # 16 ASVs were removed, 2488 ASVs

#### ---- Filtering and subsetting ---- ####

# Remove eukaryotes, archea, chloroplasts and mitochondria
ps_filt <- ps_clean_nc %>%
  subset_taxa( domain == "Bacteria" & family  != "mitochondria" & class   != "Chloroplast") %>%
  prune_taxa(taxa_sums(.) > 0, .)# removed 574 ASVs, 1914 ASvs

# Remove additional Escherichia, Staphy, serratia, cutibacterium etc. 
ps_filt1 = ps_filt %>%
  subset_taxa( genus !="Escherichia-Shigella" & genus !="Staphylococcus" & 
                 genus !="Serratia" & genus !="Cutibacterium" & 
                 genus !="Thermomonas" & genus !="Dermacoccus") %>%
  prune_taxa(taxa_sums(.) > 0, .)# removed 716 ASVs, 1198 ASVs

#### ---- Sample depths (should low read sample out)--- ####

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(ps_filt1), sampleID = sample_names(ps_filt1))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

#Max and min
min(sample_sums(ps_filt1))
max(sample_sums(ps_filt1)) 

ggplot(sample_sum_df, aes(x= reorder(sampleID, -sum), y= sum))+
  geom_bar(color= "black", fill = "indianred", stat = "identity") +
  ggtitle("Distriution of sample sequencing depth") +
  xlab("sample.ID") +
  ylab("Read counts") +
  theme(axis.text.x = element_text(angle=90, vjust = 0.4))+ 
  geom_hline(yintercept = 1000, linetype = "dashed", size = 1)

sample_sum_df = sample_sum_df[order(sample_sum_df$sum),]

# Exclude low read samples with reads below 10,000 (sampleID 5,7,10,35,140)
ps_filt2 = ps_filt1 %>%
  prune_samples(sample_names(.) != "140", .) %>%
  prune_samples(sample_names(.) != "5", .) %>%
  prune_samples(sample_names(.) != "7", .) %>%
  prune_samples(sample_names(.) != "10", .) %>%
  prune_samples(sample_names(.) != "35", .) %>%
  prune_samples(sample_names(.) != "94", .) %>%
  prune_taxa(taxa_sums(.) > 0, .)

#Save the cleaned data

# Save filtered and cleaned PS object  # Final PS object with 1091 ASVs and 126 samples 
#save(ps_filt2, file = "Data/ps.clean.RData")
# use this dataset for the Phaeobacter dynamics (part one of results)


# Additional cleaning #

# Filter to keep taxa that has more than one read in more than one sample 
ps.reduced <- filter_taxa(ps_filt2, function (x) {sum(x > 0) > 1}, prune=TRUE) # Final reduced dataset has 159 taxa (removed 932 ASVs)
# Save as Phyloseq object 
#save(ps.reduced, file = "Data/ps.reduced.RData")

### Filter off genera with too low abundance in total dataset
load("Data/ps.clean.RData")
ps.clean = ps_filt2

ps.genus = ps.clean %>%
  tax_glom(taxrank = "genus")
get_taxa_unique(ps.genus, "genus") #202 genera 
#save(ps.genus, file = "Data/ps.genus.RData")

#Find read cut-off for generas 
df.genus = ps.genus %>%
  psmelt()

df.genus.sum = df.genus %>%
  group_by(genus) %>%
  summarise(Total.abundance = sum(Abundance))

#Make relative abundance
df.genus.sum$RA = (100*(df.genus.sum$Total.abundance/sum(df.genus.sum$Total.abundance)))

# factor for below 0.5%
df.05 = df.genus.sum[df.genus.sum$RA < 0.5,]
df.05$cutoff ="Below"
df.06 = df.genus.sum[df.genus.sum$RA >= 0.5,]
df.06$cutoff ="Above"
df.genus.RA = rbind(df.06, df.05)

#Elbow plot
ggplot(df.genus.sum, aes(x= reorder(genus, -Total.abundance), y = (Total.abundance))) +
  geom_bar(color= "black", fill = "indianred", stat = "identity") +
  ggtitle("Distribution of reads on genera") +
  xlab("Genus") +
  ylab("Relative abundance %") +
  theme(axis.text.x = element_text(angle=90, vjust = 0.4))
#Looking at the plot we decide to remove all genera that has less that 100 reads across the total dataset.

# Take a look in table
genus.totalsum = df.genus.sum[order(df.genus.sum$Total.abundance),]
  
# Percentage genera with reads <100 account for of total reads = 0.01879076
100*sum(df.genus.sum$Total.abundance[df.genus.sum$Total.abundance<100])/sum(df.genus.sum$Total.abundance)

#Remove these genera from phyloseq object #WINNER FOR NOW#
ps.genus.reduced<- filter_taxa(ps.genus, function (x) {sum(x > 100) > 0}, prune=TRUE)
save(ps.genus.reduced, file = "Data/ps.genus.reduced.RData")
# Reduced ps.genus has 38 taxa (removed 164 generas)

#Check plot
df.genus.reduced = ps.genus.reduced %>%
  psmelt()
df.genus.reduced.sum = df.genus.reduced %>%
  group_by(genus) %>%
  summarise(Total.abundance = sum(Abundance))
ggplot(df.genus.reduced.sum, aes(x= reorder(genus, -Total.abundance), y = log10(Total.abundance))) +
  geom_bar(color= "black", fill = "indianred", stat = "identity") +
  ggtitle("Distribution of reads on genera") +
  xlab("Genus") +
  ylab("Reads in total") +
  theme(axis.text.x = element_text(angle=90, vjust = 0.4)) +
  geom_hline(yintercept = 2, linetype = "dashed", size = 1)

#Following the reduction by more than one read in more than one sample
ps.genus.reduced2<- filter_taxa(ps.genus, function (x) {sum(x > 1) > 1}, prune=TRUE) 
# Gives 109 genera (removed 93 genera)

### Final plots ####

#Rare curve
rarecurve(t(otu_table(ps_filt2)), step=50, cex=0.5)


# library size plot

df = as.data.frame(sample_data(ps_filt2))
df$lib.size = sample_sums(ps_filt2)
df = df[order(df$lib.size),]
df$Index = seq(nrow(df))
df1 = df[1:6,]
df1$Size ="Below"
df2 = df[7:132,]
df2$Size ="Above"
df3 = rbind(df1, df2)

#Plot for library size plot
ggplot(df3, aes(x=Index, y=lib.size, color = Size)) +
  geom_point() + 
  theme_bw() +
  xlab("\n Sample ID")+
  ylab("Library size (reads) \n")+
  scale_y_continuous(breaks = seq(10000,600000,by=100000)) +
  scale_x_continuous(breaks = seq(0,145,by=10)) +
  scale_color_manual(values= system_color2) +
  theme(axis.text.x = element_text(size=14))+
  theme(axis.title = element_text(size = 16))+
  theme(axis.text = element_text(size = 14, color = "Black"))+
  theme(legend.position = "none") 
