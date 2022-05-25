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
# use this dataset for the Phaeobacter dynamics (Fig5C)


# Additional cleaning #

### Filter off genera with too low abundance in total dataset
load("Data/ps.clean.RData")
ps.clean = ps_filt2

# ps.genus = ps.clean %>%
#   tax_glom(taxrank = "genus")
# get_taxa_unique(ps.genus, "genus") #202 genera 
# #save(ps.genus, file = "Data/ps.genus.RData")

#Find read cut-off for generas 
df.asv = ps.clean %>%
  psmelt()

df.asv.sum = df.asv %>%
  group_by(OTU) %>%
  summarise(Total.abundance = sum(Abundance))

#Make relative abundance
df.asv.sum$RA = (100*(df.asv.sum$Total.abundance/sum(df.asv.sum$Total.abundance)))

# factor for below 0.5%
df.05 = df.asv.sum[df.asv.sum$RA < 0.5,]
df.05$cutoff ="Below"
df.06 = df.asv.sum[df.asv.sum$RA >= 0.5,]
df.06$cutoff ="Above"
df.asv.RA = rbind(df.06, df.05)

#Elbow plot
ggplot(df.asv.sum, aes(x= reorder(OTU, -Total.abundance), y = (Total.abundance))) +
  geom_bar(color= "black", fill = "indianred", stat = "identity") +
  ggtitle("Distribution of reads on ASV") +
  xlab("ASV") +
  geom_hline(yintercept = 100, linetype = "dashed", color ="red", size = 1)+
  ylab("Total abundance") +
  theme(axis.text.x = element_text(angle=90, vjust = 0.4))
#Looking at the plot we decide to remove all ASVs that has less that 100 reads across the total dataset.

# Take a look in table
asv.totalsum = df.asv.sum[order(df.asv.sum$Total.abundance),]
  
# Percentage ASV with reads <100 account for of total reads = 0.07006834
100*sum(df.asv.sum$Total.abundance[df.asv.sum$Total.abundance<100])/sum(df.asv.sum$Total.abundance)

#Remove these genera from phyloseq object #WINNER FOR NOW#
ps.asv.reduced <- filter_taxa(ps.clean, function (x) {sum(x > 1) > 1}, prune=TRUE)
save(ps.asv.reduced, file = "Data/ps.asv.reduced.RData")
# Reduced ps.genus has 158 taxa (removed 933 asvs)

#Check plot
df.asv.reduced = ps.asv.reduced %>%
  psmelt()
df.asv.reduced.sum = df.asv.reduced %>%
  group_by(OTU) %>%
  summarise(Total.abundance = sum(Abundance))

ggplot(df.asv.reduced.sum, aes(x= reorder(OTU, -Total.abundance), y = log10(Total.abundance))) +
  geom_bar(color= "black", fill = "indianred", stat = "identity") +
  ggtitle("Distribution of reads on ASVs") +
  xlab("ASV") +
  ylab("log(Reads in total") +
  theme(axis.text.x = element_text(angle=90, vjust = 0.4))
  #geom_hline(yintercept = 1, linetype = "dashed", size = 1, color = "red")

### Final plots ####

#Rare curve



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
