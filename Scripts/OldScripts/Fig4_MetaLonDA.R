#####################################################################################

              ### Differential abundance analysis - MetaLonDA ###

#####################################################################################


#################
### Load data ###
#################

load("Data/ps.asv.reduced.wTree.RData") # PS with tree


####################
### TDA vs NoTDA ###
####################

# PS object with no control samples
ps.asv.TDA.noTDA = ps.new %>%
  subset_samples(exp.condition != "Control") %>%
  subset_taxa(genus != "Phaeobacter") %>%
  subset_samples(time != "0") %>%
  phyloseq::prune_taxa(taxa_sums(.) > 0, .) %>%
  phyloseq::filter_taxa(function (x) {sum(x > 0) > 5}, prune = TRUE) #Keep only ASVs that are present in more than 5 replicates


# Metadata
meta = as.data.frame(unclass(sample_data(ps.asv.TDA.noTDA)))

#Add sample id (from count table) to metadata
Sample = psmelt(ps.asv.TDA.noTDA)[,2]
sample.name = psmelt(ps.asv.TDA.noTDA)[,4]
sampleinfo = as.data.frame(cbind(Sample, sample.name))
sampleinfo = sampleinfo[!duplicated(sampleinfo$Sample),] 
meta.table = merge(meta, sampleinfo, by = 'sample.name')

#Make data input for MetaLonDA
meta.table$time <- as.numeric(as.character(meta.table$time))
meta.table$exp.condition<-as.factor(meta.table$exp.condition)
meta.table$Sample<-as.factor(meta.table$Sample)

#Make count table with tax names
#taxonomy
# tax.table = ps.asv.TDA.noTDA@tax_table@.Data 
# rownames(count.table) = tax.table[,6] # change rownames(ASV) to taxa names
# count.table = as.matrix(count.table) #make matrix
# rownames(count.table)<-as.character(rownames(count.table))
# test = as.numeric(asv.table)

#OTU table
#Transpose rows and columns
count.table = t(abundances(ps.asv.TDA.noTDA))

#Wrangle dataframe to get ID.name and correct sample name
Sample = rownames(count.table)
count.table.2 = cbind(Sample, count.table)
count.table.3 = merge(meta.table[,c(1,3:6,8)], count.table.2, by ='Sample') %>%
unite(ID.name, c(exp.condition, bio.rep), sep = "-", remove = FALSE) %>%
  arrange(ID.name, day,time)
rownames(count.table.3) = count.table.3$sample.name
count.table.3$id = 0
count.table.4 = count.table.3 %>%
  mutate(id = ifelse(ID.name == "NoTDA-1", "1", id)) %>%
  mutate(id = ifelse(ID.name == "NoTDA-2", "2", id)) %>%
  mutate(id = ifelse(ID.name == "NoTDA-3", "3", id)) %>%
  mutate(id = ifelse(ID.name == "NoTDA-4", "4", id)) %>%
  mutate(id = ifelse(ID.name == "TDA-1", "5", id)) %>%
  mutate(id = ifelse(ID.name == "TDA-2", "6", id)) %>%
  mutate(id = ifelse(ID.name == "TDA-3", "7", id)) %>%
  mutate(id = ifelse(ID.name == "TDA-4", "8", id))
  
#Make count table for use
asv.table = as.data.frame(t(count.table.4[,c(8:90)])) %>%
  dplyr::mutate(dplyr::across(.cols = 1:79, .fns=as.integer))
str(asv.table)

#Input for function
time.vec <-count.table.3$time
group.vec <-count.table.3$exp.condition
id.vec <-count.table.4$id

#Run MetaLonDA  
output.NoTDAvsTDA = metalondaAll(Count = asv.table, 
                          Time = time.vec, 
                          Group = group.vec,
                          ID = id.vec, 
                          n.perm =100, 
                          fit.method = "lowess",
                          num.intervals = 100,
                          parall = FALSE, 
                          pvalue.threshold=0.05,
                          adjust.method = "BH", 
                          time.unit = "Time_point",
                          norm.method = "median_ratio", 
                          prefix = "MetalonDA_TDAvsNoTDA_MR_lowess", 
                          ylabel = "Read Counts", col = c("#7294D4","#5B1A18"))


########################
### NoTDA vs Control ###
########################

# PS object with no TDA samples
ps.asv.Control.noTDA = ps.asv.reduced %>%
  subset_samples(exp.condition != "TDA") %>%
  subset_taxa(genus != "Phaeobacter") %>%
  subset_samples(time != "0") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  filter_taxa(function (x) {sum(x > 0) > 5}, prune=TRUE)


# Metadata
meta = as.data.frame(unclass(sample_data(ps.asv.Control.noTDA)))

#Add sample id to metadata
Sample = psmelt(ps.asv.Control.noTDA)[,2]
sample.name = psmelt(ps.asv.Control.noTDA)[,4]
sampleinfo = as.data.frame(cbind(Sample, sample.name))
sampleinfo = sampleinfo[!duplicated(sampleinfo$Sample),] 
meta.table = merge(meta, sampleinfo, by = 'sample.name')

#Make data input for MetaLonDA
meta.table$time<-as.numeric(as.character(meta.table$time))
meta.table$exp.condition<-as.factor(meta.table$exp.condition)
meta.table$Sample<-as.factor(meta.table$Sample)
count.table.test = as.data.frame(count.table[1:2,])

#OTU table
#Transpose rows and columns
count.table = t(abundances(ps.asv.Control.noTDA))

#Wrangle dataframe to get ID.name and correct sample name
Sample = rownames(count.table)
count.table.2 = cbind(Sample, count.table)
count.table.3 = merge(meta.table[,c(1,3:6,8)], count.table.2, by ='Sample') %>%
  unite(ID.name, c(exp.condition, bio.rep), sep = "-", remove = FALSE) %>%
  arrange(ID.name, day,time)
rownames(count.table.3) = count.table.3$sample.name
count.table.3$id = 0
count.table.4 = count.table.3 %>%
  mutate(id = ifelse(ID.name == "Control-1", "1", id)) %>%
  mutate(id = ifelse(ID.name == "Control-2", "2", id)) %>%
  mutate(id = ifelse(ID.name == "Control-3", "3", id)) %>%
  mutate(id = ifelse(ID.name == "Control-4", "4", id)) %>%
  mutate(id = ifelse(ID.name == "NoTDA-1", "5", id)) %>%
  mutate(id = ifelse(ID.name == "NoTDA-2", "6", id)) %>%
  mutate(id = ifelse(ID.name == "NoTDA-3", "7", id)) %>%
  mutate(id = ifelse(ID.name == "NoTDA-4", "8", id))

#Make count table for use
asv.table = as.data.frame(t(count.table.4[,c(8:93)])) %>%
  mutate(across(.cols = 1:78, .fns=as.integer))
str(asv.table)

time.vec <-count.table.3$time
group.vec <-count.table.3$exp.condition
id.vec <-count.table.4$id

#Run MetaLonDA
output.ControlvsNoTDA = metalondaAll(Count = asv.table, 
                          Time = time.vec, 
                          Group = group.vec,
                          ID = id.vec, 
                          n.perm =1000, 
                          fit.method = "nbinomial", 
                          num.intervals = 100,
                          parall = FALSE, 
                          pvalue.threshold=0.05,
                          adjust.method = "BH", 
                          time.unit = "Time_point",
                          norm.method = "median_ratio", 
                          prefix = "MetalonDA_NoTDAvsControl", 
                          ylabel = "Read Counts", col = c("#7294D4", "#D8A499"))



######################
### TDA vs Control ###
######################

# PS object with no NoTDA samples
ps.genus.Control.TDA = ps.genus.reduced %>%
  subset_samples(exp.condition != "NoTDA") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  filter_taxa(function (x) {sum(x > 0) > 5}, prune=TRUE)

#Metadata
meta = as.data.frame(unclass(sample_data(ps.genus.Control.TDA)))

#Add sample id to metadata
Sample = psmelt(ps.genus.Control.TDA)[,2]
sample.name = psmelt(ps.genus.Control.TDA)[,4]
sampleinfo = as.data.frame(cbind(Sample, sample.name))
sampleinfo = sampleinfo[!duplicated(sampleinfo$Sample),] 
meta.table = merge(meta, sampleinfo, by = 'sample.name')

#Make data input for MetaLonDA
meta.table$time<-as.numeric(as.character(meta.table$time))
meta.table$exp.condition<-as.factor(meta.table$exp.condition)
meta.table$Sample<-as.factor(meta.table$Sample)

#Inputs
time.vec <-meta.table$time
group.vec <-meta.table$exp.condition
id.vec <-meta.table$Sample

#Make count table with tax names
tax.table = ps.genus.Control.TDA@tax_table@.Data #taxonomy
count.table = abundances(ps.genus.Control.TDA) #OTU table
rownames(count.table) = tax.table[,6] # change rownames(ASV) to tax names
count.table = as.matrix(count.table) #make matrix
rownames(count.table)<-as.character(rownames(count.table))

points= seq(0,10, length.out=10)


#Run MetaLonDA started 
output.all.3 = metalondaAll(Count = count.table, 
                            Time = time.vec, 
                            Group = as.factor(group.vec),
                            ID = id.vec, 
                            n.perm =1000, 
                            fit.method = "nbinomial", 
                            num.intervals = 10,
                            parall = FALSE, 
                            pvalue.threshold=0.05,
                            adjust.method = "BH", 
                            time.unit = "Time_point",
                            norm.method = "median_ratio", 
                            prefix = "MetalonDA_TDAvsControl", 
                            ylabel = "Read Counts", col = c("#D8A499", "#5B1A18"))


##############
### Figure ###
##############

#Load combined data - combined outside R
df.metalonda.combined = read.csv("Data/MetaLonDA_combined.csv")

# TDA vs. NoTDA #
df.tdavsnotda = df.metalonda.combined[1:4,]

# Make matrix with numeric data
mat.tdavsnotda = as.matrix(df.tdavsnotda[,7:17])

# Change rownames
rownames(mat.tdavsnotda) = df.tdavsnotda$Genera

# Wide to long format
df.l.tdavsnotda = df.tdavsnotda %>% melt(id.vars = c("Genera", "Family", "Order", "Class", "System", "Comparison"), 
                                         variable.name = "Day")

#Set NA to 0.06 (for better plotting)
df.l.tdavsnotda$value[is.na(df.l.tdavsnotda$value)] <- 0.06

# Make levels correct
df.l.tdavsnotda$Genera <- factor(df.l.tdavsnotda$Genera,levels=unique(df.l.tdavsnotda$Genera))

# Figure
plot.metalonda.tdavsnotda = ggplot(df.l.tdavsnotda, aes(x = Day, y = Genera)) + 
  geom_point(aes(size = value, fill = System), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.007, 0.06), range = c(10,1), breaks = c(0.001,0.01,0.02,0.03,0.04,0.06)) + 
  labs( x= "", y = "", size = "p-value", fill = "Dominant system")  + 
  theme_bw(base_size = 8)+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = "black", face = "bold.italic"), 
        legend.text = element_text(face ="bold", colour ="black"), 
        legend.title = element_text(face = "bold.italic"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", size = 1.2), 
        legend.position = "top", legend.box = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#7294D4","#5B1A18")) +
  guides(fill = guide_legend(nrow = 2, ncol = 1,  byrow = TRUE))
plot.metalonda.tdavsnotda

ggsave("metalonda_TDAvsNoTDA.svg",
       plot = plot.metalonda.tdavsnotda,
       units="cm",
       width=90,
       height=90,
       dpi = 300)


# TDA vs. Control #
df.tdavscontrol = df.metalonda.combined[5:9,]

# Wide to long format
df.l.tdavscontrol = df.tdavscontrol %>% melt(id.vars = c("Genera", "Family", "Order", "Class", "System", "Comparison"), 
                                         variable.name = "Day")

#Set NA to 0.06
df.l.tdavscontrol$value[is.na(df.l.tdavscontrol$value)] <- 0.06

# Correct levels
df.l.tdavscontrol$Genera <- factor(df.l.tdavscontrol$Genera,levels=unique(df.l.tdavscontrol$Genera))

# Figure
plot.metalonda.tdavscontrol = ggplot(df.l.tdavscontrol, aes(x = Day, y = Genera)) + 
  geom_point(aes(size = value, fill = System), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.007, 0.06), range = c(10,1), breaks = c(0.001,0.01,0.02,0.03,0.04,0.06)) + 
  labs( x= "", y = "", size = "p-value", fill = "Dominant system")  + 
  theme_bw(base_size = 8)+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = "black", face = "bold.italic"), 
        legend.text = element_text(face ="bold", colour ="black"), 
        legend.title = element_text(face = "bold.italic"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", size = 1.2), 
        legend.position = "",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#D8A499", "#5B1A18")) 
plot.metalonda.tdavscontrol

ggsave("Figures/metalonda_TDAvsControl.png",
       plot = plot.metalonda.tdavscontrol,
       units="cm",
       width=90,
       height=90,
       dpi = 300)


# NoTDA vs. control #

df.notdavscontrol = df.metalonda.combined[10:11,]

# Wide to long format
df.l.notdavscontrol = df.notdavscontrol %>% melt(id.vars = c("Genera", "Family", "Order", "Class", "System", "Comparison"), 
                                             variable.name = "Day")

#Set NA to 0.06
df.l.notdavscontrol$value[is.na(df.l.notdavscontrol$value)] <- 0.06

# Correct levels
df.l.notdavscontrol$Genera <- factor(df.l.notdavscontrol$Genera,levels=unique(df.l.notdavscontrol$Genera))

# Figure
plot.metalonda.notdavscontrol = ggplot(df.l.notdavscontrol, aes(x = Day, y = Genera)) + 
  geom_point(aes(size = value, fill = System), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.007, 0.06), range = c(10,1), breaks = c(0.001,0.01,0.02,0.03,0.04,0.06)) + 
  labs( x= "", y = "", size = "p-value", fill = "Dominant system")  + 
  theme_bw(base_size = 8)+
  theme(axis.text.x = element_text(colour = "black", face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold.italic"), 
        legend.text = element_text(face ="bold", colour ="black"), 
        legend.title = element_text(face = "bold.italic"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", size = 1.2), 
        legend.position = "",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#7294D4", "#D8A499")) 
plot.metalonda.notdavscontrol

ggsave("Figures/metalonda_TDAvsControl.png",
       plot = plot.metalonda.tdavscontrol,
       units="cm",
       width=90,
       height=90,
       dpi = 300)

# Make Panel plot
fig4 = plot_grid(plot.metalonda.tdavsnotda,
             plot.metalonda.tdavscontrol,
             plot.metalonda.notdavscontrol,
             nrow = 3, labels = "AUTO")

ggsave("Fig4.svg",
       plot = fig4,
       units="cm",
       width=18,
       height=18,
       dpi = 300)


# All plotted at the same time #

# Wide to long format
df.l.combined = df.metalonda.combined %>% melt(id.vars = c("Genera", "Family", "Order", "Class", "System", "Comparison"), 
                                                 variable.name = "Day")

#Set NA to 0.06
df.l.combined$value[is.na(df.l.combined$value)] <- 0.06

# Correct levels
df.l.combined$Genera <- factor(df.l.combined$Genera,levels=unique(df.l.combined$Genera))

# Figure
plot.metalonda.combined = ggplot(df.l.combined, aes(x = Day, y = Genera)) + 
  geom_point(aes(size = value, fill = System), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.007, 0.06), range = c(10,1), breaks = c(0.001,0.01,0.02,0.03,0.04,0.06)) + 
  labs( x= "", y = "", size = "p-value", fill = "Dominant system")  + 
  theme_bw(base_size = 8)+
  theme(axis.text.x = element_text(colour = "black", face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold.italic"), 
        legend.text = element_text(face ="bold", colour ="black"), 
        legend.title = element_text(face = "bold.italic"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", size = 1.2), 
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#D8A499","#7294D4", "#5B1A18"))
plot.metalonda.combined

ggsave("legend.svg",
       plot = plot.metalonda.combined,
       units="cm",
       width=18,
       height=9,
       dpi = 300)
