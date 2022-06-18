# make periods nmds

load("Data/ps.asv.reduced.wTree.RData")

# Remove genus Phaeobacter from data
ps.new.NoP.p1 = subset_taxa(ps.new, genus != "Phaeobacter") %>%
  subset_samples(day > 14)

#Normalize with TSS
ps.norm.p1 = transform_sample_counts(ps.new.NoP.p1, function(x) 100000 * x/sum(x))

# Function to make OTU matrix
otu.func <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
# Function to make Metadata dataframe
metadata.func <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# Make OTU matrix and metadata dataframe 
veg.p1= round(otu.func(ps.norm.p1))

# Make metadata
metadata.veg = metadata.func(ps.norm.p1) %>%
  unite("Lineage", exp.condition, bio.rep, sep = "_",  remove = FALSE) #Make lineage variable

# Make an ordination to determine stress for
set.seed(41)
#get distances
bray.p1 = vegdist(wisconsin(sqrt(veg.p1)), method = "bray")


# betadispr 
b = betadisper(bray.p1, factor(paste(metadata.veg$exp.condition,
                                               metadata.veg$day)))
anova(b)
set.seed(41)
permutest = permutest(b, permutations = 1000, pairwise = TRUE)
df.permu.system = permutest$tab

# Adonis test (PERMANOVA - with all variables)
set.seed(41)
mdl = adonis(bray.p1 ~ exp.condition+day, data = metadata.veg, 
             permutations = 1000) 
SSQ <- mdl$aov.tab$SumsOfSqs
tb.anova <- data.frame(source = rownames(mdl$aov.tab), 
                       Df = mdl$aov.tab$Df,
                       SSQrel = 100*SSQ/SSQ[4],
                       R2 = mdl$aov.tab$R2,
                       F.Model = mdl$aov.tab$F.Model,
                       pv = mdl$aov.tab$`Pr(>F)`)

#Pairwise adonis
pwa <- pairwise.adonis2(bray.p1~ exp.condition * day * Lineage, metadata.veg)
pwa.data <- pwa[[3]]


#Make final NMDS with transformed distance matrix
set.seed(41)
mds<- metaMDS(bray.p1,k = 3, 
                        trymax = 1000, 
                        trace = T, autotransform = FALSE) #stress 0.1170633

# Shepard plot
stressplot(NMDS.bray.d3) 

# Collect scores for plotting nMDS
NMDS1 = mds$points[,1]
NMDS2 = mds$points[,2]
NMDS3 = mds$points[,3]
df.NMDS = data.frame(NMDS1 = NMDS1, NMDS2 = NMDS2, NMDS3 = NMDS3)

# Add metadata to scores for plotting
df.NMDS$System = metadata.veg$exp.condition
df.NMDS$Day = as.factor(metadata.veg$time*7)
df.NMDS$BioRep = metadata.veg$bio.rep
df.NMDS$SampleName = metadata.veg$sample.name
df.NMDS$Sample = rownames(metadata.veg)
df.NMDS$Lineage = metadata.veg$Lineage
df.NMDS$System = factor(df.NMDS$System, levels= c("TDA", "NoTDA", "Control"))


#############################
### nMDS plot(s) - System ###
#############################

# Color palettes
system_color2 = c("#9C3206", "#00A091", "black")
#line.colors = c("#5B1A18","#D8A499","#f7f7f7" ,"#7294D4")

# Plot: Color by SYSTEM (d1-2)
ggplot(df.NMDS, aes(x = NMDS2, y = NMDS3)) + 
  geom_point(size = 2, shape = 21, aes(fill = as.numeric(Day)), alpha = 0.9, stroke = 0)+ 
  theme_bw(base_size = 8)+
  facet_grid(System~.)+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(face ="bold", colour ="black"),
        legend.position = "top", axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))+
  #stat_ellipse(aes(group=Day), size = 1, alpha = 0.5)+
  # scale_fill_manual(values = system_color2)  +
  # scale_color_manual(values = system_color2)+
  labs(x = "NMDS1\n", colour = "System", y = "\n NMDS2")

