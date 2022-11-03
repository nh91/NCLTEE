
######################################################################

### Ordination using nMDS and Bray-Curtis distance (VEGAN edition) ###

######################################################################

#################
### Load data ###
#################

load("Data/ReadO/ps.asv.reduced.wTree.RData")

# Remove genus Phaeobacter from data
ps.new.NoP = subset_taxa(ps.new, genus != "Phaeobacter")

#Normalize with TSS
ps.norm = transform_sample_counts(ps.new.NoP, function(x) 100000 * x/sum(x))

###############################################
### Create VEGAN dataset from physeq object ###
###############################################

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
otu.veg.genus.reduced = round(otu.func(ps.norm))

# Make metadata
metadata.veg = metadata.func(ps.norm) %>%
  unite("Lineage", exp.condition, bio.rep, sep = "_",  remove = FALSE) #Make lineage variable


##################
### Scree plot ###
##################

# Pick dimensions to calculate stress for (i.e. 10) 
k.dim = 1:10 

# Make empty vector to save stress values in
stress = numeric(length(k.dim))

# Make an ordination to determine stress for
set.seed(41)
ord.veg.bray = metaMDS(otu.veg.genus.reduced) # Check what transformation metaMDS makes - makes squareroot and wisconsin.

#Get distances with the transfomations suggested above
bray.distances = vegdist(wisconsin(sqrt(otu.veg.genus.reduced)), method = "bray")

# Function to calculate stress 
set.seed(41) # needs distance in loop 
for(i in seq_along(k.dim)) {
  solutions <- metaMDSiter(bray.distances, k = i,
                     trace = FALSE)
  stress[i] <- solutions$stress
}

# Scree dataset
scree = as.data.frame(cbind(k.dim, stress))

# Scree plot
screeplot = ggplot(scree, aes( x = as.factor(k.dim), y = stress))+ 
  geom_bar(stat="identity", fill = "#5B1A18")+
  xlab("\n Dimensions")+
  ylab("Stress\n")+
  theme_bw(base_size = 12)+
  theme(axis.text.y = element_text(colour = "black", face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold"), 
        axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold", colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))
screeplot

ggsave("Fig3/S3_scree.pdf",
       plot = screeplot,
       units="cm",
       width=10,
       height=8,
       dpi = 300)



##################
### Beta dispr ###
##################

# Make betadispersion with time and exp.condition
betadisp.main.system = betadisper(bray.distances,
                                  factor(paste(metadata.veg$exp.condition,
                                               metadata.veg$day)))
anova(betadisp.main.system)
set.seed(41)
permutest1 = permutest(betadisp.main.system, permutations = 1000, pairwise = TRUE)


# Without time
betadisp.main.system2 = betadisper(bray.distances,
           factor(metadata.veg$exp.condition))
anova(betadisp.main.system2)
set.seed(41)
permutest = permutest(betadisp.main.system2, permutations = 1000, pairwise = TRUE)
## Tukey's Honest Significant Differences
betadisp.main.system.HSD <- TukeyHSD(betadisp.main.system2)
plot(betadisp.main.system.HSD)
dim(betadisp.main.system.HSD$group)


#Lineages witihin systems
betadisp.main.system3 = betadisper(bray.distances,
                                  factor(paste(metadata.veg$exp.condition,
                                               metadata.veg$Lineage)))
anova(betadisp.main.system3)
set.seed(41)
permutest = permutest(betadisp.main.system3, permutations = 1000, pairwise = TRUE) # No significant difference


df.permu.system = permutest$tab
df.permu.system$variable = "System"
#write.csv(df.permu.system, file = "Data/betadisper.mainNMDS.csv")



#################
### PERMANOVA ###
#################

# Adonis test (PERMANOVA - with all variables)
set.seed(41)
adonis2(bray.distances ~ exp.condition+time*Lineage, data = metadata.veg, 
              permutations = 1000) 

metadata.veg$day = as.factor(metadata.veg$day)

#Pairwise adonis
pairwiseADONIS <- pairwise.adonis2(bray.distances ~ exp.condition * day * Lineage, metadata.veg)
#write.table(pairwiseADONIS, file = "Fig3/pwa.data.csv", sep = "\t", dec = ",")

##################
### Ordination ###
##################

#Make final NMDS with transformed distance matrix
set.seed(41)
NMDS.bray.d3 <- metaMDS(bray.distances, 
                        k = 3, 
                        trymax = 1000, 
                        trace = T, autotransform = FALSE) #stress 0.1170633

# Shepard plot
shepardplot = stressplot(NMDS.bray.d3)  

# Collect scores for plotting nMDS
NMDS1 = NMDS.bray.d3$points[,1]
NMDS2 = NMDS.bray.d3$points[,2]
NMDS3 = NMDS.bray.d3$points[,3]
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
nmds.system.12 = ggplot(df.NMDS, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 2, color = "black", shape = 21, aes(fill = System), alpha = 0.9)+ 
  theme_bw(base_size = 10)+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(colour ="black"),
        legend.position = "none", axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))+
  stat_ellipse(aes(group=System, color = System), size = 1, alpha = 0.6)+
  scale_color_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  scale_fill_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  labs(x = "NMDS1\n", colour = "System", y = "\n NMDS2")
nmds.system.12

##### Plot: Color by SYSTEM (d1-3) #####
nmds.system.13 = ggplot(df.NMDS, aes(x = NMDS1, y = NMDS3)) + 
  geom_point(size = 2, color = "black", shape = 21, aes(fill = System), alpha = 0.9)+ 
  theme_bw(base_size = 12)+
  #facet_grid(.~System)+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(face ="bold", colour ="black"),
        legend.position = "none", axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.key=element_blank())+
  scale_fill_manual(values = system_color2)  +
  stat_ellipse(aes(group=System, color = System), size = 1, alpha = 0.6)+
  scale_color_manual(values = system_color2)+
  labs(x = "NMDS1\n", colour = "System", y = "\n NMDS3")

#### Plot: Color by SYSTEM (d2-3)####
nmds.system.23 = ggplot(df.NMDS, aes(x = NMDS2, y = NMDS3)) + 
  geom_point(size = 2, shape = 21, color = "black", aes(fill = System), alpha = 0.9)+ 
  theme_bw(base_size = 12)+
  #facet_grid(.~System)+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(face ="bold", colour ="black"),
        legend.position = "none", axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.key=element_blank())+
  stat_ellipse(aes(group=System, color = System), size = 1, alpha = 0.6)+
  scale_fill_manual(values = system_color2)  +
  scale_color_manual(values = system_color2)+
  labs(x = "NMDS2\n", colour = "System", y = "\n NMDS3")


###########################
### nMDS plot(s) - Time ###
###########################

# Plot: Color by TIME (d1-2)

nmds.time.12 = ggplot(df.NMDS, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 2, shape = 21,  colour = "black", aes(fill = Day), alpha = 0.9)+ 
  theme_bw(base_size = 10)+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(face ="bold", colour ="black"),
        legend.position = "right",
        axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))+
  #stat_ellipse(aes(group=Day), size = 1, alpha = 0.6)+
 scale_fill_brewer(palette = "PiYG", guide = "colourbar")+
  labs(x = "NMDS1\n", colour = "Day", y = "\n NMDS2")+
  guides(fill = guide_legend(ncol = 2, byrow = TRUE))
nmds.time.12

#### Plot: Color by TIME (d1-3) ####
nmds.time.13 = ggplot(df.NMDS, aes(x = NMDS1, y = NMDS3)) + 
  geom_point(size = 2, shape = 21, colour = "black", aes(fill =Day), alpha = 0.9)+ 
  theme_bw(base_size = 12)+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(face ="bold", colour ="black"),
        legend.position = "none", axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.key=element_blank())+
  scale_fill_brewer(palette = "PiYG")+
  labs(x = "NMDS1\n", colour = "Day", y = "\n NMDS3")

##### Plot: Color by TIME (d2-3) ####
nmds.time.23 = ggplot(df.NMDS, aes(x = NMDS2, y = NMDS3)) + 
  geom_point(size = 2, shape = 21, colour = "black", aes(fill = Day), alpha = 0.9)+ 
  theme_bw(base_size = 12)+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(face ="bold", colour ="black"),
        legend.position = "none", axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.key=element_blank())+
  scale_fill_brewer(palette = "PiYG")+
  labs(x = "NMDS2\n", colour = "Day", y = "\n NMDS3")

plot.main.nmds = plot_grid(nmds.time.12,nmds.time.13,nmds.time.23,
                           nmds.system.12, nmds.system.13, nmds.system.23,
                           ncol =3, nrow = 2, labels= c("A", "", "", "B", "", ""))
plot.main.nmds

######
Fig3BC = plot_grid(nmds.time.12, nmds.system.12, ncol =2, nrow = 1)

legend.time<- as_ggplot(ggpubr::get_legend(nmds.time.12))
legend.system<- as_ggplot(ggpubr::get_legend(nmds.system.12))


ggsave("Fig3/Fig3BC_nMDS.svg",
       plot = Fig3BC,
       units="cm",
       width=15,
       height= 7,
       dpi = 300)

ggsave("Fig3/legend.time.svg",
       plot = legend.time,
       units="cm",
       width=4,
       height= 6,
       dpi = 300)

ggsave("Fig3/legend.system.svg",
       plot = legend.system,
       units="cm",
       width=4,
       height= 3,
       dpi = 300)


##################################
### nMDS plot - individual day ###
##################################

# Plot for individual days (insert day)
# nmds.system.day14= ggplot(df.NMDS, aes(x = NMDS1, y = NMDS2)) + 
#   geom_point(size = 2, shape = 21, aes(fill = System), alpha = 0.9)+ 
#   theme_bw(base_size = 10)+
#   facet_wrap(vars(Day), nrow =4)+
#   theme(axis.text.y = element_text(colour = "black", face = "bold"), 
#         axis.text.x = element_text(colour = "black", face = "bold"), 
#         legend.text = element_text(face ="bold", colour ="black"),
#         legend.position = "none", axis.title.y = element_text(face = "bold"), 
#         plot.title = element_text(face = "bold", hjust = 0.5),
#         axis.title.x = element_text(face = "bold", colour = "black"), 
#         legend.title = element_text(colour = "black", face = "bold"), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1),
#         legend.key=element_blank())+
#   stat_ellipse(aes(group=System, color = System), size = 1, alpha = 0.5)+
#   scale_fill_manual(values = system_color2)  +
#   scale_color_manual(values = system_color2)+
#   labs(x = "NMDS1", title = "Day 14", colour = "System", y = "NMDS2")
# nmds.system.day14

##################################
### Envfit w. taxa to lineage ####
##################################

# Fit environmental factors/Taxa in the ordination 
set.seed(41)
envfit.genera = envfit(NMDS.bray.d3, 
                       env = data.frame(t(otu_table(ps.norm))), 
                       permutations = 1000,
                       choices = c(1,2,3))

# Look at R values for all asvs, (Goodness of fit /squared correlation coeff)
plot(sort(envfit.genera$vectors$r, decreasing = TRUE))
 
# Pick out only vectors with an r value above 0.5
sigIndx=which(envfit.genera$vectors$r>0.5 & envfit.genera$vectors$pvals<0.05)

# Find vector coordinates
test = envfit.genera$vectors$arrows[sigIndx,]
tesp = as.data.frame(scores(envfit.genera, "vectors"))[sigIndx,] * 0.6
#en_coord_cont = as.data.frame(scores(ASVfit, "vectors")) * ordiArrowMul(ASVfit) # ordiArrowMul gives the factor for meaningfull scaling
colnames(tesp) =c("NMDS1", "NMDS2", "NMDS3")
 
# Add taxa names
tax_forfit = as.data.frame(tax_table(ps.norm))
tax_forfit2 <- as.data.frame(tax_forfit[row.names(tax_forfit) %in% row.names(tesp), ])
rownames_genera = rownames(tesp)
tes = cbind(tesp, tax_forfit2, rownames_genera) %>%
  unite("ASV_genus", c("genus", "rownames_genera"), sep ="_")
 
#rownames_loadings <- with(tes, paste(family, ASV_genus, sep = "_"))
rownames_loadings <- with(tes, ASV_genus)
row.names(tesp) <- make.unique((rownames_loadings))
 
# nmds with tax (d1-2)
nmds.time.12.tax = ggplot(df.NMDS, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 2, shape = 21, colour = "black", aes(fill = Day), alpha = 0.9)+
  theme_bw(base_size = 10)+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(face ="bold", colour ="black"),
        legend.position = "none", axis.title.y = element_text(face = "bold"), 
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.key=element_blank())+
  scale_fill_brewer(palette = "PiYG")+
  xlim(-0.60, 0.5)+
  labs(x = "NMDS1\n", y = "\n NMDS2")+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = tesp, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = tesp, aes(x = NMDS1, y = NMDS2), colour = "grey10", 
            fontface = "italic", size = 2, label = row.names(tesp))
nmds.time.12.tax

# ggsave("nmds.time.12.tax.svg",
#        plot = nmds.time.12.tax,
#        units="cm",
#        width=18,
#        height= 15,
#        dpi = 300)

# Get p-values and r for plotted taxa

# Select data
sigIndx2=as.data.frame(which(envfit.genera$vectors$r>0.5 & envfit.genera$vectors$pvals<0.05))
r = as.data.frame(envfit.genera$vectors$r)
pvals = as.data.frame(envfit.genera$vectors$pvals)
r2 = as.data.frame(r[row.names(r) %in% row.names(sigIndx2),])
pvals2 = as.data.frame(pvals[row.names(pvals) %in% row.names(sigIndx2),])

#Combine data and clean table
df.tax.scores = cbind(sigIndx2,tax_forfit2, r2, pvals2)
colnames(df.tax.scores) = c("x", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "r", "pvals")
table.tax.scores = df.tax.scores[,-c(1,8)]

#Write csv file
write.csv(table.tax.scores, file = "Fig3/df.scores.envfit.csv")


# # Make panel plot of main (supp figure)
# plot.supp.nmds = plot_grid(nmds.time.12.tax, nmds.time.13.tax, nmds.time.23.tax,
#                            ncol = 3, nrow = 1)

##############################
### nMDS plot(s) - Lineage ###
##############################

# Make path of time for plotting lines
path<-df.NMDS
path$Day<-as.numeric(as.character(path$Day))
path<-path[order(path$Day),]

# TDA LINEAGE
TDAlin = c("#B85C48", "#80221E", "#AD7C59", "#4D392E")
nmds.lineage.12.TDA = ggplot(df.NMDS[df.NMDS$System == "TDA",], 
                             aes(x = NMDS1, y = NMDS2)) +
  scale_shape_manual(values = c(21, 22, 23, 24))+
  geom_path(data=path[path$System == "TDA",],
            aes(x=NMDS1,y=NMDS2,group= Lineage, color = Lineage),size=1, alpha = 0.7, show.legend = FALSE)+
  geom_text(data=path[path$System == "TDA",],
            aes(x=NMDS1,y=NMDS2,label=Day,vjust=-0.7, color = Lineage), size = 2, fontface = "bold", show.legend = FALSE) +
  geom_point(size = 2, colour = "black", aes(shape = Lineage, fill = Lineage), alpha = 0.9)+
  theme_bw(base_size = 10)+
  theme(axis.text.y = element_text(colour = "black", face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold"),
        legend.text = element_text(face ="bold", colour ="black"), 
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", colour = "black"),
        legend.title = element_text(colour = "black", face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.key=element_blank(), 
        legend.position = "right", 
        plot.title = element_text(hjust = 0, size = 10, face = "bold"))+
  scale_fill_manual(values = TDAlin)+
  scale_color_manual(values = TDAlin)+
  xlim(-0.65, 0.5)+
  ylim(-0.5, 0.5)+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = tesp, size =0.5, alpha = 0.7, colour = "black", linetype = "dashed") +
  geom_text(data = tesp, aes(x = NMDS1, y = NMDS2), colour = "grey10",
            fontface = "bold.italic", size = 3, label = row.names(tesp))+
  labs(x = "NMDS1\n", y = "\n NMDS2")
nmds.lineage.12.TDA

legend.envfit.lin<- as_ggplot(ggpubr::get_legend(nmds.lineage.12.TDA))

ggsave("Fig3/legend.envfit.lin.svg",
       plot = legend.envfit.lin,
       units="cm",
       width=2,
       height= 3,
       dpi = 300)

ggsave("Fig3/nmds.lineage.12.TDA.svg",
       plot = nmds.lineage.12.TDA,
       units="cm",
       width=8,
       height= 8,
       dpi = 300)

#NoTDA
nTDAlin = c("#8EBCB1", "#549895", "#4A756E", "#245254")
nmds.lineage.12.NoTDA = ggplot(df.NMDS[df.NMDS$System == "NoTDA",], 
                               aes(x = NMDS1, y = NMDS2)) +
  scale_shape_manual(values = c(21, 22, 23, 24))+
  geom_path(data=path[path$System == "NoTDA",],
            aes(x=NMDS1,y=NMDS2,group= Lineage, color = Lineage),size=1, alpha = 0.7, show.legend = FALSE)+
  geom_text(data=path[path$System == "NoTDA",],
            aes(x=NMDS1,y=NMDS2,label=Day,vjust=-0.7, color = Lineage), size=2, fontface = "bold", show.legend = FALSE) +
  geom_point(size = 2, colour = "black", aes(shape = Lineage, fill = Lineage), alpha = 0.9)+
  theme_bw(base_size = 10)+
  theme(axis.text.y = element_text(colour = "black", face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold"),
        legend.text = element_text(face ="bold", colour ="black"), 
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", colour = "black"),
        legend.title = element_text(colour = "black", face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.key=element_blank(), 
        legend.position = "right",
        plot.title = element_text(hjust = 0, size = 10, face = "bold"))+
  scale_fill_manual(values = nTDAlin)+
  scale_color_manual(values = nTDAlin)+
  xlim(-0.65, 0.5)+
  ylim(-0.5, 0.5)+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = tesp, size =0.5, alpha = 0.7, colour = "black", linetype = "dashed") +
  geom_text(data = tesp, aes(x = NMDS1, y = NMDS2), colour = "grey10",
            fontface = "bold.italic", size = 3, label = row.names(tesp))+
  labs(x = "NMDS1\n", y = "\n NMDS2")
nmds.lineage.12.NoTDA

#Control
cTDAlin = c("#8C8C8C", "#555555", "#333333", "#000000")
nmds.lineage.12.Control = ggplot(df.NMDS[df.NMDS$System == "Control",], 
                                 aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 2, colour = "black", aes(shape = Lineage, fill = Lineage), alpha = 0.9)+
  scale_shape_manual(values = c(21, 22, 23, 24))+
  geom_path(data=path[path$System == "Control",],
            aes(x=NMDS1,y=NMDS2,group= Lineage, color = Lineage),size=1, alpha = 0.7, show.legend = FALSE)+
  geom_text(data=path[path$System == "Control",],
            aes(x=NMDS1,y=NMDS2,label=Day,vjust=-0.7, color = Lineage), size =2, fontface = "bold", show.legend = FALSE) +
  theme_bw(base_size = 10)+
  theme(axis.text.y = element_text(colour = "black", face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold"),
        legend.text = element_text(face ="bold", colour ="black"), 
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", colour = "black"),
        legend.title = element_text(colour = "black", face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.key=element_blank(), 
        legend.position = "right",
        plot.title = element_text(hjust = 0, size = 10, face = "bold"))+
  scale_fill_manual(values = cTDAlin)+
  scale_color_manual(values = cTDAlin)+
  xlim(-0.65, 0.5)+
  ylim(-0.5, 0.5)+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = tesp, size =0.5, alpha = 0.7, colour = "black", linetype = "dashed") +
  geom_text(data = tesp, aes(x = NMDS1, y = NMDS2), colour = "grey10",
            fontface = "bold.italic", size = 3, label = row.names(tesp))+
  labs(x = "NMDS1\n", y = "\n NMDS2")
nmds.lineage.12.Control


plot.supp.nmds.lin = plot_grid(nmds.lineage.12.TDA, nmds.lineage.12.NoTDA, nmds.lineage.12.Control,
                           ncol = 1, nrow = 3)

ggsave("Fig3/plot.supp.nmds.lin.pdf",
       plot = plot.supp.nmds.lin,
       units="cm",
       width=16,
       height= 33,
       dpi = 300)

ggsave("Fig3/plot.supp.nmds.lin.svg",
       plot = plot.supp.nmds.lin,
       units="cm",
       width=16,
       height= 33,
       dpi = 300)



##### Bray-curtis distances comparison within systems with between system

hist(bray.distances)

#Get distances
distmat <- data.frame(as.matrix(bray.distances))
distmat$SampleID = rownames(distmat)
metadata.veg$SampleID = rownames(metadata.veg)
brayC.df = merge(distmat,metadata.veg, by = "SampleID")

#Take out and make seperate df
write_csv(brayC.df, file = "brayC.df.csv")


#Control
test = as.data.frame(readxl::read_excel("brayC.df.xlsx", sheet = 7)[,1:42])
rownames(test) = test$SampleID
test.control = as.matrix(test[,-1])

metadata.con = as.data.frame(readxl::read_excel("brayC.df.xlsx", sheet = 7)[,c(1,43:45)])

## Get df 
distmat <- test.control
sub_dist <- list()
groups_all1 <- metadata.con$Lineage
groups_all <- as.factor(metadata.con$Lineage)


for (group in levels(groups_all)) { 
  row_group <- which(groups_all == group)
  sample_group <- metadata.con$SampleID[row_group]
  sub_dist[[group]] <- distmat[ sample_group, sample_group]
  sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}

braygroups <- melt(sub_dist)
df.bray.con <- braygroups[complete.cases(braygroups), ]
df.bray.con$L1 <- factor(df.bray.con$L1, levels=names(sub_dist))
head(df.bray.con)
df.bray.con$Condition = "Control"

my_comparisons1 <- list( c("Control_1", "Control_2"), 
                         c("Control_1", "Control_3"), 
                         c("Control_1", "Control_4"),
                         c("Control_2", "Control_3"),
                         c("Control_2", "Control_4"),
                         c("Control_3", "Control_4"))

ggplot(df.bray, aes(x=L1, y=value, colour=L1)) +
  geom_jitter() + 
  geom_boxplot(alpha=0.6) +  
  theme(legend.position="none") +
  ylab("Bray-Curtis diversity") +
  stat_compare_means(comparisons = my_comparisons1)+
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), axis.text.y=element_text(size=12))


#TDA
test = as.data.frame(readxl::read_excel("brayC.df.xlsx", sheet = 5)[,1:44])
rownames(test) = test$SampleID
test.tda = as.matrix(test[,-1])

metadata.tda = as.data.frame(readxl::read_excel("brayC.df.xlsx", sheet = 5)[,c(1,45:47)])

## get df
distmat <- test.tda
sub_dist <- list()
groups_all <- as.factor(metadata.tda$Lineage)


for (group in levels(groups_all)) { 
  row_group <- which(groups_all == group)
  sample_group <- metadata.tda$SampleID[row_group]
  sub_dist[[group]] <- distmat[ sample_group, sample_group]
  sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}

braygroups <- melt(sub_dist)
df.bray.tda <- braygroups[complete.cases(braygroups), ]
df.bray.tda$L1 <- factor(df.bray.tda$L1, levels=names(sub_dist))
head(df.bray.tda)
df.bray.tda$Condition = "TDA"

my_comparisons2 <- list( c("TDA_1", "TDA_2"), 
                        c("TDA_1", "TDA_3"), 
                        c("TDA_1", "TDA_4"),
                        c("TDA_2", "TDA_3"),
                        c("TDA_2", "TDA_4"),
                        c("TDA_3", "TDA_4"))

ggplot(df.bray.tda, aes(x=L1, y=value, colour=L1)) +
  geom_jitter() + 
  geom_boxplot(alpha=0.6) +  
  theme(legend.position="none") +
  ylab("Bray-Curtis dissimilarity") +
  stat_compare_means(comparisons = my_comparisons2)+
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), axis.text.y=element_text(size=12))


#NoTDA
test = as.data.frame(readxl::read_excel("brayC.df.xlsx", sheet = 6)[,1:43])
rownames(test) = test$SampleID
test.notda = as.matrix(test[,-1])

metadata.notda = as.data.frame(readxl::read_excel("brayC.df.xlsx", sheet = 6)[,c(1,44:46)])

## get df
distmat <- test.notda
sub_dist <- list()
groups_all <- as.factor(metadata.notda$Lineage)


for (group in levels(groups_all)) { 
  row_group <- which(groups_all == group)
  sample_group <- metadata.notda$SampleID[row_group]
  sub_dist[[group]] <- distmat[ sample_group, sample_group]
  sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
}

braygroups <- melt(sub_dist)
df.bray.notda <- braygroups[complete.cases(braygroups), ]
df.bray.notda$L1 <- factor(df.bray.notda$L1, levels=names(sub_dist))
head(df.bray.notda)
df.bray.notda$Condition = "noTDA"

my_comparisons3 <- list( c("NoTDA_1", "NoTDA_2"), 
                         c("NoTDA_1", "NoTDA_3"), 
                         c("NoTDA_1", "NoTDA_4"),
                         c("NoTDA_2", "NoTDA_3"),
                         c("NoTDA_2", "NoTDA_4"),
                         c("NoTDA_3", "NoTDA_4"))

ggplot(df.bray.notda, aes(x=L1, y=value, colour=L1)) +
  geom_jitter() + 
  geom_boxplot(alpha=0.6) +  
  theme(legend.position="none") +
  ylab("Bray-Curtis diversity") +
  stat_compare_means(comparisons = my_comparisons3)+
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), axis.text.y=element_text(size=12))


#Combine all
testall = rbind(df.bray.notda, df.bray.tda, df.bray.con)

my_comparisons <- list( c("TDA", "noTDA"), c("TDA", "Control"), c("Control", "noTDA") )

ggplot(testall, aes(x=Condition, y=value, colour=L1)) +
  geom_jitter(show.legend = FALSE) + 
  geom_boxplot(alpha=0.6, show.legend = FALSE) +  
  theme(legend.position="none") +
  ylab("Bray-Curtis dissimilarity") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value   # Add global p-value
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12), axis.text.y=element_text(size=12))

compare_means(value ~ Condition,  data = testall, p.adjust.method = "BH")
compare_means(value ~ Condition,  data = testall, method = "anova")
