
######################################################################

### Ordination using nMDS and Bray-Curtis distance (VEGAN edition) ###

######################################################################

#################
### Load data ###
#################

load("Data/ps.asv.reduced.wTree.RData")

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

#Get distances
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
  theme_bw(base_size = 10)+
  theme(axis.text.y = element_text(colour = "black", face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold"), 
        axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold", colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2))
screeplot

ggsave("scree.pdf",
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
permutest = permutest(betadisp.main.system, permutations = 1000, pairwise = TRUE)

df.permu.system = permutest$tab
#df.permu.system$variable = "System"
write.csv(df.permu.system, file = "Data/betadisper.mainNMDS.csv")

## Tukey's Honest Significant Differences
(betadisp.main.system.HSD <- TukeyHSD(betadisp.main.system))
plot(betadisp.main.system.HSD)
dim(betadisp.main.system.HSD$group)


#################
### PERMANOVA ###
#################

# Adonis test (PERMANOVA - with all variables)
set.seed(41)
mdl = adonis(bray.distances ~ exp.condition+day, data = metadata.veg, 
              permutations = 1000) 
SSQ <- mdl$aov.tab$SumsOfSqs
tb.anova <- data.frame(source = rownames(mdl$aov.tab), 
                       Df = mdl$aov.tab$Df,
                       SSQrel = 100*SSQ/SSQ[4],
                       R2 = mdl$aov.tab$R2,
                       F.Model = mdl$aov.tab$F.Model,
                       pv = mdl$aov.tab$`Pr(>F)`)
write.csv(tb.anova, file = "Data/adonis.mainNMDS.csv")

#Pairwise adonis
pwa <- pairwise.adonis2(bray.distances ~ exp.condition * day * Lineage, metadata.veg)
pwa.data <- pwa[[3]]
write.table(pwa.data, file = "pwa.data.csv", sep = "\t", dec = ",")

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
  #facet_grid(.~System)+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(face ="bold", colour ="black"),
        legend.position = "none", axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))+
  stat_ellipse(aes(group=System, color = System), size = 1, alpha = 0.6)+
  scale_fill_manual(values = system_color2)  +
  scale_color_manual(values = system_color2)+
  labs(x = "NMDS1\n", colour = "System", y = "\n NMDS2")
nmds.system.12

# Plot: Color by SYSTEM (d1-3)
nmds.system.13 = ggplot(df.NMDS, aes(x = NMDS1, y = NMDS3)) + 
  geom_point(size = 2, color = "black", shape = 21, aes(fill = System), alpha = 0.9)+ 
  theme_bw(base_size = 10)+
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

# Plot: Color by SYSTEM (d2-3)
nmds.system.23 = ggplot(df.NMDS, aes(x = NMDS2, y = NMDS3)) + 
  geom_point(size = 2, shape = 21, color = "black", aes(fill = System), alpha = 0.9)+ 
  theme_bw(base_size = 10)+
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
        legend.position = "none", legend.box = "horizontal", axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))+
  #stat_ellipse(aes(group=Day), size = 1, alpha = 0.6)+
  scale_fill_brewer(palette = "PiYG")+
  labs(x = "NMDS1\n", colour = "Day", y = "\n NMDS2") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
nmds.time.12

# Plot: Color by TIME (d1-3)
nmds.time.13 = ggplot(df.NMDS, aes(x = NMDS1, y = NMDS3)) + 
  geom_point(size = 2, shape = 21, colour = "black", aes(fill =Day), alpha = 0.9)+ 
  theme_bw(base_size = 10)+
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

# Plot: Color by TIME (d2-3)
nmds.time.23 = ggplot(df.NMDS, aes(x = NMDS2, y = NMDS3)) + 
  geom_point(size = 2, shape = 21, colour = "black", aes(fill = Day), alpha = 0.9)+ 
  theme_bw(base_size = 10)+
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


ggsave("Fig3.pdf",
       plot = plot.main.nmds,
       units="cm",
       width=24,
       height= 15,
       dpi = 300)


##############################
### nMDS plot(s) - Lineage ###
##############################

# Make path of time for plotting lines
path<-df.NMDS
path$Day<-as.numeric(as.character(path$Day))
path<-path[order(path$Day),]

# Plot: Color by LINEAGE (D1-2)
nmds.lineage.12 = ggplot(df.NMDS, aes(x = NMDS1, y = NMDS2, color = System)) +
  geom_point(size = 2, colour = "black", aes(shape = as.factor(BioRep), fill = System), alpha = 0.9)+
  scale_shape_manual(values = c(21, 22, 23, 24))+
  geom_path(data=path,aes(x=NMDS1,y=NMDS2,group=Lineage,colour=System),size=1, alpha = 0.7)+
  geom_text(data=path,aes(x=NMDS1,y=NMDS2,label=Day,vjust=-0.7)) +
  theme_bw(base_size = 12)+
  facet_grid(~System)+
  theme(axis.text.y = element_text(colour = "black", face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold"),
        legend.text = element_text(face ="bold", colour ="black"), axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", colour = "black"),
        legend.title = element_text(colour = "black", face = "bold"),
        strip.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        legend.key=element_blank(), 
        legend.position = "none",
        plot.title = element_text(hjust = 0, size = 10, face = "bold"))+
  scale_fill_manual(values = system_color2)  +
  scale_color_manual(values = system_color2)+
  labs(x = "NMDS1\n", y = "\n NMDS2")
nmds.lineage.12

ggsave("nmds.lineage.12.pdf",
       plot = nmds.lineage.12,
       units="cm",
       width=37,
       height= 13,
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

#######################
### Envfit w. taxa ####
#######################

# Fit environmental factors/Taxa in the ordination 
set.seed(41)
envfit.genera = envfit(NMDS.bray.d3, 
                       env = data.frame(t(otu_table(ps.norm))), 
                       permutations = 1000,
                       choices = c(1,2,3))

# Look at R values for all asvs, (Goodness of fit /squared correlation coeff)
plot(sort(envfit.genera$vectors$r, decreasing = TRUE))
 
# Pick out only vectors with an r value above 0.2
sigIndx=which(envfit.genera$vectors$r>0.55)

# Find vector coordinates
test = envfit.genera$vectors$arrows[sigIndx,]
tesp = as.data.frame(scores(envfit.genera, "vectors"))[sigIndx,] * 0.55
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
  xlim(-0.50, 0.3)+
  labs(x = "NMDS1\n", y = "\n NMDS2")+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = tesp, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = tesp, aes(x = NMDS1, y = NMDS2), colour = "grey10", 
            fontface = "italic", size = 3, label = row.names(tesp))
nmds.time.12.tax

ggsave("nmds.time.12.tax.svg",
       plot = nmds.time.12.tax,
       units="cm",
       width=18,
       height= 15,
       dpi = 300)

# nmds with tax (d1-3)
nmds.time.13.tax = ggplot(df.NMDS, aes(x = NMDS1, y = NMDS3)) +
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
  xlim(-0.50, 0.5)+
  labs(x = "NMDS1\n", y = "\n NMDS3")+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS3), 
               data = tesp, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = tesp, aes(x = NMDS1, y = NMDS3), colour = "grey10", 
            fontface = "italic", size = 3, label = row.names(tesp))
nmds.time.13.tax

# nmds with tax (d2-3)
nmds.time.23.tax = ggplot(df.NMDS, aes(x = NMDS2, y = NMDS3)) +
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
  xlim(-0.50, 0.5)+
  labs(x = "NMDS2\n", y = "\n NMDS3")+
  geom_segment(aes(x = 0, y = 0, xend = NMDS2, yend = NMDS3), 
               data = tesp, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = tesp, aes(x = NMDS2, y = NMDS3), colour = "grey10", 
            fontface = "italic", size = 3, label = row.names(tesp))
nmds.time.23.tax

# Get p-values and r for plotted taxa

# Select data
sigIndx2=as.data.frame(which(envfit.genera$vectors$r>0.2 & envfit.genera$vectors$pvals<0.05))
r = as.data.frame(envfit.genera$vectors$r)
pvals = as.data.frame(envfit.genera$vectors$pvals)
r2 = as.data.frame(r[row.names(r) %in% row.names(sigIndx2),])
pvals2 = as.data.frame(pvals[row.names(pvals) %in% row.names(sigIndx2),])

#Combine data and clean table
df.tax.scores = cbind(sigIndx2,tax_forfit2, r2, pvals2)
colnames(df.tax.scores) = c("x", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "r", "pvals")
table.tax.scores = df.tax.scores[,-c(1,8)]

#Write csv file
write.csv(table.tax.scores, file = "df.scores.envfit.csv")


# Make panel plot of main (supp figure)
plot.supp.nmds = plot_grid(nmds.time.12.tax, nmds.time.13.tax, nmds.time.23.tax,
                           ncol = 3, nrow = 1)

# ggsave("nmds_supp.svg",
#        plot = plot.supp.nmds,
#        units="cm",
#        width=35,
#        height=10,
#        dpi = 300)
