#trajectories of each ASV

load("Data/ReadO/ps.asv.reduced.wTree.RData")

###################################
### Get dataframes for each ASV ###
###################################

# dataframe
df.genus = ps.new %>% 
  phyloseq::filter_taxa(function (x) {sum(x > 0) > 4}, prune=TRUE) %>%
  phyloseq::transform_sample_counts(function(x) {x/sum(x)}*100 ) %>%
  psmelt() %>%
  dplyr::rename(System = exp.condition)

df.genus$OTU = as.factor(df.genus$OTU)

#Summarize

sum.genus = df.genus %>%
  group_by(day, System, OTU, genus) %>%
  summarise(mean.abundance = mean(Abundance),
            sd.abundance = sd(Abundance),
            count = n())
  
############
### Plot ###
############


#Overview plot of ASVs abundance

sum.genus$System = factor(sum.genus$System , levels= c("TDA", "NoTDA", "Control"))
system_color2 = c("#9C3206", "#00A091", "black")
labels<- c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control")
names(labels) = c("TDA", "NoTDA", "Control")
sum.genus = unite(sum.genus, "TaxID", genus, OTU, sep = "_",  remove = FALSE) #Make lineage variable

Fig_trajectories = ggplot(sum.genus, 
         aes(x = day, y = mean.abundance,fill=System, color = System, shape = System)) +
  geom_line(size = 1)+
  geom_point(size = 2)+
  geom_ribbon(aes(ymax = mean.abundance + sd.abundance, 
                  ymin = mean.abundance - sd.abundance, fill = System), alpha = 0.2, color = NA) +
  facet_wrap(.~TaxID, ncol = 5, nrow = 10, scales = "free_y") + 
  theme_bw(base_size = 12) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom",
        strip.text = element_text(face = "bold", colour ="black"), 
        legend.title = element_text(face = "bold", colour ="black"),
        legend.text = element_text(face = "bold", colour ="black"),
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_color_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  scale_fill_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  labs(x = "\n Time (days)", y = " Relative abundance (%)\n") + 
  scale_x_continuous(breaks = c(0,7,14,21,28,35,42,49,56,63,70))
Fig_trajectories

ggsave("FigS6_Traject.pdf",
       plot = Fig_trajectories,
       units="cm",
       width=36,
       height=46,
       dpi = 300)

legend.trajectory<- as_ggplot(ggpubr::get_legend(Fig_trajectories))

ggsave("FigS6_Traject.legend.pdf",
       plot = legend.trajectory,
       units="cm",
       width=36,
       height=5,
       dpi = 300)

#GLM for all ASVs. Insert ASV name in first line
ASV = df.genus[df.genus$OTU =="ASV_43",]
ASV$day = as.factor(ASV$day)
ASV$System= as.factor(ASV$System)
ASV$lineage = paste(ASV$System,ASV$bio.rep, sep = "")

# Linear mixed model (initial model with interaction)
lmm.test <- lmer(Abundance ~ System *day +(1|lineage), data=ASV)
Anova(lmm.test)
qqnorm(resid(lmm.test))
hist(rstudent(lmm.test))
highResIndx=which(abs(rstudent(lmm.test))>2.5)

# Linear mixed model with interaction and removed outliers
lmm.trajectory <- lmer(Abundance ~ System*day +( 1|lineage), data=ASV[-highResIndx,])
summary(lmm.trajectory) # lineage does have an effect so should be a random effect 
qqnorm(resid(lmm.trajectory))
plot(lmm.trajectory)
Anova(lmm.trajectory)
emm.trajectory = emmeans(lmm.trajectory, ~ System | day)
test = contrast(emm.trajectory, interaction = "pairwise", adjust = "Bonferroni")

ASV_43 = data.frame(contrast(emm.trajectory, interaction = "pairwise", adjust = "Bonferroni"))
ASV_43$TaxID = "ASV_43"


trajectories.emmeans = rbind(alt20, alt28, ASV_2, asv7,ASV_10, ASV_43,ASV_15, ASV_47, ASV_26, ASV_3)
write.csv(trajectories.emmeans, file = "Fig4/trajectories.emmeans.csv")

# Alteromonas 
# PS object with Phaeobacter ASVs
ps.alteromonas = ps.new %>% 
  phyloseq::filter_taxa(function (x) {sum(x > 0) > 4}, prune=TRUE) %>%
  phyloseq::transform_sample_counts(function(x) {x/sum(x)}*100 ) %>%
  phyloseq::tax_glom('genus') %>%
  phyloseq::subset_taxa(genus == "Alteromonas")

# Get reference sequences to inspect
refseq.altermonas = as.data.frame(refseq(ps.alteromonas))
#write.csv(refseq.altermonas, file = "refseq.alteromonas.csv") # one nucleotide differences, keep seperate

# Get abundance trajectory for main figure
significant = c("ASV_20","ASV_28","ASV_10","ASV_43","ASV_15","ASV_2","ASV_47","ASV_26","ASV_3","ASV_7")

ps.genus = ps.new %>% 
  phyloseq::filter_taxa(function (x) {sum(x > 0) > 4}, prune=TRUE) %>%
  phyloseq::tax_glom('genus') %>%
  phyloseq::transform_sample_counts(function(x) {x/sum(x)}*100 )

# plot the tree,

#get bootstrap 
plot_tree(ps.genus, color="family", label.tips="genus",
          nodelabf=nodeplotboot(), base.spacing=0.02, ladderize = "right")

ggtree(ps.genus) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab() # Node number

fam.color.2 = c("#C51B7D", "#8E0152", "#DE77AE","#F1B6DA","#80CDC1","#CB4335", "#35978F",
                "#BF812D","#DFC27D","#16A085","#01665E", "#003C30","#543005","#8C510A")

tr = ggtree(ps.genus, color="black", linetype=1, 
            branch.length='none', layout='fan', open.angle = 300) + 
  # geom_text2(aes(subset=!isTip, label=label), hjust=-.2, size=4) +
  geom_tiplab(aes(label=genus),fontface="bold.italic", size=3, hjust=-.1) +
  geom_tippoint(aes(color=family), size=4, alpha=1)+
  theme_tree(legend.position= "bottom") +
  # geom_text(aes(label=node), hjust=-.2) +
  scale_fill_manual(values = fam.color.2) +
  scale_color_manual(values = fam.color.2) +
  # scale_fill_brewer(palette="BrBG")+
  # scale_color_brewer(palette="BrBG")+
  theme(legend.text = element_text(face = "italic", size = 6))+
  xlim(-20, 20)
  tr

  ggsave(filename = "Tranjectoriestree.bold.pdf", plot = tr , width=30, height=30, units="cm", dpi = 400)
  legend.treetraj<- as_ggplot(ggpubr::get_legend(tr))
  
  ggsave("legend.treetraj.pdf",
         plot = legend.treetraj,
         units="cm",
         width=20,
         height=10,
         dpi = 300)
  
  
ggplot(sum.genus[sum.genus$OTU %in% significant,],aes(x = day, y = mean.abundance, fill=System, color = System, shape = System)) +
  geom_line(size = 1)+
  geom_point(size = 2)+
  geom_ribbon(aes(ymax = mean.abundance + sd.abundance, 
                  ymin = mean.abundance - sd.abundance, fill = System), alpha = 0.2, color = NA) +
  facet_wrap(.~TaxID, ncol = 5, nrow = 10, scales = "free_y") + 
  theme_bw(base_size = 12) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom",
        strip.text = element_text(face = "bold", colour ="black"), 
        legend.title = element_text(face = "bold", colour ="black"),
        legend.text = element_text(face = "bold", colour ="black"),
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_color_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  scale_fill_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  labs(x = "\n Time (days)", y = " Relative abundance (%)\n") + 
  scale_x_continuous(breaks = c(0,7,14,21,28,35,42,49,56,63,70))

#Seperate
ASV10 = ggplot(sum.genus[sum.genus$OTU == "ASV_10",], aes(x = day, y = mean.abundance, fill=System, color = System, shape = System)) +
  geom_line(size = 1)+
  geom_point(size = 2)+
  geom_ribbon(aes(ymax = mean.abundance + sd.abundance, 
                  ymin = mean.abundance - sd.abundance, fill = System), alpha = 0.2, color = NA) +
  facet_wrap(.~TaxID, ncol = 5, nrow = 10, scales = "free_y") + 
  theme_bw(base_size = 8) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = "white"),
        legend.position = "none",
        strip.text = element_text(face = "bold", colour ="black", size = 6), 
        legend.title = element_text(face = "bold", colour ="black"),
        legend.text = element_text(face = "bold", colour ="black"),
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black", size = 6),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_color_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  scale_fill_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  labs(x = "Time (days)", y = " Relative abundance (%)") + 
  scale_x_continuous(breaks = c(0,7,14,21,28,35,42,49,56,63,70))

ASV10
CD
Rhodo
ros
dong
rosei
hen
alt
Flavo

ggsave("ASV10.svg",
       plot = ASV10,
       units="cm",
       width=5,
       height=4,
       dpi = 300)

