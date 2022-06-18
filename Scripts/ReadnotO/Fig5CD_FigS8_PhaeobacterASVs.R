######################################################################################

              ### Relative abundance of Phaeobacter ASVs ###

######################################################################################


#################
### Load data ###
#################

#PS object with ASVs
load("Data/ps.clean.RData")
load("Data/ps.asv.reduced.wTree.RData")

########################
### Phaeobacter ASVs ###
########################

# PS object with Phaeobacter ASVs
ps.phaeobacter = ps.new %>% 
  phyloseq::transform_sample_counts(function(x) {x/sum(x)}*100 ) %>%
  phyloseq::filter_taxa(function(x) sum(x) > .5, TRUE) %>% # Remove ASV with abundance below 0.5%
  phyloseq::subset_taxa(genus == "Phaeobacter") %>%
  phyloseq::prune_taxa(taxa_sums(.) >0, .)

# Get reference sequences for BLAST 
refseq.phaeobacter = as.data.frame(refseq(ps.phaeobacter))
write.csv(refseq.phaeobacter, file = "Data/refseq.phaeobacter.csv")


####################
### Inspect ASVs ###
####################

# Make df for plotting
df.phaeobacter = ps.phaeobacter %>%
  psmelt() %>%
  dplyr::rename(System = exp.condition)
df.phaeobacter$OTU = as.factor(df.phaeobacter$OTU)

# Get reference sequences out for inspection manually
refseq.phaeobacter = as.data.frame(refseq(ps.phaeobacter))
write.csv(refseq.phaeobacter, file = "Data/refseq.phaeobacter.csv")

#Overview plot of ASVs abundance
FigS5 = ggplot(df.phaeobacter, aes(day, Abundance, fill = OTU)) + 
  geom_bar(alpha = 0.7, stat ="identity", position=position_dodge(), width = 5) + 
  facet_grid(bio.rep~System) + 
  labs(x = "\n Time (day)", y = " Relative abundance (%)\n") +
  guides(fill=guide_legend(title="ASV", ncol = 2))+
  theme_bw(base_size = 10) + 
  theme(strip.text.y = element_text(angle = 0), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom",
        strip.text = element_text(face = "bold", colour ="black"), 
        legend.title = element_text(face = "bold", colour ="black"),
        legend.text = element_text(face = "bold", colour ="black"),
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black"))
FigS5 
ggsave("Fig5S.pdf",
       plot = FigS5,
       units="cm",
       width=18,
       height=18,
       dpi = 300)



# Conclusion: ASV_35/36 is identified to be one population and ASV_56/57 another


####################################
### Make df with ASV35 and ASV36 ###
####################################

# Make df with ASV35 and 36
df.phaeobacter.35 = df.phaeobacter[df.phaeobacter$OTU == "ASV_35",]
df.phaeobacter.36 = df.phaeobacter[df.phaeobacter$OTU == "ASV_36",]
df.phaeobacter.3536 = rbind(df.phaeobacter.36, df.phaeobacter.35)

# Sum up to total abundance 
sum.3536 = df.phaeobacter.3536 %>%
  group_by(day, bio.rep, System) %>%
  summarise(sum.abundance = sum(Abundance),
            count = n())
sum.3536$ASV = "ASV_35/36"


####################################
### Make df with ASV56 and ASV57 ###
####################################

# Make df with ASV56 and 57
df.phaeobacter.56 = df.phaeobacter[df.phaeobacter$OTU == "ASV_56",]
df.phaeobacter.57 = df.phaeobacter[df.phaeobacter$OTU == "ASV_57",]
df.phaeobacter.5657 = rbind(df.phaeobacter.56, df.phaeobacter.57)

# Sum up to total abundance
sum.5657 = df.phaeobacter.5657 %>%
  group_by(day, bio.rep, System) %>%
  summarise(sum.abundance = sum(Abundance),
            count = n())
sum.5657$ASV = "ASV_56/57"


###########################
### Get total abundance ###
###########################

# df.phaeobacter.total = df.phaeobacter
# sum.phaeobacter.total = df.phaeobacter.total %>%
#   group_by(day, bio.rep, System) %>%
#   summarise(sum.abundance = sum(Abundance),
#             count = n())
# sum.phaeobacter.total$ASV = "Total"

########################
### Make combined df ###
########################

# Combine dfs
df.phaeobacter.subset = rbind(sum.3536,sum.5657)

#Calculate mean for bio.rep
sum.phaeobacter.subset = df.phaeobacter.subset %>%
  group_by(day, System, ASV) %>%
  summarise(mean.abundance = mean(sum.abundance),
            sd.abundance = sd(sum.abundance),
            count = n())


##############
### Figure ###
##############

# Change levels
sum.phaeobacter.subset$System = factor(sum.phaeobacter.subset$System , levels= c("TDA", "NoTDA", "Control"))
system_color2 = c("#9C3206", "#00A091", "black")
labels<- c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control")
names(labels) = c("TDA", "NoTDA", "Control")
# Line plot
Fig6C = ggplot(sum.phaeobacter.subset, 
               aes(x = day, y = mean.abundance, 
                   fill=ASV, color = ASV, shape = System)) +
  geom_line(size = 1)+
  geom_point(size = 2)+
  geom_ribbon(aes(ymax = mean.abundance + sd.abundance, ymin = mean.abundance - sd.abundance, fill = ASV), alpha = 0.2, color = NA) +
  facet_wrap(vars(System), nrow = 3) + 
  theme_bw(base_size = 10) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = "white"),
        legend.position = "none",
        strip.text = element_text(face = "bold", colour ="black"), 
        legend.title = element_text(face = "bold", colour ="black"),
        legend.text = element_text(face = "bold", colour ="black"),
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = c("#CE738C", "Black"))+
  scale_color_manual(values = c( "#CE738C", "Black"))+
  labs(x = "\n Time (days)", y = " Relative abundance (%)\n") + 
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14))+ 
  scale_x_continuous(breaks = c(0,7,14,21,28,35,42,49,56,63,70))
Fig6C

#Save figure
ggsave("Fig6C.pdf",
       plot = Fig6C,
       units="cm",
       width=10,
       height=18,
       dpi = 300)

legend.Fig6C <- cowplot::get_legend(Fig6C)
grid.newpage()
grid.draw(legend.Fig6C)

# zoom in on day 42-70
df.phaeobacter.subset$System = factor(df.phaeobacter.subset$System , levels= c("TDA", "NoTDA", "Control"))
Fig6C.TDA = ggplot(sum.phaeobacter.subset %>% filter(35<day) %>% filter(System == "TDA"),
                aes(x = day, y = mean.abundance, fill=ASV, color = ASV)) +
  geom_point(data = df.phaeobacter.subset %>% filter(35<day) %>% filter(System == "TDA"),
             aes(x = day, y = sum.abundance,fill=ASV, stroke = 0),shape = 16,
               alpha = 0.6, stat ="identity", position=position_dodge(3))+
  geom_point(alpha = 1, stat ="identity", position=position_dodge(3), size = 1.5, shape = 16)+
  geom_errorbar(aes(x = day, ymax = mean.abundance + sd.abundance, ymin = mean.abundance - sd.abundance, 
                    color = ASV, group = ASV), size = 0.5, width = 0, position=position_dodge(3), alpha = 0.6)+
  theme_bw(base_size = 8) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_blank(), 
        legend.title = element_text(face = "bold", colour ="black"),
        legend.text = element_text(face = "bold", colour ="black"), 
        legend.position = "bottom",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = c("#CE738C", "Black"))+
  scale_color_manual(values = c( "#CE738C", "Black")) +
  labs(x = "", y = "") + 
  scale_y_continuous(limits = c(-0.1,1.6), breaks = c(0,0.5,1,1.5))+
  scale_x_continuous(breaks = c(42,49,56,63,70))
Fig6C.TDA

ggsave("Fig6C.TDA.pdf",
       plot = Fig6C.TDA,
       units="cm",
       width=4,
       height=4,
       dpi = 300)

#Get legend
legend.Fig6C <- as_ggplot(ggpubr::get_legend(Fig6C))

ggsave("legend.Fig6C.pdf",
       plot = legend.Fig6C,
       units="cm",
       width=3,
       height=7,
       dpi = 300)


#############
### Stats ###
#############

#Make df for ASV56/57
df.p.5657 = df.phaeobacter.5657 %>% select(Abundance, System, bio.rep, day)
df.p.5657$ASV = "ASV_56/57"
df.p.5657$lineage = paste(df.p.5657$System,df.p.5657$bio.rep, sep = "")
df.p.5657$day = as.factor(df.p.5657$day)

#Make df for ASV35/36
df.p.3536 = df.phaeobacter.3536 %>% select(Abundance, System, bio.rep, day)
df.p.3536$ASV = "ASV_35/36"
df.p.3536$lineage = paste(df.p.3536$System,df.p.3536$bio.rep, sep = "")
df.p.3536$day = as.factor(df.p.3536$day)

#Linear mixed model for ASV56/57 - difference in abundance between system 
lmm.5657 <- lmer(Abundance ~ System*day +( 1|lineage), data=df.p.5657)
summary(lmm.5657) # lineage does have an effect so should be a random effect 
Anova(lmm.5657)
emm.5657 = emmeans(lmm.5657, ~ System | day)
contrast(emm.5657, interaction = "pairwise", adjust = "Bonferroni")

#Linear mixed model ASV35/36 - difference in abundance between system 
lmm.3536 <- lmer(Abundance ~ System*day +( 1|lineage), data=df.p.3536)
summary(lmm.3536) # lineage does have an effect so should be a random effect 
Anova(lmm.3536)
emm.3536 = emmeans(lmm.3536, ~ System | day)
contrast(emm.3536, interaction = "pairwise", adjust = "Bonferroni")

# Get 42 to 70 for each asv
Fig6D.ASVs = ggplot(sum.phaeobacter.subset %>% filter(35<day),
                   aes(x = day, y = mean.abundance, fill=System, color = System)) +
  geom_point(data = df.phaeobacter.subset %>% filter(35<day),
             aes(x = day, y = sum.abundance,fill=System, stroke = 0, shape = System),
             alpha = 0.6, stat ="identity", position=position_dodge(3))+
  geom_point(aes(shape = System),alpha = 1, stat ="identity", position=position_dodge(3), size = 2.5)+
  geom_errorbar(aes(x = day, ymax = mean.abundance + sd.abundance, ymin = mean.abundance - sd.abundance, 
                    color = System, group = System), size = 0.5, width = 0, position=position_dodge(3), alpha = 0.6)+
  facet_wrap(vars(ASV), nrow = 2)+
  theme_bw(base_size = 10) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.title = element_text(face = "bold", colour ="black"),
        legend.text = element_text(face = "bold", colour ="black"), 
        strip.background = element_rect(fill = "white"),
        legend.position = c(0.8,0.9),
        strip.text = element_text(face = "bold", colour ="black"), 
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = system_color2 )+
  scale_color_manual(values = system_color2 ) +
   labs(x = "\nTime (days)", y = "Relative abundance (%)\n") + 
  scale_y_continuous(limits = c(-0.1,3.0), breaks = c(0,0.5,1,1.5,2.0,2.5,3.0))+
  scale_x_continuous(breaks = c(42,49,56,63,70))
Fig6D.ASVs

ggsave("Fig6D.ASVs.pdf",
       plot = Fig6D.ASVs,
       units="cm",
       width=8,
       height=16,
       dpi = 300)
