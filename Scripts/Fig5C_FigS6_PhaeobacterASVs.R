######################################################################################

              ### Relative abundance of Phaeobacter ASVs ###

######################################################################################


#################
### Load data ###
#################

#PS object with ASVs
load("Data/ps.clean.RData")


########################
### Phaeobacter ASVs ###
########################

# PS object with Phaeobacter ASVs
ps.phaeobacter = ps_filt2 %>% 
  transform_sample_counts(function(x) {x/sum(x)}*100 ) %>%
  filter_taxa(function(x) sum(x) > .5, TRUE) %>% # Remove ASV with abundance below 0.5%
  subset_taxa(genus == "Phaeobacter") %>%
  prune_taxa(taxa_sums(.) >0, .)

# Get reference sequences for BLAST 
refseq.phaeobacter = as.data.frame(refseq(ps.phaeobacter))
write.csv(refseq.phaeobacter, file = "Data/refseq.phaeobacter.csv")


####################
### Inspect ASVs ###
####################

# Make df for plotting
df.phaeobacter = ps.phaeobacter %>%
  psmelt() %>%
  rename(System = exp.condition)
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
  theme_bw(base_size = 8) + 
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

ggsave("Fig5S.svg",
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

# Line plot
Fig6C = ggplot(sum.phaeobacter.subset, aes(x = day, y = mean.abundance, fill=ASV, color = ASV)) +
  geom_line(size = 1)+
  geom_point(size = 2)+
  geom_ribbon(aes(ymax = mean.abundance + sd.abundance, ymin = mean.abundance - sd.abundance, fill = ASV), alpha = 0.2, color = NA) +
  facet_grid(.~System) + 
  theme_bw(base_size = 8) + 
  theme(strip.text.y = element_text(angle = 0), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = "white"),
        legend.position = c(0.92,0.85), legend.background = element_rect(colour = "grey"),
        strip.text = element_text(face = "bold", colour ="black"), 
        legend.title = element_text(face = "bold", colour ="black"),
        legend.text = element_text(face = "bold", colour ="black"),
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = c("#4E2A1E",  "#D8B70A"))+
  scale_color_manual(values = c( "#4E2A1E",  "#D8B70A")) +
  labs(x = "\n Time (day)", y = " Relative abundance (%)\n") + 
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14))+ scale_x_continuous(breaks = c(0,7,14,21,28,35,42,49,56,63,70))
Fig6C

#Save figure
ggsave("Fig6C.svg",
       plot = Fig6C,
       units="cm",
       width=18,
       height=18,
       dpi = 300)

# Dodge barplot (not fixed geom_bar) zoomed in on day 35-70
Fig6C2 = ggplot(sum.phaeobacter.subset %>% filter(28<day) %>% filter("TDA" == System), aes(x = day, y = mean.abundance, fill=ASV, color = ASV)) +
  geom_bar(alpha = 1, stat ="identity", position=position_dodge(), width = 5)+
  geom_errorbar(aes(x = day, ymax = mean.abundance + sd.abundance, ymin = mean.abundance - sd.abundance, 
                    color = ASV, group = ASV), size = 1, width = 0, position=position_dodge(5), alpha = 1)+
  #facet_grid(.~System) + 
  theme_bw(base_size = 8) + 
  theme(strip.text.y = element_text(angle = 0), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", colour ="black"), 
        legend.title = element_text(face = "bold", colour ="black"),
        legend.text = element_text(face = "bold", colour ="black"), legend.position = "none",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = c("#4E2A1E",  "#D8B70A"))+
  scale_color_manual(values = c( "#4E2A1E",  "#D8B70A")) +
  labs(x = "\n Time (day)", y = " Relative abundance (%)\n") + 
  scale_x_continuous(breaks = c(35,42,49,56,63,70))
Fig6C2


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

