#####################################################################################

                     ### Relative abundance of Phaeobacter genus ###

#####################################################################################


#################
### Load data ###
#################

load("Data/ps.genus.reduced.RData")


############################
### Polish PS and get df ###
############################

# Relative abundance
genus_relative= ps.genus.reduced %>%
  transform_sample_counts(function(x) {100*x/sum(x)} )

# Subset to Phaeobacter and melt PS
phaeobacter = subset_taxa(genus_relative, genus == "Phaeobacter")%>%
  psmelt()

# Polish metadata
phaeobacter$bio.rep <- as.factor(phaeobacter$bio.rep)

# Summarize across bio.rep
sum.phaeobacter = phaeobacter %>%
  group_by(day, exp.condition) %>%
  summarise( mean.RA = mean(Abundance),
             sd.RA = sd(Abundance))

#Rename factor and order level
sum.phaeobacter$exp.condition = factor(sum.phaeobacter$exp.condition, levels= c("TDA", "NoTDA", "Control"))
sum.phaeobacter = rename(sum.phaeobacter, System = exp.condition)

##############
### Figure ###
##############

# Color palette
system_color2 = c("#5B1A18", "#7294D4", "#D8A499")

# Line and ribbon plot
pRA = ggplot(as.data.frame(sum.phaeobacter), aes( x = day, y=mean.RA, color= System, shape = System)) + 
  geom_line(size = 1)+
  geom_point(size = 2)+
  geom_ribbon(aes(ymax = mean.RA + sd.RA, ymin = mean.RA - sd.RA, fill = System), alpha = 0.2, color = NA) +
  labs(x="\n Time (day) \n",y=" \n Relative abundance (%) \n") +
  theme_bw(base_size = 8) + 
  scale_color_manual(values=system_color2) + scale_fill_manual(values=system_color2)+
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14))+ scale_x_continuous(breaks = c(0,7,14,21,28,35,42,49,56,63,70))+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(face ="bold", colour ="black"),
        legend.position = c(0.85,0.8), legend.background = element_rect(colour = "grey"),
        axis.title.y = element_text(face = "bold"), 
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
pRA


#############
### Stats ###
#############

#Make factors ready
phaeobacter$day = as.factor(phaeobacter$day)
phaeobacter$lineage = paste(phaeobacter$exp.condition,phaeobacter$bio.rep, sep = "")

#Linear mixed model
lmm.RAPhaeobacter <- lmer(Abundance ~ exp.condition*day +( 1|lineage), data=phaeobacter)
summary(lmm.RAPhaeobacter) # lineage does have an effect so should be a random effect 
Anova(lmm.RAPhaeobacter)
emm.RAPhaeobacter = emmeans(lmm.RAPhaeobacter, ~ exp.condition | day)
contrast(emm.RAPhaeobacter, interaction = "pairwise", adjust = "Bonferroni")
