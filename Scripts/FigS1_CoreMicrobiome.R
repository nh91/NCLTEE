######################################################################################

                                ### Core Microbiome ###

#####################################################################################

#################
### Load data ###
#################

load("Data/ps.genus.reduced.RData")

#######################
### Core microbiome ###
#######################

# Taxa present in more than five replicates 

ps.genus.core0 = ps.genus.reduced %>%
  subset_taxa( genus !="Phaeobacter") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) {x/sum(x)}*100 ) %>%
  filter_taxa(function (x) {sum(x > 0) > 5}, prune=TRUE)


# Summarize across entire study
df.genus.core.sum = psmelt(ps.genus.core0) %>%
  group_by(genus) %>%
  summarise(mean.abundance = mean(Abundance),
            sd.abundance = sd(Abundance),
            count = n())

# Figure
S2 = ggplot(df.genus.core.sum, aes(x = reorder(genus, -mean.abundance), y = mean.abundance, color = genus, fill = genus)) + 
  geom_bar(stat = "identity", alpha = 0.5) +
  theme_bw() +
  geom_errorbar(aes(x = genus, ymax = mean.abundance + sd.abundance, 
                    ymin = mean.abundance - sd.abundance, 
                    color = genus, group = genus), size = 0.5, width = 0, 
                position=position_dodge(5), alpha = 1)+
  ylab("Relative abundance (%) \n") +
  xlab("")+
  ggtitle("") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(face = "bold.italic", angle = 270,hjust=0.05, vjust= 0.1),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.position = "none")


ggsave("FigS2-coremicrobiome.pdf",
       plot = S2,
       units="cm",
       width=18,
       height=16,
       dpi = 300)
