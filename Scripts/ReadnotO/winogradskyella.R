# winogradskyella

load("Data/ps.asv.reduced.wTree.RData")

########################
### Phaeobacter ASVs ###
########################

# PS object with Phaeobacter ASVs
ps.wino = ps.new %>% 
  phyloseq::transform_sample_counts(function(x) {x/sum(x)}*100 ) %>%
  phyloseq::filter_taxa(function(x) sum(x) > .5, TRUE) %>% # Remove ASV with abundance below 0.5%
  phyloseq::subset_taxa(genus == "Winogradskyella") %>%
  phyloseq::prune_taxa(taxa_sums(.) >0, .)

####################
### Inspect ASVs ###
####################

# Make df for plotting
df.wino = ps.wino %>%
  psmelt() %>%
  dplyr::rename(System = exp.condition)
df.wino$OTU = as.factor(df.wino$OTU)


#Overview plot of ASVs abundance
ggplot(df.wino, aes(day, Abundance, fill = OTU)) + 
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
