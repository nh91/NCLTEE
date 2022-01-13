######################################################################################

                                    ### Composition ###

#####################################################################################

#################
### Load data ###
#################

load("Data/ps.genus.reduced.RData")


#######################
### Stacked barplot ###
#######################

# Make df
df.genus = ps.genus.reduced %>%
  transform_sample_counts(function(x) 100*{x/sum(x)} ) %>%
  psmelt() %>%                                          
  arrange(genus)                                      

# Tabel
#write.csv(df.genus, file = "genus_composition.csv")

# Color palette
n <- 40
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Plot with no division between lineages
ggplot(df.genus, aes(x = bio.rep, y = Abundance, fill = genus)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance (Genus > 10%) \n") +
  xlab("\n Lineage")+
  scale_fill_manual(values = col_vector)+
  facet_grid(exp.condition~time) +
  ggtitle("genus")

# Plot with division into lineages
comp.lineage.RA = ggplot(df.genus, aes(x = day, y = Abundance, fill = genus)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("Relative abundance (%)\n") +
  xlab("\n Time (day)")+
  scale_fill_manual(values = col_vector)+
  facet_grid(bio.rep~exp.condition) +
  guides(fill=guide_legend(title="Genus", ncol = 4, byrow = TRUE))+
  #scale_x_continuous(limits=c(0,70), breaks=seq(0,70,by =7))+
  theme(strip.background = element_rect(color="black", fill="#FFFFFF"),
        strip.text = element_text(face = "bold"),
        strip.text.y = element_text(angle=0)) +
  theme(legend.text = element_text(face = "bold.italic"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")+
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"))+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
comp.lineage.RA


ggsave("FigS4.svg",
       plot = comp.lineage.RA,
       units="cm",
       width=18,
       height=25,
       dpi = 300)
