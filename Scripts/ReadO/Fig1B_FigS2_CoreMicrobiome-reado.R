######################################################################################

                                ### Core Microbiome ###

#####################################################################################

#################
### Load data ###
#################
load("Data/ReadO/ps.asv.reduced.wTree.RData")

#######################
### Core microbiome ###
#######################

ps.RA = ps.new %>%
  phyloseq::subset_taxa(genus !="Phaeobacter") %>%
  phyloseq::prune_taxa(taxa_sums(.) > 0, .) %>%
  phyloseq::transform_sample_counts(function(x) {x/sum(x)}*100 )
  

#Make only top 20 across total dataset 
ps.RA.genus <-tax_glom(ps.RA, taxrank="genus")
TopASV20 = names(sort(taxa_sums(ps.RA.genus), TRUE)[1:20])
ps.top20 = prune_taxa(TopASV20, ps.RA.genus)
test = merge_samples(ps.top20, "exp.condition")

#Check tree
plot_tree(test, color="family", label.tips="genus",
          nodelabf=nodeplotboot(), base.spacing=0.02, ladderize = "right")
  
# Make tree with ggtree

#Color scheme for family
fam.color = c("#e488a0", "#6cb643","#a35bcd","#bdad46","#596ec6","#d05335","#4ea9cf","#d9436c","#55b483",
              "#d352a9","#a44d51","#647a35","#b392d9","#c48242","#9a508b","#a44d51")

fam.color.2 = c("#C51B7D", "#8E0152", "#DE77AE","#F1B6DA","#80CDC1", "#35978F",
                "#BF812D","#DFC27D","#01665E", "#003C30","#543005","#8C510A")

# tree 
tree.top20 = ggtree(test, color="black", size=1.5, linetype=1, right=TRUE) + 
  #geom_text2(aes(subset=!isTip, label=label), hjust=-.2, size=4) +
  geom_tiplab(aes(label=genus, color = family, 
                  fontface="bold", size=4), hjust=-.1) +
  geom_tippoint(aes(color=family), size=5, alpha=1)+
  theme_tree(legend.position= "none") +
  #geom_text(aes(label=node), hjust=-.2) +
  scale_fill_manual(values = fam.color.2) +
  scale_color_manual(values = fam.color.2) +
  # scale_fill_brewer(palette="BrBG")+
  # scale_color_brewer(palette="BrBG")+
  theme(legend.text = element_text(face = "italic", size = 8))+
  xlim(0, 0.9)
tree.top20

ggsave(filename = "top20tree.pdf", plot = tree.top20 , width=14, height=12, units="cm", dpi = 400)

legend.top20tree<- as_ggplot(ggpubr::get_legend(tree.top20))

ggsave("legend.top20tree.pdf",
       plot = legend.top20tree,
       units="cm",
       width=10,
       height=10,
       dpi = 300)

#Export tree for Itol 
#ape::write.tree(phy_tree(test), "Data/tree_top20.tree")

# Get abundance data
# Summarize across entire study
df.genus.core.sum = psmelt(ps.top20) %>%
  group_by(genus, family) %>%
  summarise(mean.abundance = mean(Abundance),
            sd.abundance = sd(Abundance),
            count = n())
# Figure
S2 = ggplot(df.genus.core.sum, aes(x = reorder(genus, -mean.abundance), y = mean.abundance, color = family, fill = family)) + 
  geom_bar(stat = "identity", alpha = 0.7) +
  theme_bw(base_size = 10) +
  geom_errorbar(aes(x = genus, ymax = mean.abundance + sd.abundance, 
                    ymin = mean.abundance - sd.abundance, 
                    color = family, group = genus), size = 0.5, width = 0, 
                position=position_dodge(5), alpha = 1)+
  ylab("Relative abundance (%) \n") +
  scale_fill_manual(values = fam.color.2) +
  scale_color_manual(values = fam.color.2) +
  xlab("")+
  ggtitle("") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(face = "bold", angle = 270,hjust=0.05, vjust= 0.1),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        legend.position = "none")
S2

ggsave("FigS2-coremicrobiome.pdf",
       plot = S2,
       units="cm",
       width=12,
       height=10,
       dpi = 300)


