# Rhodobactereceae tree

# Taxa present in more than five replicates 

ps.rhodobac = ps.new %>%
  phyloseq::subset_taxa(family =="Rhodobacteraceae") %>%
  phyloseq::transform_sample_counts(function(x) {x/sum(x)}*100 ) %>%
  phyloseq::filter_taxa(function(x) sum(x) > .5, TRUE)


ps.RA.genus <-tax_glom(ps.rhodobac, taxrank="genus")
TopASV20 = names(sort(taxa_sums(ps.RA.genus), TRUE)[1:20])
ps.top20 = prune_taxa(TopASV20, ps.RA.genus)

#test tree look
plot_tree(tip_glom(ps.rhodobac, h = 0.005),
          color = "OTU") + coord_polar(theta="y")

#get Tree with tip_gloom
tree.Rhodobacteraceae = phy_tree(tip_glom(ps.rhodobac, h = 0.005))
#For export
ape::write.tree(tree.Rhodobacteraceae, "tree_Rhodobacteraceae.tree")

#make tree
tree.rhodo= ggtree(ps.rhodobac, color="black", size=1, linetype=1,
                   layout="fan", open.angle=120) +
  #geom_text2(aes(subset=!isTip, label=label), hjust=-.2, size=4) +
  geom_tiplab(aes(label=label, 
                  fontface="bold", size=4), hjust=-.2) +
  geom_tippoint(aes(color=genus), size=4, alpha=0.6, stroke = 0)+
  theme_tree(legend.position= "bottom") +
  #geom_text(aes(label=node), hjust=-.2) +
  scale_fill_manual(values = fam.color) +
  scale_color_manual(values = fam.color) +
  theme(legend.text = element_text(face = "italic", size = 8))
tree.rhodo


ggsave(filename = "tree.rhodo.svg", plot = tree.rhodo , width=18, height=18, units="cm", dpi = 400)

legend.top20tree<- as_ggplot(ggpubr::get_legend(tree.top20))

ggsave("legend.top20tree.pdf",
       plot = legend.top20tree,
       units="cm",
       width=20,
       height=4,
       dpi = 300)
