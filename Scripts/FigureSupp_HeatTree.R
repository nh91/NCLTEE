load("Data/ReadO/ps.asv.reduced.wTree.RData")

ps.TDA.NoTDA.day42 <- ps.new %>%
  subset_samples(day == "42") %>%
  subset_taxa(genus != "Phaeobacter") %>%
  transform_sample_counts(function(x) x/sum(x))

#ps.TDA.NoTDA.day42 <- tax_glom(ps.TDA.NoTDA.day42, "genus")



tax_table(ps.TDA.NoTDA.day42) <- tax_table(ps.TDA.NoTDA.day42)[,1:6]

# Create metacoder environment
#TopASV50 = names(sort(taxa_sums(ps.TDA.NoTDA.day42), TRUE)[1:30])
#ps.day42.top50 = prune_taxa(TopASV50, ps.TDA.NoTDA.day42)

#Convert to metacoder object
metacoder <- parse_phyloseq(ps.TDA.NoTDA.day42)
#metacoder <- parse_phyloseq(physeq_100)
metacoder$data$tax_abund <- calc_taxon_abund(metacoder, data = "otu_table")
metacoder$data$tax_occ <- calc_n_samples(metacoder, "tax_abund", groups = "exp.condition")

#Compare groups
metacoder$data$diff_table <- compare_groups(metacoder, 
                                            data = "tax_abund", 
                                            cols = metacoder$data$sample_data$sample_id,
                                            groups = metacoder$data$sample_data$exp.condition)
#Adjust p-value
metacoder$data$diff_table$adjusted_p_value <- p.adjust(metacoder$data$diff_table$wilcox_p_value,
                                                       method = "fdr")

#metacoder$data$diff_table$log2_median_ratio[metacoder$data$diff_table$adjusted_p_value > 0.05] <- 0

set.seed(41)
metacoder.plot42 <- heat_tree_matrix(metacoder,
                                   data = "diff_table",
                                   node_size = n_obs,
                                   node_label = taxon_names,
                                   node_color = log2_median_ratio,
                                   node_color_range = diverging_palette(),
                                   node_color_trans = "linear",
                                   node_color_interval = c(-5, 5),
                                   edge_color_interval = c(-5, 5),
                                   key_size = 0.7,
                                   node_label_size_range = c(0.025, 0.03),
                                   overlap_avoidance = 5,
                                   layout = "davidson-harel", # The primary layout algorithm
                                   initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                                   node_size_axis_label = "# of ASVs",
                                   node_color_axis_label = "Log2 ratio median proportions",
                                   output_file = "day42-heattree.svg")
metacoder.plot42

