# Assign taxonomy - September 2021

ls()
#"seqtab.nochim"
tt = assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", tryRC = TRUE, verbose = TRUE, multithread = TRUE)
ls()
#"seqtab.nochim" "tt"
tt.plus = addSpecies(tt, "silva_species_assignment_v138.1.fa.gz", verbose = TRUE, tryRC = TRUE, multithread = TRUE)
ls()
#"seqtab.nochim" "tt"
tt.plus = addSpecies(tt, "silva_species_assignment_v138.1.fa.gz", verbose = TRUE, tryRC = TRUE)
#303 out of 25390 were assigned to the species level.
#Of which 254 had genera consistent with the input table.>