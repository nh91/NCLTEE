### Phylogenetic tree and Unifrac distances ###

# Phyloseq object to use (cleaned) 
load("Data/ReadO/ps.asv.reduced.reado.RData")
ps.asv.reduced

#Make DNAbin
seqs = as.DNAbin(refseq(ps.asv.reduced))

#alignen with mafft tool and putting on gap penalty of 1
seqs_align <- mafft(seqs, exec = "/usr/local/bin/mafft", 
                    options = c("--adjustdirection"), ep = 1) # ep = gap extension penalty.
#View(as.character(seqs_align))

#visualize the alginment
image(seqs_align, cex.lab = 0.2)
table(nchar(as.character(refseq(ps.asv.reduced)))) # check length 

#collapsing the alignment into a format that we can make a fasta file of
seqs_alignX <- apply(as.character(seqs_align), 1, function(x) paste0(x, collapse = ""))
cat(file="asv_aligned.fasta", paste(paste0(">",names(seqs_alignX)), seqs_alignX, sep="\n"), sep="\n");

#correct ASV names
asv_aligned = seqinr::read.fasta("asv_aligned.fasta")
df.asv_aligned = data.frame(ASV=names(asv_aligned), Seqs=unlist(getSequence(asv_aligned, as.string=T)))
df.asv_aligned$ASV <- gsub("_R_","", as.character(df.asv_aligned$ASV))
cat(file="asv_aligned_correctnames.fasta", paste(paste0(">",df.asv_aligned$ASV), df.asv_aligned$Seqs, sep="\n"), sep="\n");

# Fast tree infers approximately-maximum-likelihood phylogenetic, GTR generalized time-reversible 
#runing Fasttree through the console by using the system2 function
system2("/Users/nnseh/FastTree", args = c("-gtr", "-nt", "-out", "asv_aligned.tree", "asv_aligned_correctnames.fasta"))

#importing the new tree 
treeNew <- read.tree("asv_aligned.tree")
plot(treeNew, cex = 0.5)
#Picking an outgroup automatically as the longest branch
pick_new_outgroup = function(tree.unrooted) {
  require("magrittr")
  require("data.table")
  require("ape") #ape::Ntip
  #tablify parts of tree that we need
  treeDT = 
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup = treeDT[which.max(length)]$id
  return(new.outgroup)
}

new_root = pick_new_outgroup(treeNew)

#rooting the new tree
treeNewR <- root(treeNew, outgroup = new_root, resolve.root = TRUE)

#plotting
plot(treeNewR, cex = 0.5)

#putting tree back in
ps.new = ps.asv.reduced
phy_tree(ps.new) <- phy_tree(treeNewR)
ps.new

#Save PS with tree
save(ps.new, file = "Data/ReadO/ps.asv.reduced.wTree.RData")
