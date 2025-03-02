## kkyle@KlassenLab-server:~/kek16SallITStrachyDecons/fastq/p1p2_16S/
## R v4.4.1

##########################################################################################################################################
######### SET THESE VARIABLES TO YOUR SPECIFICS ##########################################################################################
##########################################################################################################################################

n.threads = 20
libs = c("phyloseq", "tidyverse", "phytools")

work.dir = "."
pdf.dir = file.path(work.dir,"pdf")
rds.dir = file.path(work.dir,"rds")
csv.dir = file.path(work.dir,"csv")
fna.dir = file.path(work.dir,"fasta")
rdat.dir = file.path(work.dir,"rdata")
tree.dir = file.path(work.dir,"trees")

##########################################################################################################################################


## ----load libraries---------------------------------------------------------------------------------------------------------------------

options("Ncpus" = n.threads)

if(!require(pacman, quietly=T)) {
	install.packages("pacman")
}
library(pacman)
p_load(char=libs)

cat(">>> Libraries loaded:\n")
for( i in seq_along(libs) ) {
	cat("\t", libs[i], "version", as.character(packageVersion(libs[i])),"\n")
}



## ---- add phylogenetic trees ---------------------------------------------------------------------------- ##


decons.ps.p.filt = readRDS("./decons/decons.ps.filt.RDS")
tsinf.ps.p.filt = readRDS("./tsinf/tsinf.ps.filt.RDS")

# get ASVs as fasta files
tsinf.ps.p.filt %>% 
	refseq() %>% 
	Biostrings::writeXStringSet( 
		"./tsinf/tsinf.ps.p.filt_ASVs.fasta", 
		format="fasta",
		append=F,
		compress=F,
		compression_level=NA
	)



sn.l = grepl(sample_names(decons.ps.p.filt), pattern="^d3")
decon3.ps.p.filt = subset_samples(decons.ps.p.filt, sn.l)
# get ASVs as fasta files
decon3.ps.p.filt %>% 
	refseq() %>% 
	Biostrings::writeXStringSet( 
		"./decons/decon3/decon3.ps.p.filt_ASVs.fasta", 
		format="fasta",
		append=F,
		compress=F,
		compression_level=NA
	)

sn.l = grepl(sample_names(decons.ps.p.filt), pattern="^d4")
decon4.ps.p.filt = subset_samples(decons.ps.p.filt, sn.l)
decon4.ps.p.filt %>% 
	refseq() %>% 
	Biostrings::writeXStringSet( 
		"./decons/decon4/decon4.ps.p.filt_ASVs.fasta", 
		format="fasta",
		append=F,
		compress=F,
		compression_level=NA
	)

sn.l = grepl(sample_names(decons.ps.p.filt), pattern="^d5")
decon5.ps.p.filt = subset_samples(decons.ps.p.filt, sn.l)
decon5.ps.p.filt %>% 
	refseq() %>% 
	Biostrings::writeXStringSet( 
		"./decons/decon5/decon5.ps.p.filt_ASVs.fasta", 
		format="fasta",
		append=F,
		compress=F,
		compression_level=NA
	)


## ---- BASH ---------------------------------------------------------------------------------------------- ##

### kkyle@KlassenLab-server:~/kek16SallITStrachyDecons/fastq/p1p2_16S
# conda activate mafft
# mafft --auto --thread 20 ps.nonc_ASVs.fasta > ps.nonc_ASVs.fasta.fft-ns-2_align 
# mafft --thread 20 --maxiterate 1000 ps.nonc_ASVs.fasta > ps.nonc_ASVs.fasta.fft-ns-i1k_align
# conda deactivate; conda activate fasttree
# fasttreeMP -gtr -gamma -nome -boot 500 -seed 1103 -log fasttree-gtr-gamma-boot500.log -nt ps.nonc_ASVs.fasta.fft-ns-i1k_align > ps.nonc_ASVs.fasta.fft-ns-i1k_align.ft-gtr-gamma-boot500.tree
# conda deactivate

## ---- \BASH --------------------------------------------------------------------------------------------- ##

treeLs = list(decon3.tree = read.newick(file="decons/decon3/decon3.ps.p.filt_ASVs.fasta.fft-ns-i1k_align.ft-gtr-gamma-nome-boot500.tree"),
	decon4.tree = read.newick(file="decons/decon4/decon4.ps.p.filt_ASVs.fasta.fft-ns-i1k_align.ft-gtr-gamma-nome-boot500.tree"),
	decon5.tree = read.newick(file="decons/decon5/decon5.ps.p.filt_ASVs.fasta.fft-ns-i1k_align.ft-gtr-gamma-nome-boot500.tree"),
	tsinf.tree = read.newick(file="tsinf/tsinf.ps.p.filt_ASVs.fasta.fft-ns-i1k_align.ft-gtr-gamma-nome-boot500.tree") 
	)
# midpoint root and write to file
for( i in names(treeLs) ) {
	rtree = midpoint_root(treeLs[[i]])
	f = paste0(i, "midptree.nwk")
	write.tree(rtree, f)
	treeLs[[i]] = rtree
}
# add rooted trees to phyloseq objects
phy_tree(decon3.ps.p.filt) = treeLs$decon3.tree; saveRDS(decon3.ps.p.filt, "decon3.ps.p.filt.rtree.RDS")
phy_tree(decon4.ps.p.filt) = treeLs$decon4.tree; saveRDS(decon4.ps.p.filt, "decon4.ps.p.filt.rtree.RDS")
phy_tree(decon5.ps.p.filt) = treeLs$decon5.tree; saveRDS(decon5.ps.p.filt, "decon5.ps.p.filt.rtree.RDS")
phy_tree(tsinf.ps.p.filt) = treeLs$tsinf.tree; saveRDS(tsinf.ps.p.filt, "tsinf.ps.p.filt.rtree.RDS")

psList = list(decon3 = decon3.ps.p.filt,
				decon4 = decon4.ps.p.filt,
				decon5 = decon5.ps.p.filt,
				tsinf = tsinf.ps.p.filt 
			)

for( i in seq_along(psList) ) {
	ps = psList[[i]]
	ps.top = prune_taxa(taxa_names(ps)[1:100], ps)
	f = paste0(names(psList)[i],"top100.tree.pdf")
	pdf(f)
	print( plot_tree(ps.top, nodelabf=nodeplotboot(), ladderize="left", color="Phylum") )
	dev.off()
}