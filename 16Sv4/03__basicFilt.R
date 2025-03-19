## kkyle@KlassenLab-server:~/kek16SallITStrachyDecons/fastq/p1p2_16S/dada2
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


## ---- filtering taxa ------------------------------------------------------------------- ##
f = file.path(work.dir,"decontam","ps.p.decontam.prev0.5.RDS")
ps = readRDS(f)
ps.nonc = prune_samples(
	!(sample_names(ps) %>% grepl(pattern = "NC|NFW")), 
	ps
)

# subset_taxa() removes instances of NA if part of match
ps.nonc = na_to_unclassified_taxa(ps.nonc)

# starting with 2165 taxa and 166 samples
ps.filt = subset_taxa(ps.nonc, Kingdom == "Bacteria")					# 2158 taxa
ps.filt = subset_taxa(ps.filt, Phylum != "unclass. Bacteria")	# 2150 taxa
ps.filt = subset_taxa(ps.filt, Order != "Chloroplast")				# 2105 taxa 
ps.filt = subset_taxa(ps.filt, Family != "Mitochondria")			# 2068 taxa
saveRDS(ps.filt, "./rds/ps.p.filt.RDS")

# make a table of reads lost
filt.df <- as.data.frame( 
	list(input = 
			 	ps.nonc %>% sample_sums(),
			 k.bact = 
			 	ps.nonc %>% subset_taxa(Kingdom == "Bacteria") %>% 
			 	sample_sums(),
			 p.nna = 
			 	ps.nonc %>% subset_taxa(Kingdom == "Bacteria") %>% 
			 	subset_taxa(Phylum != "unclass. Bacteria") %>% 
			 	sample_sums(),
			 o.nchl = 
			 	ps.nonc %>% subset_taxa(Kingdom == "Bacteria") %>% 
			 	subset_taxa(Phylum != "unclass. Bacteria") %>% 
			 	subset_taxa(Order != "Chloroplast") %>%
			 	sample_sums(),
			 f.nmit = 
			 	ps.nonc %>% subset_taxa(Kingdom == "Bacteria") %>% 
			 	subset_taxa(Phylum != "unclass. Bacteria") %>% 
			 	subset_taxa(Order != "Chloroplast") %>%
			 	subset_taxa(Family != "Mitochondria") %>%
			 	sample_sums()
) )
write.csv(filt.df, "./csv/filt.df.csv")


## ---- track ASVs through pipeline  ---------------------------------------------------------------------- ##
# UNpooled
	seqtab.p = as.data.frame(readRDS("./rds/seqtab.p.RDS"))
	colnames(seqtab.p) = seq(1,ncol(seqtab.p))
# pooled
	seqtab.nochim.p = as.data.frame(readRDS("./rds/seqtab.p.nochim.RDS"))
	colnames(seqtab.nochim.p) = seq(1,ncol(seqtab.nochim.p))

# count number of ASVs per sample
asv.track = list()
asv.track$merged= rowSums(seqtab.p != 0)
asv.track$nochim= rowSums(seqtab.nochim.p != 0)
asv.track$decontam.prev0.1 = rowSums( otu_table(ps.prev0.1) != 0 )
asv.track$decontam.prev0.25 = rowSums( otu_table(ps.prev0.25) != 0 )
asv.track$decontam.prev0.5 = rowSums( otu_table(ps) != 0 )			# ps = ps.prev0.5
ps = na_to_unclassified_taxa(ps)
asv.track$kng.bact = 
	rowSums(otu_table(ps %>% 
										subset_taxa(Kingdom == "Bacteria") ) != 0 )
asv.track$phy.nona = 
	rowSums(otu_table( ps %>% 
										subset_taxa(Kingdom == "Bacteria") %>% 
										subset_taxa(Phylum != "unclass. Bacteria") ) != 0 )
asv.track$ord.nchl = 
	rowSums(otu_table( ps %>% 
										subset_taxa(Kingdom == "Bacteria") %>% 
										subset_taxa(Phylum != "unclass. Bacteria") %>%
										subset_taxa(Order != "Chloroplast") ) != 0 )
asv.track$fam.nmit = 
	rowSums(otu_table( ps %>% 
										subset_taxa(Kingdom == "Bacteria") %>% 
										subset_taxa(Phylum != "unclass. Bacteria") %>%
										subset_taxa(Order != "Chloroplast") %>%
										subset_taxa(Family != "Mitochondria") ) != 0 )

write.csv( as.data.frame(asv.track), 
					"./csv/asv.track.csv",
					quote=F,row.names=T )
