## kkyle@KlassenLab-server:~/kek16SallITStrachyDecons/fastq/p1p2_16S/
## R v4.4.1

##########################################################################################################################################
######### SET THESE VARIABLES TO YOUR SPECIFICS ##########################################################################################
##########################################################################################################################################

n.threads = 20
libs = c("phyloseq", "tidyverse", "phytools")

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

## ---- load filtered phyloseq object ----------------------------------------------------------------------------------- ##

ps.filt = readRDS("dada2/rds/ps.p.filt.RDS")



## ---- add sample data ------------------------------------------------------------------------------------------------- ##

decons.nonc.samdat = read.csv("../decons_master_sampledata.csv")  %>% filter(!grepl(pattern="NC|NFW|zymo",seq.id))
tsinf.nonc.samdat = read.csv("../TSinf_master_sampledata.csv")  %>% filter(!grepl(pattern="NC|NFW|zymo",seq.id))
zymo.nonc.samdat = read.csv("../decons_master_sampledata.csv")  %>% filter(grepl(pattern="zymo",seq.id))



## ---- make experiment specific objects

decons.sn.l = grepl(sample_names(ps.filt), pattern="^d3|^d4|^d5")
decons.ps.filt = subset_samples(ps.filt, decons.sn.l)
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2068 taxa and 121 samples ]
# sample_data() Sample Data:       [ 121 samples by 3 sample variables ]
# tax_table()   Taxonomy Table:    [ 2068 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 2068 tips and 2067 internal nodes ]
# refseq()      DNAStringSet:      [ 2068 reference sequences ]

tsinf.sn.l = grepl(sample_names(ps.filt), pattern="^TS")
tsinf.ps.filt = subset_samples(ps.filt, tsinf.sn.l)
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2068 taxa and 42 samples ]
# sample_data() Sample Data:       [ 42 samples by 3 sample variables ]
# tax_table()   Taxonomy Table:    [ 2068 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 2068 tips and 2067 internal nodes ]
# refseq()      DNAStringSet:      [ 2068 reference sequences ]

zymo.sn.l = grepl(sample_names(ps.filt), pattern="^zymo")
zymo.ps.filt = subset_samples(ps.filt, zymo.sn.l)
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2068 taxa and 3 samples ]
# sample_data() Sample Data:       [ 3 samples by 3 sample variables ]
# tax_table()   Taxonomy Table:    [ 2068 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 2068 tips and 2067 internal nodes ]
# refseq()      DNAStringSet:      [ 2068 reference sequences ]

psList = list(
		decons = decons.ps.filt,
		tsinf = tsinf.ps.filt,
		zymo = zymo.ps.filt
	)
sdList = list(
		decons = decons.nonc.samdat,
		tsinf = tsinf.nonc.samdat,
		zymo = zymo.nonc.samdat
	)

# check if sample names match between ps and samdat
for( i in seq_along(psList) ){
	if(  all( sdList[[i]]$seq.id %in% sample_names(psList[[i]]) ) ) {
		cat(">>> sample names for", 
			names(psList)[i], 
			"phyloseq object are the same as seq.id in", 
			names(sdList)[i], 
			"samdat.\n"
		)
		
		rownames(sdList[[i]]) = sdList[[i]]$seq.id
		sample_data(psList[[i]]) = sdList[[i]]
		
		cat(">>> sample_data slot filled for phyloseq object named", 
			names(psList)[i], 
			".\n"
		)
		
		f = paste0(names(psList)[i], ".ps.filt.RDS")
		saveRDS(psList[[i]], f)
		cat(">>> Saving file", f,"... ... ...\n")

	} else {
		cat(">>> sample names for", 
			names(psList)[i],
			"phyloseq object do not match seq.id's for", 
			names(sdList)[i],
			"samdat! Please reorder.\n"
		)
	}
}

