#!/usr/bin/Rscript

## automating the standard DADA2 16S pipeline from:
## 		https://benjjneb.github.io/dada2/tutorial.html

## suggested command to run: 
##		./DADA2_16Sv4.R > DADA2_16Sv4.R.out 2>&1



##########################################################################################################################################
######### SET THESE VARIABLES TO YOUR SPECIFICS ##########################################################################################
##########################################################################################################################################

# rm(list=ls())					                        # clear environment variables if needed

# add required libraries here
libs = c( "dada2", 
          "phyloseq", 
          "ShortRead", 
          "Biostrings", 
          "ggplot2", 
          "tidyverse" )	

# for getting input fastq files
fqpath="./fastq"	# folder containing fastq files
	## ^don't put trailing '/' for this path!
R1.pattern="_R1.fastq.gz"		                    # unique pattern to match your forward reads
R2.pattern="_R2.fastq.gz"		                    # unique pattern to match your forward reads

# argument to set threads for cmds that support it
n.threads=20

# path to database for taxa calling
# SILVA dada2 formatted databases
togen.silva.db = "./tax_dbs/silva/v138.2/silva_nr99_v138.2_toGenus_trainset.fa.gz"
tospe.silva.db = "./tax_dbs/silva/v138.2/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
asspe.silva.db = "./tax_dbs/silva/v138.2/silva_v138.2_assignSpecies.fa.gz"

# directories for output files
work.dir = "."
pdf.dir = file.path(work.dir,"pdf")
rds.dir = file.path(work.dir,"rds")
csv.dir = file.path(work.dir,"csv")
fna.dir = file.path(work.dir,"fasta")
rdat.dir = file.path(work.dir,"rdata")
tree.dir = file.path(work.dir,"trees")
for( d in c(pdf.dir,rds.dir,csv.dir,fna.dir,rdat.dir,tree.dir) ){
  if(!dir.exists(d))
    dir.create(d)
}

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



## ---- get input files setup --------------------------------------------------------------------------------------------------------------

fnFs = sort(list.files(fqpath, pattern = R1.pattern, full.names = TRUE))
fnRs = sort(list.files(fqpath, pattern = R2.pattern, full.names = TRUE))



## ---- get sample names -------------------------------------------------------------------------------------------------------------------

# Extract sample names, assuming filenames are formatted with sample name at very beginning followed by "_"
# 		ex: 	samplename_blah_blah_R1.fastq.gz
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))

# make sure looks hunky dory
cat(">>> Sample Names, Input R1, Input R2:\n")
head(cbind(sample.names, fnFs, fnRs))



## ---- raw qual plots -----------------------------------------------------------------------------------------------------------

# forward reads
f = file.path(pdf.dir, "fnFs.plotQualProf.pdf")
if( !file.exists(f) ) {
	cat(">>> Creating", f, "with plotQualityProfile() ... ... ...\n")
	pdf(f)
	lapply(fnFs, function(x) { print(plotQualityProfile(x)) } )
	dev.off()
}
# reverse reads
f = file.path(pdf.dir,"fnRs.plotQualProf.pdf")
if( !file.exists(f) ) {
	cat(">>> Creating", f, "with plotQualityProfile() ... ... ...\n")
	pdf(f)
	lapply(fnRs, function(x) { print(plotQualityProfile(x)) } )
	dev.off()
}



## ---- filter and trim --------------------------------------------------------------------------------------------------------

fqpath.filt <- file.path(fqpath, "filtered")
cat(">>> Path for filterAndTrim() output:", fqpath.filt, "... ... ...\n")

filtFs <- file.path(fqpath.filt, gsub("\\.fastq\\.gz$", "_filt.fastq.gz", basename(fnFs))) 
filtRs <- file.path(fqpath.filt, gsub("\\.fastq\\.gz$", "_filt.fastq.gz", basename(fnRs))) 
names(filtFs) = sample.names
names(filtRs) = sample.names

fo = "filterAndTrim.out"; ff = file.path(rds.dir,"filtFs.RDS"); fr = file.path(rds.dir,"filtRs.RDS")
if( file.exists(fo) & file.exists(ff) & file.exists(fr) ) {
	cat(">>> Reading table", fo, "... ... ...\n")
	out = read.table(fo, header = TRUE)
	cat(">>> Reading RDS files", ff, "and", fr,"... ... ...\n")
	filtFs = readRDS(ff)
	filtRs = readRDS(fr)
} else {
	cat(">>> Running filterAndTrim() ... ... ...\n")
	out <- 
		filterAndTrim(	
			fnFs, filtFs, 
			fnRs, filtRs,
			truncLen = c(225, 200), 
			maxN = 0, 
			maxEE = c(2, 2),
			truncQ = 2,
			rm.phix = TRUE,
			multithread = TRUE,
			compress = TRUE
		)
	cat(">>> Writing out to", fo, "... ... ...\n")
	write.table(out, fo)

	# remove any filenames that no longer have any reads
	fe = file.exists(filtFs)	# keep these files with >0 reads
	cat(">>> Samples with 0 reads after filterAndTrim removed:\n", 
		names(filtFs[which(!fe)]), 
		"\n"
	)
	cat(">>> Writing RDS files", ff, "and", fr, "... ... ...\n")
	filtFs = filtFs[fe]; saveRDS(filtFs, ff)
	filtRs = filtRs[fe]; saveRDS(filtRs, fr)
	sample.names = names(filtFs)	# update sample.names to remove any that had 0 reads
}

## ---- filtered qual plots -----------------------------------------------------------------------------------------------------------

# forward reads
f = file.path(pdf.dir, "filtFs.plotQualProf.pdf")
if( !file.exists(f) ) {
	cat(">>> Creating", f, "with plotQualityProfile() ... ... ...\n")
	pdf(f)
	lapply(filtFs, function(x) { print(plotQualityProfile(x)) } )
	dev.off()
}
# reverse reads
f = file.path(pdf.dir,"filtRs.plotQualProf.pdf")
if( !file.exists(f) ) {
	cat(">>> Creating", f, "with plotQualityProfile() ... ... ...\n")
	pdf(f)
	lapply(filtRs, function(x) { print(plotQualityProfile(x)) } )
	dev.off()
}

## ----learn error rates and plot---------------------------------------------------------------------------------------------

# setting randomize=T in learnErrors() bc i dont feel like the first 10 samples in order necessarily represent all the data

# forward reads
f = file.path(rds.dir,"errF.RDS")
if( file.exists(f) ) {
	cat(">>> Reading RDS file", f, "... ... ...\n")
	errF = readRDS(f)
} else {
	cat(">>> Creating object", f, "with learnErrors() ... ... ...\n")
	errF <- learnErrors(filtFs, multithread = TRUE, randomize = TRUE)
	saveRDS(errF, f)
}

# reverse reads
f = file.path(rds.dir,"errR.RDS")
if( file.exists(f) ) {
	cat(">>> Reading RDS file", f, "... ... ...\n")
	errR = readRDS(f)
} else {
	cat(">>> Creating object", f, "with learnErrors() ... ... ...\n")
	errR <- learnErrors(filtRs, multithread = TRUE, randomize = TRUE)
	saveRDS(errR, f)
}

# can look at PDf if don't want to rerun errF/R
f = file.path(pdf.dir,"filt.plotErrors.pdf")
if( !file.exists(f) ) {
	cat(">>> Creating file", f, "with plotErrors() ... ... ...\n")
	pdf(f)
	print(plotErrors(errF, nominalQ = TRUE))
	print(plotErrors(errR, nominalQ = TRUE))
	dev.off()
}



## ----derep identical reads--------------------------------------------------------------------------------------------------

# this step is not included in the v1.16 dada2 16S tutorial... see this: https://github.com/benjjneb/dada2/issues/1095

# forward reads
f = file.path(rds.dir,"derepFs.RDS")
if( file.exists(f) ) {
	cat(">>> Reading RDS file", f, "... ... ...\n")
	derepFs = readRDS(f)
} else {
	cat(">>> Creating object", f, "with derepFastq() ... ... ...\n")
	derepFs <- derepFastq(filtFs, verbose = FALSE) 
	saveRDS(derepFs, f)
}

# reverse
f = file.path(rds.dir,"derepRs.RDS")
if( file.exists(f) ) {
	cat(">>> Reading RDS file", f, "... ... ...\n")
	derepRs = readRDS(f)
} else {
	cat(">>> Creating object", f, "with derepFastq() ... ... ...\n")
	derepRs <- derepFastq(filtRs, verbose = FALSE) 
	saveRDS(derepRs, f)
}



## ----sample inference (dada function!)----------------------------------------------------------------------------------------

## unpooled -------

# forward reads
f = file.path(rds.dir,"dadaFs.RDS")
if( file.exists(f) ){
	cat(">>> Reading RDS file", f, "... ... ...\n")
	dadaFs = readRDS(f)
} else {
	cat(">>> Creating object", f, "with dada() ... ... ...\n")
	dadaFs <- dada(derepFs, err = errF, multithread = TRUE, verbose = FALSE)
	saveRDS(dadaFs, f)
}

# reverse reads
f = file.path(rds.dir,"dadaRs.RDS")
if( file.exists(f) ){
	cat(">>> Reading RDS file", f, "... ... ...\n")
	dadaRs = readRDS(f)
} else {
	cat(">>> Creating object", f, "with dada() ... ... ...\n")
	dadaRs <- dada(derepRs, err = errR, multithread = TRUE, verbose = FALSE)
	saveRDS(dadaRs, f)
}

## pooled -------

# forward reads
f = file.path(rds.dir,"dadaFs.p.RDS")
if( file.exists(f) ){
	cat(">>> Reading RDS file", f, "... ... ...\n")
	dadaFs.p = readRDS(f)
} else {
	cat(">>> Creating object", f, "with dada() ... ... ...\n")
	dadaFs.p <- dada(derepFs, err = errF, pool = TRUE, multithread = TRUE, verbose = FALSE)
	saveRDS(dadaFs.p, f)
}

# reverse reads
f = file.path(rds.dir,"dadaRs.p.RDS")
if( file.exists(f) ){
	cat(">>> Reading RDS file", f, "... ... ...\n")
	dadaRs.p = readRDS(f)
} else {
	cat(">>> Creating object", f, "with dada() ... ... ...\n")
	dadaRs.p <- dada(derepRs, err = errR, pool = TRUE, multithread = TRUE, verbose = FALSE)
	saveRDS(dadaRs.p, f)
}



## ----merge paired reads -----------------------------------------------------------------------------------------------------

## unpooled -------
f = file.path(rds.dir,"mergers.RDS")
if( file.exists(f) ) {
	cat(">>> Reading RDS file", f, "... ... ...\n")
	mergers = readRDS(f)
} else {
	cat(">>> Creating object", f, "with mergePairs() ... ... ...\n")
	mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = FALSE)
	saveRDS(mergers, f)
}
# pooled
f = file.path(rds.dir,"mergers.p.RDS")
if( file.exists(f) ) {
	cat(">>> Reading RDS file", f, "... ... ...\n")
	mergers.p = readRDS(f)
} else {
	cat(">>> Creating object", f, "with mergePairs() ... ... ...\n")
	mergers.p <- mergePairs(dadaFs.p, derepFs, dadaRs.p, derepRs, verbose = FALSE)
	saveRDS(mergers.p, f)
}



## ----construct sequence table-----------------------------------------------------------------------------------------------------------
min.length=250; max.length=256	

## unpooled -------
f = file.path(rds.dir,"seqtab.RDS")
f.untrim = file.path(rds.dir,"seqtab.untrim.RDS")

if( file.exists(f) & file.exists(f.untrim) ) {
	cat(">>> Reading RDS files", f, "and", f.untrim, "... ... ...\n")
	seqtab = readRDS(f)
	seqtab.untrim = readRDS(f.untrim)
} else {
	cat(">>> Making objects", f, "and", f.untrim, "with makeSequenceTable() ... ... ...\n")
	seqtab <- makeSequenceTable(mergers)
	rownames(seqtab) = sample.names
	
	seqtab.untrim <- seqtab
	saveRDS(seqtab.untrim, f.untrim)

	# prevent vector return when only 1 sample
	if(length(sample.names) == 1) {	
		seqtab <- t(as.matrix(seqtab[ , nchar(colnames(seqtab)) %in% min.length:max.length]))
	} else {
		seqtab <- seqtab[ , nchar(colnames(seqtab)) %in% min.length:max.length]
	}
	saveRDS(seqtab, f)
}

cat(">>> seqtab.untrim dimensions:\n")
dim(seqtab.untrim)
cat(">>> seqtab.untrim seq lengths:\n")
table(nchar(getSequences(seqtab.untrim)))

cat(">>> seqtab trimmed to seq lengths ", min.length, "-", max.length, "... ... ...\n")
cat(">>> seqtab dimensions:\n")
dim(seqtab)
cat(">>> seqtab seq lengths:\n")
table(nchar(getSequences(seqtab)))

## pooled -------
f = file.path(rds.dir,"seqtab.p.RDS")
f.untrim = file.path(rds.dir,"seqtab.p.untrim.RDS")

if(  file.exists(f) & file.exists(f.untrim) ) {
	cat(">>> Reading RDS files", f, "and", f.untrim, "... ... ...\n")
	seqtab.p = readRDS(f)
	seqtab.p.untrim = readRDS(f.untrim)
} else {
	cat(">>> Making objects", f, "and", f.untrim, "with makeSequenceTable() ... ... ...\n")
	seqtab.p <- makeSequenceTable(mergers.p)
	rownames(seqtab.p) = sample.names
	
	seqtab.p.untrim <- seqtab.p
	saveRDS(seqtab.p.untrim, f.untrim)

	# prevent vector return when only 1 sample
	if(length(sample.names) == 1) {	
		seqtab.p <- t(as.matrix(seqtab.p[ , nchar(colnames(seqtab.p)) %in% min.length:max.length]))
	} else {
		seqtab.p <- seqtab.p[ , nchar(colnames(seqtab.p)) %in% min.length:max.length]
	}
	saveRDS(seqtab.p, f)
}

cat(">>> seqtab.p.untrim dimensions:\n")
dim(seqtab.p.untrim)
cat(">>> seqtab.p.untrim seq lengths:\n")
table(nchar(getSequences(seqtab.p.untrim)))

cat(">>> seqtab.p trimmed to seq lengths", min.length, "-", max.length, "... ... ...\n")
cat(">>> seqtab.p dimensions:\n")
dim(seqtab.p)
cat(">>> seqtab.p seq lengths:\n")
table(nchar(getSequences(seqtab.p)))



## ----remove chimeras--------------------------------------------------------------------------------------------------------------------

## unpooled -------

f = file.path(rds.dir,"seqtab.nochim.RDS")
if( file.exists(f) ){
	cat(">>> Reading RDS file", f, "... ... ...\n")
	seqtab.nochim <- readRDS(f)
} else {
	cat(">>> Creating object", f, "with removeBimeraDenovo() ... ... ...\n")
	seqtab.nochim <- 
		removeBimeraDenovo(
			seqtab, 
			method = "consensus", 
			multithread = TRUE, 
			verbose = FALSE
		)
		saveRDS(seqtab.nochim, f)
}

cat(">>> Percent reads retained as nonchimeric:\n")
sum(seqtab.nochim)/sum(seqtab)*100

cat(	">>> ASVs: ",				dim(seqtab)[2], 		"\t",
			"nonchimeric ASVs: ",	dim(seqtab.nochim)[2], 	"\n"
)
cat(">>> Percent ASVs retained as nonchimeric:\n")
dim(seqtab.nochim)[2]/dim(seqtab)[2]*100

cat(">>> Non-chimeric ASV seq lengths:\n")
table(nchar(getSequences(seqtab.nochim)))


## pooled -------

f = file.path(rds.dir,"seqtab.p.nochim.RDS")
if( file.exists(f) ) {
	cat(">>> Reading RDS file", f, "... ... ...\n")
	seqtab.p.nochim <- readRDS(f)
} else {
	cat(">>> Creating object", f, "with removeBimeraDenovo() ... ... ...\n")
	seqtab.p.nochim <- 
		removeBimeraDenovo(
			seqtab.p, 
			method = "consensus", 
			multithread = TRUE, 
			verbose = FALSE
		)
	saveRDS(seqtab.p.nochim, f)
}

cat(">>> Percent reads retained as .p.nonchimeric:\n")
sum(seqtab.p.nochim)/sum(seqtab.p)*100

cat(	">>> .p ASVs: ",			dim(seqtab.p)[2], 			"\t",
		".p.nonchimeric ASVs: ",	dim(seqtab.p.nochim)[2], 	"\n"
)

cat(">>> Percent .p ASVs retained as .p.nonchimeric:\n")
dim(seqtab.p.nochim)[2]/dim(seqtab.p)[2]*100

cat(">>> .p Non-chimeric ASV seq lengths:\n")
table(nchar(getSequences(seqtab.p.nochim)))



## ----track reads through pipeline and plot---------------------------------------------------------------------------------------

getN <- function(x) sum(getUniques(x))

## unpooled -------

f = file.path(csv.dir,"track.csv")
if( file.exists(f) ) {
	cat(">>> Reading file", f, "... ... ...\n")
	track = read.csv(f)
} else {
	cat(">>> Creating file", f, "... ... ...\n")
	out = as.data.frame(out)
	# check for samples that failed filterAndTrim() 
	if( nrow(out) != length(sample.names) ) {
		fat.fail = filter(out, reads.out == 0)		# save rows that failed
		out = filter(out, reads.out > 0)			# remove rows that failed
	}
	# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
	if(nrow(out) == 1) {
  		track <-cbind( out, 
                 getN(dadaFs), 
                 getN(dadaRs),
                 getN(mergers),
                 rowSums(seqtab.nochim)
  		)
	} else {
  		track <-cbind( out, 
                 sapply(dadaFs, getN), 
                 sapply(dadaRs, getN),
                 sapply(mergers, getN),
                 rowSums(seqtab.nochim)
  		)
	}
	colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
	rownames(track) <- sample.names
	# add back the samples that failed filterAndTrim()
	colnames(fat.fail) = colnames(track)[1:2]
	rownames(fat.fail) = get.sample.name(rownames(fat.fail))
	fat.fail <- 
		tibble::add_column( 
			fat.fail, 
			"denoisedF" = 0, 
			"denoisedR" = 0, 
			"merged" = 0, 
			"nonchim" = 0 
		)
	track = rbind(track, fat.fail)

	# output csv
	cat(">>> Writing table", f, "... ... ...\n")
	write.csv(track, file = f)
}

# plot to pdf
f = file.path(pdf.dir,"track.barplot.pdf")
if( !file.exists(f) ) {
	cat(">>> Creating file", f, "with barplot()... ... ...\n")
	pdf(f, width = 12, height = 6)
	barplot( 
		t(track),
		beside = TRUE, legend.text = TRUE, las = 2, cex.names = 1,
		col = c("black", "red", "blue", "yellow", "purple", "green"),
		space = c(0.2, 2), ylim = c(0, max( t(track) )),
		args.legend = list(cex = 0.6)
	)
	dev.off()
}

## pooled -------

f = file.path(csv.dir,"track.p.csv")
if( file.exists(f) ) {
	cat(">>> Reading file", f, "... ... ...\n")
	track.p = read.csv(f)
} else {
	cat(">>> Creating file", f, "... ... ...\n")
	out.p = as.data.frame(out)

	# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
	if(nrow(out.p) == 1) {
	  track.p <- 
	  	cbind( 
	  		out.p, 
	      getN(dadaFs.p), 
	      getN(dadaRs.p),
	      getN(mergers.p),
	      rowSums(seqtab.p.nochim)
	  )
	} else {
	  track.p <-
	  	cbind( 
	  		out.p, 
	      sapply(dadaFs.p, getN), 
	      sapply(dadaRs.p, getN),
	      sapply(mergers.p, getN),
	      rowSums(seqtab.p.nochim)
	  )
	}
	colnames(track.p) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
	rownames(track.p) <- sample.names
	# add samples that failed filterAndTrim() ## should be same as for unpooled
	track.p = rbind(track.p, fat.fail)

	# output csv 
	cat(">>> Writing table", f, "... ... ...\n")
	write.csv(track.p, file = f)
}

# plot to pdf
f = file.path(pdf.dir,"track.p.barplot.pdf")
if( !file.exists(f) ) {
	cat(">>> Creating file", f, "with barplot() ... ... ...\n")
	pdf("track.p.barplot.pdf", width = 12, height = 6)
	barplot( 
		t(track.p),
		beside = TRUE, legend.text = TRUE, las = 2, cex.names = 1,
		col=c("black","red","blue","yellow","purple","green"),
		space = c(0.2, 2), ylim = c(0, max( t(track.p) )),
		args.legend = list(cex = 0.6)
	)
	dev.off()
}



## ---- dada2::assignTaxonomy() --------------------------------------------------------------------------------------------------------

## UNpooled ------- 
f = file.path(rds.dir,"taxa.RDS")
if( file.exists(f) ) {
	cat(">>> Reading file", f, "... ... ...\n")
	taxa = readRDS(f)
} else {
	cat(">>> Creating file", f, "with dada2::assignTaxonomy() and addSpecies() ... ... ...\n")
	taxa = assignTaxonomy(seqtab.nochim, togen.silva.db, multithread = TRUE)
	taxa = addSpecies(taxa, asspe.silva.db)
	saveRDS(taxa, f)
}

## pooled ------- 
f = file.path(rds.dir,"taxa.p.RDS")
if( file.exists(f) ) {
	cat(">>> Reading file", f, "... ... ...\n")
	taxa.p = readRDS(f)
} else {
	cat(">>> Creating file", f, "with dada2::assignTaxonomy() and addSpecies() ... ... ...\n")
	taxa.p = assignTaxonomy(seqtab.p.nochim, togen.silva.db, multithread = TRUE)
	taxa.p = addSpecies(taxa.p, asspe.silva.db)
	saveRDS(taxa.p, f)
}

## inspect UNpooled -------
taxa.print <- taxa				# Removing sequence rownames for display only
rownames(taxa.print) <- NULL
cat(">>> Inspect UNpooled taxa -- head taxa.print:\n")
head(taxa.print)  						# looks like bacteria

## inspect pooled -------
taxa.p.print = taxa.p
rownames(taxa.p.print) = NULL
cat(">>> Inspect pooled taxa -- head taxa.p.print:\n")
head(taxa.p.print)



## ---- phyloseq -------------------------------------------------------------------------------------------------------------------------------------------- ##

## UNpooled -------

f = file.path(rds.dir,"ps.RDS")
if( file.exists(f) ){
	cat(">>> Reading file", f, "... ... ...\n")
	ps = readRDS(f)
} else {
	cat(">>> Creating object", f, "with phyloseq() ... ... ...\n")
	ps = phyloseq( 
	  		otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
	  		tax_table(taxa) 
	  	)
	cat(">>> Adding refseq slot to ps and renaming ASVs by number ... ... ...\n")
	sequences <- Biostrings::DNAStringSet(taxa_names(ps))
	names(sequences) <- taxa_names(ps)
	ps <- merge_phyloseq(ps, sequences)
	taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
	
	cat(">>> Writing file", f, "... ... ...\n")
	saveRDS(ps, f)
}

# if you want to write these tables to read into excel
f = file.path(csv.dir,"ps.otutable.csv")
if( !file.exists(f) ) {
	cat(">>> Writing count table", f, " ... ... ...\n")
	write.csv(seqtab.nochim, file = f)
}

f = file.path(csv.dir,"ps.taxtable.csv")
if( !file.exists(f) ) {
	cat(">>> Writing taxa table", f, " ... ... ...\n")
	write.csv(taxa, file = f)
}

f = file.path(fna.dir,"ps.ASV.fasta")
if( !file.exists(f) ) {
	cat(">>> Writing ASV fasta", f, " ... ... ...\n")
	Biostrings::writeXStringSet(refseq(ps), f)
}



## pooled -------

f = file.path(rds.dir,"ps.p.RDS")
if( file.exists(f) ) {
	cat(">>> Reading file", f, "... ... ...\n")
	ps.p = readRDS(f)
} else {
		cat(">>> Creating object", f, "with phyloseq() ... ... ...\n")
		ps.p = phyloseq(
					otu_table(seqtab.p.nochim, taxa_are_rows = FALSE), 
					tax_table(taxa.p)
  				)
		cat(">>> Adding refseq slot to ps.p and renaming ASVs by number ... ... ...\n")
		sequences <- Biostrings::DNAStringSet(taxa_names(ps.p))
		names(sequences) <- taxa_names(ps.p)
		ps.p <- merge_phyloseq(ps.p, sequences)
		taxa_names(ps.p) <- paste0("ASV", seq(ntaxa(ps.p)))
		
		cat(">>> Writing file", f, "... ... ...\n")
		saveRDS(ps.p, f)
}

# write non-phyloseq object tables to read into other stuff
f = file.path(csv.dir,"ps.p.otutable.csv")
if( !file.exists(f) ) {
	cat(">>> Writing count table", f, " ... ... ...\n")
	write.csv(seqtab.p.nochim, file = f)
}

f = file.path(csv.dir,"ps.p.taxtable.csv")
if( !file.exists(f) ) {
	cat(">>> Writing taxa table", f, " ... ... ...\n")
	write.csv(taxa.p, file = f)
}

f = file.path(fna.dir,"ps.p.ASV.fasta")
if( !file.exists(f) ) {
	cat(">>> Writing .p ASV fasta", f, " ... ... ...\n")
	Biostrings::writeXStringSet(refseq(ps.p), f)
}



## ----plot-----------------------------------------------------------------------------------------------------------------------------

theme_set(theme_bw())

## UNpooled -------

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

# plot just top 20 ASVs
top20 = names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]
ps.prop.top20 = prune_taxa(top20,ps.prop)

f = file.path(pdf.dir,"ps.top20__gen_stackedbar.pdf")
if( !file.exists(f) ) {
	cat(">>> Creating file", f, "with plot_bar() ... ... ...\n")
	pdf(f, width = 12, height = 6)
	print(	plot_bar(	
				ps.prop.top20, 
				fill = "Genus", 
				title = "Top 20 taxa"
			) + 
			theme( axis.text = element_text(size = 14))
		)
	dev.off()
}

# color by family
f = file.path(pdf.dir,"ps.top20__fam_stackedbar.pdf")
if( !file.exists(f) ) {
	cat(">>> Creating file", f, "with plot_bar() ... ... ...\n")
	pdf(f, width = 12, height = 6)
	print(	plot_bar(	
				ps.prop.top20, 
				fill = "Family", 
				title = "Top 20 taxa"
			) + 
			theme(axis.text = element_text(size = 14))
		)
	dev.off()
}

## pooled -------

ps.p.prop <- transform_sample_counts(ps.p, function(otu) otu/sum(otu))

# plot just top 20 ASVs
top20.p = names(sort(taxa_sums(ps.p), decreasing = TRUE))[1:20]
ps.p.prop.top20 = prune_taxa(top20.p, ps.p.prop)

# color by genus
f = file.path(pdf.dir,"ps.ptop20__gen_stackedbar.pdf")
if( !file.exists(f) ) {
	cat(">>> Creating file", f, "with plot_bar() ... ... ...\n")
	pdf(f, width = 12, height = 6)
	print(	plot_bar(	
				ps.p.prop.top20, 
				fill = "Genus", 
				title = "Top 20 taxa"
			) + 
			theme(axis.text = element_text(size = 14))
		)
	dev.off()
}

# color by family
f = file.path(pdf.dir,"ps.p.top20__fam_stackedbar.pdf")
if( !file.exists(f) ) {
	cat(">>> Creating file", f, "with plot_bar() ... ... ...\n")
	pdf(f, width = 12, height = 6)
	print(	plot_bar(	
				ps.p.prop.top20, 
				fill = "Family", 
				title = "Top 20 taxa"
			) + 
			theme(axis.text = element_text(size = 14))
		)
	dev.off()
}



## ----save image-----------------------------------------------------------------------------------------------------------------------------
cat(">>> sessionInfo() ... ... ...\n")
sessionInfo()

image.file = file.path(rdat.dir, format(Sys.time(), "%b%d%Y_%X%Z"),".RData")
cat(">>> Saving image file", image.file, "... ... ...\n")
save.image(image.file)
cat(">>> COMPLETED\n")

