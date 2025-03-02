#!/usr/bin/Rscript

## automating the DADA2 ITS2 pipeline from:
## 		https://benjjneb.github.io/dada2/ITS_workflow.html
## includes first trimming the reads and running ITSx on ASVs

## suggested command(s) to run: 
##    ./00_fastp.bash
##		./DADA2_ITS2.R > DADA2_ITS2.R.out 2>&1



##########################################################################################################################################
######### SET THESE VARIABLES TO YOUR SPECIFICS ##########################################################################################
##########################################################################################################################################

# add required libraries here
libs = c( "dada2", 
          "phyloseq", 
          "ShortRead", 
          "Biostrings", 
          "ggplot2", 
          "tidyverse",
          "devtools" )	

# argument to set threads for cmds that support it
n.threads=20

# folder containing fastq files
fqpath="../fastp/fastp.fq"
      ## ^don't put trailing '/' for this path!

# my files are named like: 
#     prefix__sid_R1.fastq.gz  prefix__sid_R2.fastq.gz
# zipped or unzipped seems to both work fine
R1.pattern="_R1.fastq.gz"		                    # unique pattern to match your forward reads
R2.pattern="_R2.fastq.gz"		                    # unique pattern to match your reverse reads

# Read in input files
fnFs = sort(list.files(fqpath, pattern=R1.pattern, full.names=T))	
fnRs = sort(list.files(fqpath, pattern=R2.pattern, full.names=T))
# process file names
sample.names = sapply(strsplit(basename(fnFs),"_"),`[`,3)			# changed to 3rd substring for my files
preprocc.str = sapply(strsplit(basename(fnFs),"_"),`[`,1)     # keep the prefix of the files for record keeping

# primer sequences for cutadapt 
  # fITS7: GTGARTCATCGAATCTTTG
  # ITS4: TCCTCCGCTTATTGATATGC
# i think these MARS versions are 10bp pad, 2bp linker, then fPrim
FWD <- "GCCGGCTGCGACGTGARTCATCGAATCTTTG"  		# MARS version
REV <- "AGGCAGTCAGCCTCCTCCGCTTATTGATATGC"     # MARS version

# conda environment names for external programs
# will be used with conda run like so:
#     conda run -n $conda.n ...
cutadapt.conda.n = "cutadapt" 
itsx.conda.n = "itsx"

# path to database for taxa calling
base.dir = "/home/kkyle/kek16SallITStrachyDecons/fastq/its2/dada2/"
unite.ref = file.path(base.dir, "sh_general_release_dynamic_04.04.2024.fasta")

# directories for output files
work.dir = "."
pdf.dir = file.path(work.dir,"pdf")
rds.dir = file.path(work.dir,"rds")
csv.dir = file.path(work.dir,"csv")
fna.dir = file.path(work.dir,"fasta")
rdat.dir = file.path(work.dir,"rdata")
itsx.dir = file.path(work.dir, "itsx")
for( d in c(pdf.dir,rds.dir,csv.dir,fna.dir,rdat.dir,tree.dir,itsx.dir) ){
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

# load some custom functions
source_url("https://raw.githubusercontent.com/kek12e/my_r_functions/refs/heads/main/my_r_functions.R")



## ---- get sample names -------------------------------------------------------------------------------------------------------------------

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



## ---- (ITS specific) check for and remove primer readthrough -----------------------------------------------------------------

# ----funcs/----
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  # The Biostrings works w/ DNAString objects rather than character vectors
  dna <- DNAString(primer)  
  orients <- c(Forward = dna, 
               Complement = Biostrings::complement(dna), 
               Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna) )
  return(sapply(orients, toString))  # Convert back to character vector
}

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
# ---\funcs----

# primer seq variations
FWD.orients = allOrients(FWD)
REV.orients = allOrients(REV)
FWD.RC = dada2:::rc(FWD) # reverse-complement
REV.RC = dada2:::rc(REV) # reverse-complement

# Put N-filtered files in filtN/ subdirectory
fnFs.filtN = 
  file.path( fqpath, "filtN", 
             paste0( preprocc.str, "_", sample.names, "_R1_filtN.fastq.gz") )
fnRs.filtN =
  file.path( fqpath, "filtN", 
             paste0( preprocc.str, "_", sample.names, "_R2_filtN.fastq.gz") )

# filter N's
if( !all(file.exists(fnFs.filtN,fnRs.filtN)) ) {
  cat(">>> Creating N-filtered files ... ... ...\n")
  filterAndTrim( fnFs, fnFs.filtN, 
                 fnRs, fnRs.filtN, 
                 maxN = 0, multithread = TRUE )
}

# check
rbind( FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]) )

# cutadapt prep
# output directory
path.cut = file.path(fqpath, "cutadapt")
cat(">>> Path for cutadapt output:", path.cut, "... ... ...\n")
if( !dir.exists(path.cut) )
  dir.create(path.cut)
# output files
fnFs.cut = 
  file.path( path.cut, 
             paste0(preprocc.str, "_", sample.names, "_R1_filtN_cut.fastq.gz"))
fnRs.cut = 
  file.path( path.cut, 
             paste0(preprocc.str, "_", sample.names, "_R2_filtN_cut.fastq.gz"))
# cutadapt flags
R1.flags <- paste("-g", FWD, "-a", REV.RC) # Trim FWD and REV.RC off of forward reads
R2.flags <- paste("-G", REV, "-A", FWD.RC) # Trim REV and FWD.RC off of reverse reads

if( !all(file.exists(fnFs.cut,fnRs.cut)) ) {
  cat(">>> Running cutadapt ... ... ...\n")
  # --- Run Cutadapt --- #
  for( i in seq_along(fnFs) ) {
    # make cutadapt command string
    cutadapt.cmd = 
      paste("cutadapt", R1.flags, R2.flags, 
            "-n", 2,                              # -n 2 required to remove FWD and REV from reads
            "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
            fnFs.filtN[i], fnRs.filtN[i]          # input files
      )
    # use conda run to run cutadapt
    system2("conda", 
            args = c("run -n", cutadapt.conda.n,
                     cutadapt.cmd
            ) )
  }
}

# sanity check
cat(">>> There should NOT be any primers now ... ... ...\n")
rbind( FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]) )



## ---- filter and trim --------------------------------------------------------------------------------------------------------

fqpath.filt = file.path(fqpath, "filtered")
cat(">>> Path for filterAndTrim() output:", fqpath.filt, "... ... ...\n")

filtFs = file.path(fqpath.filt, gsub("\\.fastq\\.gz$", "_filt.fastq.gz", basename(fnFs))) 
filtRs = file.path(fqpath.filt, gsub("\\.fastq\\.gz$", "_filt.fastq.gz", basename(fnRs))) 
names(filtFs) = sample.names
names(filtRs) = sample.names

fo = file.path(csv.dir,"filterAndTrim.csv")
ff = file.path(rds.dir,"filtFs.RDS")
fr = file.path(rds.dir,"filtRs.RDS")
if( file.exists(fo) & file.exists(ff) & file.exists(fr) ) {
  cat(">>> Reading table", fo, "... ... ...\n")
  out = read.table(fo, sep=",", row.names=1, header=T)
  cat(">>> Reading RDS files", ff, "and", fr,"... ... ...\n")
  filtFs = readRDS(ff)
  filtRs = readRDS(fr)
} else {
  cat(">>> Running filterAndTrim() ... ... ...\n")
  out <- 
    filterAndTrim(	
      fnFs, filtFs, 
      fnRs, filtRs,
      # truncLen = c(225, 200), 
      maxN = 0, 
      maxEE = c(2, 2),
      truncQ = 2,
      rm.phix = TRUE,
      multithread = TRUE,
      compress = TRUE
    )
  cat(">>> Writing out to", fo, "... ... ...\n")
  write.csv(out, fo)
  
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
# setting randomize=T in learnErrors() bc i dont feel like the first 10 samples 
# in order necessarily represent all the data

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
# this step is not included in the v1.16 dada2 16S tutorial... 
# see this: https://github.com/benjjneb/dada2/issues/1095

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

# reverse reads
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

## ---- unpooled -------
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

## ---- pooled -------
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

## ---- unpooled -------
  f = file.path(rds.dir,"mergers.RDS")
  if( file.exists(f) ) {
    cat(">>> Reading RDS file", f, "... ... ...\n")
    mergers = readRDS(f)
  } else {
    cat(">>> Creating object", f, "with mergePairs() ... ... ...\n")
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = FALSE)
    saveRDS(mergers, f)
  }
  
# ---- pooled ----------
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

## ---- unpooled -------
  f = file.path(rds.dir,"seqtab.RDS")
  if( file.exists(f) ) {
    cat(">>> Reading RDS file", f, "... ... ...\n")
    seqtab = readRDS(f)
  } else {
    cat(">>> Making object", f, "with makeSequenceTable() ... ... ...\n")
    seqtab <- makeSequenceTable(mergers)
    rownames(seqtab) = sample.names
    saveRDS(seqtab, f)
  }
  cat(">>> seqtab dimensions:\n")
  dim(seqtab)
  cat(">>> seqtab seq lengths:\n")
  table(nchar(getSequences(seqtab)))

## ---- pooled -------
  f = file.path(rds.dir,"seqtab.p.RDS")
  if(  file.exists(f) ) {
    cat(">>> Reading RDS file", f, "... ... ...\n")
    seqtab.p = readRDS(f)
  } else {
    cat(">>> Making object", f, "with makeSequenceTable() ... ... ...\n")
    seqtab.p <- makeSequenceTable(mergers.p)
    rownames(seqtab.p) = sample.names
    saveRDS(seqtab.p, f)
  }
  cat(">>> seqtab.p dimensions:\n")
  dim(seqtab.p)
  cat(">>> seqtab.p seq lengths:\n")
  table(nchar(getSequences(seqtab.p)))



## ----remove chimeras--------------------------------------------------------------------------------------------------------------------

## ---- unpooled -------
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

## ---- pooled -------
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

## ---- unpooled -------
  fat.fail = data.frame()
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
    if(length(fat.fail) > 0) {
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
    }

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

## ---- pooled -------
  f = file.path(csv.dir,"track.p.csv")
  if( file.exists(f) ) {
    cat(">>> Reading file", f, "... ... ...\n")
    track.p = read.csv(f)
  } else {
    cat(">>> Creating file", f, "... ... ...\n")
    out.p = as.data.frame(out)
    
    # If processing a single sample, remove the sapply calls: 
    #   e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
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
    colnames(track.p) <- c("input", "filtered", "denoisedF", 
                           "denoisedR", "merged", "nonchim")
    rownames(track.p) <- sample.names
    # add samples that failed filterAndTrim()
    if(length(fat.fail) > 0) {
      track.p = rbind(track.p, fat.fail)  ## should be same as for unpooled
    }
    
    # output csv 
    cat(">>> Writing table", f, "... ... ...\n")
    write.csv(track.p, file = f)
  }
  
  # plot to pdf
  f = file.path(pdf.dir,"track.p.barplot.pdf")
  if( !file.exists(f) ) {
    cat(">>> Creating file", f, "with barplot() ... ... ...\n")
    pdf(f, width = 12, height = 6)
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

## ---- UNpooled ------- 
  f = file.path(rds.dir,"taxa.RDS")
  if( file.exists(f) ) {
    cat(">>> Reading file", f, "... ... ...\n")
    taxa = readRDS(f)
  } else {
    cat(">>> Creating file", f, "with dada2::assignTaxonomy() ... ... ...\n")
    taxa = assignTaxonomy(seqtab.nochim, unite.ref, multithread=T, tryRC=T)
    saveRDS(taxa, f)
  }


## ---- pooled ------- 
  f = file.path(rds.dir,"taxa.p.RDS")
  if( file.exists(f) ) {
    cat(">>> Reading file", f, "... ... ...\n")
    taxa.p = readRDS(f)
  } else {
    cat(">>> Creating file", f, "with dada2::assignTaxonomy() ... ... ...\n")
    taxa.p = assignTaxonomy(seqtab.p.nochim, unite.ref, multithread=T, tryRC=T)
    saveRDS(taxa.p, f)
  }

## ---- inspect UNpooled -------
  taxa.print <- taxa				# Removing sequence rownames for display only
  rownames(taxa.print) <- NULL
  cat(">>> Inspect UNpooled taxa -- head taxa.print:\n")
  head(taxa.print)  						## Fungi o_O ??

## ---- inspect pooled -------
  taxa.p.print = taxa.p
  rownames(taxa.p.print) = NULL
  cat(">>> Inspect pooled taxa -- head taxa.p.print:\n")
  head(taxa.p.print)            ## Fungi o_O ??



## ---- phyloseq -------------------------------------------------------------------------------------------------------------------------------------------- ##

## ---- UNpooled -------
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

## ---- pooled -------
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

## ---- UNpooled -------
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

## ---- pooled -------
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



## ----ITSx-----------------------------------------------------------------------------------------------------------------------------

system2( "conda",
         args = c("run -n", itsx.conda.n, # set conda env
                  "ITSx --bugs" ))			  # version 1.1.3

## ---- UNpooled -------
  f = file.path(fna.dir,"ps.ASV.fasta")
  fo = file.path(itsx.dir, "ITSx.out")
  fe = file.path(itsx.dir, "ITSx.err")
  o.pfx = file.path(itsx.dir, "ITSx_ps.ASV")
  if( !file.exists(fo) ) {
    cat(">>> Running ITSx ... ... ...\n")
    system2( "conda", 
             args = c("run -n",itsx.conda.n,  # set conda env
                      "ITSx --cpu", n.threads,
                      "-i", f,		            # input fasta
                      "-o", o.pfx	),	        # output prefix
             stdout = fo,
             stderr = fe 
    )
  }

## ---- pooled -------
  f = file.path(fna.dir,"ps.p.ASV.fasta")
  fo = file.path(itsx.dir, "ITSx.p.out")
  fe = file.path(itsx.dir, "ITSx.p.err")
  o.pfx = file.path(itsx.dir, "ITSx_ps.p.ASV")
  if( !file.exists(fo) ) {
    cat(">>> Running ITSx ... ... ...\n")
    system2( "conda", 
             args = c("run -n", itsx.conda.n,     #  set conda env
                      "ITSx --cpu", n.threads, 
                      "-i", f,                    # input fasta
                      "-o", o.pfx ),			        # output prefix
             stdout = fo,
             stderr = fe 
    )
  }



## ----ITSx processing------------------------------------------------------------------------------------------------------------------

## ---- UNpooled -------
  
  # update phyloseq obj to have ITSx-treated ASVs
  fna = file.path(itsx.dir, "ITSx_ps.ASV.ITS2.fasta")
  f = file.path(itsx.dir,"ps.itsx.RDS")
  if( !file.exists(f) ) {
    cat(">>> Creating file", f, " ... ... ...\n")
    
    itsx.fna = readDNAStringSet(fna)
    # change names to just be ASV#
    names(itsx.fna) =
      names(itsx.fna) %>% 
      sub(pattern = "\\|.*$",replacement ="" )
    # get seqs as character vector
    itsx.asv.char = as.character(itsx.fna)
    
    # Check if any ASVs did not pass ITSx and filter those out of phyloseq object
    tn.filt = taxa_names(ps)[!(taxa_names(ps) %in% names(itsx.fna))]
    if(length(tn.filt) > 0) {
      cat(">>> Removing ITSx-failed ASVs:", tn.filt," ... ... ...\n")
    }
    ps.itsx = prune_taxa(taxa_names(ps) %in% names(itsx.fna), ps)

    # duplicated() determines which elements of a vector or data frame 
    #   are duplicates of elements with smaller subscripts
    w.dup = which(duplicated(itsx.asv.char))        # these indices are duplicates of smaller indices
    if( length(w.dup) > 0 ) {
      dup.m = match( itsx.asv.char[w.dup], itsx.asv.char ) # these are the indices that the dups match to
      cat(">>> Identical ASVs after ITSx:\n")
      print( cbind( "orig_asv" = names(itsx.asv.char)[dup.m],
                    "dup_asv" = names(itsx.asv.char)[w.dup]) )
      cat(">>> ... ... ...\n")
      
      # merged taxa/otu tables
      merge.ps.itsx = ps.itsx
      for( i in seq_along(w.dup) ) {
        # match results (dup.m) must come first, so that duplicates (w.dup)
        #   are merged into the smaller indexed ASV (dup.m)
        # multiple matches to same ASV are merged one by one during loop
        merge.ps.itsx = 
          merge_taxa( merge.ps.itsx, 
                      taxa_names(ps.itsx)[c(dup.m[i],w.dup[i])] )
      }
      new.otutab = otu_table(merge.ps.itsx)
      
      # remove duplicate ASVs from ITSx refseq table
      new.seqs = Biostrings::DNAStringSet(itsx.asv.char[-w.dup])
      
      # recall taxonomy for ITSx treated ASVs
      new.taxtab = assignTaxonomy(new.seqs, unite.ref, multithread = T)
      rownames(new.taxtab) = colnames(new.otutab)
      
      # remake phyloseq object with new ITSx data
      ps.itsx <- phyloseq( tax_table(new.taxtab), 
                           otu_table(new.otutab, taxa_are_rows = F),
                           refseq(new.seqs) )
    }
    
    cat(">>> Writing file", f, "... ... ...\n")
    saveRDS(ps.itsx, f)
  }


## ---- pooled -------
  
  # update phyloseq obj to have ITSx-treated ASVs
  fna = file.path(itsx.dir, "ITSx_ps.p.ASV.ITS2.fasta")
  f = file.path(itsx.dir,"ps.p.itsx.RDS")
  if( !file.exists(f) ) {
    cat(">>> Creating file", f, " ... ... ...\n")
    
    itsx.fna = readDNAStringSet(fna)
    # change names to just be ASV#
    names(itsx.fna) =
      names(itsx.fna) %>% 
      sub(pattern = "\\|.*$",replacement ="" )
    # get seqs as character vector
    itsx.asv.char = as.character(itsx.fna)
    
    # Check if any ASVs did not pass ITSx and filter those out of phyloseq object
    tn.filt = taxa_names(ps.p)[!(taxa_names(ps.p) %in% names(itsx.fna))]
    if(length(tn.filt) > 0) {
      cat(">>> Removing ITSx-failed ASVs:", tn.filt," ... ... ...\n")
    }
    ps.p.itsx = prune_taxa(taxa_names(ps.p) %in% names(itsx.fna), ps.p)
    
    # duplicated() determines which elements of a vector or data frame 
    #   are duplicates of elements with smaller subscripts
    w.dup = which(duplicated(itsx.asv.char))        # these indices are duplicates of smaller indices
    if( length(w.dup) > 0 ) {
      dup.m = match( itsx.asv.char[w.dup], itsx.asv.char ) # these are the indices that the dups match to
      cat(">>> Identical ASVs after ITSx:\n")
      print( cbind( "orig_asv" = names(itsx.asv.char)[dup.m],
                    "dup_asv" = names(itsx.asv.char)[w.dup]) )
      cat(">>> ... ... ...\n")
      
      # merged taxa/otu tables
      merge.ps.p.itsx = ps.p.itsx
      for( i in seq_along(w.dup) ) {
        # match results (dup.m) must come first, so that duplicates (w.dup)
        #   are merged into the smaller indexed ASV (dup.m)
        # multiple matches to same ASV are merged one by one during loop
        merge.ps.p.itsx = 
          merge_taxa( merge.ps.p.itsx, 
                      taxa_names(ps.p.itsx)[c(dup.m[i],w.dup[i])] )
      }
      new.otutab = otu_table(merge.ps.p.itsx)
      
      # remove duplicate ASVs from ITSx refseq table
      new.seqs = Biostrings::DNAStringSet(itsx.asv.char[-w.dup])
      
      # recall taxonomy for ITSx treated ASVs
      new.taxtab = assignTaxonomy(new.seqs, unite.ref, multithread = T)
      rownames(new.taxtab) = colnames(new.otutab)
      
      # remake phyloseq object with new ITSx data
      ps.p.itsx <- phyloseq( tax_table(new.taxtab), 
                             otu_table(new.otutab, taxa_are_rows = F),
                             refseq(new.seqs) )
    }
    
    cat(">>> Writing file", f, "... ... ...\n")
    saveRDS(ps.p.itsx, f)
  }



## ----save image-----------------------------------------------------------------------------------------------------------------------------

cat(">>> sessionInfo() ... ... ...\n")
sessionInfo()
f = paste0(format(Sys.time(), "%b%d%Y_%X%Z"),".RData") %>% gsub(pattern=" ",replacement="")
image.file = file.path(rdat.dir, f)
cat(">>> Saving image file", image.file, "... ... ...\n")
save.image(image.file)
cat(">>> COMPLETED\n")

