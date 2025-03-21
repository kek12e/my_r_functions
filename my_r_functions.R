# library(dada2)
# library(vegan)
# library(tidyverse)
# library(ggplot2)
# library(phyloseq)
# 
# set.seed(14)

my_load_libs <- function(libs) {
  if(!require(pacman, quietly=T)) {
	  install.packages("pacman")
  }
  library(pacman)
  p_load(char=libs)

  cat(">>> Libraries loaded:\n")
  for( i in seq_along(libs) ) {
	  cat("\t", libs[i], "version", as.character(packageVersion(libs[i])),"\n")
  }
}

my_rarefy <- function(x, sample) {
  x <- x[x>0]
  sum(1-exp(lchoose(sum(x) - x, sample) - lchoose(sum(x), sample)))
}

make_tidy <- function(ps_otutab) {
  # takes in phyloseq otutab, returns tidy tibble with cols: sid, asv, count
  ps_otutab  <-  as.data.frame(ps_otutab) 
  ps_otutab$sid = rownames(ps_otutab)
  ps_otutab %>% #    make tidy 
    relocate(sid) %>%          #    make sid the first col
    pivot_longer(!sid, names_to = "asv", values_to="count")
}

coll_curve <- function(data, group) {
  data %>% 
    filter(sid == group) %>% 
    uncount(count) %>% 
    sample_n(n()) %>% 
    mutate(observation = row_number()) %>% 
    arrange(asv, sid) %>% 
    group_by(asv) %>% 
    mutate(distinct = row_number() == 1) %>%
    ungroup() %>% 
    arrange(observation) %>% 
    mutate(s = cumsum(distinct)) %>% 
    select(observation, s)
}

taxnames_lt_threshold <- function(physeq, threshold=0.01) {
  ## returns the taxa names that are < threshold in every sample
  ## physeq object must be relative abundance data
  
  ps = physeq
  tn = names(which(filter_taxa(ps, function(x) sum(x >= threshold, na.rm=T) == 0)))
  return(tn)
}

per_sample_taxafilt <- function(ps, threshold=0.01, glom="none") {
  ## ps is a phyloseq object with a tax_table
  ## taxa less than threshold are set to "other $threshold" at all taxonomy ranks
  ## taxa are not glommed by default; set glom to tax_rank you want to glom at
  
  ps.rab = ps
  # 0.) check if rel abund
  if( !all( sample_sums(ps) == 1 ) )
    ps.rab = make_rel_abund(ps)
  # 1.) ASVs to set to "other"
  tn = taxnames_lt_threshold(ps.rab, threshold)
  # 2.) set all tax ranks to "other"
  if(length(tn) > 0) {
    n.taxranks = ncol(tax_table(ps))
    tax_table(ps)[tn,1:n.taxranks] = paste("other",threshold) # ex: "other 0.01"
  }
  # 3.) glom
  if( glom %in% rank_names(ps) ) {
    cat("\n>>> Glomming by", glom, "... ... ...\n")
    ps = tax_glom(ps, glom)
  }
  # save
  return(ps)
}

make_rel_abund <- function(psList) {
  ## input list of phyloseq objects that have otu_table defined
  ## returns same list after transforming counts to relative abundance
  
  if(is.list(psList)) {
    for( i in seq_along(psList) ) {
      
      ps.prop = transform_sample_counts(psList[[i]], function(otu) otu/sum(otu) )
      
      # check for NaNs and set to 0
      w.nan = which(is.nan(sample_sums(ps.prop)))
      n.t = ntaxa(ps.prop)
      if(length(w.nan) > 0) {
        otu_table(ps.prop)[w.nan, 1:n.t] = 0
      }
      
      # save
      psList[[i]] = ps.prop
    }
  }
  else {
    ps.prop = transform_sample_counts(psList, function(otu) otu/sum(otu) )
    psList = ps.prop
  }
  return(psList)
}

na_to_unclassified_taxa <- function(physeq, rank.start=1, glom=0) {
  ## set "NA"s in tax_table to "unclass. prior_tax_rank"
  ## glom = 0 means do not glom. glom = 1...n means glom to rank_names(tax_table(ps))[1...n]
  ## rank.start is the highest (left most) taxonomic rank to start at
  ##  typical phyloseq objects have taxranks: 
  ##  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
  ## Example: 
  ##  
  ##
  ##
  ps = physeq
  rank.last = length(rank_names(ps))
  if ( rank.start > rank.last ) {
    print("Warning: chosen rank.start is outside upper bound of rank_names(tax_table(physeq))!")
    print(paste("Setting rank.start to rank.last =", rank.last))
    rank.start = rank.last
  } else if (rank.start < 1) {
    print(paste("Warning: rank.start cannot be < 1!"))
    print("Setting rank.start to 1")
    rank.start = 1
  }
  # Have to start at "Phylum" bc there is no prior tax.rank for "Kingdom"
  # i.e. rank_names(ps)[0] is nothing
  if (rank.start == 1) {
    rank.start = rank.start+1
  }
  for ( i in rank.start:rank.last ) {
    w.na = which(is.na(tax_table(ps)[,i]))
    if(length(w.na) > 0) {
      unclass.str = paste("unclass.",tax_table(ps)[w.na,rank_names(ps)[i-1]])
      tax_table(ps)[w.na,i:rank.last] = unclass.str
    }
  }
  # check if need to tax_glom and for improper input
  if(glom < 0) {
    print("Cannot glom to rank < 0. No glomming performed.")
  }
  else if(glom > 0) {
    if(glom > rank.last) {
      print(
        paste(
          "Cannot glom to tax_rank > rank.last. Setting glom to rank.last =", 
          rank.last, 
          "..."
        )
      )
      glom = rank.last
    }
    
    print(
      paste(
        "Glomming taxa to level:",
        rank_names(ps)[glom]
      )
    )
    ps = 
      tax_glom( 
        ps,
        rank_names(ps)[glom]
      )
  }
  
  return(ps)
}

make_cultivar_fungus = function(ps, pfx = "") {
  ## defines "cultivar fungus" genus
  ## input: phyloseq object, optional filename prefix
  ## return: phyloseq object
  
  # check for NA's in Genus names
  if( all( is.na(tax_table(ps)[,"Genus"]) ) ) {
    cat(">>> NA's detected in Genus names, running na_to_unclassified_taxa() ... ... ...")
    ps = na_to_unclassified_taxa(ps)
  }
  
  # get "cultivar fungus" indices
  w.cf = which( tax_table(ps)[,"Genus"] == "g__Leucocoprinus" |
                tax_table(ps)[,"Genus"] == "g__Leucoagaricus" |
                tax_table(ps)[,"Genus"] == "unclassified f__Agaricaceae" )
  
  # create "cultivar fungus" fasta for BLASTn confirm
  if( length(w.cf) > 0 ) {
    cat(">>> Potential cultivar fungus ASVs detected. Creating fasta ... ... ...\n")
    
    # make fasta of "cultivar fungus" ASVs
    cf.dnass = Biostrings::DNAStringSet(unlist(lapply(refseq(ps)[w.cf],toString)))
    fa.names = paste0(names(cf.dnass),"_"," ",tax_table(ps)[w.cf,"Genus"])
    names(cf.dnass) = fa.names
    f = file.path(fna.dir, paste0(pfx,"cultivar_fungus_ASVs.fasta"))
    cat(">>> Writing file", f," ... ... ...\n")
    cf.dnass %>% Biostrings::writeXStringSet(f, format="fasta")

    # set new genus name
    tax_table(ps)[w.cf,"Genus"] = "Cultivar Fungus"

    # glom by genus
    ps = tax_glom(ps, "Genus", NArm=T)
  }
  
  return(ps)
  
} #\end of func
