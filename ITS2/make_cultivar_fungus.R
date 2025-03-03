make_cultivar_fungus = function(ps) {
  
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
    
    # write to file
    f = file.path(fna.dir,"cultivar_fungus_ASVs.fasta")
    cat(">>> Writing file", f," ... ... ...\n")
    cf.dnass %>% Biostrings::writeXStringSet(f, format="fasta")

    # set new genus name
    tax_table(ps)[w.cf,"Genus"] = "Cultivar Fungus"

    # glom by genus
    ps = tax_glom(ps, "Genus", NArm=T)
  }
  
  return(ps)
  
} #\end of func
