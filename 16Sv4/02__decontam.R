## kkyle@KlassenLab-server:~/kek16SallITStrachyDecons/fastq/p1p2_16S/
## R v4.4.1

##########################################################################################################################################
######  SET THESE VARIABLES TO YOUR SPECIFICS   ##########################################################################################
##########################################################################################################################################

# how many threads when multithreading?
n.threads = 20

# required packages
libs = c( "decontam",
          "phyloseq", 
          "rlist",
          "tidyverse", 
          "phytools"#, 
          #"ShortRead",
          #"Biostrings"
        )

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



## ---- DECONTAM ---------------------------------------------------------------------------- ##

ps.p = readRDS("")
samdat = 
  data.frame( 
    "seq.id" = sample_names(ps.p), 
    "Sample_or_Control"="True Sample",
    stringsAsFactors=FALSE
  )
controls.i = grep("NFW|NC",samdat$seq.id)
samdat$Sample_or_Control[controls.i] = "Control Sample"
rownames(samdat) = samdat$seq.id

sample_data(ps.p) = samdat

df <- as.data.frame(sample_data(ps.p)) 
df$LibrarySize <- sample_sums(ps.p)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
pdf("libsize_trueVScontrol.pdf")
ggplot( data=df, 
        aes(x=Index, y=LibrarySize, color=Sample_or_Control) 
      ) + geom_point()
dev.off()

sample_data(ps.p)$is.neg <- sample_data(ps.p)$Sample_or_Control == "Control Sample"



## ---- Call contaminants ---------------------------------------------------------------------------- ##

thresh = c(0.1, 0.25, 0.5)
contamdf.ls = list()
for( i in seq_along(thresh) ){
  ln = paste0( "prev", as.character(thresh[i]) )
  contamdf.ls[[ln]] <- 
    isContaminant(  ps.p, 
                    method="prevalence", 
                    neg="is.neg", 
                    threshold=thresh[i]
    )
  table(contamdf.ls[[ln]]$contaminant)
}
# 0.01
#  FALSE  TRUE 
#   2274     8 
# 0.25
#   FALSE  TRUE 
#    2231    51 
# 0.5
#   FALSE  TRUE 
#    2165    117 

# ---- Make phyloseq object of presence-absence in negative controls and true samples
ps.p.pa <- transform_sample_counts(ps.p, function(abund) 1*(abund>0))
ps.p.pa.neg <- 
  prune_samples(  sample_data(ps.p.pa)$Sample_or_Control == "Control Sample", 
                  ps.p.pa 
  )
ps.p.pa.pos <- 
  prune_samples(  sample_data(ps.p.pa)$Sample_or_Control == "True Sample", 
                  ps.p.pa
  )

# ---- Make data.frame of prevalence in positive and negative samples
ps.nocontam.ls = list()
for( ln in names(contamdf.ls) ) {
  #ln = names(contamdf.ls)[[i]]
  df.pa <- 
    data.frame( pa.pos = taxa_sums(ps.p.pa.pos), 
                pa.neg = taxa_sums(ps.p.pa.neg),
                contaminant = contamdf.ls[[ln]]$contaminant
  )
  f = paste0("decontam.", ln, ".pdf")
  pdf(f)
  print(
    ggplot( data=df.pa, 
          aes(x=pa.neg, y=pa.pos, color=contaminant)
    ) + 
    geom_point() +
    xlab("Prevalence (Negative Controls)") + 
    ylab("Prevalence (True Samples)")
  )
  dev.off()

  # ---- prune contaminants
  asv.keep <- 
    contamdf.ls[[ln]] %>% filter(contaminant == "FALSE") %>% rownames()
  ps.nocontam.ls[[ln]] <- 
    prune_taxa(asv.keep, ps.p)
  f = paste0("ps.p.decontam.", ln, ".RDS")
  saveRDS(ps.nocontam.ls[[ln]], f)
}

## ---- make relabund
rab.ps.nocontam.ls = make_rel_abund(ps.nocontam.ls)

## ---- plot by phylum
for( ln in names(rab.ps.nocontam.ls) ) {
  f = paste0("plotbar.phy_", ln, ".pdf")
  pdf(f, width=12, height=6)
  plot_bar(rab.ps.nocontam.ls[[ln]], fill="Phylum")
  dev.off()
}

## ---- decontam stats 
decontam.totsize.df <- 
  data.frame( ps.p.nreads = sum(sample_sums(ps.p)),
              prev0.1.nreads = sum(sample_sums(ps.nocontam.ls$prev0.1)),
              prev0.25.nreads = sum(sample_sums(ps.nocontam.ls$prev0.25)),
              prev0.5.nreads = sum(sample_sums(ps.nocontam.ls$prev0.5))
  )
decontam.totsize.df <-
  decontam.totsize.df %>% 
  mutate(prev0.1.percL = (ps.p.nreads-prev0.1.nreads)/ps.p.nreads*100) %>%
  mutate(prev0.25.percL = (ps.p.nreads-prev0.25.nreads)/ps.p.nreads*100) %>%
  mutate(prev0.5.percL = (ps.p.nreads-prev0.5.nreads)/ps.p.nreads*100)

samsums.decontam.df <- 
  data.frame(
    ps.p.samsums = sample_sums(ps.p),
    prev0.1.samsums = sample_sums(ps.nocontam.ls$prev0.1),
    prev0.25.samsums = sample_sums(ps.nocontam.ls$prev0.25),
    prev0.5.samsums = sample_sums(ps.nocontam.ls$prev0.5)
  )

decontam.stats.df <- 
  samsums.decontam.df %>% 
  mutate(prev0.1.nreads.lost = ps.p.samsums-prev0.1.samsums) %>% 
  mutate(prev0.25.nreads.lost = ps.p.samsums-prev0.25.samsums) %>%
  mutate(prev0.5.nreads.lost = ps.p.samsums-prev0.5.samsums) %>%
  mutate(prev0.1.percL = prev0.1.nreads.lost/ps.p.samsums*100) %>%
  mutate(prev0.25.percL = prev0.25.nreads.lost/ps.p.samsums*100) %>%
  mutate(prev0.5.percL = prev0.5.nreads.lost/ps.p.samsums*100)
write.csv(decontam.stats.df, "decontam.stats.df.csv")



## ----save image-----------------------------------------------------------------------------------------------------------------------------
cat(">>> sessionInfo() ... ... ...\n")
sessionInfo()

image.file = paste0("./rdata/",format(Sys.time(), "%b%d%Y_%X%Z"),".RData")
cat(">>> Saving image file", image.file, "... ... ...\n")
save.image(image.file)
cat(">>> COMPLETED\n")
