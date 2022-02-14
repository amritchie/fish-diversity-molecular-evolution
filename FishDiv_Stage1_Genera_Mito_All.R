
require(ape)
require(tidyverse)
require(phytools)

source("R\\my.drop.tip.R")
source("R\\bigfish_find_xtets.R")
source("R\\get_cherries.R")
source("R\\phylo_halves.R")

fishgenes_genera_accessions <- read.csv("fishgenes_genera_accessions_all.csv")
fishtree <- read.tree("actinopt_12k_treePL.tre")


## Remove rogue or blacklisted taxa
fish2.with.mito.all <- fishgenes_genera_accessions %>% 
  filter(TaxonRemoved == "In Tree", !grepl("coi|12s|16s|cytb", blacklisted))


## Filter for SPECIES with available sequence
fish2.spp.with.mito.all <- fish2.with.mito.all %>%
  group_by(Species) %>% 
  summarise(MITO = any(!is.na(X12s.acc) & !is.na(X16s.acc) & !is.na(COI.acc) & !is.na(CYTB.acc) & 
                         X12s.acc != "" & X16s.acc != "" &  COI.acc != "" &  CYTB.acc != "" &  
                         X12s.cov > 0 &  X16s.cov > 0 &  COI.cov > 0 &  CYTB.cov !=0)) %>% 
  filter(MITO) %>% 
  pull(Species)

## Filter for GENERA with available sequence
fish2.genera.with.mito.all <- fish2.with.mito.all %>% 
  group_by(genus) %>% 
  summarise(MITO = any(!is.na(X12s.acc) & !is.na(X16s.acc) & !is.na(COI.acc) & !is.na(CYTB.acc) & 
                         X12s.acc != "" & X16s.acc != "" &  COI.acc != "" &  CYTB.acc != "" &  
                         X12s.cov > 0 &  X16s.cov > 0 &  COI.cov > 0 &  CYTB.cov !=0)) %>% 
  filter(MITO) %>% 
  pull(genus)


## Create a tree with one tip per GENUS with sequence


fish2.genera.with.mito.all.tree <- keep.tip(fishtree, 
                                         fishgenes_genera_accessions %>% filter(TaxonRemoved == "In Tree", 
                                                                                family %in% fish2.genera.with.mito.all) %>% 
                                           mutate(our_name = str_replace_all(Species, " ", "_")) %>% group_by(genus) %>% 
                                           summarise(our_name = first(our_name)) %>% 
                                           pull(our_name))

fish2.genera.with.mito.all.tree$tip.label <- sapply(fish2.genera.with.mito.all.tree$tip.label, 
                                                 function (x) (fishgenes_genera_accessions %>% 
                                                                 mutate(sp = str_replace_all(Species, " ", "_")) %>% 
                                                                 filter(sp == x) %>% 
                                                                 pull(genus))[1]
)

## Select all immediately adjacent sister pairs

fish2.genera.with.mito.all.pairs <- get.cherries(fish2.genera.with.mito.all.tree)

## Extract all members of sister pairs with available sequence

fish2.genera.with.mito.all.pairclades.seq <- sapply(
  fish2.genera.with.mito.all.pairs, 
  function (x) keep.tip(fishtree, 
                        str_replace_all(fishgenes_genera_accessions %>% 
                                          filter(genus %in% x, 
                                                 Species %in% fish2.spp.with.mito.all) %>% 
                                          pull (Species), " ", "_")
  ), 
  simplify = F)
  
## Extract all available members of sister pairs

fish2.genera.with.mito.all.pairclades.all <- sapply(
  fish2.genera.with.mito.all.pairs, 
  function (x) extract.clade(fishtree, 
                             getMRCA(fishtree, 
                                     str_replace_all(fishgenes_genera_accessions %>% 
                                                       filter(TaxonRemoved == "In Tree", 
                                                              genus %in% x) %>% 
                                                       pull (Species), " ", "_")
                             )
  ), 
  simplify = F)

## Node density effect compensation - randomly sample down larger sister clade so that sequence numbers are even

fish2.genera.with.mito.all.pairclades.sampled <- sapply(
  fish2.genera.with.mito.all.pairclades.seq, 
  phylo.sample.even, 
  simplify = F
)

## Make a tree showing relationships between pairs

fish2.genera.with.mito.all.pairtree <- keep.tip(fishtree, 
                                             sapply(fish2.genera.with.mito.all.pairclades.sampled, 
                                                    function (x) x$tip.label[1]))

## Initial arrangement of sister pairs into quartets as outgroup (actual quartets later adjusted manually)
                                             
fish2.genera.with.mito.all.pairtree$tip.label <- 1:Ntip(fish2.genera.with.mito.all.pairtree)

## Make up final table of pair, accession and quartet assignments

fish2.genera.with.mito.all.xtets <- find.xtets.sisterpairs(fish2.genera.with.mito.all.pairtree)

fish2.genera.with.mito.all.accspairs.t <- bind_rows(
  sapply(1:length(fish2.genera.with.mito.all.pairclades.sampled),  
         function (x) fishgenes_genera_accessions %>% 
           filter(Species %in% str_replace_all(
             fish2.genera.with.mito.all.pairclades.sampled[[x]]$tip.label, "_", " ")
           ) %>% 
           select(Species, X12s.acc, X16s.acc, COI.acc, CYTB.acc) %>% 
           filter(X12s.acc != "" & X16s.acc != "" & COI.acc != "" & CYTB.acc != "") %>% 
           mutate(pair_no = x) %>%
           
           distinct(Species, .keep_all=T), simplify = F)
)

halves <- sapply(fish2.genera.with.mito.all.pairclades.sampled, 
                 function (x) sapply(list(phylo.firsthalf(x), phylo.secondhalf(x)), 
                                     function (y) y$tip.label, simplify = F), 
                 simplify=F)

clade.names <- data.frame(Species = gsub("_", " ", unlist(halves)), 
                          clade_name = unlist(sapply(
                            1:length(halves), 
                            function (x) paste(x, 
                                               c(rep("1", length(halves[[x]][[1]])), 
                                                 rep("2", length(halves[[x]][[2]]))), 
                                               sep="_"))
                            )
                          )

fish2.genera.with.mito.all.accspairs.t2 <- fish2.genera.with.mito.all.accspairs.t %>% 
  left_join(clade.names, by="Species")

fish2.genera.with.mito.all.accspairs.t2$quartet_str <- sapply(
  fish2.genera.with.mito.all.accspairs.t$pair_no, 
  function (x) paste(which(sapply(fish2.genera.with.mito.all.xtets, 
                                  function (y) x %in% y)
  ), collapse=" ")
)

## Write accessions/pairs table

write.csv(fish2.genera.with.mito.all.accspairs.t2, "outputs\\Accessions_pairs_tables\\Genera_with_Mito_All_AccsPairsTable.csv")

## Gather and write sequences


fish2.genera.with.mito.all.12S.seqs <- read.GenBank(fish2.genera.with.mito.all.accspairs.t$X12s.acc)
fish2.genera.with.mito.all.16S.seqs <- read.GenBank(fish2.genera.with.mito.all.accspairs.t$X16s.acc)
fish2.genera.with.mito.all.COI.seqs <- read.GenBank(fish2.genera.with.mito.all.accspairs.t$COI.acc)
fish2.genera.with.mito.all.CYTB.seqs <- read.GenBank(fish2.genera.with.mito.all.accspairs.t$CYTB.acc)

write.FASTA(fish2.genera.with.mito.all.12S.seqs, 
            "alignments\\raw_seqs\\fish2_genera_with_mito_all_12s.fasta")
write.FASTA(fish2.genera.with.mito.all.16S.seqs, 
            "alignments\\raw_seqs\\fish2_genera_with_mito_all_16s.fasta")
write.FASTA(fish2.genera.with.mito.all.COI.seqs, 
            "alignments\\fish2_genera_with_mito_all_COI.fasta")
write.FASTA(fish2.genera.with.mito.all.CYTB.seqs, 
            "alignments\\fish2_genera_with_mito_all_CYTB.fasta")

## Steps for Stage 2:

#### Manually convert sequences into quartet-wise alignments, edit, check and trim
#### Manually adjust quartets to minimise long inter-pair branches and triplets
#### Extract quartet trees as constraints for PAML and label each sister clade with a different branch model ($1, $2, etc) for codeml
#### Run baseml and codeml on each table