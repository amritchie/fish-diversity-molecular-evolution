
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

fish2.with.mito.coding <- fishgenes_genera_accessions %>% 
  filter(TaxonRemoved == "In Tree", !grepl("coi|cytb", blacklisted))

## Filter for SPECIES with available sequence

fish2.spp.with.nuc.all <- fish2.with.nuc.all %>%
  group_by(Species) %>% 
  summarise(NUC = any(!is.na(RAG1.acc) & RAG1.acc != "" & RAG1.cov != 0)) %>% 
  filter(NUC) %>% 
  pull(Species)

## Filter for GENERA with available sequence

fish2.genera.with.mito.coding <- fish2.with.mito.coding %>% 
  group_by(genus) %>% 
  summarise(MITO = any(!is.na(COI.acc) & !is.na(CYTB.acc) & COI.acc != "" &  
                         CYTB.acc != "" & COI.cov > 0 & CYTB.cov !=0)) %>% 
  filter(MITO) %>% 
  pull(genus)

## Create a tree with one tip per GENUS with sequence

fish2.genera.with.mito.coding.tree <- keep.tip(fishtree, 
                                            fishgenes_genera_accessions %>% filter(TaxonRemoved == "In Tree", 
                                                                                   family %in% fish2.genera.with.mito.coding) %>% 
                                              mutate(our_name = str_replace_all(Species, " ", "_")) %>% group_by(genus) %>% 
                                              summarise(our_name = first(our_name)) %>% 
                                              pull(our_name))

fish2.genera.with.mito.coding.tree$tip.label <- sapply(fish2.genera.with.mito.coding.tree$tip.label, 
                                                    function (x) (fishgenes_genera_accessions %>% 
                                                                    mutate(sp = str_replace_all(Species, " ", "_")) %>% 
                                                                    filter(sp == x) %>% 
                                                                    pull(genus))[1]
)

## Select all immediately adjacent sister pairs

fish2.genera.with.mito.coding.pairs <- get.cherries(fish2.genera.with.mito.coding.tree)

## Extract all members of sister pairs with available sequence

fish2.genera.with.mito.coding.pairclades.seq <- sapply(
  fish2.genera.with.mito.coding.pairs, 
  function (x) keep.tip(fishtree, 
                        str_replace_all(fishgenes_genera_accessions %>% 
                                          filter(genus %in% x, 
                                                 Species %in% fish2.spp.with.mito.coding) %>% 
                                          pull (Species), " ", "_")
  ), 
  simplify = F)

## Extract all available members of sister pairs

fish2.genera.with.mito.coding.pairclades.all <- sapply(
  fish2.genera.with.mito.coding.pairs, 
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

fish2.genera.with.mito.coding.pairclades.sampled <- sapply(
  fish2.genera.with.mito.coding.pairclades.seq, 
  phylo.sample.even, 
  simplify = F
)

## Make a tree showing relationships between pairs

fish2.genera.with.mito.coding.pairtree <- keep.tip(fishtree, 
                                             sapply(fish2.genera.with.mito.coding.pairclades.sampled, 
                                                    function (x) x$tip.label[1]))

## Initial arrangement of sister pairs into quartets as outgroup (actual quartets later adjusted manually)
                                             
fish2.genera.with.mito.coding.pairtree$tip.label <- 1:Ntip(fish2.genera.with.mito.coding.pairtree)

## Make up final table of pair, accession and quartet assignments

fish2.genera.with.mito.coding.xtets <- find.xtets.sisterpairs(fish2.genera.with.mito.coding.pairtree)


fish2.genera.with.mito.coding.accspairs.t <- bind_rows(
  sapply(1:length(fish2.genera.with.mito.coding.pairclades.sampled),  
         function (x) fishgenes_genera_accessions %>% 
           filter(Species %in% str_replace_all(
             fish2.genera.with.mito.coding.pairclades.sampled[[x]]$tip.label, "_", " ")
           ) %>% 
           select(Species, COI.acc, CYTB.acc) %>% 
           filter(COI.acc != "" & CYTB.acc != "") %>% 
           mutate(pair_no = x) %>%
           
           distinct(Species, .keep_all=T), simplify = F)
)

halves <- sapply(fish2.genera.with.mito.coding.pairclades.sampled, 
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

fish2.genera.with.mito.coding.accspairs.t2 <- fish2.genera.with.mito.coding.accspairs.t %>% 
  left_join(clade.names, by="Species")

fish2.genera.with.mito.coding.accspairs.t2$quartet_str <- sapply(
  fish2.genera.with.mito.coding.accspairs.t$pair_no, 
  function (x) paste(which(sapply(fish2.genera.with.mito.coding.xtets, 
                                  function (y) x %in% y)
  ), collapse=" ")
)

## Write accessions/pairs table

write.csv(fish2.genera.with.mito.coding.accspairs.t2, "outputs\\Accessions_pairs_tables\\Genera_with_Mito_Coding_AccsPairsTable.csv")

## Gather and write sequences

fish2.genera.with.mito.coding.COI.seqs <- read.GenBank(fish2.genera.with.mito.goding.accspairs.t$COI.acc)
fish2.genera.with.mito.coding.CYTB.seqs <- read.GenBank(fish2.genera.with.mito.coding.accspairs.t$CYTB.acc)

write.FASTA(fish2.genera.with.mito.coding.COI.seqs, 
            "alignments\\fish2_genera_with_mito_coding_COI.fasta")
write.FASTA(fish2.genera.with.mito.coding.CYTB.seqs, 
            "alignments\\fish2_genera_with_mito_coding_CYTB.fasta")

## Steps for Stage 2:

#### Manually convert sequences into quartet-wise alignments, edit, check and trim
#### Manually adjust quartets to minimise long inter-pair branches and triplets
#### Extract quartet trees as constraints for PAML and label each sister clade with a different branch model ($1, $2, etc) for codeml
#### Run baseml and codeml on each table