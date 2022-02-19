
require(ape)
require(tidyverse)
require(phytools)

source("R/my.drop.tip.R")
source("R/bigfish_find_xtets.R")
source("R/get_cherries.R")
source("R/phylo_halves.R")
source("R/my_read_genbank.R")

fishgenes_genera_accessions <- read.csv("fishgenes_genera_accessions_all.csv")
fishtree <- read.tree("trees/Rabosky_actinopt_12k_treePL.tre")


## Remove rogue or blacklisted taxa

fish2.with.nuc.all <- fishgenes_genera_accessions %>% 
  filter(TaxonRemoved == "In Tree", !grepl("rag1", blacklisted))


## Filter for SPECIES with available sequence

fish2.spp.with.nuc.all <- fish2.with.nuc.all %>%
  group_by(Species) %>% 
  summarise(NUC = any(!is.na(RAG1.acc) & RAG1.acc != "" & RAG1.cov != 0)) %>% 
  filter(NUC) %>% 
  pull(Species)


## Filter for FAMILIES with available sequence

fish2.fam.with.nuc.all <- fish2.with.nuc.all %>% 
  group_by(family) %>% 
  summarise(NUC = any(!is.na(RAG1.acc) & RAG1.acc != "")) %>% 
  filter(NUC) %>% 
  pull(family)

## Create a tree with one tip per FAMILY with sequence

fish2.fam.with.nuc.all.tree <- keep.tip(fishtree, 
                                            fishgenes_genera_accessions %>% filter(TaxonRemoved == "In Tree", 
                                                                                   family %in% fish2.fam.with.nuc.all) %>% 
                                              mutate(our_name = str_replace_all(Species, " ", "_")) %>% group_by(family) %>% 
                                              summarise(our_name = first(our_name)) %>% 
                                              pull(our_name))

fish2.fam.with.nuc.all.tree$tip.label <- sapply(fish2.fam.with.nuc.all.tree$tip.label, 
                                                    function (x) (fishgenes_genera_accessions %>% 
                                                                    mutate(sp = str_replace_all(Species, " ", "_")) %>% 
                                                                    filter(sp == x) %>% 
                                                                    pull(family))[1]
)


## Select all immediately adjacent sister pairs

fish2.fam.with.nuc.all.pairs <- get.cherries(fish2.fam.with.nuc.all.tree)

## Extract all members of sister pairs with available sequence

fish2.fam.with.nuc.all.pairclades.seq <- sapply(
  fish2.fam.with.nuc.all.pairs, 
  function (x) keep.tip(fishtree, 
                        str_replace_all(fishgenes_genera_accessions %>% 
                                          filter(family %in% x, 
                                                 Species %in% fish2.spp.with.nuc.all) %>% 
                                          pull (Species), " ", "_")
  ), 
  simplify = F)

## Extract all available members of sister pairs

fish2.fam.with.nuc.all.pairclades.all <- sapply(
  fish2.fam.with.nuc.all.pairs, 
  function (x) extract.clade(fishtree, 
                             getMRCA(fishtree, 
                                     str_replace_all(fishgenes_genera_accessions %>% 
                                                       filter(TaxonRemoved == "In Tree", 
                                                              family %in% x) %>% 
                                                       pull (Species), " ", "_")
                             )
  ), 
  simplify = F)

## Node density effect compensation - randomly sample down larger sister clade so that sequence numbers are even

fish2.fam.with.nuc.all.pairclades.sampled <- sapply(
  fish2.fam.with.nuc.all.pairclades.seq, 
  phylo.sample.even, 
  simplify = F
)

## Make a tree showing relationships between pairs

fish2.fam.with.nuc.all.pairtree <- keep.tip(fishtree, 
                                             sapply(fish2.fam.with.nuc.all.pairclades.sampled, 
                                                    function (x) x$tip.label[1]))                        

## Initial arrangement of sister pairs into quartets as outgroup (actual quartets later adjusted manually)
                                             
fish2.fam.with.nuc.all.pairtree$tip.label <- 1:Ntip(fish2.fam.with.nuc.all.pairtree)

## Make up final table of pair, accession and quartet assignments

fish2.fam.with.nuc.all.xtets <- find.xtets.sisterpairs(fish2.fam.with.nuc.all.pairtree)

fish2.fam.with.nuc.all.accspairs.t <- bind_rows(
  sapply(1:length(fish2.fam.with.nuc.all.pairclades.sampled),  
         function (x) fishgenes_genera_accessions %>% 
           filter(Species %in% str_replace_all(
             fish2.fam.with.nuc.all.pairclades.sampled[[x]]$tip.label, "_", " ")
           ) %>% 
           select(Species, RAG1.acc) %>% 
           filter(RAG1.acc != "") %>% 
           mutate(pair_no = x) %>%
           
           distinct(Species, .keep_all=T), simplify = F)
)

halves <- sapply(fish2.fam.with.nuc.all.pairclades.sampled, 
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

fish2.fam.with.nuc.all.accspairs.t2 <- fish2.fam.with.nuc.all.accspairs.t %>% 
  left_join(clade.names, by="Species")

fish2.fam.with.nuc.all.accspairs.t2$quartet_str <- sapply(
  fish2.fam.with.nuc.all.accspairs.t$pair_no, 
  function (x) paste(which(sapply(fish2.fam.with.nuc.all.xtets, 
                                  function (y) x %in% y)
  ), collapse=" ")
)

## Write accessions/pairs table

write.csv(fish2.fam.with.nuc.all.accspairs.t2, "outputs/Stage1_Accessions_pairs_tables/Fam_with_Nuc_All_AccsPairsTable.csv")

## Gather and write sequences.


## 'New Data' entries from Rabosky et al. 2018)
fish2.fam.with.nuc.all.newdata <- read.FASTA("alignments/newdata/fish2_fam_with_nuc_accspairs_t_newdata.fasta")

fish2.fam.with.nuc.all.RAG1.gb.seqs <- my.read.GenBank(
  fish2.fam.with.nuc.all.accspairs.t2$RAG1.acc[sapply(fish2.fam.with.nuc.all.accspairs.t2$RAG1.acc,
                                                    function (x) ifelse(str_detect(x, "New data"), F, T))])

fish2.fam.with.nuc.all.RAG1.gb.seqnames <- fish2.fam.with.nuc.all.accspairs.t2$Species[sapply(
  fish2.fam.with.nuc.all.accspairs.t2$RAG1.acc, 
  function (x) ifelse(str_detect(x, "New data"), F, T))]

fish2.fam.with.nuc.all.RAG1.seqs <- sapply(
  sapply(
    1:length(fish2.fam.with.nuc.all.accspairs.t2$RAG1.acc), 
    function (x) ifelse(str_detect(fish2.fam.with.nuc.all.accspairs.t2$RAG1.acc[x], "New data"), 
                        list(fish2.fam.with.nuc.all.newdata[str_replace_all(fish2.fam.with.nuc.all.accspairs.t2[x, "Species"], " ", "_")]),
                        list(fish2.fam.with.nuc.all.RAG1.gb.seqs[fish2.fam.with.nuc.all.RAG1.gb.seqnames == fish2.fam.with.nuc.all.accspairs.t2[x, "Species"]])), 
    simplify = F), function (x) x[1])

fish2.fam.with.nuc.all.RAG1.DNAbin <- fish2.fam.with.nuc.all.RAG1.seqs[[1]]

for (i in 2:length(fish2.fam.with.nuc.all.RAG1.seqs)) 
{
  if (length(fish2.fam.with.nuc.all.RAG1.seqs[[i]][[1]]) > 0) 
  {
    fish2.fam.with.nuc.all.RAG1.DNAbin <- c(fish2.fam.with.nuc.all.RAG1.DNAbin, 
                                            fish2.fam.with.nuc.all.RAG1.seqs[[i]])
  }
}

write.FASTA(fish2.fam.with.nuc.all.RAG1.DNAbin, 
            "alignments/raw_sequence/fish2_fam_with_nuc_all_RAG1.fasta")

## Steps for Stage 2:

#### Manually convert sequences into quartet-wise alignments, edit, check and trim
#### Manually adjust quartets to minimise long inter-pair branches and triplets
#### Extract quartet trees as constraints for PAML and label each sister clade with a different branch model ($1, $2, etc) for codeml
#### Run baseml and codeml on each table