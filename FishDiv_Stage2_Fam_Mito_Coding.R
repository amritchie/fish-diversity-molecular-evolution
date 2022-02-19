library(ape)
library(tidyverse)
library(fs)
library(foreach)
library(doParallel)

source("R/FishDiv_Analyse_Pair_Data.R")
source("R/baseml.extract.tree.R")
source("R/phylo.average.brlen.R")
source("R/ordered.ls.dir.R")
source("R/welch_test.R")

fishgenes_genera_accessions <- read.csv("fishgenes_genera_accessions_all.csv")
rabosky.ata.trees <- read.tree("trees/Rabosky_actinopt_full.trees")
genera.mito.all.accspairs.t2 <- read.csv("outputs/Stage2_Accessions_pairs_tables/fish2_genera_with_mito_accspairs_t_2.csv", stringsAsFactors = F)
genera.mito.all.pairclades.all.nodup <- sapply(
  ordered.ls.dir("trees/Stage2_pair_trees/pairclades_all_species/Genera_with_mito/", 
                 regexp=".*tre$"), 
  read.tree, 
  simplify=F
)

genera.mito.all.pairclades.sampled <- sapply(
  ordered.ls.dir("trees/Stage2_pair_trees/pairclades_sampled/Genera_with_mito/", 
                 regexp=".*tre$"), 
  read.tree, 
  simplify=F
)

genera.mito.all.baseml.trees <- sapply(ordered.ls.dir("outputs/PAML_output/Genera_Mito_All/baseml", regexp = ".*txt"), 
										baseml.extract.tree, 
										simplify = F
										)

genera.mito.all.out <- fishdiv.analyse.pair.data(genera.mito.all.accspairs.t2, 
												fishgenes_genera_accessions, 
												genera.mito.all.baseml.trees, 
												genera.mito.all.pairclades.all.nodup, 
												"Genus"
												)

##############################

## Rectifying N_spp with ata trees
## Commented out as this part can take several hours. Pre-baked sets of ATA tree
## pairs are saved in the repository.
#
# ncores = detectCores()
# registerDoParallel(ncores)
#
# genera.mito.all.pairs.ata <- foreach(x=1:(dim(genera.mito.all.out)[1]), .packages="ape") %dopar% {sapply(rabosky.ata.trees, function(tr) sapply(c(phylo.firsthalf, phylo.secondhalf), function (f) f(my.extract.clade(phy=tr, node=my.getMRCA(tr, genera.mito.all.pairclades.sampled[[genera.mito.all.out$Pair[x]]]$tip.label)))$tip.label, simplify = F), simplify = F)}

load("ata_tree_pairs_prebaked/tmp_ata_pairs_welch/fish2_genera_with_mito_pairs_ata")
# If using monophyly filter, load("ata_tree_pairs_prebaked/tmp_ata_pairs_mono_welch/fish2_genera_with_mito_pairs_ata")

## Orientation of sister clades in ATA trees is not guaranteed
genera.mito.all.pairs.ata.sorted <- sapply(1:(dim(genera.mito.all.out)[1]), 
											function (x) sapply(fish2.genera.with.mito.pairs.ata[[x]], 
																function (y) if(length(intersect(y[[1]], phylo.firsthalf(genera.mito.all.pairclades.sampled[[genera.mito.all.out$Pair[x]]])$tip.label)) > 0) y else (rev(y)), 
																simplify = F), 
											simplify = F
											)
 
## Get clade sizes as average species numbers across ATA trees
genera.mito.all.nspp.ata.sorted <- sapply(genera.mito.all.pairs.ata.sorted, 
											function (x) apply(sapply(x, function (y) c(length(y[[1]]), length(y[[2]]))), 1, mean))

## Format final data sets with substitutions and clade sizes

## NO MONOPHYLY FILTER
genera.mito.all.ata.out <- genera.mito.all.out %>% 
							ungroup() %>% 
							mutate(N_spp_ata_first = log(genera.mito.all.nspp.ata.sorted[1,]), N_spp_ata_last = log(genera.mito.all.nspp.ata.sorted[2,])) %>% 
							mutate(swap_ata = sign(N_spp_ata_first - N_spp_ata_last)) %>% 
							mutate(N_spp_ata_diff = abs(N_spp_ata_first - N_spp_ata_last), Blen_ata_diff = swap_ata*(Blen_first - Blen_last)) %>%
							mutate(N_spp_ata_std = N_spp_ata_diff/sqrt(Age_first), Blen_ata_std = Blen_ata_diff/sqrt(Age_first)) %>% 
							filter(N_spp_ata_std > 0)

## MONOPHYLY FILTER - run this code instead of the above
# genera.mito.all.nonmono.pairs.ata <- sapply(genera.mito.all.pairs.ata.sorted, function (x) sum(sapply(x, function (y) length(intersect(str_extract(y[[1]], "^.*_"), str_extract(y[[2]], "^.*_"))) > 0)) > 20)
# 
# genera.mito.all.ata.out <- genera.mito.all.out %>% ungroup() %>% mutate(N_spp_ata_first = log(genera.mito.all.nspp.ata.sorted[1,]), N_spp_ata_last = log(genera.mito.all.nspp.ata.sorted[2,])) %>% filter(!genera.mito.all.nonmono.pairs.ata) %>% mutate(swap_ata = sign(N_spp_ata_first - N_spp_ata_last)) %>% mutate(N_spp_ata_diff = abs(N_spp_ata_first - N_spp_ata_last), Blen_ata_diff = swap_ata*(Blen_first - Blen_last)) %>% mutate(N_spp_ata_std = N_spp_ata_diff/sqrt(Age_first), Blen_ata_std = Blen_ata_diff/sqrt(Age_first)) %>% filter(N_spp_ata_std > 0)


write.csv(genera.mito.all.ata.out,"outputs/Genera_Mito_All_Output.csv")
