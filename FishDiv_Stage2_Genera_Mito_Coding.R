library(ape)
library(tidyverse)
library(fs)
library(foreach)
library(doParallel)

source("R\\FishDiv_Analyse_Pair_Data.R")
source("R\\baseml_extract_tree.R")
source("R\\phylo_average_blen.R")

fishgenes_genera_accessions <- read.csv("fishgenes_genera_accessions_all.csv")
rabosky.ata.trees <- read.tree("trees\\Rabosky_actinopt_full.trees")
genera.mito.coding.accspairs.t2 <- read.csv("outputs\\Accessions_pairs_tables\\Genera_with_Mito_Coding_AccsPairsTable.csv")

genera.mito.coding.codeml.trees <- sapply(ordered.ls.dir("outputs\\PAML_output\\Genera_Mito_Coding\\codeml", regexp = ".*txt"), 
									codeml.extract.tree, simplify = F)

genera.mito.coding.dN.out <- fishdiv.analyse.pair.data(genera.mito.coding.accspairs.t3, 
														fishgenes_genera_accessions, 
														sapply(genera.mito.all.coding.codeml_trees, function (x) x$dN, simplify = F),
														genera.mito.coding.pairclades.all.nodup, "Genus")
genera.mito.coding.dS.out <- fishdiv.analyse.pair.data(genera.mito.coding.accspairs.t3, 
														fishgenes_genera_accessions, 
														sapply(genera.mito.all.coding.codeml_trees, function (x) x$dS, simplify = F), 
														genera.mito.coding.pairclades.all.nodup, "Genus")

##############################

## Rectifying N_spp with ata trees
## Commented out as this part can take several hours. Pre-baked sets of ATA tree
## pairs are saved in the repository.
#
# ncores = detectCores()
# registerDoParallel(ncores)

# genera.mito.coding.dN.pairs.ata <- foreach(x=1:(dim(genera.mito.coding.dN.out)[1]), .packages="ape") %dopar% {sapply(rabosky.ata.trees, function(tr) sapply(c(phylo.firsthalf, phylo.secondhalf), function (f) f(my.extract.clade(phy=tr, node=my.getMRCA(tr, genera.mito.coding.pairclades.sampled[[genera.mito.coding.dN.out$Pair[x]]]$tip.label)))$tip.label, simplify = F), simplify = F)}
# genera.mito.coding.dS.pairs.ata <- foreach(x=1:(dim(genera.mito.coding.dS.out)[1]), .packages="ape") %dopar% {sapply(rabosky.ata.trees, function(tr) sapply(c(phylo.firsthalf, phylo.secondhalf), function (f) f(my.extract.clade(phy=tr, node=my.getMRCA(tr, genera.mito.coding.pairclades.sampled[[genera.mito.coding.dS.out$Pair[x]]]$tip.label)))$tip.label, simplify = F), simplify = F)}

## Orientation of sister clades in ATA trees is not guaranteed
genera.mito.coding.dS.pairs.ata.sorted <- sapply(1:(dim(genera.mito.coding.dS.out)[1]), 
													function (x) sapply(genera.mito.coding.dS.pairs.ata[[x]], 
																		function (y) if(length(intersect(y[[1]], phylo.firsthalf(genera.mito.coding.pairclades.sampled[[genera.mito.coding.dS.out$Pair[x]]])$tip.label)) > 0) y else (rev(y)), 
																		simplify = F), 
												simplify = F
												)
genera.mito.coding.dN.pairs.ata.sorted <- sapply(1:(dim(genera.mito.coding.dN.out)[1]), 
													function (x) sapply(genera.mito.coding.dN.pairs.ata[[x]], function (y) if(length(intersect(y[[1]], phylo.firsthalf(genera.mito.coding.pairclades.sampled[[genera.mito.coding.dN.out$Pair[x]]])$tip.label)) > 0) y else (rev(y)), simplify = F), 
													simplify = F
												)

## Get clade sizes as average species numbers across ATA trees
genera.mito.coding.dS.nspp.ata.sorted <- sapply(genera.mito.coding.dS.pairs.ata.sorted, 
												function (x) apply(sapply(x, function (y) c(length(y[[1]]), length(y[[2]]))), 1, mean))
genera.mito.coding.dN.nspp.ata.sorted <- sapply(genera.mito.coding.dN.pairs.ata.sorted, 
												function (x) apply(sapply(x, function (y) c(length(y[[1]]), length(y[[2]]))), 1, mean))


## Format final data sets with substitutions and clade sizes
genera.mito.coding.dS.ata.out <- genera.mito.coding.dS.out %>% 
									ungroup() %>% 
									mutate(N_spp_ata_first = transfunc.genera(genera.mito.coding.dS.nspp.ata.sorted[1,]), N_spp_ata_last = transfunc.genera(genera.mito.coding.dS.nspp.ata.sorted[2,])) %>% 
									mutate(swap_ata = sign(N_spp_ata_first - N_spp_ata_last)) %>% 
									mutate(N_spp_ata_diff = abs(N_spp_ata_first - N_spp_ata_last), Blen_ata_diff = swap_ata*(Blen_first - Blen_last)) %>% 
									mutate(N_spp_ata_std = N_spp_ata_diff/sqrt(Age_first), Blen_ata_std = Blen_ata_diff/sqrt(Age_first))  %>% 
									filter(N_spp_ata_std > 0)
genera.mito.coding.dN.ata.out <- genera.mito.coding.dN.out %>% 
									ungroup() %>% 
									mutate(N_spp_ata_first = transfunc.genera(genera.mito.coding.dN.nspp.ata.sorted[1,]), N_spp_ata_last = transfunc.genera(genera.mito.coding.dN.nspp.ata.sorted[2,])) %>% 
									mutate(swap_ata = sign(N_spp_ata_first - N_spp_ata_last)) %>% 
									mutate(N_spp_ata_diff = abs(N_spp_ata_first - N_spp_ata_last), Blen_ata_diff = swap_ata*(Blen_first - Blen_last)) %>% 
									mutate(N_spp_ata_std = N_spp_ata_diff/sqrt(Age_first), Blen_ata_std = Blen_ata_diff/sqrt(Age_first))  %>% 
									filter(N_spp_ata_std > 0)

write.csv(genera.mito.coding.dS.ata.out,"outputs\\Genera_Mito_Coding_dS_Output.csv")
write.csv(genera.mito.coding.dN.ata.out,"outputs\\Genera_Mito_Coding_dN_Output.csv")
