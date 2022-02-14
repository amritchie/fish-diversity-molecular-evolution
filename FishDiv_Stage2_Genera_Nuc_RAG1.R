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
genera.nuc.all.accspairs.t2 <- read.csv("outputs\\Accessions_pairs_tables\\Genera_with_Nuc_RAG1_AccsPairsTable.csv")

genera.nuc.all.baseml.trees <- sapply(ordered.ls.dir("outputs\\PAML_output\\Genera_Nuc_RAG1\\baseml", regexp = ".*txt"), baseml.extract.tree, simplify = F)
genera.nuc.all.codeml.trees <- sapply(ordered.ls.dir("outputs\\PAML_output\\Genera_Nuc_RAG1\\codeml", regexp = ".*txt"), codeml.extract.tree, simplify = F)

genera.nuc.all.out <- fishdiv.analyse.pair.data(genera.nuc.all.accspairs.t2, 
												fishgenes_genera_accessions, 
												genera.nuc.all.baseml.trees, 
												genera.nuc.all.pairclades.all.nodup, "Genus")
genera.nuc.all.dS.out <- fishdiv.analyse.pair.data(genera.nuc.all.accspairs.t2, 
													fishgenes_genera_accessions, 
													sapply(genera.nuc.all.codeml.trees, function (x) x$dS, simplify = F), 
													genera.nuc.all.pairclades.all.nodup, 
													"Genus")
genera.nuc.all.dN.out <- fishdiv.analyse.pair.data(genera.nuc.all.accspairs.t2, 
													fishgenes_genera_accessions, 
													sapply(genera.nuc.all.codeml.trees, function (x) x$dN, simplify = F), 
													genera.nuc.all.pairclades.all.nodup, 
													"Genus")

##############################

## Rectifying N_spp with ata trees
## Commented out as this part can take several hours. Pre-baked sets of ATA tree
## pairs are saved in the repository.
#
# ncores = detectCores()
# registerDoParallel(ncores)
#
# genera.nuc.all.pairs.ata <- foreach(x=1:(dim(genera.nuc.all.out)[1]), .packages="ape") %dopar% {sapply(rabosky.ata.trees, function(tr) sapply(c(phylo.firsthalf, phylo.secondhalf), function (f) f(my.extract.clade(phy=tr, node=my.getMRCA(tr, genera.nuc.all.pairclades.sampled[[genera.nuc.all.out$Pair[x]]]$tip.label)))$tip.label, simplify = F), simplify = F)}
# genera.nuc.all.dS.pairs.ata <- foreach(x=1:(dim(genera.nuc.all.dS.out)[1]), .packages="ape") %dopar% {sapply(rabosky.ata.trees, function(tr) sapply(c(phylo.firsthalf, phylo.secondhalf), function (f) f(my.extract.clade(phy=tr, node=my.getMRCA(tr, genera.nuc.all.pairclades.sampled[[genera.nuc.all.dS.out$Pair[x]]]$tip.label)))$tip.label, simplify = F), simplify = F)}
# genera.nuc.all.dN.pairs.ata <- foreach(x=1:(dim(genera.nuc.all.dN.out)[1]), .packages="ape") %dopar% {sapply(rabosky.ata.trees, function(tr) sapply(c(phylo.firsthalf, phylo.secondhalf), function (f) f(my.extract.clade(phy=tr, node=my.getMRCA(tr, genera.nuc.all.pairclades.sampled[[genera.nuc.all.dN.out$Pair[x]]]$tip.label)))$tip.label, simplify = F), simplify = F)}

## Orientation of sister clades in ATA trees is not guaranteed
genera.nuc.all.pairs.ata.sorted <- sapply(1:(dim(genera.nuc.all.out)[1]), 
											function (x) sapply(genera.nuc.all.pairs.ata[[x]], 
																function (y) if(length(intersect(y[[1]], phylo.firsthalf(genera.nuc.all.pairclades.sampled[[genera.nuc.all.out$Pair[x]]])$tip.label)) > 0) y else (rev(y)), 
																simplify = F), 
											simplify = F
										)
genera.nuc.all.dS.pairs.ata.sorted <- sapply(1:(dim(genera.nuc.all.dS.out)[1]), 
												function (x) sapply(genera.nuc.all.dS.pairs.ata[[x]], 
																	function (y) if(length(intersect(y[[1]], phylo.firsthalf(genera.nuc.all.pairclades.sampled[[genera.nuc.all.dS.out$Pair[x]]])$tip.label)) > 0) y else (rev(y)),
																	simplify = F), 
												simplify = F
											)
genera.nuc.all.dN.pairs.ata.sorted <- sapply(1:(dim(genera.nuc.all.dN.out)[1]), 
												function (x) sapply(genera.nuc.all.dN.pairs.ata[[x]], 
																	function (y) if(length(intersect(y[[1]], phylo.firsthalf(genera.nuc.all.pairclades.sampled[[genera.nuc.all.dN.out$Pair[x]]])$tip.label)) > 0) y else (rev(y)),
																	simplify = F), 
												simplify = F
											)

## Get clade sizes as average species numbers across ATA trees
genera.nuc.all.nspp.ata.sorted <- sapply(genera.nuc.all.pairs.ata.sorted, 
										function (x) apply(sapply(x, function (y) c(length(y[[1]]), length(y[[2]]))), 1, mean))
genera.nuc.all.dS.nspp.ata.sorted <- sapply(genera.nuc.all.dS.pairs.ata.sorted, 
											function (x) apply(sapply(x, function (y) c(length(y[[1]]), length(y[[2]]))), 1, mean))
genera.nuc.all.dN.nspp.ata.sorted <- sapply(genera.nuc.all.dN.pairs.ata.sorted, 
											function (x) apply(sapply(x, function (y) c(length(y[[1]]), length(y[[2]]))), 1, mean))

## Format final data sets with substitutions and clade sizes
genera.nuc.all.ata.out <- genera.nuc.all.out %>% 
							ungroup() %>% 
							mutate(N_spp_ata_first = transfunc.genera(genera.nuc.all.nspp.ata.sorted[1,]), N_spp_ata_last = transfunc.genera(genera.nuc.all.nspp.ata.sorted[2,])) %>% 
							mutate(swap_ata = sign(N_spp_ata_first - N_spp_ata_last)) %>% 
							mutate(N_spp_ata_diff = abs(N_spp_ata_first - N_spp_ata_last), Blen_ata_diff = swap_ata*(Blen_first - Blen_last)) %>%
							mutate(N_spp_ata_std = N_spp_ata_diff/sqrt(Age_first), Blen_ata_std = Blen_ata_diff/sqrt(Age_first)) %>% 
							filter(N_spp_ata_std > 0)
genera.nuc.all.dS.ata.out <- genera.nuc.all.dS.out %>% 
								ungroup() %>% 
								mutate(N_spp_ata_first = transfunc.genera(genera.nuc.all.dS.nspp.ata.sorted[1,]), N_spp_ata_last = transfunc.genera(genera.nuc.all.dS.nspp.ata.sorted[2,])) %>% 
								mutate(swap_ata = sign(N_spp_ata_first - N_spp_ata_last)) %>% 
								mutate(N_spp_ata_diff = abs(N_spp_ata_first - N_spp_ata_last), Blen_ata_diff = swap_ata*(Blen_first - Blen_last)) %>%
								mutate(N_spp_ata_std = N_spp_ata_diff/sqrt(Age_first), Blen_ata_std = Blen_ata_diff/sqrt(Age_first))  %>% 
								filter(N_spp_ata_std > 0)
genera.nuc.all.dN.ata.out <- genera.nuc.all.dN.out %>% 
								ungroup() %>% 
								mutate(N_spp_ata_first = transfunc.genera(genera.nuc.all.dN.nspp.ata.sorted[1,]), N_spp_ata_last = transfunc.genera(genera.nuc.all.dN.nspp.ata.sorted[2,])) %>% 
								mutate(swap_ata = sign(N_spp_ata_first - N_spp_ata_last)) %>% 
								mutate(N_spp_ata_diff = abs(N_spp_ata_first - N_spp_ata_last), Blen_ata_diff = swap_ata*(Blen_first - Blen_last)) %>%
								mutate(N_spp_ata_std = N_spp_ata_diff/sqrt(Age_first), Blen_ata_std = Blen_ata_diff/sqrt(Age_first))  %>% 
								filter(N_spp_ata_std > 0)

write.csv(genera.nuc.all.ata.out,"outputs\\Genera_Nuc_RAG1_Total_Output.csv")
write.csv(genera.nuc.all.dS.ata.out,"outputs\\Genera_Nuc_RAG1_dS_Output.csv")
write.csv(genera.nuc.all.dN.ata.out,"outputs\\Genera_Nuc_RAG1_dN_Output.csv")
