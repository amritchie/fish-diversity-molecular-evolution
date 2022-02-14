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
fam.mito.all.accspairs.t2 <- read.csv("outputs\\Accessions_pairs_tables\\Fam_with_Mito_All_AccsPairsTable.csv")

fam.mito.all.baseml.trees <- sapply(ordered.ls.dir("outputs\\PAML_output\\Fam_Mito_All\\baseml", regexp = ".*txt"), baseml.extract.tree, simplify = F)

fam.mito.all.out <- fishdiv.analyse.pair.data(fam.mito.all.accspairs.t2, fishgenes_genera_accessions, fam.mito.all.baseml.trees, fam.mito.all.pairclades.all, "Genus")


##############################

## Rectifying N_spp with ata trees
## Commented out as this part can take several hours. Pre-baked sets of ATA tree
## pairs are saved in the repository.
#
# ncores = detectCores()
# registerDoParallel(ncores)
#
# fam.mito.all.pairs.ata <- foreach(x=1:(dim(fam.mito.all.out)[1]), .packages="ape") %dopar% {sapply(rabosky.ata.trees, function(tr) sapply(c(phylo.firsthalf, phylo.secondhalf), function (f) f(my.extract.clade(phy=tr, node=my.getMRCA(tr, fam.mito.all.pairclades.sampled[[fam.mito.all.out$Pair[x]]]$tip.label)))$tip.label, simplify = F), simplify = F)}


## Orientation of sister clades in ATA trees is not guaranteed
fam.mito.all.pairs.ata.sorted <- sapply(1:(dim(fam.mito.all.out)[1]), 
										function (x) sapply(fam.mito.all.pairs.ata[[x]], 
															function (y) if(length(intersect(y[[1]], phylo.firsthalf(fam.mito.all.pairclades.sampled[[fam.mito.all.out$Pair[x]]])$tip.label)) > 0) y else (rev(y)), simplify = F), simplify = F)

## Get clade sizes as average species numbers across ATA trees
fam.mito.all.nspp.ata.sorted <- sapply(fam.mito.all.pairs.ata.sorted, 
										function (x) apply(sapply(x, 
																	function (y) c(length(y[[1]]), length(y[[2]]))), 1, mean))

## Format final data sets with substitutions and clade sizes
fam.mito.all.ata.out <- fam.mito.all.out %>% 
							ungroup() %>% 
							mutate(N_spp_ata_first = transfunc.fam(fam.mito.all.nspp.ata.sorted[1,]), N_spp_ata_last = transfunc.fam(fam.mito.all.nspp.ata.sorted[2,])) %>% 
							mutate(swap_ata = sign(N_spp_ata_first - N_spp_ata_last)) %>% 
							mutate(N_spp_ata_diff = abs(N_spp_ata_first - N_spp_ata_last), Blen_ata_diff = swap_ata*(Blen_first - Blen_last)) %>% mutate(N_spp_ata_std = N_spp_ata_diff/sqrt(Age_first), Blen_ata_std = Blen_ata_diff/sqrt(Age_first)) %>% 
							filter(N_spp_ata_std > 0)

write.csv(fam.mito.all.ata.out,"outputs\\Fam_Mito_All_Output.csv")
