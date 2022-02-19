library(ape)
library(tidyverse)
library(fs)
library(foreach)
library(doParallel)

source("R/FishDiv_Analyse_Pair_Data.R")
source("R/baseml.extract.tree.R")
source("R/phylo.average.brlen.R")
source("R/ordered.ls.dir.R")
source("R/ordered.ls.dir.R")
source("R/welch_test.R")

fishgenes_genera_accessions <- read.csv("fishgenes_genera_accessions_all.csv")
rabosky.ata.trees <- read.tree("trees/Rabosky_actinopt_full.trees")
fam.mito.coding.accspairs.t2 <- read.csv("outputs/Stage2_Accessions_pairs_tables/fish2_fam_with_mito_coding_accspairs_t_2.csv", stringsAsFactors = F)
fam.mito.coding.pairclades.all <- sapply(
  ordered.ls.dir("trees/Stage2_pair_trees/pairclades_all_species/Fam_with_mito_coding/", 
                 regexp=".*tre$"), 
  read.tree, 
  simplify=F
)

fam.mito.coding.pairclades.sampled<- sapply(
  ordered.ls.dir("trees/Stage2_pair_trees/pairclades_sampled/Fam_with_mito_coding/", 
                 regexp=".*tre$"), 
  read.tree, 
  simplify=F
)


fam.mito.coding.codeml.trees <- sapply(ordered.ls.dir("outputs/PAML_output/Fam_Mito_Coding/codeml", regexp = ".*txt"), 
										codeml.extract.trees, simplify = F)

fam.mito.coding.dS.out <- fishdiv.analyse.pair.data(fam.mito.coding.accspairs.t2, 
													fishgenes_genera_accessions, 
													sapply(fam.mito.coding.codeml.trees, 
															function (x) x$dS, simplify = F), fam.mito.coding.pairclades.all, "Genus")
fam.mito.coding.dN.out <- fishdiv.analyse.pair.data(fam.mito.coding.accspairs.t2, 
													fishgenes_genera_accessions, 
													sapply(fam.mito.coding.codeml.trees, 
															function (x) x$dN, simplify = F), 
													fam.mito.coding.pairclades.all, "Genus")

##############################

## Rectifying N_spp with ata trees
## Commented out as this part can take several hours. Pre-baked sets of ATA tree
## pairs are saved in the repository.
#
# ncores = detectCores()
# registerDoParallel(ncores)


# fam.mito.coding.dS.pairs.ata <- foreach(x=1:(dim(fam.mito.coding.dS.out)[1]), .packages="ape") %dopar% {sapply(rabosky.ata.trees, function(tr) sapply(c(phylo.firsthalf, phylo.secondhalf), function (f) f(my.extract.clade(phy=tr, node=my.getMRCA(tr, fam.mito.coding.pairclades.sampled[[fam.mito.coding.dS.out$Pair[x]]]$tip.label)))$tip.label, simplify = F), simplify = F)}
# fam.mito.coding.dN.pairs.ata <- foreach(x=1:(dim(fam.mito.coding.dN.out)[1]), .packages="ape") %dopar% {sapply(rabosky.ata.trees, function(tr) sapply(c(phylo.firsthalf, phylo.secondhalf), function (f) f(my.extract.clade(phy=tr, node=my.getMRCA(tr, fam.mito.coding.pairclades.sampled[[fam.mito.coding.dN.out$Pair[x]]]$tip.label)))$tip.label, simplify = F), simplify = F)}

load("ata_tree_pairs_prebaked/tmp_ata_pairs_welch/fish2_fam_with_mito_coding_dS_pairs_ata")
load("ata_tree_pairs_prebaked/tmp_ata_pairs_welch/fish2_fam_with_mito_coding_dN_pairs_ata")
# If using monophyly filter:
# load("ata_tree_pairs_prebaked/tmp_ata_pairs_mono_welch/fish2_fam_with_mito_coding_dS_pairs_ata")
# load("ata_tree_pairs_prebaked/tmp_ata_pairs_mono_welch/fish2_fam_with_mito_coding_dN_pairs_ata")

## Orientation of sister clades in ATA trees is not guaranteed
fam.mito.coding.dS.pairs.ata.sorted <- sapply(1:(dim(fam.mito.coding.dS.out)[1]), 
											function (x) sapply(fish2.fam.with.mito.coding.dS.pairs.ata[[x]], 
																function (y) if(length(intersect(y[[1]], phylo.firsthalf(fam.mito.coding.pairclades.sampled[[fam.mito.coding.dS.out$Pair[x]]])$tip.label)) > 0) y else (rev(y)), simplify = F), 
																simplify = F)
fam.mito.coding.dN.pairs.ata.sorted <- sapply(1:(dim(fam.mito.coding.dN.out)[1]), 
												function (x) sapply(fish2.fam.with.mito.coding.dN.pairs.ata[[x]], 
																function (y) if(length(intersect(y[[1]], phylo.firsthalf(fam.mito.coding.pairclades.sampled[[fam.mito.coding.dN.out$Pair[x]]])$tip.label)) > 0) y else (rev(y)), simplify = F), simplify = F)

## Get clade sizes as average species numbers across ATA trees
fam.mito.coding.dS.nspp.ata.sorted <- sapply(fam.mito.coding.dS.pairs.ata.sorted, 
												function (x) apply(sapply(x, function (y) c(length(y[[1]]), length(y[[2]]))), 1, mean))
fam.mito.coding.dN.nspp.ata.sorted <- sapply(fam.mito.coding.dN.pairs.ata.sorted, 
												function (x) apply(sapply(x, function (y) c(length(y[[1]]), length(y[[2]]))), 1, mean))

## Format final data sets with substitutions and clade sizes

## NO MONOPHYLY FILTER
fam.mito.coding.dS.ata.out <- fam.mito.coding.dS.out %>% 
								ungroup() %>% mutate(N_spp_ata_first = log(fam.mito.coding.dS.nspp.ata.sorted[1,]), N_spp_ata_last = log(fam.mito.coding.dS.nspp.ata.sorted[2,])) %>% 
								mutate(swap_ata = sign(N_spp_ata_first - N_spp_ata_last)) %>% 
								mutate(N_spp_ata_diff = abs(N_spp_ata_first - N_spp_ata_last), Blen_ata_diff = swap_ata*(Blen_first - Blen_last)) %>% 
								mutate(N_spp_ata_std = N_spp_ata_diff/sqrt(Age_first), Blen_ata_std = Blen_ata_diff/sqrt(Age_first))  %>% 
								filter(N_spp_ata_std > 0)
fam.mito.coding.dN.ata.out <- fam.mito.coding.dN.out %>% 
								ungroup() %>% 
								mutate(N_spp_ata_first = log(fam.mito.coding.dN.nspp.ata.sorted[1,]), N_spp_ata_last = log(fam.mito.coding.dN.nspp.ata.sorted[2,])) %>% 
								mutate(swap_ata = sign(N_spp_ata_first - N_spp_ata_last)) %>% 
								mutate(N_spp_ata_diff = abs(N_spp_ata_first - N_spp_ata_last), Blen_ata_diff = swap_ata*(Blen_first - Blen_last)) %>%
								mutate(N_spp_ata_std = N_spp_ata_diff/sqrt(Age_first), Blen_ata_std = Blen_ata_diff/sqrt(Age_first))  %>% 
								filter(N_spp_ata_std > 0)

## MONOPHYLY FILTER - run this code instead of the above
# fam.mito.coding.nonmono.pairs.ata <- sapply(fam.mito.coding.pairs.ata.sorted, function (x) sum(sapply(x, function (y) length(intersect(str_extract(y[[1]], "^.*_"), str_extract(y[[2]], "^.*_"))) > 0)) > 20)
# fam.mito.coding.dS.nonmono.pairs.ata <- sapply(fam.mito.coding.dS.pairs.ata.sorted, function (x) sum(sapply(x, function (y) length(intersect(str_extract(y[[1]], "^.*_"), str_extract(y[[2]], "^.*_"))) > 0)) > 20)
# fam.mito.coding.dN.nonmono.pairs.ata <- sapply(fam.mito.coding.dN.pairs.ata.sorted, function (x) sum(sapply(x, function (y) length(intersect(str_extract(y[[1]], "^.*_"), str_extract(y[[2]], "^.*_"))) > 0)) > 20)
# 
# fam.mito.coding.ata.out <- fam.mito.coding.out %>% ungroup() %>% mutate(N_spp_ata_first = log(fam.mito.coding.nspp.ata.sorted[1,]), N_spp_ata_last = log(fam.mito.coding.nspp.ata.sorted[2,])) %>% filter(!fam.mito.coding.nonmono.pairs.ata) %>% mutate(swap_ata = sign(N_spp_ata_first - N_spp_ata_last)) %>% mutate(N_spp_ata_diff = abs(N_spp_ata_first - N_spp_ata_last), Blen_ata_diff = swap_ata*(Blen_first - Blen_last)) %>% mutate(N_spp_ata_std = N_spp_ata_diff/sqrt(Age_first), Blen_ata_std = Blen_ata_diff/sqrt(Age_first)) %>% filter(N_spp_ata_std > 0)
# fam.mito.coding.dS.ata.out <- fam.mito.coding.dS.out %>% ungroup() %>% mutate(N_spp_ata_first = log(fam.mito.coding.dS.nspp.ata.sorted[1,]), N_spp_ata_last = log(fam.mito.coding.dS.nspp.ata.sorted[2,])) %>% filter(!fam.mito.coding.dS.nonmono.pairs.ata) %>% mutate(swap_ata = sign(N_spp_ata_first - N_spp_ata_last)) %>% mutate(N_spp_ata_diff = abs(N_spp_ata_first - N_spp_ata_last), Blen_ata_diff = swap_ata*(Blen_first - Blen_last)) %>% mutate(N_spp_ata_std = N_spp_ata_diff/sqrt(Age_first), Blen_ata_std = Blen_ata_diff/sqrt(Age_first))  %>% filter(N_spp_ata_std > 0);
# fam.mito.coding.dN.ata.out <- fam.mito.coding.dN.out %>% ungroup() %>% mutate(N_spp_ata_first = log(fam.mito.coding.dN.nspp.ata.sorted[1,]), N_spp_ata_last = log(fam.mito.coding.dN.nspp.ata.sorted[2,])) %>% filter(!fam.mito.coding.dN.nonmono.pairs.ata) %>% mutate(swap_ata = sign(N_spp_ata_first - N_spp_ata_last)) %>% mutate(N_spp_ata_diff = abs(N_spp_ata_first - N_spp_ata_last), Blen_ata_diff = swap_ata*(Blen_first - Blen_last)) %>% mutate(N_spp_ata_std = N_spp_ata_diff/sqrt(Age_first), Blen_ata_std = Blen_ata_diff/sqrt(Age_first))  %>% filter(N_spp_ata_std > 0);


write.csv(fam.mito.coding.dS.ata.out,"outputs/Fam_Mito_Coding_dS_Output.csv")
write.csv(fam.mito.coding.dN.ata.out,"outputs/Fam_Mito_Coding_dN_Output.csv")


