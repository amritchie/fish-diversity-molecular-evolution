library(ape)
library(fs)
library(tidyverse)
library(stringr)
library(phytools)
library(Matrix)
library(matrixcalc)
library(rfishbase)

source("R\\get_cherries.R")

outgroup.fams <- c("Rajidae",
                   "Callorhinchidae",
                   "Latimeriidae",
                   "Neoceratodontidae",
                   "Protopteridae",
                   "Lepidosirenidae",
                   "Pipidae",
                   "Didelphidae",
                   "Muridae",
                   "Hominidae"
)
noncoding.markers <- c("16S", "hoxc6a")

get.cherries <- function (phy) sapply(which(node.depth(phy) == 2), function (x) phy$tip.label[getDescendants(phy, x)], simplify = F)

etol.tree <- read.tree("trees\\RAxMLThree_Plus_24_part.tre")
etol.aln <- read.nexus.data("alignments\\ETOL_Concat_Three_Plus.nex")
etol.parts <- read.csv("etol_partitions.csv", stringsAsFactors = F)


# ALIGNMENT AND COVERAGE -------------------------------------------------------

etol.aln.m <- t(sapply(etol.aln, identity))
etol.gaps <- matrix(etol.aln.m %in% c('A','C','T','G','a','c','g','t'), 
                    nrow = dim(etol.aln.m)[1])
rownames(etol.gaps) = rownames(etol.aln.m)
etol.covs <- mapply(function (x,y) etol.gaps[,x:y] %*% rep(1,y-x+1), 
                    etol.parts$Start, etol.parts$Stop)
colnames(etol.covs) <- etol.parts$Marker
rownames(etol.covs) <- rownames(etol.gaps)
etol.counts <- etol.covs %*% rep(1, dim(etol.covs)[2])
rownames(etol.counts) <- rownames(etol.covs)

# PAIR SELECTION ---------------------------------------------------------------

etol.accessions <- read.csv("ETOL_TABLE_S1.csv", stringsAsFactors = F)
etol.accessions.tips <- tibble(Label=etol.tree$tip.label) %>% 
  left_join(etol.accessions, by="Label")

## Check and filter for family monophyly

etol.accessions.tips$FamIsMono <- sapply(etol.accessions.tips$Family, 
                                         function (x) is.monophyletic(etol.tree, 
                                                                      etol.accessions.tips %>% 
                                                                        filter(Family == x) %>% 
                                                                        pull(Label)
                                                                      )
                                         )

etol.tips.tokeep <- etol.accessions.tips %>% 
  filter(FamIsMono & !(Family %in% outgroup.fams)) %>% 
  group_by(Family) %>% 
  summarise(Label=first(Label))

## Make a tree with one tip per family

etol.family.tree <- keep.tip(etol.tree, etol.tips.tokeep$Label)
etol.family.tree$tip.label <- sapply(etol.family.tree$tip.label, 
                                     function (x) etol.tips.tokeep$Family[etol.tips.tokeep$Label==x])

## Select clades with size 2 as family pairs

etol.family.pairs <- get.cherries(etol.family.tree)
etol.tip.pairs <- sapply(etol.family.pairs, 
                         function (x) sapply(x, 
                                             function (y) etol.accessions.tips %>% 
                                               filter(Family == y) %>% 
                                               pull(Label), 
                                             simplify=F), 
                         simplify=F
                         )

for (i in 1:length(etol.tip.pairs)) 
{
  names(etol.tip.pairs[[i]]) <- etol.family.pairs[[i]]
}

## Node density effect - sample larger clade to even tip numbers

etol.tip.pairs.sampled <- list()

for (pair in etol.tip.pairs)
{
  sislengths  <- sapply(pair, length)
  bigsister <- which.max(sislengths)
  littlesister <- c(1,2)[-bigsister]
  
  bigsister.counts <- etol.counts[pair[[bigsister]],1]
  bigsister.sampled <- pair[[bigsister]][order(bigsister.counts, decreasing = T)][1:min(sislengths)]

  new.pair <- list(list(bigsister.sampled, pair[[littlesister]]))
  
  names(new.pair[[1]]) <- c(names(pair)[[bigsister]], names(pair)[[littlesister]])
  
  etol.tip.pairs.sampled <- c(etol.tip.pairs.sampled, new.pair)
}


# OUTGROUPS --------------------------------------------------------------------

## Select nearest tip to pair root with most available sequence as outgroup

outgroup.ids <- character(length(etol.tip.pairs.sampled))
for (i in 1:length(etol.tip.pairs.sampled))
{
  etol.tree.seq <- keep.tip(etol.tree, intersect(etol.tree$tip.label, rownames(etol.counts)))
  pairmrca <- getMRCA(etol.tree.seq, unlist(etol.tip.pairs.sampled[[i]]))
  parent.node <- getParent(etol.tree.seq, pairmrca)
  outgroup.parent.node <- etol.tree.seq$edge[etol.tree.seq$edge[,1] == parent.node, 2]
  outgroup.descs <- getDescendants(etol.tree.seq, 
                                   outgroup.parent.node[outgroup.parent.node != pairmrca]
  )
  poss.outgroups <- etol.tree.seq$tip.label[
    outgroup.descs[
      outgroup.descs <= Ntip(etol.tree.seq)
    ]
  ]
  outgroup.ids[i] <- names(which.max(etol.counts[poss.outgroups,1]))
  
}

etol.tip.pairs.with.outgroups <- mapply(function (x, y) c(unlist(x), y), etol.tip.pairs.sampled, outgroup.ids, SIMPLIFY=F)

## Create constraint trees for PAML and label sister clades 
## with separate branch models for codeml

etol.pair.trees <- sapply(etol.tip.pairs.with.outgroups, function (x) keep.tip(etol.tree, x), simplify=F)
for (i in seq_along(etol.pair.trees)){etol.pair.trees[[i]]$node.label<-NULL}

etol.pair.trees.lab <- list()
for ( i in seq_along(etol.tip.pairs.sampled))
{
  labtree <- etol.pair.trees[[i]]
  currpair <- etol.tip.pairs.sampled[[i]]
  if (Ntip(labtree) > 3)
  {
    mrcas <- sapply(currpair, function (x) getMRCA(labtree, x))
    nodelabs <- character(Nedge(labtree) + 1)
    nodelabs[mrcas] <- c('$1', '$2')
    labtree$node.label <- nodelabs[(Ntip(labtree) + 1):length(nodelabs)]
  } else {
    tipnames <- unlist(currpair)
    whichtips <- which(labtree$tip.label %in% tipnames)
    labtree$tip.label[whichtips] <- paste(labtree$tip.label[whichtips], c('#1', '#2'))
  }
  etol.pair.trees.lab[[i]] <- labtree
}

## Make sister pair alignments
etol.good.poses <- unlist(apply(etol.parts[1:19,], 1, function (x) x[2]:x[3]))
etol.pair.alns <- sapply(etol.tip.pairs.with.outgroups, 
                         function (x) etol.aln.m[x, unlist(etol.good.poses)], 
                         simplify = F)

# UNGAP ALIGNMENT --------------------------------------------------------------

## Trim alignments so only codons that are complete in the outgroup and at least
## one tip in both sister clades remain

etol.pair.alns.ungap <- list()

for (aln in etol.pair.alns) {
  ncodon <- dim(aln)[2] %/% 3
  seq.p.pair <- dim(aln)[1] %/% 2
  print(paste("Trimming alignment with", ncodon, "codons and", seq.p.pair*2, "sequences plus outgroup") )
  aln.m <- matrix(as.numeric(aln %in% c('A', 'T', 'C', 'G', 'a', 't', 'c', 'g')), nrow = dim(aln)[1])
  good.codons <- !((!aln.m) %&% (diag(ncodon) %x% rep(1,3)))
  good.sites <- !(t(rep(1,3)) %&% !(t(direct.sum(diag(2) %x% rep(1, seq.p.pair), 1)) %&% good.codons)) %x% rep(1,3)
  etol.pair.alns.ungap <- c(etol.pair.alns.ungap, list(aln[,as.logical(good.sites)]))
  
}

## Filter out sequences that have been reduced to nothing, 
## write alignments and trees to file

etol.tip.pairs.ungap <- list()

for (i in 1:length(etol.pair.alns.ungap)) {
  
  if (dim(etol.pair.alns.ungap[[i]])[2] == 0)
  {
    print(paste("pair", i, "has no sequence remaining, skipping"))
    next
  }
  
  etol.tip.pairs.ungap <- c(etol.tip.pairs.ungap, list(etol.tip.pairs.sampled[[i]]))
  
  print(paste0("alignments//raw_sequence//ETOL//ETOL_", 
               i, 
               ".fasta"))
  write.FASTA(as.DNAbin(etol.pair.alns.ungap[[i]]), 
              paste0("alignments//raw_sequence//ETOL//ETOL_", 
                     i, 
                     ".fasta")
  )
  
  write.tree(etol.pair.trees[[i]], 
             paste0("trees//ETOL//ETOL_",
                    i,
                    ".tre")
  )
  
  trstr <- gsub("_\\$", " \\$", write.tree(etol.pair.trees.lab[[i]]))
  write(trstr, file = paste0("trees//ETOL//ETOL_",
                             i,
                             "_lab.tre")
  )
  
}
## Source to here before STAGE 2

# ------------------------------------------------------------------------------
# For STAGE 2

## Run codeml and baseml on trees and alignments above
## Then source the following

etol.baseml.trees <- sapply(ordered.ls.dir("outputs\\PAML_output\\Fam_Nuc_ETOL", regexp = ".*baseml.txt"), baseml.extract.tree, simplify = F)
etol.codeml.trees <- sapply(ordered.ls.dir("outputs\\PAML_output\\Fam_Nuc_ETOL", regexp = ".*codeml.txt"), codeml.extract.trees, simplify = F)

etol.blens <- mapply(function (x,y,q) sapply(c(list(q),y[2:3]), function (z) sapply(x, function (a) phylo.average.brlen(my.drop.tip(z, z$tip.label[str_detect(paste(a, collapse=' '), z$tip.label, negate = T)], root.edge = 1)))), etol.tip.pairs.ungap, etol.codeml.trees, etol.baseml.trees)
etol.famnames <- sapply(etol.tip.pairs.ungap, names)

etol.valid.spp <- apply(etol.famnames, c(1,2),function (x) length(species_list(Family=x)))
etol.pair.blensums <- sapply(etol.tip.pairs.ungap, function (x) sum(keep.tip(etol.tree, unlist(x))$edge.length))

etol.out <- t(rbind(1:(length(etol.tip.pairs.ungap)), etol.famnames, etol.blens, etol.valid.spp, etol.pair.blensums))
colnames(etol.out) <- c("Pair_no", "Family_first", "Family_last", "dS_first", "dS_last", "dN_first", "dN_last", "blen_first", "blen_last","N_spp_first", "N_spp_last", "Blensum")
etol.out <- as_tibble(etol.out) %>% mutate_at(vars(blen_first, blen_last, dN_first, dN_last, dS_first, dS_last, N_spp_first, N_spp_last, Blensum), as.numeric, .keep="unused")
etol.out.t <- as_tibble(etol.out)

etol.analysis.blen <- etol.out.t %>% 
  filter(N_spp_first > 0, N_spp_last > 0) %>%
  mutate(diff = log(blen_first) - log(blen_last)) %>% 
  welch.test(diff, Blensum) %>% 
  select(-diff) %>% 
  mutate_if(is.numeric, log) %>% 
  mutate(blen_std = abs(blen_first - blen_last), 
         N_spp_std = sign(blen_first - blen_last) * (N_spp_first - N_spp_last)
  )

etol.analysis.dS <- etol.out.t %>% 
  filter(N_spp_first > 0, N_spp_last > 0) %>%
  mutate(diff = log(dS_first) - log(dS_last)) %>% 
  welch.test(diff, Blensum) %>% 
  select(-diff) %>% 
  mutate_if(is.numeric, log) %>% 
  mutate(dS_std = abs(dS_first - dS_last), 
         N_spp_std = sign(dS_first - dS_last) * (N_spp_first - N_spp_last)
  )

etol.analysis.dN <- etol.out.t %>% 
  filter(N_spp_first > 0, N_spp_last > 0) %>%
  mutate(diff = log(dN_first) - log(dN_last)) %>% 
  welch.test(diff, Blensum) %>% 
  select(-diff) %>% 
  mutate_if(is.numeric, log) %>% 
  mutate(dN_std = abs(dN_first - dN_last), 
         N_spp_std = sign(dN_first - dN_last) * (N_spp_first - N_spp_last)
  )

etol.blen.lm <- summary(lm(data=etol.analysis.blen, formula=N_spp_std ~ blen_std-1))$coeff
etol.dS.lm <- summary(lm(data=etol.analysis.dS, formula=N_spp_std ~ dS_std-1))$coeff
etol.dN.lm <- summary(lm(data=etol.analysis.dN, formula=N_spp_std ~ dN_std-1))$coeff

etol.blen.pl <- ggplot(etol.analysis.blen) + geom_point(aes(x = blen_std, y = N_spp_std)) + xlab("Family Contrast in Average Total Substitutions/Site") +  ylab("Family contrast in clade size")+ theme(panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.5), text = element_text(family="serif")) + geom_abline(linetype = "dashed", size=1, colour="red", slope=etol.blen.lm[1], intercept = 0)
etol.dS.pl <- ggplot(etol.analysis.dS) + geom_point(aes(x = dS_std, y = N_spp_std)) + xlab("Family Contrast in Average Synonymous Substitutions/Site") +  ylab("Family contrast in clade size")+ theme(panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.5), text = element_text(family="serif")) + geom_abline(linetype = "dashed", size=1, colour="red", slope=etol.dS.lm[1], intercept = 0)
etol.dN.pl <- ggplot(etol.analysis.dN) + geom_point(aes(x = dN_std, y = N_spp_std)) + xlab("Family Contrast in Average Nonsynonymous Substitutions/Site") +  ylab("Family contrast in clade size")+ theme(panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.5), text = element_text(family="serif")) + geom_abline(linetype = "dashed", size=1, colour="red", slope=etol.dN.lm[1], intercept = 0)
