# File: phylo_halves.R

phylo.secondhalf <- function (phy) 
{
  res <- my.drop.tip(phy, 
                     phy$edge[1:Nedge(phy) < which(phy$edge[,1] == Ntip(phy) + 1)[2] & 
                                phy$edge[,2] <= Ntip(phy), 2], 
                     root.edge = 1)
  return (res)
}

phylo.firsthalf <- function (phy)
{
  res <- my.drop.tip(phy, 
                     phy$edge[1:Nedge(phy) >= which(phy$edge[,1] == Ntip(phy) + 1)[2] & 
                                phy$edge[,2] <= Ntip(phy), 2], 
                     root.edge = 1)
  return (res)
}

phylo.sample.even <- function (phy) {
  fst <- phylo.firsthalf(phy)$tip.label
  snd <- phylo.secondhalf(phy)$tip.label
  return(keep.tip(phy, c(sample(fst, min(length(fst),length(snd))), 
                         sample(snd, min(length(fst),length(snd)))
  )
  )
  )
}
