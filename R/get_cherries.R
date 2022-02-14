#File:get_cherries.R

get.cherries <- function (phy) 
{
  sapply(which(node.depth(phy) == 2), 
  function (x) phy$tip.label[getDescendants(phy, x)], 
  simplify = F
  )
}