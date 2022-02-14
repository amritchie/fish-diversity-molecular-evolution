require(phytools)

find.xtets.sisterpairs <- function (pairtree)
{
  
  init.xtets <- pairup.tips(pairtree)
  xtet.mrcas <- sapply(init.xtets, function (x) getMRCA(pairtree, x))
  
  if (any(sapply(init.xtets, length) < 2)) {
    lone.pair.idx <- which(sapply(init.xtets, length) < 2)
    lone.pair <- init.xtets[[lone.pair.idx]]
    xtet.mrcas <- xtet.mrcas[-lone.pair.idx]
  } else {
    lone.pair <- integer(0)
  }

  if (length(lone.pair) > 0)
  {
    nnbour <- which.min(
      sapply(xtet.mrcas, function (x) node.pair.dist(pairtree, x, lone.pair))
    )
    init.xtets[[nnbour]] <- c(init.xtets[[nnbour]], lone.pair)

    init.xtets <- init.xtets[1:(length(init.xtets) - 1)]
    xtet.mrcas <- sapply(init.xtets, function (x) getMRCA(pairtree, x))
  }
  
  xtet.root.dists <- sapply(init.xtets, function (x) node.pair.dist(pairtree, x[1], x[2]))

  outliers <- sapply(init.xtets, 
                     function (x) length(x) < 3 & 
                       abs(node.pair.dist(
                         pairtree, x[1], x[2]) > quantile(xtet.root.dists, 0.75)) 
                        #x[1], x[2]) - median(xtet.root.dists)) > 1.5 * IQR(xtet.root.dists)
                     )
  qtree <- keep.tip(pairtree, c(sapply(init.xtets[!outliers], 
                                       function (x) x[1]),
                                unlist(init.xtets[outliers]))
                    )
  final.quartets <- init.xtets

  for (i in init.xtets[outliers]) {
    dists <- sapply(sapply(init.xtets[!outliers], 
                           function (x) x[1]), 
                    function (x) node.pair.dist(qtree,
                                                getMRCA(qtree, which(qtree$tip.label %in% i)),
                                                which(qtree$tip.label == x)))
    nnbor <- unlist(init.xtets[!outliers])[which.min(dists)]

    final.quartets[[which(sapply(final.quartets, function (x) all(x %in% i)))]] <- unlist(union(i, nnbor[1]))
  }

  out.tree <- pairtree
  out.tree$tip.label <-  sapply(pairtree$tip.label, 
                               function(x) which(sapply(final.quartets, 
                                                        function (y) x %in% y)))
  plot(out.tree)
  
  return(final.quartets)
  
}

pairup.tips <- function (phy) 
{
  newpairs <- list()
  tips <- phy$tip.label

  while(length(tips) > 1) {
    
    phy <- keep.tip(phy, which(phy$tip.label %in% tips))
    newpairs <- c(newpairs, get.cherries(phy))
    if (length(tips) > 0) tips <- tips[-which(tips %in% unlist(newpairs))]
    
  }
  
  return(c(newpairs, as.list(tips)))
}


node.pair.dist <- function (phy, i, j) 
{
  sum(phy$edge.length[which(phy$edge[,2] %in% nodepath(phy, i, j))])
} 
