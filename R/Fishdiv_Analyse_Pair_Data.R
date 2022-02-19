
fishdiv.analyse.pair.data <- function (qtable, lifetable, est.trees, lifetrees, taxon){

  fish2.life <- fishgenes_genera_accessions %>% 
    group_by(Species) %>% 
    summarise(Family = first(family), 
              Genus = first(genus), 
              Length = mean(log(Length), na.rm = T), 
              Depth = mean(log(DepthRangeDeep), na.rm = T), 
              Lat = mean(log(abs(lat_centroid)), na.rm = T)
    )

  clades <- unique(qtable$clade_name)
  pairs <- as.integer(sapply(str_split(clades, "_"), function (x) x[[1]]))

  fams <- character(0) 
  genera <- character(0)
  lifeout <- array(0, c(0,3))
  blen <- numeric(0)
  validspp <- integer(0)
  pairages <- integer(0)
  
  for (sister in clades)
  {
    spp <- gsub(" ", "_", qtable %>% 
                  filter(clade_name == sister) %>%
                  pull(Species))

    pair <- as.integer(str_split(as.character(sister), "_")[[1]][1])
    pairages <- c(pairages, node.depth.edgelength(lifetrees[[pair]])[1])

    halves <- list(phylo.firsthalf(lifetrees[[pair]]), 
                   phylo.secondhalf(lifetrees[[pair]]))

    myhalf <- halves[sapply(halves, 
                             function (x) setequal(intersect(x$tip.label, spp), spp))][[1]]

    lifestats <- fish2.life %>% 
      filter(Species %in% gsub("_", " ", myhalf$tip.label), !is.infinite(Depth)) %>% 
      group_by(Genus) %>% 
      summarise(Family = first(Family),
                Length = mean(Length, na.rm = T), 
                Depth = mean(Depth, na.rm = T), 
                Lat = mean(Lat, na.rm = T))

    lifeout <- rbind(lifeout, 
                     c(mean(lifestats$Length, na.rm =T), 
                       mean(lifestats$Depth, na.rm = T), 
                       mean(lifestats$Lat, na.rm = T)))

    genera <- c(genera, paste(unique(lifestats$Genus), collapse = ","))
    fams <- c(fams, paste(unique(lifestats$Family), collapse = ","))
    if (taxon == "Genus") {
      validspp <- c(validspp, dim(fish2.life %>% 
                                     filter(Genus %in% unique(lifestats$Genus)))[1])
    } else {
      validspp <- c(validspp, dim(fish2.life %>% 
                                           filter(Family %in% unique(lifestats$Family)))[1])
    }

    quartets <- strsplit((qtable %>% filter(clade_name == sister))[1,] %>%
      pull(quartet_str), ",")[[1]]

    qtrees <- sapply(quartets, 
                     function (x) est.trees[[as.integer(x)]], 
                     simplify= F)

    qtrees <- qtrees[sapply(qtrees,
                            function (x) !is.null(x) & class(x) != "numeric")]

    halves <- sapply(qtrees,
                     function (qtree) my.drop.tip(qtree,
                                                  qtree$tip.label[!(qtree$tip.label %in% spp)],
                                                  root.edge = 1), 
                     simplify = F)

    blen <- c(blen, mean(sapply(halves, phylo.average.brlen)))

  }
  print(length(pairs)/2)
  out.t <- tibble(
    Family = fams,
    Genus = genera,
    Sister = clades,
    Pair = pairs,
    Blen = log(blen),
    Length = lifeout[,1],
    Depth = lifeout[,2],
    Lat = lifeout[,3],
    N_spp = log(validspp),
    Age = pairages
  ) %>% arrange(Sister)%>% 
    group_by(Pair) %>% 
    summarise_all(c("first" = first, "last" = last)) %>%
    mutate(swap = ifelse(N_spp_first > N_spp_last, 1, -1),
           N_spp_diff = (N_spp_first - N_spp_last) * swap,
           Blen_diff = (Blen_first - Blen_last) * swap,
           Length_diff = (Length_first - Length_last) * swap,
           Depth_diff = (Depth_first - Depth_last) * swap,
           Lat_diff = (Lat_first - Lat_last) * swap) %>% 
    mutate(Blen_std = Blen_diff/sqrt(Age_first), 
           N_spp_std = N_spp_diff/sqrt(Age_first),
           Length_std = Length_diff/sqrt(Age_first),
           Depth_std = Depth_diff/sqrt(Age_first),
           Lat_std = Lat_diff/sqrt(Age_first)) %>%
    filter(Blen_first < log(2) & Blen_last < log(2) & !is.infinite(Blen_diff))
  print(dim(out.t)[1])

## MONOPHYLY FILTER - COMMENT OUT FOR 'FULL' DATASET, RESTORE FOR 'FILTERED' DATASET##
  # if (taxon == "Genus") {
  #   out.t <- out.t[sapply(1:length(out.t$Sister_first),
  #                         function (x) length(intersect(str_split(out.t[x,"Genus_first"],",")[[1]],
  #                                                str_split(out.t[x,"Genus_last"],",")[[1]])) <= 0),
  #                  ]
  # } else {
  #   out.t <- out.t[sapply(1:length(out.t$Sister_first),
  #                         function (x) length(intersect(str_split(out.t[x,"Family_first"], ",")[[1]],
  #                                                 str_split(out.t[x,"Family_last"],",")[[1]])) <= 0),
  #                  ]
  # 
  # }
## END COMMENT OUT ##
  print(dim(out.t)[1])
  out.t <- out.t %>% welch.test(Blen_diff, Age_first) %>% means.test(Blen_first, Blen_last)
  print(dim(out.t)[1])

  return(out.t)

}

