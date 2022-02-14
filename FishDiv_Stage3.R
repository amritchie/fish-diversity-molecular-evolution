
fam.mito.all.ata.lm <- read.csv(fam.mito.all.ata.out,"outputs\\Fam_Mito_All_Output.csv")
fam.nuc.all.ata.lm <- read.csv(fam.nuc.all.ata.out,"outputs\\Fam_Nuc_RAG1_Total_Output.csv")
genera.mito.all.ata.lm <- read.csv(genera.mito.all.ata.out,"outputs\\Genera_Mito_All_Output.csv")
genera.nuc.all.ata.lm <- read.csv(genera.mito.all.ata.out,"outputs\\Genera_Mito_All_Output.csv")

fam.mito.coding.dS.ata.lm <- read.csv(fam.mito.coding.dS.ata.out,"outputs\\Fam_Mito_Coding_dS_Output.csv")
fam.mito.coding.dN.ata.lm <- read.csv(fam.mito.coding.dN.ata.out,"outputs\\Fam_Mito_Coding_dN_Output.csv")
fam.nuc.all.dS.ata.lm <- read.csv(fam.nuc.all.dS.ata.out,"outputs\\Fam_Nuc_RAG1_dS_Output.csv")
fam.nuc.all.dN.ata.lm <- read.csv(fam.nuc.all.dN.ata.out,"outputs\\Fam_Nuc_RAG1_dN_Output.csv")
genera.mito.coding.dS.ata.out <- read.csv(genera.mito.coding.dS.ata.out,"outputs\\Genera_Mito_Coding_dS_Output.csv")
genera.mito.coding.dN.ata.out <- read.csv(genera.mito.coding.dN.ata.out,"outputs\\Genera_Mito_Coding_dN_Output.csv")
genera.nuc.all.dS.ata.out <- read.csv(genera.nuc.all.dS.ata.out,"outputs\\Genera_Nuc_RAG1_dS_Output.csv")
genera.nuc.all.dN.ata.out <- read.csv(genera.nuc.all.dN.ata.out,"outputs\\Genera_Nuc_RAG1_dN_Output.csv")

# Conduct regressions
fam.mito.all.ata.lm <- summary(lm(fam.mito.all.ata.out, formula = N_spp_ata_std ~ Blen_ata_std-1))$coefficient
fam.nuc.all.ata.lm <- summary(lm(fam.nuc.all.ata.out, formula = N_spp_ata_std ~ Blen_ata_std-1))$coefficient
genera.mito.all.ata.lm <- summary(lm(genera.mito.all.ata.out, formula = N_spp_ata_std ~ Blen_ata_std-1))$coefficient
genera.nuc.all.ata.lm <- summary(lm(genera.nuc.all.ata.out, formula = N_spp_ata_std ~ Blen_ata_std-1))$coefficient

fam.mito.coding.dS.ata.lm <- summary(lm(fam.mito.coding.dS.ata.out, formula = N_spp_ata_std ~ Blen_ata_std-1))$coefficient
fam.mito.coding.dN.ata.lm <- summary(lm(fam.mito.coding.dN.ata.out, formula = N_spp_ata_std ~ Blen_ata_std-1))$coefficient
fam.nuc.all.dS.ata.lm <- summary(lm(fam.nuc.all.dS.ata.out, formula = N_spp_ata_std ~ Blen_ata_std-1))$coefficient;
fam.nuc.all.dN.ata.lm <- summary(lm(fam.nuc.all.dN.ata.out, formula = N_spp_ata_std ~ Blen_ata_std-1))$coefficient;
genera.mito.coding.dS.ata.lm <- summary(lm(genera.mito.coding.dS.ata.out, formula = N_spp_ata_std ~ Blen_ata_std-1))$coefficient
genera.mito.coding.dN.ata.lm <- summary(lm(genera.mito.coding.dN.ata.out, formula = N_spp_ata_std ~ Blen_ata_std-1))$coefficient
genera.nuc.all.dS.ata.lm <- summary(lm(genera.nuc.all.dS.ata.out, formula = N_spp_ata_std ~ Blen_ata_std-1))$coefficient
genera.nuc.all.dN.ata.lm <- summary(lm(genera.nuc.all.dN.ata.out, formula = N_spp_ata_std ~ Blen_ata_std-1))$coefficient


# Format output
fish2.ata.pairnums <- c(dim(fam.mito.all.ata.out)[1],
                        dim(fam.nuc.all.ata.out)[1],
                        dim(genera.mito.all.ata.out)[1],
                        dim(genera.nuc.all.ata.out)[1],
                        dim(fam.mito.coding.dS.ata.out)[1],
                        dim(fam.mito.coding.dN.ata.out)[1],
                        dim(fam.nuc.all.dS.ata.out)[1],
                        dim(fam.nuc.all.dN.ata.out)[1],
                        dim(genera.mito.coding.dS.ata.out)[1],
                        dim(genera.mito.coding.dN.ata.out)[1],
                        dim(genera.nuc.all.dS.ata.out)[1],
                        dim(genera.nuc.all.dN.ata.out)[1]
)

fish2.ata.lms.out <- rbind(fam.mito.all.ata.lm,
                           fam.nuc.all.ata.lm,
                           genera.mito.all.ata.lm ,
                           genera.nuc.all.ata.lm,
                           fam.mito.coding.dS.ata.lm,
                           fam.mito.coding.dN.ata.lm,
                           fam.nuc.all.dS.ata.lm,
                           fam.nuc.all.dN.ata.lm,
                           genera.mito.coding.dS.ata.lm,
                           genera.mito.coding.dN.ata.lm,
                           genera.nuc.all.dS.ata.lm,
                           genera.nuc.all.dN.ata.lm
)
fish2.ata.lms.out.t <- as_tibble(fish2.ata.lms.out)
fish2.ata.lms.out.t$Rank <- c(rep("Families", 2), 
                              rep("Genera", 2),
                              rep("Families", 4),
                              rep("Genera", 4)
)
fish2.ata.lms.out.t$Sequence <- c(rep(c("mito.all.All", "nuc.all"), 2),
                                  rep(c(rep("mito.coding",2), rep("nuc.all",2)),2)
)
fish2.ata.lms.out.t$SubstitutionType <- c(rep("Total", 4),
                                          rep(c("dS", "dN"), 4)
)
fish2.ata.lms.out.t$NumPairs <- fish2.ata.pairnums

fish2.ata.lms.out.t <- fish2.ata.lms.out.t  %>% select(Rank,Sequence,SubstitutionType,Estimate,"Std. Error",
                                                       "t value", "Pr(>|t|)", NumPairs)