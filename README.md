# fish-diversity-molecular-evolution
Data and code for Ritchie et al. (2022) Diversification rate is associated with rate of molecular evolution in ray-finned fish (Actinopterygii)

The final outputs used for analysis in the paper are in outputs/Final_sister_pair_tables.

The main R code is run in three stages, which are to be sourced in order. Stage 1 produces sister pair selections and downloads molecular sequences. This results in a table of taxa linked to sequence accessions and pair assignments. Sister pairs are then manually grouped together into quartets, associated quartet-wise alignments are produced, and PAML programs are run to infer branch lengths. 

In Stage 2, the PAML output is read, clade sizes are calculated and the final tables of sister pair comparisons are produced. 

In Stage 3, data transformations are performed and linear models are run.
