#for creating heatmap using phytools
library(phytools)
packageVersion("phytools")
X <- read.csv("sml.csv", row.names=(1))
X <- X[-13,]
tree <- read.newick("final_newick")
#replacing 0s with very small number to avoid NaN errors
tree$edge.length[tree$edge.length == 0.000] <- 0.00001
sexualMaturity <- as.matrix(read.csv("sml.csv", row.names=1))[,1]
sexualMaturity <- sexualMaturity[-length(sexualMaturity)]
longevity <- as.matrix(read.csv("sml.csv", row.names=1))[,2]
longevity <- longevity[-length(longevity)]
ratio <- as.matrix(read.csv("sml.csv", row.names=1))[,3]
ratio <- ratio[-length(ratio)]
phylo.heatmap(tree, X, standardize=TRUE)