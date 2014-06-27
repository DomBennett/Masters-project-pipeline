## MRes Project 2013
## Stage 5: Create figures of example trees
## In: 0_data | Out: 2_figures
## 01/09/2013

##Libraries
require(reshape)
source(file.path('functions','EcoDataTools.R'))

## Dirs
input.dir <- "0_data"
output.dir <- "2_figures"

## Input
all.data <- read.csv(file.path(input.dir, "urban_data.csv"), stringsAsFactors = FALSE)
studies <-paste0(all.data$Source_ID, all.data$Study_number)
ustudies <- unique(studies)
phyla <- all.data$Phylum[match(ustudies, studies)]

## Example trees from distribution
phylo.dir <- file.path(input.dir, "mass_phylos")
phylo.files <- list.files(path = phylo.dir, pattern = "^.*\\.tre$")
phylo.studies <- sub("_phylo\\.tre$", "", phylo.files)
choice <- c(16, 15, 14, 14)
tree.i <- c(1,1,1,15)
phylo.list <- list()
for (i in 1:length(choice)) {
  phylo <- read.tree(file.path(phylo.dir, phylo.files[choice[i]]))
  phylo.list <- c(phylo.list, list(phylo[[tree.i[i]]])) # take first in distribution
}

## Plot some examples
pdf(file = file.path(output.dir, "exampletrees.pdf"), height = 8, width = 10.5)
par(mfrow = (c(2,2)), mar = c(1, 1, 2, 1) + 0.1)
 # choose examples
groups <- c("Mammals", "Birds", "Ants", "Ants")
genes <- c("COI", "COI", "28S", "COI")
for (i in 1:length(phylo.list)) {
  phylo <- phylo.list[[i]]
  plot(phylo, cex = 0.6, no.margin = FALSE, edge.width = 0.3)
  mtext(paste(phylo.studies[choice[i]], "|", groups[i], "|", genes[i]), cex = 0.5)
}
dev.off()