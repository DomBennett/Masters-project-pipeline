## MRes Project 2013
## Stage 1: Compare mass phylos to published phylos
## In: 0_data | Out: 1_results
## 28/08/2013

## Print stage
cat("\n\nThis is stage 1: compare\n")

## Parameters
nrands <- 1000 # number of randomisations

## Functions
fuzzyRFDist <- function(tree1, tree2, max.dist.genus = 2) {
  # Take tree1 and tree2, shrink to same size and calculate RF dist. Uses fuzzy-matching
  #  to match tip labels. If no matches -- returns NA.
  run.dist = FALSE
  matches <- tree2$tip.label %in% tree1$tip.label
  if (sum(matches) < 4) {
    # if too few matches, reduce to genus
    new.1labels <- sub("_.*$", "", tree1$tip.label)
    new.2labels <- sub("_.*$", "", tree2$tip.label)
    matches <- amatch(new.1labels, new.2labels, method = "lv", maxDist = max.dist.genus)
    matches <- unique(matches[!is.na(matches)])
    nmatches <- length(matches)
    if (nmatches > 4) {
      to.drop <- tree2$tip.label[-matches]
      test.tree2 <- drop.tip(phy = tree2, tip = to.drop)
      test.tree2$tip.label <- new.2labels[matches]
      # make names match
      pull <- amatch(test.tree2$tip.label, new.1labels, method = "lv",
                     maxDist = max.dist.genus)
      new.1labels[pull] <- test.tree2$tip.label
      # remove species of same genus
      test.tree1 <- drop.tip(tree1, tree1$tip.label[duplicated(new.1labels)])
      test.tree1$tip.label <- new.1labels[!duplicated(new.1labels)]
      to.drop <- test.tree1$tip.label[!test.tree1$tip.label %in% test.tree2$tip.label]
      test.tree1 <- drop.tip(phy = test.tree1, tip = to.drop)
      if (length(test.tree1$tip.label) == length(test.tree2$tip.label)) {
        if (test.tree1$Nnode > 1 & test.tree2$Nnode > 1) {
          run.dist = TRUE
        }
      }
    }
  } else {
    to.drop <- tree2$tip.label[!matches]
    test.tree2 <- drop.tip(phy = tree2, tip = to.drop)
    matches <- tree1$tip.label %in% test.tree2$tip.label
    to.drop <- tree1$tip.label[!matches]
    test.tree1 <- drop.tip(phy = tree1, tip = to.drop)
    if (length(test.tree1$tip.label) == length(test.tree2$tip.label)) {
      if (test.tree1$Nnode > 1 & test.tree2$Nnode > 1) {
        run.dist = TRUE
      }
    }
  }
  if (run.dist) {
    test.tree1$edge.length <- test.tree1$edge.length/sum(test.tree1$edge.length)
    test.tree2$edge.length <- test.tree2$edge.length/sum(test.tree2$edge.length)
    dist <- dist.topo(test.tree1, test.tree2, "score")
    return(dist)
  }
  else {
    return(NA)
  }
}

## Libraries
source(file.path('functions','EcoDataTools.R'))

## Dirs
input.dir <- "0_data"
output.dir <- "1_results"
if(!file.exists(output.dir))
{
  dir.create(output.dir)
}

## Treecomp.csv -- for identifying which pub phylo for which mass phylo
cat('\nReading in treecomp.csv ...')
treecomp <- read.csv(file.path(input.dir, 'pub_phylos', 'treecomp.csv'),
                     stringsAsFactors = FALSE)
treecomp <- treecomp[!is.na(treecomp$pubphylo), ]
cat(paste0('\nDone. Read in data for [', nrow(treecomp), '] comparisons.'))

## Phylogenies
cat('\nIdentifying studies with phylogenies ...')
phylo.dir <- file.path(input.dir, "mass_phylos")
phylo.files <- list.files(path = phylo.dir, pattern = "^.*\\.tre$")
phylo.studies <- unique(sub("_phylo\\.tre$", "", phylo.files))
cat(paste0('\nDone. Found [', length(phylo.studies), '] studies with phylogenies.'))

## Create phylophylo object -- list of pub phylo and mass phylos
cat("\nMatching published phylogeny to study phylogenies ...")
phylophylo <- list()
for (i in 1:length(phylo.studies)) {
  cat(paste0("\nWorking on [", phylo.studies[i], "] ..."))
  pubphylo.file <- treecomp$pubphylo[treecomp$study == phylo.studies[i]]
  if (length(pubphylo.file) > 0) {
    pubphylo <- read.tree(file.path(input.dir, 'pub_phylos', pubphylo.file))
    massphylos <- read.tree(file.path(phylo.dir, phylo.files[i]))
    std.obj <- list(phylo.studies[i], pubphylo.file, pubphylo, massphylos)
    phylophylo <- c(phylophylo, list(std.obj))
  } else {
    cat("\nNo published phylogenies specified.")
  }
}
cat("\nDone. Data for [", length(phylophylo), "] comparisons.", sep = "")

## Compare!
cat("\n\nComparing with nrands = [", nrands, "] ...", sep = "")
# out vectors
study <- pub <- ntips.pub <- ntips.pub.sd <- ntips.mass <- ntips.mass.sd <-
  dist.mean <- dist.sd <- p.value <- vector()
counter <- 0
for (i in 1:length(phylophylo)) {
  # study info
  cstudy <- phylophylo[[i]][[1]]
  cat("\nWorking on [", cstudy, "] ...", sep = "")
  pubphylo.name <- phylophylo[[i]][[2]]
  pubphylo <- phylophylo[[i]][[3]]
  massphylos <- phylophylo[[i]][[4]]
  # calc dists + ntips etc...
  mntips <- pntips <- dists <- rep(NA, length(massphylos))
  for (j in 1:length(massphylos)) {
    massphylo <- massphylos[[j]]
    mntips[j] <- length(massphylo$tip.label)
    pntips[j] <- sum(pubphylo$tip.label %in% massphylo$tip.label)
    dists[j] <- fuzzyRFDist(massphylo, pubphylo)
  }
  # randomisation test -- scrap randomisation test. it doesnt work for all studies
  # and i'm ot convinced it makes them comprable between studies.
  #null.dist <- rep(NA, nrands)
  #for (j in 1:nrands) {
  #  rn <- sample(1:length(massphylos), 1)
  #  massphylo <- massphylos[[rn]]
  #  null.dist[j] <- fuzzyRFDist(massphylo, pubphylo)
  #}
  # writing out
  study <- c(study, cstudy)
  pub <- c(pub, pubphylo.name)
  ntips.pub <- c(ntips.pub, mean(pntips))
  ntips.pub.sd <- c(ntips.pub.sd, sd(pntips))
  ntips.mass <- c(ntips.mass, mean(mntips))
  ntips.mass.sd <- c(ntips.mass.sd, sd(mntips))
  dist.mean <- c(dist.mean, mean(dists, na.rm = TRUE))
  dist.sd <- c(dist.sd, sd(dists, na.rm = TRUE))
  #p.value <- c(p.value, sum(null.dist <= mean(dists, na.rm = TRUE))/
  #               sum(!is.na(null.dist)))
  counter <- counter + 1
}
res <- data.frame(study, pub, ntips.pub, ntips.pub.sd, ntips.mass, ntips.mass.sd, 
                  dist.mean, dist.sd)#, p.value)
cat("\nDone. Compared [", counter, "] mass phylogeny distributions to published phylogenies",
    sep = "")
# Note: dist generates a warning "trees are not binary". This is ok as it simply,
#  affects the upper bound (i.e. the maximum RF distance given, completely different trees)
#  see -- https://stat.ethz.ch/pipermail/r-sig-phylo/2011-December/001794.html

## Output
cat("\n\nOutputting ...")
write.csv(x = res, file = file.path(output.dir, "treecomp_res.csv"), row.names = FALSE)
cat("\nDone.")
cat("\n\nStage finished.\n")