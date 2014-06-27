## MRes Project 2013
## Stage 3: Run nulls
## In: 0_data | Out: 1_results
## 25/08/2013

## Print stage
cat("\n\nThis is stage 3: nulls\n")

## Parameters
nrands <- 1000
comps <- list(c("4", "3"), c("4", "2"), c("3", "2"))
# where 4 is min, 3 is light and 2 is intense

## Libraries
require(reshape)
source(file.path('functions','EcoDataTools.R'))

## Dirs
input.dir <- "0_data"
input.file <- "urban_data.csv"
output.dir <- "1_results"

cat("\nImporting name modifiers ...")
load(file = file.path(input.dir, "name_modifiers.RData"))
unresolved.names <- name.modifiers[[1]]
resolved.names <- name.modifiers[[2]]
cat("\nDone")

cat('\nImporting PREDICTS data ...')
data <- read.csv(file.path(input.dir, input.file), stringsAsFactors = FALSE)
data.studies <-paste0(data$Source_ID, data$Study_number)
for (i in 1:length(unresolved.names)) {
  pull <- data$Taxon %in% unresolved.names[i]
  data$Taxon[pull] <- rep(resolved.names[i], sum(pull))
}
# account for added underscore in phylos
data$Taxon <- sub(" ", "_", data$Taxon)
# specify order -- so that i know whether there's an INCREASE or DECREASE!!
use.intensities <- c("Cannot decide", "Intense use", "Light use", "Minimal use")
for (i in 1:length(use.intensities)) {
  pull <- data$Use_intensity %in% use.intensities[i]
  data$Use_intensity[pull] <- rep(i, sum(pull)) # i.e. so that they numbered instead
}
print("\nDone.")

## Phylogenies
cat('\nImporting phylogenies ...')
phylo.dir <- file.path(input.dir, "mass_phylos")
phylo.files <- list.files(path = phylo.dir, pattern = "^.*\\.tre$")
phylo.studies <- unique(sub("_phylo\\.tre$", "", phylo.files))
phylo.list <- list()
for (phylo.file in phylo.files) {
  temp.phylos <- read.tree(file.path(phylo.dir, phylo.file))
  temp.phylos <- list(temp.phylos)
  phylo.list <- c(phylo.list, temp.phylos)
}
cat(paste0('\nDone. Found [', length(phylo.studies), '] studies with phylogenies.'))

## Generate phylodata (a list consisting of both phylogeny and community data)
cat("\nGenerating phylodata ...")
meta.phylodata <- list()
for (i in 1:length(phylo.studies)) {
  cat(paste0("\n\nWorking on [", phylo.studies[i], "] ..."))
  study.phylodata <- list()
  study.data <- data[data.studies %in% phylo.studies[i], ]
  total.meas <- sum(study.data$Measurement)
  unit <- study.data$Diversity_metric_unit[1]
  unames <- sort(unique(study.data$Taxon))
  temp.ntaxa <- temp.meas <- 0
  for (j in 1:length(phylo.list[[i]])){
    temp.phylo <- phylo.list[[i]][[j]]
    temp.data <- study.data[study.data$Taxon %in% temp.phylo$tip.label, ]
    temp.data <- temp.data[temp.data$Predominant_habitat == "Urban", ] # limit to urban
    molten.data <- melt.data.frame(temp.data,  measure.vars = 'Measurement',
                                   na.rm = TRUE)
    pull <- names(molten.data) %in% c('Site_number','Use_intensity',
                                      'Taxon', 'variable', 'value')
    molten.data <- molten.data[ ,pull]
    temp.data <- cast(molten.data, Site_number + Use_intensity ~ Taxon, mean)
    temp.data[is.na(temp.data)] <- 0 #convert nas to 0
    temp.ntaxa <- temp.ntaxa + ncol(temp.data) - 2
    temp.meas <- temp.meas + sum(temp.data[,-c(1:3)])
    temp.phylodata <- list(temp.data, temp.phylo)
    study.phylodata <- c(study.phylodata, list(temp.phylodata))
  }
  study.phylodata <- list(phylo.studies[i], unique(study.data$Site_number), study.phylodata)
  meta.phylodata <- c(meta.phylodata, list(study.phylodata))
  prop.taxa <- signif((temp.ntaxa/length(phylo.list[[i]]))*100/length(unames), 3)
  prop.meas <- signif((temp.meas/length(phylo.list[[i]]))*100/total.meas, 3)
  cat(paste0('\n... found [', prop.taxa, '%] of taxa across phylogenies ....'))
  cat(paste0('\n... representing [', prop.meas, '%] of [', unit, '] ...' ))
}
cat("\nDone.")

## Nulls!
cat("\nRunning nulls ...")
null.res <- list()
for(i in 1:length(meta.phylodata)) {
  study <- meta.phylodata[[i]][[1]]
  study.phylodata <- meta.phylodata[[i]][[3]]
  pd.null <- mpd.null <- rep(list(matrix(nrow = 6, ncol = nrands)), 3)
  for (j in 1:3) { # loop through comparisons
    comp <- comps[[j]]
    rownames(pd.null[[j]]) <- rownames(pd.null[[j]]) <- paste0("n", rep(0:5))
    touse <- vector()
    for (k in 1:length(study.phylodata)) {
      # identify phylodata that have the comp of interest
      # this means data where each comparison has 2 or more sites
      # otherwise -- how can you shuffle the rows and columns?
      temp.uses <- study.phylodata[[k]][[1]][,2]
      if (sum(temp.uses %in% comp[1]) > 1) {
        if (sum(temp.uses %in% comp[2]) > 1) {
          touse <- c(touse, k)
        }
      }
    }
    if (is.logical(touse)) {
      next
    }
    for (k in 0:5) { #loop for each null
      rand.phylodata <- sample(touse, nrands, replace = TRUE)
      for (l in 1:length(rand.phylodata)) { # loop for each randomisation
        temp.phylo <- study.phylodata[[rand.phylodata[l]]][[2]]
        temp.data <- study.phylodata[[rand.phylodata[l]]][[1]][,-c(1:3)]
        htypes <- study.phylodata[[rand.phylodata[l]]][[1]][,2]
        # pull for current comparison
        pull <- htypes %in% comp
        temp.data <- temp.data[pull, ]
        htypes <- htypes[pull]
        ord <- order(htypes)
        temp.data <- temp.data[ord, ]
        htypes <- htypes[ord]
        # and then the null....
        pd.null[[j]][k + 1, l] <- genNullDist(temp.phylo, temp.data, htypes, null = k,
                                              nrand = 1, metric = "pd", pd.type = 2,
                                              pd.min = 0)
        mpd.null[[j]][k + 1, l] <- genNullDist(temp.phylo, temp.data, htypes, null = k,
                                               nrand = 1, metric = "mpd")
      }
    }
  }
  study.res <- list(study, list(pd.null = pd.null, mpd.null = mpd.null))
  null.res <- c(null.res, list(study.res))
}
# explanation of null.res structure:
# -- contains a series of lists for each study
# -- each study has a list consisting of two lists: pd and mpd
# -- each of these lists contain three matrices containing results for each comparisons:
#     minimal --> light, minimal --> intense, light --> intense (see comps at beginning)
# -- each matrix contains the null distribution for each null model
cat("Done.")
save(null.res, file = file.path(output.dir, "null_res.RData"))