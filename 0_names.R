## MRes Project 2013
## Stage 0: Extracting taxon names from predicts data
## In: 0_data | Out: 1_names
## 23/07/2013

## Print stage
cat("\n\nThis is stage 0: names")

## Directories
input.dir <- "0_data"
input.file <- "urban_data.csv"

## Input
cat('\nImporting PREDICTS data ...')
data <- read.csv(file.path(input.dir, input.file), stringsAsFactors = FALSE)
studies <-paste0(data$Source_ID, data$Study_number)
ustudies <- unique(studies)
nstudies_all <- length(ustudies)
cat("\nDone. Imported data for [", nstudies_all, "] studies.", sep = "")

## Parse names to optimise taxon name resolution
#print ("Parsing names ...")
#parsed <- sub("\\.", " ", data$Taxon_name_entered) # get rid of dots as spaces
#for (i in 1:length(parsed)) { # remove characters beyond 1 space
#  temp.split <- strsplit(parsed[i], "\\s+")[[1]]
#  if (length(temp.split) > 2) {
#    parsed[i] <- paste(temp.split[1:2], collapse = " ")
#  }
#}
#print ("Done.")
#data$Parsed_name <- parsed
#print(paste0('[', length(unique(parsed)), '] parsed names.'))

## Extract names by study
cat("\nExtracting names by study ...")
output.file.names <- list()
nstudies <- nsites <- nurbansites <- 0
nhabitats <- ntaxa <- norders <- nbiomes <- vector()
kept.studies <- htgs <- vector()
for (i in 1:nstudies_all) {
  study.taxnames <- paste0(ustudies[i], "_taxnames")
  temp.data <- data[studies %in% ustudies[i], ]
  study.taxa <- temp.data$Taxon[temp.data$Resolution_entered %in%
                                  c("Scientific", "Uncertain", "Family", "Genus")]
  study.taxa <- unique(study.taxa)
  keep.study <- c(length(unique(temp.data$Use_intensity)) > 1,
                  length(study.taxa) > 2)
  if (all(keep.study)) {
    assign(x = study.taxnames, value = study.taxa)
    output.file.names <- c(output.file.names, list(c(study.taxnames)))
    nstudies <- nstudies + 1
    sites <- unique(paste0(temp.data$Source_ID, temp.data$Site_number))
    nsites <- nsites + length(sites)
    habitats <- data$Predominant_habitat[match(unique(sites), sites)]
    nurbansites <- nurbansites + sum(habitats == "Urban")
    nhabitats <- c(nhabitats, length(habitats))
    ntaxa <- c(ntaxa, unique(temp.data$Taxon))
    orders <- unique(temp.data$Order)
    norders <- c(norders, orders)
    nbiomes <- c(nbiomes, unique(temp.data$Biome))
    kept.studies <- c(kept.studies, ustudies[i])
    # htg - highest shared taxonomic group
    if (temp.data$Genus[1] != "" & length(unique(temp.data$Genus)) == 1) {
      htg <- unique(temp.data$Genus)
    } else if (temp.data$Family[1] != "" & length(unique(temp.data$Family)) == 1) {
      htg <- unique(temp.data$Family)
    } else if (temp.data$Order[1] != "" & length(unique(temp.data$Order)) == 1) {
      htg <- unique(temp.data$Order)
    } else if (temp.data$Class[1] != "" & length(unique(temp.data$Class)) == 1) {
      htg <- unique(temp.data$Class)
    } else if (temp.data$Phylum[1] != "" & length(unique(temp.data$Phylum)) == 1) {
      htg <- unique(temp.data$Phylum)
    } else {
      htg <- paste(unique(temp.data$Kingdom), collapse = "|")
    }
    htgs <- c(htgs, htg)
  }
}
cat('Done. Loaded: [', nrow(data), '] records containing:', sep = "")
cat(' -- [', nstudies, '] studies', sep = "")
cat(' -- [', nsites, '] sites (of which [', nurbansites,'] are urban)', sep = "")
cat(' -- [', length(unique(ntaxa)), '] taxa', sep = "")
cat(' -- [', length(unique(norders)), '] orders of life', sep = "")
cat(' -- [', length(unique(nbiomes)), '] biomes', sep = "")

## Outputting
cat("\nOutputting...")
for (n in output.file.names) {
  write.table(x = get(n[[1]]), sep = "\n", quote = FALSE, col.names = FALSE,
              file = file.path(output.dir, paste0(n[[1]], '.txt')), row.names = FALSE)
}
write.csv(x = data.frame(study = kept.studies, htg = htgs),
          file.path(input.dir, 'prem_taxadata.csv'), row.names = FALSE)
cat("\nDone.")
cat("\nStage finished.")