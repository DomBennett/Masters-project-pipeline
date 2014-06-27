## MRes Project 2013
## Stage 2: Calculate PD and MPD by site + run nulls
## In: 0_data | Out: 1_results, 2_figures
## 25/08/2013

## Print stage
cat("\n\nThis is stage 2: pmetrics\n")

## Libraries
require(reshape)
source(file.path('functions','EcoDataTools.R'))

## Dirs
input.dir <- "0_data"
input.file <- "urban_data.csv"
output.dir <- "1_results"
figures.dir <- "2_figures"
if(!file.exists(figures.dir))
{
  dir.create(figures.dir)
}

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
cat("\nDone.")

## Phylogenies
cat('\nImporting phylogenies ...')
phylo.dir <- file.path(input.dir, "mass_phylos")
phylo.files <- list.files(path = phylo.dir, pattern = "^.*\\.tre$")
phylo.studies <- unique(sub("_phylo\\.tre$", "", phylo.files))
phylo.list <- list()
for (phylo.file in phylo.files) {
  temp.phylos <- read.tree(file.path(phylo.dir,phylo.file))
  temp.phylos <- list(temp.phylos)
  phylo.list <- c(phylo.list, temp.phylos)
}
cat(paste0('\nDone. Found [', length(phylo.studies), '] studies with phylogenies.'))

## Generate phylodata (a list consisting of both phylogeny and community data)
cat("\nGenerating phylodata ...")
meta.phylodata <- list()
props.studies <- props.taxa <- props.meas <- mean.ntaxa <- mean.meas <- total.taxa <-
  vector()
for (i in 1:length(phylo.studies)) {
  cat(paste0("\n\nWorking on [", phylo.studies[i], "] ..."))
  study.phylodata <- list()
  study.data <- data[data.studies %in% phylo.studies[i], ]
  total.meas <- sum(study.data$Measurement)
  unit <- study.data$Diversity_metric_unit[1]
  unames <- sort(unique(study.data$Taxon))
  temp.ntaxa <- temp.taxa <- temp.meas <- vector()
  for (j in 1:length(phylo.list[[i]])){
    temp.phylo <- phylo.list[[i]][[j]]
    temp.data <- study.data[study.data$Taxon %in% temp.phylo$tip.label, ]
    molten.data <- melt.data.frame(temp.data,  measure.vars = 'Measurement',
                                   na.rm = TRUE)
    pull <- names(molten.data) %in% c('Site_number','Use_intensity', 'Predominant_habitat',
                                      'Taxon', 'variable', 'value')
    molten.data <- molten.data[ ,pull]
    temp.data <- cast(molten.data, Site_number + Use_intensity + Predominant_habitat ~
                      Taxon, mean)
    temp.data[is.na(temp.data)] <- 0 #convert nas to 0
    temp.ntaxa <- c(temp.ntaxa, ncol(temp.data) - 3)
    temp.taxa <- c(temp.taxa, colnames(temp.data[,-c(1:3)]))
    temp.meas <- c(temp.meas, sum(temp.data[,-c(1:3)]))
    temp.phylodata <- list(temp.data, temp.phylo)
    study.phylodata <- c(study.phylodata, list(temp.phylodata))
  }
  study.phylodata <- list(phylo.studies[i], unique(study.data$Site_number), study.phylodata)
  meta.phylodata <- c(meta.phylodata, list(study.phylodata))
  prop.taxa <- signif(mean(temp.ntaxa)*100/length(unames), 3)
  prop.meas <- signif(mean(temp.meas)*100/total.meas, 3)
  cat(paste0('\n... found [', prop.taxa, '%] of taxa across phylogenies ....'))
  cat(paste0('\n... representing [', prop.meas, '%] of [', unit, '] ...' ))
  props.studies <- c(props.studies, phylo.studies[i])
  total.taxa <- c(total.taxa, length(unique(temp.taxa)))
  mean.ntaxa <- c(mean.ntaxa, mean(temp.ntaxa))
  mean.meas <- c(mean.meas, mean(temp.meas))
  props.taxa <- c(props.taxa, prop.taxa)
  props.meas <- c(props.meas, prop.meas)
}
props <- data.frame(study = props.studies, mean.ntaxa, mean.meas, p.taxa = props.taxa,
                    p.meas = props.meas, t.taxa = total.taxa)
cat("\nDone.")

## Generate some example community plots
par(mar = c(1, 1, 2, 1) + 0.1)
pdf(file = file.path(figures.dir, "comm_plots.pdf"), width = 10.5, height = 7)
cat("Generating example community plots ...")
for (i in 1:length(meta.phylodata)) {
  study <- meta.phylodata[[i]][[1]]
  cat(paste0("\n\nWorking on [", study, "] ..."))
  comm.data <- meta.phylodata[[i]][[3]][[1]][[1]][,-c(1:3)]
  use.intensities <- meta.phylodata[[i]][[3]][[1]][[1]][,2]
  phylo <- meta.phylodata[[i]][[3]][[1]][[2]]
  # re-order use.intensities to guarantee order is the same
  use.levels <- c("Cannot decide", "Minimal use", "Light use", "Intense use")
  for (j in 1:length(use.levels)) {
    pull <- use.intensities %in% use.levels[j]
    use.intensities[pull] <- rep(j, sum(pull)) # i.e. so that they're numbered instead
  }
  ord <- order(use.intensities)
  plotComm(comm.data[ord,], phylo, groups = use.intensities[ord], no.margin = FALSE)
  uses <- paste0(use.levels[unique(sort(as.numeric(use.intensities)))], collapse = ", ")
  mtext(paste0(study, " | ", uses))
  pa <- (comm.data > 0) * 1
  plotComm(pa[ord,], phylo, groups = use.intensities[ord], no.margin = FALSE)
  uses <- paste0(use.levels[unique(sort(as.numeric(use.intensities)))], collapse = ", ")
  mtext(paste0(study, " | ", uses))
}
dev.off()
par(mar = c(5, 4, 4, 2) + 0.1)

## Calculate PD and MPD by site
cat("Calculating PD and MPD ...")
studies <- sites <- PD.mean <- PD.stderr <- MPD.mean <- MPD.stderr <- vector()
for(i in 1:length(meta.phylodata)) {
  study <- meta.phylodata[[i]][[1]]
  cat(paste0("\n\nWorking on [", study, "] ..."))
  study.sites <- meta.phylodata[[i]][[2]]
  study.phylodata <- meta.phylodata[[i]][[3]]
  pd.study <- mpd.study <- matrix(ncol = length(study.phylodata), nrow = length(study.sites))
  rownames(pd.study) <- rownames(mpd.study) <- study.sites
  for (j in 1:length(study.phylodata)) {
    temp.data <- study.phylodata[[j]][[1]][,-c(1:3)]
    temp.sites <- as.character(study.phylodata[[j]][[1]][,c(1)])
    temp.phylo <- study.phylodata[[j]][[2]]
    pd <- commPD(temp.phylo, temp.data, type = 2, min.spp = 0)[,1]
    mpd <- mpd(temp.data, cophenetic(temp.phylo), abundance.weighted = TRUE)
    pd.study[match(temp.sites, rownames(pd.study)),j] <- pd
    mpd.study[match(temp.sites, rownames(pd.study)),j] <- mpd
  }
  studies <- c(studies, rep(study, length(study.sites)))
  sites <- c(sites, study.sites)
  PD.mean <- c(PD.mean, apply(pd.study, 1, mean, na.rm = TRUE))
  PD.stderr <- c(PD.stderr, apply(pd.study, 1, function(x) sd(x)/sqrt(length(x))))
  MPD.mean <- c(MPD.mean, apply(mpd.study, 1, mean, na.rm = TRUE))
  MPD.stderr <- c(MPD.stderr, apply(mpd.study, 1, function(x) sd(x)/sqrt(length(x))))
}
cat("\nDone.")

cat("\nOutputting results ...")
pdmpd <- data.frame(studies, sites, PD.mean, PD.stderr, MPD.mean, MPD.stderr)
write.csv(x = pdmpd, file = file.path(output.dir, "pdmpd.csv"))
write.csv(x = props, file = file.path(output.dir, "props.csv"))
cat("\nDone.")
cat("\n\nStage finished.\n")