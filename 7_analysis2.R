## MRes Project 2013
## Stage 7: Is there phylogenetic signal in the species gains and losses?
## In: 0_data, 1_results | Out: 1_results, 2_figures
## 01/09/2013

## Parameters
comps <- list(c("4", "3"), c("4", "2"), c("3", "2"))
# where 4 is min, 3 is light and 2 is intense

## Functions
extractFromNullObj <- function(null.obj, comp.obj, pdmpd, metric = "pd", null = 0) {
  if (metric == "pd") {
    metric.null <- 'pd.null'
  } else {
    metric.null <- 'mpd.null'
  }
  null <- null + 1
  output <- data.frame()
  studies <- rep(NA, length(null.obj))
  for (i in 1:length(null.obj)){
    cstudy <- null.obj[[i]][[1]]
    studies[i] <- cstudy
    pd <- pdmpd$PD.mean[pdmpd$studies %in% cstudy]
    uses <- pdmpd$use.intensities[pdmpd$studies %in% cstudy]
    means <- tapply(pd, uses, mean)
    diffs <- ps <- zs <- rep(NA, length(comp.obj))
    for (j in 1:length(comp.obj)) {
      if (sum(names(means) %in% comp.obj[[j]]) > 1) {
        pull1 <- names(means) == comp.obj[[j]][2]
        pull2 <- names(means) == comp.obj[[j]][1]
        difference <- means[pull1] - means[pull2]
        null.dist <- null.obj[[i]][[2]][metric.null][[1]][[j]][null,]
        nrands <- sum(!is.na(null.dist))
        if (nrands == 0) {
          diffs[j] <- NA
          ps[j] <- NA
          zs[j] <- NA
        } else {
          diffs[j] <- difference
          ps[j] <- sum(null.dist <= difference)/nrands
            # what proportion have a SMALLER value because genNullDist calculates the
            # difference as the htype 1 - htype 2. I regulated the order by which
            # genNullDist does by changing the use_intensities to numbers. Because we 
            # expect there to be a reduction in PD as a result of urbanisation, we expect
            # the difference to be negative. So we need to test if anything in the null
            # dist has values less than or equal to it.
          zs[j] <- (difference - mean(null.dist))/sd(null.dist)
            # again, because of the comparisons we expect the difference to be smaller than
            # the nulls.
        }
      } else {
        diffs[j] <- NA
        ps[j] <- NA
        zs[j] <- NA
      }
    }
    output <- rbind(output, c(diffs, ps, zs))
  }
  rownames(output) <- studies
  colnames(output) <- c("ml.diff", "mi.diff", "li.diff", "ml.p", "mi.p", "li.p", "ml.z",
                        "mi.z", "li.z")
  return(output)
}

genNullFingerPrints <- function(transition = 'ml') {
  pvalue <- paste0(transition, '.p')
  obs <- n0 <- n1 <- n2 <- n3 <- n4 <- n5 <- rep(NA, length(rownames(null0.res)))
  for (i in 1:length(rownames(null0.res))){
    obs[i] <- null0.res[paste0(transition, '.diff')][i,]
    temp.res <- c(null0.res[pvalue][i,], null1.res[pvalue][i,], null2.res[pvalue][i,],
                  null3.res[pvalue][i,], null4.res[pvalue][i,], null5.res[pvalue][i,])
    codes <- rep(NA, length(temp.res))
    for (j in 1:length(temp.res)) {
      if (!is.na(temp.res[j])) {
        if (temp.res[j] < 0.025) {
          codes[j] <- -1
        } else if (temp.res[j] > 0.975) {
          codes[j] <- 1
        } else {
          codes[j] <- 0
        }
      }
    }
    n0[i] <- codes[1]
    n1[i] <- codes[2]
    n2[i] <- codes[3]
    n3[i] <- codes[4]
    n4[i] <- codes[5]
    n5[i] <- codes[6]
  }
  return(data.frame(study = rownames(null0.res), obs, n0, n1, n2, n3, n4, n5))
}


## Dirs
data.dir <- "0_data"
results.dir <- "1_results"
figures.dir <- "2_figures"

## 1 Input + data manipulation
load(file.path(results.dir, "null_res.RData"))
props <- read.csv(file.path(results.dir, "props.csv"), stringsAsFactors = FALSE)
all.data <- read.csv(file.path(data.dir, "urban_data.csv"), stringsAsFactors = FALSE)
# specify order -- so that i know whether there's an INCREASE or DECREASE!!
use.intensities <- c("Cannot decide", "Intense use", "Light use", "Minimal use")
for (i in 1:length(use.intensities)) {
  pull <- all.data$Use_intensity %in% use.intensities[i]
  all.data$Use_intensity[pull] <- rep(i, sum(pull)) # i.e. so that they numbered instead
}
pdmpd <- read.csv(file.path(results.dir, "pdmpd.csv"))
#pdmpd <- pdmpd[pdmpd$PD.mean !=0,] # I can't drop 0s, because
# add site IDs
all.data$siteID <- paste(all.data$Site_number, paste0(all.data$Source_ID,
                                                      all.data$Study_number), sep = "_")
pdmpd$siteID <- paste(pdmpd$sites, pdmpd$studies, sep = "_")
# add hab info to pdmpd
pull <- match(pdmpd$siteID, all.data$siteID)
pdmpd$use.intensities <- all.data$Use_intensity[pull]
pdmpd$htype <- all.data$Predominant_habitat[pull]

## 2. Nulls for PD
# 2.1 Generate p and z values
null0.res <- extractFromNullObj(null.res, comps, pdmpd, null = 0)
null1.res <- extractFromNullObj(null.res, comps, pdmpd, null = 1)
null2.res <- extractFromNullObj(null.res, comps, pdmpd, null = 2)
null3.res <- extractFromNullObj(null.res, comps, pdmpd, null = 3)
null4.res <- extractFromNullObj(null.res, comps, pdmpd, null = 4)
null5.res <- extractFromNullObj(null.res, comps, pdmpd, null = 5)
# 2.2 Generate fignerprints
mtol <- na.omit(genNullFingerPrints('ml'))
mtoi <- na.omit(genNullFingerPrints('mi'))
ltoi <- na.omit(genNullFingerPrints('li'))
# 2.3 write out
write.csv(mtol, file.path(results.dir, "mtol_pd.csv"), row.names = FALSE)
write.csv(mtoi, file.path(results.dir, "mtoi_pd.csv"), row.names = FALSE)
write.csv(ltoi, file.path(results.dir, "ltoi_pd.csv"), row.names = FALSE)

## 2. Nulls for MPD
# 2.1 Generate p and z values
null0.res <- extractFromNullObj(null.res, comps, pdmpd, null = 0, metric = "mpd")
null1.res <- extractFromNullObj(null.res, comps, pdmpd, null = 1, metric = "mpd")
null2.res <- extractFromNullObj(null.res, comps, pdmpd, null = 2, metric = "mpd")
null3.res <- extractFromNullObj(null.res, comps, pdmpd, null = 3, metric = "mpd")
null4.res <- extractFromNullObj(null.res, comps, pdmpd, null = 4, metric = "mpd")
null5.res <- extractFromNullObj(null.res, comps, pdmpd, null = 5, metric = "mpd")
# 2.2 Generate fignerprints
mtol <- na.omit(genNullFingerPrints('ml'))
mtoi <- na.omit(genNullFingerPrints('mi'))
ltoi <- na.omit(genNullFingerPrints('li'))
# 2.3 write out
write.csv(mtol, file.path(results.dir, "mtol_mpd.csv"), row.names = FALSE)
write.csv(mtoi, file.path(results.dir, "mtoi_mpd.csv"), row.names = FALSE)
write.csv(ltoi, file.path(results.dir, "ltoi_mpd.csv"), row.names = FALSE)
#ignore helden3 for mpd