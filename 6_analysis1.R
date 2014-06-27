## MRes Project 2013
## Stage 6: How does phylogenetic diversity change along the urban gradient?
## In: 0_data, 1_results | Out: 1_results, 2_figures
## 21/08/2013

## Libraries
source(file.path('functions','EcoDataTools.R'))
source(file.path('functions','plot.fixed.effects.R'))
require(vioplot)

## Functions
# http://lamages.blogspot.co.uk/2013/04/how-to-change-alpha-value-of-colours-in.html
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

## Dirs
data.dir <- "0_data"
results.dir <- "1_results"
figures.dir <- "2_figures"
if(!file.exists(figures.dir))
{
  dir.create(figures.dir)
}

## 1 Input + data manipulation
props <- read.csv(file.path(results.dir, "props.csv"), stringsAsFactors = FALSE)
all.data <- read.csv(file.path(data.dir, "urban_data.csv"), stringsAsFactors = FALSE)
pdmpd <- read.csv(file.path(results.dir, "pdmpd.csv"))
# add site IDs
all.data$siteID <- paste(all.data$Site_number, paste0(all.data$Source_ID,
                                                      all.data$Study_number), sep = "_")
pdmpd$siteID <- paste(pdmpd$sites, pdmpd$studies, sep = "_")
# add hab info to pdmpd
pull <- match(pdmpd$siteID, all.data$siteID)
pdmpd$use.intensities <- all.data$Use_intensity[pull]
pdmpd$htype <- all.data$Predominant_habitat[pull]
pdmpd$phylum <- all.data$Phylum[pull]
pdmpd$biome <- all.data$Biome[pull]

## 2. Data stats
mean(props$mean.ntaxa) # a mean mean of 13 taxa
sum(props$t.taxa) # 949
length(unique(pdmpd$studies)) # 28 studies
nrow(pdmpd) # 649 sites
length(unique(pdmpd$phylum)) # 6 phyla
length(unique(pdmpd$biome)) # 9 biomes
mean(props$p.taxa) # 36% of taxa represented across studies
sd(props$p.taxa) # 24%
mean(props$p.meas) # 35% of measurements
sd(props$p.meas) # 25%
site.sums.all.data <- tapply(all.data$Measurement, all.data$Site_name, sum)
sum(site.sums.all.data == 0)/length(site.sums.all.data) # only 10%!
sum(pdmpd$PD.mean == 0)/nrow(pdmpd) # 45%!
# drop all null observations
pdmpd <- pdmpd[pdmpd$PD.mean != 0,]

## 3. Create urbanisation level
urbanisation <- rep(NA, nrow(pdmpd))
for (i in 1:nrow(pdmpd)) {
  if (pdmpd$htype[i] == "Urban"){
    if (pdmpd$use.intensities[i] == "Minimal use") {
      urbanisation[i] <- 1
    } else if (pdmpd$use.intensities[i] == "Light use") {
      urbanisation[i] <- 2
    } else if (pdmpd$use.intensities[i] == "Intense use") {
      urbanisation[i] <- 3
    }
  } else {
    urbanisation[i] <- 0
  }
}
pdmpd$urbanisation <- urbanisation
remove(urbanisation)

## 4 PD Analysis

# 4.1. basic stats
mean(pdmpd$PD.mean) # mean PD 3.9
sd(pdmpd$PD.mean) # sd PD 5.29
var(pdmpd$PD.mean) # 28

# 4.2. normalisation
hist(pdmpd$PD.mean) # not normal
hist(log(pdmpd$PD.mean)) # more normal.... eeeee?
pdmpd.log <- pdmpd
pdmpd.log$PD.mean <- log(pdmpd.log$PD.mean) # log
mean(pdmpd.log$PD.mean) # mean PD 0.33
sd(pdmpd.log$PD.mean) # sd PD 1.64
var(pdmpd.log$PD.mean) # 2.69

# 4.3. explore by different factors:
tapply(pdmpd.log$PD.mean, as.factor(pdmpd.log$urbanisation != 0), mean)
tapply(pdmpd.log$PD.mean, as.factor(pdmpd.log$urbanisation != 0), sd)
  # much high PD for non-urban sites
# let's remove non-urban
pdmpd.urb <- pdmpd.log[pdmpd.log$urbanisation != 0,]
# now explore...
plot(pdmpd.urb$PD.mean ~ as.factor(pdmpd.urb$urbanisation))
tapply(pdmpd.urb$PD.mean, as.factor(pdmpd.urb$urbanisation), mean) # large diffs
tapply(pdmpd.urb$PD.mean, as.factor(pdmpd.urb$urbanisation), sd)
table(pdmpd.urb$urbanisation)
plot(pdmpd.urb$PD.mean ~ as.factor(pdmpd.urb$biome))
tapply(pdmpd.urb$PD.mean, as.factor(pdmpd.urb$biome), mean) # mangroves?
tapply(pdmpd.urb$PD.mean, as.factor(pdmpd.urb$biome), sd)
plot(pdmpd.urb$PD.mean ~ as.factor(pdmpd.urb$phylum))
tapply(pdmpd.urb$PD.mean, as.factor(pdmpd.urb$phylum), mean)
tapply(pdmpd.urb$PD.mean, as.factor(pdmpd.urb$phylum), sd)
plot(pdmpd.urb$PD.mean ~ as.factor(pdmpd.urb$urbanisation))
tapply(pdmpd.urb$PD.mean, as.factor(pdmpd.urb$urbanisation), mean)
tapply(pdmpd.urb$PD.mean, as.factor(pdmpd.urb$urbanisation), sd)
# look at violin plot
rural.pd <- pdmpd.urb$PD.mean[pdmpd.urb$urbanisation == 1]
subur.pd <- pdmpd.urb$PD.mean[pdmpd.urb$urbanisation == 2]
urban.pd <- pdmpd.urb$PD.mean[pdmpd.urb$urbanisation == 3]
graphics.off()
pdf(file.path(figures.dir, "pd_vioplot.pdf"), width = 10.5, height = 8)
vioplot(rural.pd, subur.pd, urban.pd, names = c("Rural", "Suburban", "Urban"),
        col = "cornflowerblue")
title(ylab = "Log mean PD by site", xlab = "Level of Urbanisation")
dev.off()
# if anything there is a slight increase with urbanisation, ignoring non-urban habs

# 4.4. linear mixed effects model
# simple model:
m1 <- lmer(as.factor(PD.mean)~ 1 + (1|studies), data = pdmpd.urb, REML = FALSE)
summary(m1)
qqnorm(resid(m1))
qqline(resid(m1))
plot(resid(m1))
plot(resid(m1) ~ fitted(m1))
abline(h=0)
plot(pdmpd.urb$PD.mean ~ fitted(m1))
abline(a=0,b=1)
# as a function of use intensity:
m2 <- lmer(as.factor(PD.mean) ~ as.factor(urbanisation) +
              (1|studies), REML = FALSE, data = pdmpd.urb)
summary(m2)
qqnorm(resid(m2))
qqline(resid(m2))
plot(resid(m2))
plot(resid(m2) ~ fitted(m2))
abline(h=0)
plot(pdmpd.urb$PD.mean ~ fitted(m2))
abline(a=0,b=1)
dotplot(ranef(m2))
anova(m1,m2)
plot.fixed.effects(m2, main = "Urbanisation")
# yes there is a decrease -- but it is not significant

# 4.5. what about for a subset of the data?
# pull data with high numbers of taxa
pull <- pdmpd.urb$studies %in% props$study[props$mean.ntaxa > 20]
pdmpd.urb.bigstds <- pdmpd.urb[pull, ]
plot(pdmpd.urb.bigstds$PD.mean ~
       as.factor(pdmpd.urb.bigstds$urbanisation))
# simple model:
m1 <- lmer(as.factor(PD.mean)~ 1 + (1|studies), data = pdmpd.urb.bigstds,
           REML = FALSE)
summary(m1)
qqnorm(resid(m1))
qqline(resid(m1))
plot(resid(m1))
plot(resid(m1) ~ fitted(m1))
abline(h=0)
plot(pdmpd.urb.bigstds$PD.mean ~ fitted(m1))
abline(a=0,b=1)
# as a function of use intensity:
m2 <- lmer(as.factor(PD.mean) ~ as.factor(urbanisation) +
             (1|studies), REML = FALSE, data = pdmpd.urb.bigstds)
summary(m2)
qqnorm(resid(m2))
qqline(resid(m2))
plot(resid(m2))
plot(resid(m2) ~ fitted(m2))
abline(h=0)
plot(pdmpd.urb.bigstds$PD.mean ~ fitted(m2))
abline(a=0,b=1)
dotplot(ranef(m2))
anova(m1,m2)
plot.fixed.effects(m2, main = "Urbanisation")
# shows the same result... pd doesn't change with urbanisation.

## 5. MPD Analysis

# 5.1. basic stats
mpd <- na.omit(pdmpd) # remove NAs
mean(mpd$MPD.mean) # mean MPD 0.89
sd(mpd$MPD.mean) # sd MPD 1.39
var(mpd$MPD.mean) # 1.938

# 5.2. normalisation
hist(mpd$MPD.mean) # not normal
hist(log(mpd$MPD.mean)) # logged much better
mpd.log <- mpd
mpd.log$MPD.mean <- log(mpd$MPD.mean) # log
mean(mpd.log$MPD.mean) # mean MPD 0.89
sd(mpd.log$MPD.mean) # sd MPD 1.39
var(mpd.log$MPD.mean) # 1.938

# 5.3. explore by different factors:
tapply(mpd.log$PD.mean, as.factor(mpd.log$urbanisation != 0), mean)
tapply(mpd.log$PD.mean, as.factor(mpd.log$urbanisation != 0), sd)
# much high PD for non-urban sites
# remove non urban again and model
mpd.urb <- mpd.log[mpd.log$urbanisation != 0,]
plot(mpd.urb$MPD.mean ~ as.factor(mpd.urb$urbanisation))
tapply(mpd.urb$MPD.mean, as.factor(mpd.urb$urbanisation), mean)
tapply(mpd.urb$MPD.mean, as.factor(mpd.urb$urbanisation), sd)
plot(mpd.urb$MPD.mean ~ as.factor(mpd.urb$biome))
plot(mpd.urb$MPD.mean ~ as.factor(mpd.urb$phylum))
table(mpd.urb$phylum)
tapply(mpd.urb$MPD.mean, as.factor(mpd.urb$phylum), mean)
tapply(mpd.urb$MPD.mean, as.factor(mpd.urb$phylum), sd)
plot(mpd.urb$MPD.mean ~ as.factor(mpd.urb$urbanisation))

# look at violin plot
rural.mpd <- mpd.urb$MPD.mean[mpd.urb$urbanisation == 1]
subur.mpd <- mpd.urb$MPD.mean[mpd.urb$urbanisation == 2]
urban.mpd <- mpd.urb$MPD.mean[mpd.urb$urbanisation == 3]
graphics.off()
pdf(file.path(figures.dir, "mpd_vioplot.pdf"), width = 10.5, height = 8)
vioplot(rural.mpd, subur.mpd, urban.mpd, names = c("Rural", "Suburban", "Urban"),
        col = "cornflowerblue")
title(xlab = "Level of urbanisation", ylab = "Log mean MPD by site")
dev.off()
# hmmm... looks like a decrease with an inversed suburban hump

# 5.4 linear mixed effects model
# simple model:
m1 <- lmer(as.factor(MPD.mean)~ 1 + (1|studies), data = mpd.urb, REML = FALSE)
summary(m1)
qqnorm(resid(m1))
qqline(resid(m1))
plot(resid(m1))
plot(resid(m1) ~ fitted(m1))
abline(h=0)
plot(mpd.urb$MPD.mean ~ fitted(m1))
abline(a=0,b=1)
# as a function of use intensity:
m2 <- lmer(as.factor(MPD.mean) ~ as.factor(urbanisation) +
             (1|studies), REML = FALSE, data = mpd.urb)
summary(m2)
qqnorm(resid(m2))
qqline(resid(m2))
plot(resid(m2))
plot(resid(m2) ~ fitted(m2))
abline(h=0)
plot(mpd.urb$MPD.mean ~ fitted(m2))
abline(a=0,b=1)
dotplot(ranef(m2))
anova(m1,m2)
plot.fixed.effects(m2, main = "Urbanisation")
# again yes there is a decrease -- but it is not significant

## 6. What do MPD and PD look like together?
correlation <- signif(cor(mpd.urb$MPD.mean,log(mpd.urb$PD.mean)), 3) # high correlation 0.82
plot(mpd.urb$MPD.mean ~ log(mpd.urb$PD.mean))
colour <- add.alpha("cornflowerblue", 0.7)
line.model <- lm(mpd.urb$MPD.mean ~ log(mpd.urb$PD.mean))
graphics.off()
pdf(file.path(figures.dir, "pd_mpd_corr.pdf"), width = 10.5, height = 8)
plot(mpd.urb$MPD.mean ~ log(mpd.urb$PD.mean), col = colour, pch = 19, xlab = "Log PD",
     ylab = "Log MPD")
abline(line.model, lty = 3)
text(-1, 1, paste0("Pearson's R:    ", correlation))
text(-0.75, 0.5, paste0("m:    ", signif(line.model$coefficient[2], 2)))
dev.off()