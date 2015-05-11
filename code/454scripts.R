#!/usr/bin/env R

# Installation of ArrayExpress package:
#
# source("http://bioconductor.org/biocLite.R")
# biocLite("ArrayExpress")

library("ArrayExpress")
library("waveform")

# Downloading the dataset separately from analysis is sensible
# Navigate to an appropriate folder with setwd() or select with the RStudio GUI
# From this working directory, download with getAE from the ArrayExpress package

# getAE("E-MTAB-454")

# It is then possible to run analysis without waiting to download a dataset each time
# Note the setting of getAE(local=TRUE)

emtab454 = getAE(accession = 'E-MTAB-454', path = getwd(), type='raw', extract=TRUE, local=TRUE, sourcedir = getwd())
emtab454raw = ae2bioc(mageFiles = emtab454)

# Attempt to build ExpressionSet object fails with anonymous error handling

# cnames <- getcolproc(emtab454)
# emtab454proc = procset(emtab454, cnames[2])

# ==> Error in eval(expr, envir, enclos) : Unable to read processed data

# Developer suggests instead use processed data as only difference is RMA and log(2),
# which I would apply before clustering anyway. See protocol notes for E-MTAB-454

fmA <- read.table('FinalMatrixA.txt', sep='\t', header = TRUE)
fmB <- read.table('FinalMatrixB.txt', sep='\t', header = TRUE)
# The first row contains CompositeElementRef: "log2/RMA normalized" for all columns, remove these
fmAmatrix <- fmA[-c(1),]
fmBmatrix <- fmB[-c(1),]

# To take the mean across replicates and thus give a single value per sample, group with dplyr
library(dplyr)
# library(reshape2)
# reshape2 has been succeeded by tidyr, as of Oct 2014 melt performance similar to gather
# http://stackoverflow.com/a/24884418/2668831
library(tidyr)

#fmAmelted <- melt(fmAmatrix, id.vars=c("Hybridization.REF"))
fmAmelted <- gather(fmAmatrix, variable, value, -Hybridization.REF)
#fmBmelted <- melt(fmBmatrix, id.vars=c("Hybridization.REF"))
fmBmelted <- gather(fmBmatrix, variable, value, -Hybridization.REF)

# The values are now 'grouped by' sample name ("Hybridization.REF") and can be averaged
# as soon as their column names (indicating time point and replicate number) are summarised in dplyr

fmAmelted['variable'] = substring(fmAmelted$variable,5)
fmBmelted['variable'] = substring(fmBmelted$variable,5)

fmAgrouped <- group_by(fmAmelted, variable, Hybridization.REF) %>%
    summarise(mean=mean(as.numeric(value)), sd=sd(value))

fmBgrouped <- group_by(fmBmelted, variable, Hybridization.REF) %>%
    summarise(mean=mean(as.numeric(value)), sd=sd(value))

fmAmean <- dcast(fmAgrouped, Hybridization.REF ~ variable, value.var = "mean")

fmAsd <- dcast(fmAgrouped, Hybridization.REF ~ variable, value.var = "sd")

fmBmean <- dcast(fmBgrouped, Hybridization.REF ~ variable, value.var = "mean")

fmBsd <- dcast(fmBgrouped, Hybridization.REF ~ variable, value.var = "sd")

write.table(fmAmean, file = 'fmAmean.tsv', quote = FALSE, sep = "\t", row.names = FALSE)
write.table(fmBmean, file = 'fmBmean.tsv', quote = FALSE, sep = "\t", row.names = FALSE)
write.table(fmAsd, file = 'fmAsd.tsv', quote = FALSE, sep = "\t", row.names = FALSE)
write.table(fmBsd, file = 'fmBsd.tsv', quote = FALSE, sep = "\t", row.names = FALSE)


# Mapping Affymetrix IDs to Entrez gene IDs

probesA <- as.character(fmAmean$
U133A.probes <- as.character(fmAmean$Hybridization.REF)
U133B.probes <- as.character(fmBmean$Hybridization.REF)

# Install the chip database from BioConductor:

# source("http://bioconductor.org/biocLite.R")
# biocLite("hgu133a.db")

library(hgu133a.db)

U133A.mapped <- select(hgu133a.db, U133A.probes, c("SYMBOL", "ENTREZID", "GENENAME"))
U133B.mapped <- select(hgu133a.db, U133B.probes, c("SYMBOL", "ENTREZID", "GENENAME"))

# Merge later

fmAmean.mat <- fmAmean[,2:13]
rownames(fmAmean.mat) <- fmAmean[,1]

fmBmean.mat <- fmBmean[,2:13]
rownames(fmBmean.mat) <- fmBmean[,1]

fmAmean.mat <- as.matrix(fmAmean.mat)
fmBmean.mat <- as.matrix(fmBmean.mat)

# NB: default cutoff method where unspecified is p value, and the cutoff is p < 0.05

dsetA = waveform(fmAmean.mat)
dsetB = waveform(fmBmean.mat)

# These objects were saved as .Rdata for ease of reproducibility, but may be generated as described above:

# save(dsetA, dsetB, file="dsetA_and_dsetB.Rdata")

# Subset above R squared of 0.6 as recommended in Conesa 2006 for gene selection in microarray time series expression data

rsq = 0.6

good.fit.A <- dsetA$waveform$fit$r.squared > rsq
good.fit.B <- dsetB$waveform$fit$r.squared > rsq

# Visualisation of heat maps

library(RColorBrewer)
library(scales)

# Order genes by phase angle for visual representation of periodicity

angle.order.A <- order(dsetA$angle[good.fit.A,4])
angle.order.B <- order(dsetB$angle[good.fit.B,4])

# NB: "raw" as in not yet [frequency] statistically modelled (these are normalised transcript levels, not raw)

raw.data.A = dsetA$raw[good.fit.A,]
model.A = dsetA$waveform$filtered.data[good.fit.A,]
raw.data.B = dsetB$raw[good.fit.B,]
model.B = dsetB$waveform$filtered.data[good.fit.B,]

# showWave <- function() {
#   poss_osci <- colSums(2^model.A)
#   regbase <- unlist(lapply(1:ncol(model.A), function(x) {min(model.A[,x])}))
#   regscale <- regbase/max(regbase) - (min(regbase/max(regbase)))
#   regplot <- 100 * regscale/max(regscale)
#   ggplot(data.frame(poss_osci),
#          aes(x=sampling_f * 1:length(poss_osci),
#              y=regplot)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
# }

showWave <- function() {
  poss_osci <- colSums(2^model.A)
  regplot <- poss_osci/max(poss_osci)
  regscale <- 100 * regplot - min(regplot)
  regscale2 <- (regscale - min(regscale))
  regscale3 <- regscale2 * 100/max(regscale2)
  regscaled <- 100 - regscale3
  ggplot(data.frame(poss_osci),
         aes(x=sampling_f * 1:length(poss_osci),
             y=regscaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
}

png(paste0("../heat_maps/454-wave-overall.png"), width=1100, height=300)
showWave()
dev.off()

# Scale timepoints by sampling frequency, 3 (hours)

sampling_f <- 3

x.A=sampling_f * 1:ncol(raw.data.A)
y.A=1:nrow(raw.data.A)
x.B=sampling_f * 1:ncol(raw.data.B)
y.B=1:nrow(raw.data.B)

# Rescale these lines for heat map visualisation

z.raw.data.A=apply(raw.data.A, 1, function(point) rescale(point, to=c(0,1)))
z.model.A=apply(model.A, 1, function(point) rescale(point, to=c(0,1)))
z.raw.data.B=apply(raw.data.B, 1, function(point) rescale(point, to=c(0,1)))
z.model.B=apply(model.B, 1, function(point) rescale(point, to=c(0,1)))



# Generate heat maps and save to file

png("heat_maps/PDaUntargRawHM.png", width=1100, height=1600)
image(x.A,y.A,z.raw.data.A[,angle.order.A], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Array probe", main="Untargeted waveform transcript levels, raw data, U-133A Pena-Diaz 2013")
dev.off()

png("heat_maps/PDbUntargRawHM.png", width=1100, height=1600)
image(x.B,y.B,z.raw.data.B[,angle.order.B], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Array probe", main="Untargeted waveform transcript levels, raw data, U-133B Pena-Diaz 2013")
dev.off()

png("heat_maps/PDaUntargModelHM.png", width=1100, height=1600)
image(x.A,y.A,z.model.A[,angle.order.A], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Array probe", main="Untargeted waveform transcript levels, model data, U-133A Pena-Diaz 2013")
dev.off()

png("heat_maps/PDbUntargModelHM.png", width=1100, height=1600)
image(x.B,y.B,z.model.B[,angle.order.B], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Array probe", main="Untargeted waveform transcript levels, model data, U-133B Pena-Diaz 2013")
dev.off()


# Targeted mode

# A function to produce the major components from a dataset

MajorComponents <- function(waveform.processed.dataset, rsq=0.6) {
  good.fit <- waveform.processed.dataset$waveform$fit$r.squared > rsq
  modelled.components = waveform.processed.dataset$waveform$component.used[good.fit,]
  gene.count = colSums(modelled.components)
  major.component = which(gene.count >  10) #use components in at least 100 genes
  groups = apply(modelled.components[,major.component], 1, function(component) paste(major.component[component], collapse=".") ) 
  sort(table(groups), decreasing=T) #lists components in decreasing order of number of containing each combination of major powers
}

# The frequency of oscillation, f, is calculated from period, t
# A generic function is used to calculate the above summaries for a range of ultradian periods
# Default R squared value is set to 0.6 in the function declaration, heatmap filenames are generated from the input frequency

TargettedWaveform <- function(transcript.matrix, period, array.letter=c("A","B"), rsq = 0.6, sampling_f="3") {
  dset <- waveform(transcript.matrix, )
  good.fit <- dset$waveform$fit$r.squared > rsq
  angle.order <- order(dset$angle[good.fit,4])
  raw.data = dset$raw[good.fit,]
  model = dset$waveform$filtered.data[good.fit,]
  x=sampling_f * 1:ncol(raw.data)
  y=1:nrow(raw.data)
  z.raw.data=apply(raw.data, 1, function(point) rescale(point, to=c(0,1)))
  z.model=apply(model, 1, function(point) rescale(point, to=c(0,1)))
  
  png(paste0("heat_maps/PD",tolower(array.letter),"Targ",freq,"RawHM.png"), width=1100, height=1600)
  image(x,y,z.raw.data[,angle.order], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Array probe",
        main=paste0("Untargeted waveform transcript levels, raw data, U-133",array.letter," Pena-Diaz 2013"))
  dev.off()
  
  png(paste0("heat_maps/PD",tolower(array.letter),"Targ",freq,"ModelHM.png"), width=1100, height=1600)
  image(x,y,z.raw.data[,angle.order], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Array probe",
        main=paste0("Untargeted waveform transcript levels, model data, U-133",array.letter," Pena-Diaz 2013"))
  dev.off()
  
}

# Targeted calculations for a given frequency
#   f = Number of cycles per 37 hour time course (as used as parameter to waveform function)
#   t = 36/f, expression oscillation period (approximate, in hours, given for reference)

# f = 5, t = 7.2

dset.A5 <- waveform(fmAmean.mat, target=5)
dset.B5 <- waveform(fmBmean.mat, target=5)

# f = 6, (t = 6)

dset.A6 <- waveform(fmAmean.mat, target=6) # Incalculable
dset.B6 <- waveform(fmBmean.mat, target=6) # Incalculable

# f = 7, (t = 5.1)

dset.A7 <- waveform(fmAmean.mat, target=7) # Incalculable
dset.B7 <- waveform(fmBmean.mat, target=7) # Incalculable

# f = 8, (t = 4.5)

dset.A8 <- waveform(fmAmean.mat, target=8) # Incalculable
dset.B8 <- waveform(fmBmean.mat, target=8) # Incalculable

# f = 9, (t = 4)

dset.A9 <- waveform(fmAmean.mat, target=9) # Incalculable
dset.B9 <- waveform(fmBmean.mat, target=9) # Incalculable

# f = 10, (t = 3.6)

dset.A5 <- waveform(fmAmean.mat, target=10) # Incalculable
dset.B5 <- waveform(fmBmean.mat, target=10) # Incalculable

# f = 11, (t = 3.2)

dset.A11 <- waveform(fmAmean.mat, target=11) # Incalculable
dset.B11 <- waveform(fmBmean.mat, target=11) # Incalculable

# f = 12, (t = 3)

dset.A12 <- waveform(fmAmean.mat, target=12) # Incalculable
dset.B12 <- waveform(fmBmean.mat, target=12) # Incalculable