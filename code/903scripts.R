#!/usr/bin/env R

library("waveform")

# The command line produced version gives a warning about quoted strings
# whit.data <- read.table('E-SMDB-903_data.txt', sep='\t', header = TRUE, fill=TRUE)

# Install the Entrez-GenBank database from BioConductor:

# source("http://bioconductor.org/biocLite.R")
# biocLite("org.Hs.eg.db")

##### To retrieve the original table run:
##### load('whitfield_matrix.Rdata')

# I used LibreOffice export instead to make dataPlusScores_E-SMDB-903.txt though not best practice...

whit.proc.data <- read.table('dataPlusScores_E-SMDB-903.txt', sep='\t', header = TRUE, fill=TRUE)

# Filter columns of Stanford format data to only the log2 normalised mean SMD:LOG_RAT2N_MEAN (numbered ending in 8)

# The following removes the Reporter.REF row (with SMD:LOG_RAT2N_MEAN on it) 
# Transcript descriptors are both row names and a column (for gather, see below)

whit.proc.data.sub <- whit.proc.data[, c(1,grep("8$", colnames(whit.proc.data)))]
whit.proc.data.sub$Scan.REF <- sub("/MAGE.Reporter", "", sub(".*:", "", whit.proc.data.sub$Scan.REF))
colnames(whit.proc.data.sub) <- sub(".*.hr", "", colnames(whit.proc.data.sub))
whit.proc.renamed <- whit.proc.data.sub
colnames(whit.proc.renamed) <- sub("^X", "", sub(".Hela.*", "", colnames(whit.proc.data.sub)))
col.numeric.sort <- function(col.n){col.n[order(as.numeric(sub(".hr","",col.n)))]}
reordered.cols <- c(col.numeric.sort(reorder.colnames)[48],col.numeric.sort(reorder.colnames)[1:47])
whit.proc.reordered <- whit.proc.renamed[reordered.cols][-c(1), ]
colnames(whit.proc.reordered) <- sub(".hr", "", colnames(whit.proc.reordered))
whitfield.matrix <- data.matrix(whit.proc.reordered)

write.table(whitfield.matrix, file = 'whitfield_matrix.tsv', quote = FALSE, sep = "\t", row.names = FALSE)

# NB: default cutoff method where unspecified is p value, and the cutoff is p < 0.05

# The values are mean normalised so can be used in waveform straight away,
# however NA values were found and had to be removed (unusable.probes)

unusable.probes.truth.table <- apply(whitfield.matrix, 1, function(x) any(is.na(x)))
SCAN.REFs <- as.character(whitfield.matrix[,1])
# save(SCAN.REFs, file='Scan_REF_accessions.Rdata')

# 12% of the full transcripts available were discarded due to contaminating NA values on at least 1 time point

# table(unusable.probes.truth.table)
# FALSE  TRUE 
# 39186  4974

# The complement of the set was therefore usable probes, wMmean (abbreviation for whitfield matrix, mean [transcript expression level] values)
wMmean <- whitfield.matrix[!unusable.probes.truth.table, ]

# dim(wMmean)
# 39186    48

# dim(whitfield.matrix)
# 44160    48

waveform.ready.Wmatrix <- wMmean[,-c(1)]
rownames(waveform.ready.Wmatrix) <- wMmean[,c(1)]

# Since this is the first step to reproduce waveform calculations, this matrix is also saved
# save(waveform.ready.Wmatrix, file="waveform_ready_Wmatrix.Rdata")

dsetW = waveform(waveform.ready.Wmatrix)
# save(dsetW, file="E-SMDB-903_analysis/dsetW.Rdata")

# This matrix may be generated as described above, but was saved as .Rdata for ease of reproducibility:

# save(whitfield.matrix, file="whitfield_matrix.Rdata")

# Subset above R squared of 0.6 as recommended in Conesa 2006 for gene selection in microarray time series expression data

rsq = 0.6

good.fit.W <- dsetW$waveform$fit$r.squared > rsq

# Visualisation of heat maps

library(RColorBrewer)
library(scales)

# Order genes by phase angle for visual representation of periodicity

angle.order.W <- order(dsetW$angle[good.fit.W,4])

# NB: "raw" as in not yet [frequency] statistically modelled (these are normalised transcript levels, not raw)

raw.data.W = dsetW$raw[good.fit.W,]
model.W = dsetW$waveform$filtered.data[good.fit.W,]

sampling_f <- 1

x.W=sampling_f * 1:ncol(raw.data.W)
y.W=1:nrow(raw.data.W)
# Rescale these lines for heat map visualisation

z.raw.data.W=apply(raw.data.W, 1, function(point) rescale(point, to=c(0,1)))
z.model.W=apply(model.W, 1, function(point) rescale(point, to=c(0,1)))

# Generate heat maps and save to file

png("heat_maps/WhitfieldUntargRawHM.png", width=1100, height=1600)
image(x.W,y.W,z.raw.data.W[,angle.order.W], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Array probe", main="Untargeted waveform transcript levels, raw data, Whitfield 2002")
dev.off()

png("heat_maps/WhitfieldUntargModelHM.png", width=1100, height=1600)
image(x.W,y.W,z.model.W[,angle.order.W], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Array probe", main="Untargeted waveform transcript levels, model data, Whitfield 2002")
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

TargettedWaveform <- function(freq, transcript.matrix = waveform.ready.Wmatrix, rsq = 0.6) {
  time.course <- 47
  period <- format(round(time.course / freq, 1), nsmall = 1)
  
  dset <- waveform(transcript.matrix, target = freq)
  good.fit <- dset$waveform$fit$r.squared > rsq
  angle.order <- order(dset$angle[good.fit,4])
  raw.data = dset$raw[good.fit,]
  model = dset$waveform$filtered.data[good.fit,]
  x=1:ncol(raw.data) # no scaling required for single-hourly sampling
  y=1:nrow(raw.data)
  z.raw.data=apply(raw.data, 1, function(point) rescale(point, to=c(0,1)))
  z.model=apply(model, 1, function(point) rescale(point, to=c(0,1)))
  
  png(paste0("heat_maps/WhitfieldTarg",freq,"RawHM.png"), width=1100, height=1600)
  image(x,y,z.raw.data[,angle.order], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Array probe",
        main=paste0("Targeted waveform transcript levels, ", period, " hour period (",freq, " cycles), raw data, Whitfield 2002"))
  dev.off()
  
  png(paste0("heat_maps/WhitfieldTarg",freq,"ModelHM.png"), width=1100, height=1600)
  image(x,y,z.raw.data[,angle.order], col=rev(brewer.pal(11, "RdBu")), xlab="Timepoint", ylab="Array probe",
        main=paste0("Targeted waveform transcript levels, ", period, " hour period (",freq, " cycles), model data, Whitfield 2002"))
  dev.off()
  
  return(dset)
}

# Targeted calculations for a given frequency
#   f = Number of cycles per 47 hour time course (as used as parameter to TargettedWaveform function)
#   t = 47/f, expression oscillation period (approximate, in hours, given for reference)

# Application note - if running the following algorithms in succession run rm(dsetW.6) etc. after each to avoid maxing out RAM

# f = 6, (t = 7.8)

dsetW.6 <- waveform(waveform.ready.Wmatrix, target = 6)
save(dsetW.6, file="dsetW_6.Rdata")

# f = 7, (t = 6.7)

dsetW.7 <- waveform(waveform.ready.Wmatrix, target = 7)
save(dsetW.7, file="dsetW_7.Rdata")

# f = 8, (t = 5.9)

dsetW.8 <- waveform(waveform.ready.Wmatrix, target = 8)
save(dsetW.8, file="dsetW_8.Rdata")

# f = 9, (t = 5.2)

dsetW.9 <- waveform(waveform.ready.Wmatrix, target = 9)
save(dsetW.9, file="dsetW_9.Rdata")

# f = 10, (t = 4.7)

dsetW.10 <- waveform(waveform.ready.Wmatrix, target = 10)
save(dsetW.10, file="dsetW_10.Rdata")

# f = 11, (t = 4.2)

dsetW.11 <- waveform(waveform.ready.Wmatrix, target = 11)
save(dsetW.11, file="dsetW_11.Rdata")
rm(dsetW.11)

# f = 12, (t = 3.9)

dsetW.12 <- waveform(waveform.ready.Wmatrix, target = 12)
save(dsetW.12, file="dsetW_12.Rdata")
rm(dsetW.12)

# f = 13, (t = 3.6)

dsetW.13 <- waveform(waveform.ready.Wmatrix, target = 13)
save(dsetW.13, file="dsetW_13.Rdata")
rm(dsetW.13)

# f = 14, (t = 3.35)

dsetW.14 <- waveform(waveform.ready.Wmatrix, target = 14)
save(dsetW.14, file="dsetW_14.Rdata")
rm(dsetW.14)

# f = 15, (t = 3.13)

dsetW.15 <- waveform(waveform.ready.Wmatrix, target = 15)
save(dsetW.15, file="dsetW_15.Rdata")
rm(dsetW.15)

# f = 16, (t = 2.9)

dsetW.16 <- waveform(waveform.ready.Wmatrix, target = 16)
save(dsetW.16, file="dsetW_16.Rdata")
rm(dsetW.16)

# f = 17, (t = 2.7)

dsetW.17 <- waveform(waveform.ready.Wmatrix, target = 17)
save(dsetW.17, file="dsetW_17.Rdata")
rm(dsetW.17)

# f = 18, (t = 2.6)

dsetW.18 <- waveform(waveform.ready.Wmatrix, target = 18)
save(dsetW.18, file="dsetW_18.Rdata")
rm(dsetW.18)

# f = 19, (t = 2.5)

dsetW.19 <- waveform(waveform.ready.Wmatrix, target = 19)
save(dsetW.19, file="dsetW_19.Rdata")
rm(dsetW.19)

# f = 20, (t = 2.3)

dsetW.20 <- waveform(waveform.ready.Wmatrix, target = 20)
save(dsetW.20, file="dsetW_20.Rdata")
rm(dsetW.20)

library(ggplot2)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

showhist <- function(n) {
  data <- eval(as.name(paste0('dsetW.',n)))$waveform$fit$r.squared
  return(
    ggplot(data.frame(data), aes(x=data))
      + geom_density(fill=NA)
      + theme_bw()
      + xlab("R^2")
      + ggtitle(paste0("R-squared histogram for frequency",n))
  )
}

showWave <- function(n) {
  checkLoaded(n)
  assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
  assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
  poss_osci <- colSums(eval(as.name(paste0('model',n))))
  osci_scaling <- poss_osci/max(poss_osci)
  osci_scaling2 <- osci_scaling - min(osci_scaling)
  osci_scaled <- osci_scaling2 * (100/max(osci_scaling2))
  unloadDset(n)
  ggplot(data.frame(poss_osci),
    aes(x=1:length(poss_osci),
    y=osci_scaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
}

saveWave <- function(n) {
  png(paste0("../heat_maps/903-wave-",n,".png"), width=1100, height=300)
  showWave(n)
  dev.off()
}

# lapply(6:20, saveWave) #doesn't work arg
n<-6
png(paste0("../heat_maps/903-wave-",n,".png"), width=1100, height=300)
checkLoaded(n)
assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
poss_osci <- colSums(eval(as.name(paste0('model',n))))
osci_scaling <- poss_osci/max(poss_osci)
osci_scaling2 <- osci_scaling - min(osci_scaling)
osci_scaled <- osci_scaling2 * (100/max(osci_scaling2))
unloadDset(n)
ggplot(data.frame(poss_osci),
       aes(x=1:length(poss_osci),
           y=osci_scaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
dev.off()

n<-7
png(paste0("../heat_maps/903-wave-",n,".png"), width=1100, height=300)
checkLoaded(n)
assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
poss_osci <- colSums(eval(as.name(paste0('model',n))))
osci_scaling <- poss_osci/max(poss_osci)
osci_scaling2 <- osci_scaling - min(osci_scaling)
osci_scaled <- osci_scaling2 * (100/max(osci_scaling2))
unloadDset(n)
ggplot(data.frame(poss_osci),
       aes(x=1:length(poss_osci),
           y=osci_scaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
dev.off()

n<-8
png(paste0("../heat_maps/903-wave-",n,".png"), width=1100, height=300)
checkLoaded(n)
assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
poss_osci <- colSums(eval(as.name(paste0('model',n))))
osci_scaling <- poss_osci/max(poss_osci)
osci_scaling2 <- osci_scaling - min(osci_scaling)
osci_scaled <- osci_scaling2 * (100/max(osci_scaling2))
unloadDset(n)
ggplot(data.frame(poss_osci),
       aes(x=1:length(poss_osci),
           y=osci_scaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
dev.off()

n<-9
png(paste0("../heat_maps/903-wave-",n,".png"), width=1100, height=300)
checkLoaded(n)
assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
poss_osci <- colSums(eval(as.name(paste0('model',n))))
osci_scaling <- poss_osci/max(poss_osci)
osci_scaling2 <- osci_scaling - min(osci_scaling)
osci_scaled <- osci_scaling2 * (100/max(osci_scaling2))
unloadDset(n)
ggplot(data.frame(poss_osci),
       aes(x=1:length(poss_osci),
           y=osci_scaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
dev.off()

n<-10
png(paste0("../heat_maps/903-wave-",n,".png"), width=1100, height=300)
checkLoaded(n)
assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
poss_osci <- colSums(eval(as.name(paste0('model',n))))
osci_scaling <- poss_osci/max(poss_osci)
osci_scaling2 <- osci_scaling - min(osci_scaling)
osci_scaled <- osci_scaling2 * (100/max(osci_scaling2))
unloadDset(n)
ggplot(data.frame(poss_osci),
       aes(x=1:length(poss_osci),
           y=osci_scaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
dev.off()

n<-11
png(paste0("../heat_maps/903-wave-",n,".png"), width=1100, height=300)
checkLoaded(n)
assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
poss_osci <- colSums(eval(as.name(paste0('model',n))))
osci_scaling <- poss_osci/max(poss_osci)
osci_scaling2 <- osci_scaling - min(osci_scaling)
osci_scaled <- osci_scaling2 * (100/max(osci_scaling2))
unloadDset(n)
ggplot(data.frame(poss_osci),
       aes(x=1:length(poss_osci),
           y=osci_scaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
dev.off()

n<-12
png(paste0("../heat_maps/903-wave-",n,".png"), width=1100, height=300)
checkLoaded(n)
assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
poss_osci <- colSums(eval(as.name(paste0('model',n))))
osci_scaling <- poss_osci/max(poss_osci)
osci_scaling2 <- osci_scaling - min(osci_scaling)
osci_scaled <- osci_scaling2 * (100/max(osci_scaling2))
unloadDset(n)
ggplot(data.frame(poss_osci),
       aes(x=1:length(poss_osci),
           y=osci_scaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
dev.off()

n<-13
png(paste0("../heat_maps/903-wave-",n,".png"), width=1100, height=300)
checkLoaded(n)
assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
poss_osci <- colSums(eval(as.name(paste0('model',n))))
osci_scaling <- poss_osci/max(poss_osci)
osci_scaling2 <- osci_scaling - min(osci_scaling)
osci_scaled <- osci_scaling2 * (100/max(osci_scaling2))
unloadDset(n)
ggplot(data.frame(poss_osci),
       aes(x=1:length(poss_osci),
           y=osci_scaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
dev.off()

n<-14
png(paste0("../heat_maps/903-wave-",n,".png"), width=1100, height=300)
checkLoaded(n)
assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
poss_osci <- colSums(eval(as.name(paste0('model',n))))
osci_scaling <- poss_osci/max(poss_osci)
osci_scaling2 <- osci_scaling - min(osci_scaling)
osci_scaled <- osci_scaling2 * (100/max(osci_scaling2))
unloadDset(n)
ggplot(data.frame(poss_osci),
       aes(x=1:length(poss_osci),
           y=osci_scaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
dev.off()

n<-15
png(paste0("../heat_maps/903-wave-",n,".png"), width=1100, height=300)
checkLoaded(n)
assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
poss_osci <- colSums(eval(as.name(paste0('model',n))))
osci_scaling <- poss_osci/max(poss_osci)
osci_scaling2 <- osci_scaling - min(osci_scaling)
osci_scaled <- osci_scaling2 * (100/max(osci_scaling2))
unloadDset(n)
ggplot(data.frame(poss_osci),
       aes(x=1:length(poss_osci),
           y=osci_scaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
dev.off()

n<-16
png(paste0("../heat_maps/903-wave-",n,".png"), width=1100, height=300)
checkLoaded(n)
assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
poss_osci <- colSums(eval(as.name(paste0('model',n))))
osci_scaling <- poss_osci/max(poss_osci)
osci_scaling2 <- osci_scaling - min(osci_scaling)
osci_scaled <- osci_scaling2 * (100/max(osci_scaling2))
unloadDset(n)
ggplot(data.frame(poss_osci),
       aes(x=1:length(poss_osci),
           y=osci_scaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
dev.off()

n<-17
png(paste0("../heat_maps/903-wave-",n,".png"), width=1100, height=300)
checkLoaded(n)
assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
poss_osci <- colSums(eval(as.name(paste0('model',n))))
osci_scaling <- poss_osci/max(poss_osci)
osci_scaling2 <- osci_scaling - min(osci_scaling)
osci_scaled <- osci_scaling2 * (100/max(osci_scaling2))
unloadDset(n)
ggplot(data.frame(poss_osci),
       aes(x=1:length(poss_osci),
           y=osci_scaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
dev.off()

n<-18
png(paste0("../heat_maps/903-wave-",n,".png"), width=1100, height=300)
checkLoaded(n)
assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
poss_osci <- colSums(eval(as.name(paste0('model',n))))
osci_scaling <- poss_osci/max(poss_osci)
osci_scaling2 <- osci_scaling - min(osci_scaling)
osci_scaled <- osci_scaling2 * (100/max(osci_scaling2))
unloadDset(n)
ggplot(data.frame(poss_osci),
       aes(x=1:length(poss_osci),
           y=osci_scaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
dev.off()

n<-19
png(paste0("../heat_maps/903-wave-",n,".png"), width=1100, height=300)
checkLoaded(n)
assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
poss_osci <- colSums(eval(as.name(paste0('model',n))))
osci_scaling <- poss_osci/max(poss_osci)
osci_scaling2 <- osci_scaling - min(osci_scaling)
osci_scaled <- osci_scaling2 * (100/max(osci_scaling2))
unloadDset(n)
ggplot(data.frame(poss_osci),
       aes(x=1:length(poss_osci),
           y=osci_scaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
dev.off()

n<-20
png(paste0("../heat_maps/903-wave-",n,".png"), width=1100, height=300)
checkLoaded(n)
assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
poss_osci <- colSums(eval(as.name(paste0('model',n))))
osci_scaling <- poss_osci/max(poss_osci)
osci_scaling2 <- osci_scaling - min(osci_scaling)
osci_scaled <- osci_scaling2 * (100/max(osci_scaling2))
unloadDset(n)
ggplot(data.frame(poss_osci),
       aes(x=1:length(poss_osci),
           y=osci_scaled)) + geom_line() + theme_bw() + xlab("Time point (hour)") + ylab("Intensity (%)")
dev.off()

checkLoaded <- function(n) {
  # if the dataset corresponding to n cycles isn't available in the environment, load it in
  if (!exists(paste0('dsetW.',n))) {
    load(paste0('dsetW_',n,'.Rdata'), envir=globalenv())
  }
}

saveWave(6)

unloadDset <- function(n) {
  dsetname <- paste0('dsetW.',n) # To use a constructed variable name like this, must supply list = character vector (1-tuple)
  rm(list = dsetname, envir=globalenv()) # otherwise would remove the (non-existent) global variable dsetname
}

# Outputs a truth table, which can be cross-referenced against the name list SCAN.REFs if desired:
listPeriodicTranscripts <- function(n, rsq=0.5) {
  checkLoaded(n)
  assign(paste0('dsetW.',n,'.conf'), eval(as.name(paste0('dsetW.',n)))$waveform$fit$r.squared > rsq)
  unloadDset(n)
  return(eval(as.name(paste0('dsetW.',n,'.conf'))))
  # dsetW.n is now a logical vector corresponding to the values in the Whitfield matrix filtered for unusable probes (wMmean)
}

#raw.whitfield.table <- read.table('Whitfield_processed_data.tsv', sep='\t', header = TRUE, fill=TRUE)
array_description <- read.table('ADF_array_probe_description.tsv', sep='\t', header=TRUE, fill=TRUE)

getGenBankId <- function(reporter.suid)  {
  as.character(array_description[which(
    sub("/MAGE.*", "", array_description$Reporter.Name) == as.character(reporter.suid)
  ),]$Reporter.Database.Entry.embl)
}

# save(array_description, SCAN.REFs, checkLoaded, getGenBankId, listPeriodicTranscripts, namePeriodicTranscripts, unloadDset, file = 'Core_functions.Rdata')

namePeriodicTranscripts <- function(n, rsq=0.5, return.gene.names=c(TRUE,FALSE)) {
  library('org.Hs.eg.db')
  if (length(return.gene.names) > 1) {return.gene.names = TRUE} # will be c(TRUE,FALSE) if unchanged from default hence length 2
  # SCAN.REFs contains the identifiers from wMmean, that is the nth value in SCAN.REFs is the array probe label for the nth item in dsetW.n
  periodic.SCAN.REF.indices <- as.numeric(names(listPeriodicTranscripts(n)[which(listPeriodicTranscripts(n))]))
  if (length(periodic.SCAN.REF.indices) == 0) { return(data.frame()) }
  if (return.gene.names == TRUE) {
    not.yet.gene.names <- lapply(SCAN.REFs[periodic.SCAN.REF.indices], getGenBankId)
    still.not.yet.gene.names <- unlist(lapply(not.yet.gene.names, function(each_row_of_transcripts) {
      each_row_of_transcripts[which(unlist(lapply(
        each_row_of_transcripts, nchar)) != 0 && !grepl('mitoch', each_row_of_transcripts)
                                                 )][1]
      # Some rows come back empty, others come back with mitochondrial misannotation (manual inspection shows no usable ID)
      # The [1] index reference will mean any zero-length (empty) lists will be filled with NA
      # and duplicates (there are some 2-tuple duplicates) will appear normal
    }))
    not.na.but.still.not.yet.gene.names <- still.not.yet.gene.names[which(!is.na(still.not.yet.gene.names))]
    UniGene.mapped <- unlist(lapply(not.na.but.still.not.yet.gene.names, GenBankToUnigene))
    Unigene.IDs <- UniGene.mapped[which(!is.na(UniGene.mapped))]
    Ensembl.genes.with.retirees <- tryCatch(
      {select(org.Hs.eg.db, keys=unique(Unigene.IDs), columns=c("ENSEMBL"), keytype="UNIGENE")},
      error=function(cond){
        return(NA)
      })
    if (is.na(Ensembl.genes.with.retirees)) { return(data.frame()) }
    Ensembl.genes <- Ensembl.genes.with.retirees[which(!is.na(Ensembl.genes.with.retirees$ENSEMBL)),]
    return(Ensembl.genes)
# not this:
# periodic.probe.labels <- SCAN.REFs[periodic.SCAN.REF.indices]
# return(periodic.probe.labels)
  } else {return(periodic.SCAN.REF.indices)}
}

library('stringr')

# Equivalent to the Stanford data table, ran `tr -d '"'` to remove quotation marks
# Only want first 2 columns
Stanford.conversion.table <- read.table('Full_Stanford_E-SMDB-903-table.tsv', sep='\t', header = TRUE, fill=TRUE)[,1:2]
# Stanford.conversion.table <- read.table('ALTFull_Stanford_E-SMDB-903-table.tsv', sep='\t', header = TRUE, fill=TRUE)[,1:2]
Stanford.conversion.table$clone.id <- sub('IMAGE:','',Stanford.conversion.table$UID)
Stanford.conversion.table$gene.name <- sub(" .*", "", sub('.*? ','',Stanford.conversion.table$NAME))
Stanford.conversion.table$unigene.id <- str_match(Stanford.conversion.table$NAME, '(Hs\\..*?) ')[,2]
Stanford.conversion.table$genbank.acc <- str_match(Stanford.conversion.table$NAME, 'Hs\\..*? (.*?) ')[,2]
Stanford.conversion.table <- Stanford.conversion.table[-c(1:2)]

GenBankToUnigene <- function(genbank.id) {
  # unigene.id <- Stanford.conversion.table[which(Stanford.conversion.table$genbank.acc == genbank.id),3]
  # Table is totally worthless, R is awful on tables and dropped about half the values
  full.table.row <- system(paste0("grep ", genbank.id, " Full_Stanford_E-SMDB-903-table.tsv;"), intern = TRUE, ignore.stderr = TRUE)
  if (length(full.table.row) > 0) {
    unigene.id <- str_match(full.table.row, '(Hs\\..*?) ')[,2]
  } else {
    unigene.id <- NA
  }
  return(unigene.id)
}

# GenBankToUnigene('R02578')
# [1] "Hs.107387"

# save(Stanford.conversion.table, file="Stanford_conversion_table.Rdata")
# write.table(Stanford.conversion.table, file = 'ALTStanford_conversion_table.tsv', quote = FALSE, sep = "\t", row.names = TRUE)

write.table(Stanford.conversion.table, file = 'Stanford_conversion_table.tsv', quote = FALSE, sep = "\t", row.names = TRUE)
# function checkline () { sed "$1"'n;d' Stanford_conversion_table.tsv | sed 's/\tNA/\t/g' | cut -f 4 | sed '/[a-z]/d'; }
# function Stanfordcheck() { checkline "$1" | (read lineinput; if [ ! $(echo "$lineinput" | wc -c) -lt 2 ]; then echo "$lineinput"; fi) | awk '{print "Hs.* "$0 }' | grep -nf - Full_Stanford_E-SMDB-903-table.tsv | cut -f 1-2; }

# There are 42,921 lines in Full_Stanford_E-SMDB-903-table.tsv
#    for num in {1..42921}; do printf "$num\t$(checkline $num)\n" >> myfile.txt; done

library('org.Hs.eg.db')

GenBank.to.EntrezGenes <- as.list(org.Hs.egACCNUM2EG)
# save(GenBank.to.EntrezGenes, 'GenBank_to_EntrezGenes.Rdata')
# doesn't work, must convert to AccNum manually

listMeanPCC <- function(n, rsq=0.5) {
  checkLoaded(n)
  assign(paste0('dsetW.',n,'.pcc'), eval(as.name(paste0('dsetW.',n)))$waveform$pearson$correlation[!is.na(eval(as.name(paste0('dsetW.',n)))$waveform$pearson$correlation)])
  unloadDset(n)
  return(eval(as.name(paste0('dsetW.',n,'.pcc'))))
}

listMeanPVal <- function(n, rsq=0.5) {
  checkLoaded(n)
  assign(paste0('dsetW.',n,'.pval'), eval(as.name(paste0('dsetW.',n)))$waveform$pearson$p.value[!is.na(eval(as.name(paste0('dsetW.',n)))$waveform$pearson$p.value)])
  unloadDset(n)
  return(eval(as.name(paste0('dsetW.',n,'.pval'))))
}

# R^2 > 0.5

osc_freqs <- 6:20

transcripts_Rsq_0_5 <- unlist(lapply(osc_freqs, function(osc_freq) {
  nrow(namePeriodicTranscripts(osc_freq))
}))

# R^2 > 0.6

high_conf_transcripts_Rsq_0_6 <- unlist(lapply(osc_freqs, function(osc_freq) {
  nrow(namePeriodicTranscripts(osc_freq, rsq=0.6))
}))

mean_PVal_Rsq_0_5 <- unlist(lapply(osc_freqs, function(osc_freq, rsq=0.5) {mean(listMeanPVal(osc_freq)) }))
mean_PVal_Rsq_0_6 <- unlist(lapply(osc_freqs, function(osc_freq, rsq=0.6) {mean(listMeanPVal(osc_freq)) }))

time.course <- 47
osc_periods <- format(round(time.course / osc_freqs, 1), nsmall = 1)

mean_PVal_Rsq_0_5_rounded <- format(round(mean_PVal_Rsq_0_5, 3), nsmall = 3)
mean_PVal_Rsq_0_6_rounded <- format(round(mean_PVal_Rsq_0_6, 3), nsmall = 3)

rsq_freqs_confidence_table <- t(data.frame(osc_freqs, osc_periods, transcripts_Rsq_0_5, high_conf_transcripts_Rsq_0_6,
                                           mean_PVal_Rsq_0_5_rounded, mean_PVal_Rsq_0_6_rounded))


# osc_freqs                 6    7    8    9   10   11   12   13   14    15    16    17    18    19    20
# high_conf_list_Rsq_0.5 3431 1269  705  564  329  611   94   47   47    47    47     0    47     0     0
# high_conf_list_Rsq_0.6  329   94   47   47    0   47    0    0    0     0     0     0     0     0     0

write.table(rsq_freqs_confidence_table, file = 'rsq_freqs_confidence_table.tsv', quote = FALSE, sep = "\t", row.names = TRUE)
save(rsq_freqs_confidence_table, file='rsq_freqs_confidence_table.Rdata')

# Confirmed results:
# tail --lines=+14 mirgate-query-results-*.csv | cut -d ',' -f 1,4,12 | sed '/,0/d'

#==> mirgate-query-results-11.csv <==
#  Input Gene,Mirna,Confirmed Predictions
#ENSG00000124785,hsa-mir-28-5p,1
#MIMAT0000085

#==> mirgate-query-results-6.csv <==
#  Input Gene,Mirna,Confirmed Predictions
#ENSG00000206503,hsa-mir-148b-3p,2
#MIMAT0000759

#==> mirgate-query-results-7.csv <==
#  Input Gene,Mirna,Confirmed Predictions
# ENSG00000197951,hsa-mir-193b-3p,1
# MIMAT0002819

# cat mirgate-query-mirna-names.txt | sort | uniq | sed '/kshv-mir/d' | sed '/ebv-mir/d' > mirgate-query-filtered-sorted-names.txt
# cat aliases.txt | python alias_process.py | grep 'hsa-' > aliases_longform.txt

mirbase.aliases <- read.table('aliases_longform.txt', sep='\t')
colnames(mirbase.aliases) <- c("Accession", "Name")
mirbase.aliases$UpperName <- toupper(mirbase.aliases$Name)

# 6
mirgate.output.6 <- read.table('mirgate-query-filtered-sorted-names-6.txt', sep='\t')
colnames(mirgate.output.6) <- "miRGate.name"
mirgate.output.6$UpperName <- toupper(mirgate.output.6$miRGate.name)
mirgate.mirbase.mapping.6 <- merge.data.frame(mirgate.output.6, mirbase.aliases)
dup.poss.6 <- do.call("rbind", lapply(mirgate.mirbase.mapping.6[which(base::duplicated(mirgate.mirbase.mapping.6$UpperName)),1], function(upper.name) {
  mirgate.mirbase.mapping.6[which(mirgate.mirbase.mapping.6$UpperName == upper.name),]
}))
duplist.per.id.6 <- lapply(mirgate.mirbase.mapping.6[which(base::duplicated(mirgate.mirbase.mapping.6$UpperName)),1], function(upper.name) {
  dups <- mirgate.mirbase.mapping.6[which(mirgate.mirbase.mapping.6$UpperName == upper.name),]
  duplist.to.keep <- dups[which(as.character(dups$miRGate.name) != as.character(dups$Name)),] # has a pre-miRNA accession, unwanted
  if (nrow(duplist.to.keep == 2)) {
    if (length(unique(as.character(dups$Accession))) == 1) {
      # would be different if one was a pre-miRNA (i.e. length 1) here they're just identical so drop the second
      duplist.to.keep <- dups[1,]
    } else {
      # they're not simply the same in mir/miR forms, so provide both (i.e. leave duplist.to.keep as dups) but warn
      warning(paste0("Check that both these accessions exist:",paste0(as.character(dups$Accession), collapse=", ")))
    }
  }
  return(duplist.to.keep)
})
duplist.6 <- do.call("rbind", duplist.per.id.6)
good.rows.6 <- rownames(duplist.6)
bad.rows.6 <- setdiff(rownames(dup.poss.6), good.rows.6)
use.rows.6 <- setdiff(rownames(mirgate.mirbase.mapping.6), bad.rows.6)
use.list.6 <- mirgate.mirbase.mapping.6[use.rows.6,]
save(mirgate.output.6, mirgate.mirbase.mapping.6, dup.poss.6, duplist.6, duplist.per.id.6, use.rows.6, use.list.6, file = "mirgate_mapping_6.Rdata")
rm(mirgate.output.6, mirgate.mirbase.mapping.6, duplist.6, duplist.per.id.6) # keep use.list

# 7
mirgate.output.7 <- read.table('mirgate-query-filtered-sorted-names-7.txt', sep='\t')
colnames(mirgate.output.7) <- "miRGate.name"
mirgate.output.7$UpperName <- toupper(mirgate.output.7$miRGate.name)
mirgate.mirbase.mapping.7 <- merge.data.frame(mirgate.output.7, mirbase.aliases)
dup.poss.7 <- do.call("rbind", lapply(mirgate.mirbase.mapping.7[which(base::duplicated(mirgate.mirbase.mapping.7$UpperName)),1], function(upper.name) {
  mirgate.mirbase.mapping.7[which(mirgate.mirbase.mapping.7$UpperName == upper.name),]
}))
duplist.per.id.7 <- lapply(mirgate.mirbase.mapping.7[which(base::duplicated(mirgate.mirbase.mapping.7$UpperName)),1], function(upper.name) {
  dups <- mirgate.mirbase.mapping.7[which(mirgate.mirbase.mapping.7$UpperName == upper.name),]
  duplist.to.keep <- dups[which(as.character(dups$miRGate.name) != as.character(dups$Name)),] # has a pre-miRNA accession, unwanted
  if (nrow(duplist.to.keep == 2)) {
    if (length(unique(as.character(dups$Accession))) == 1) {
      # would be different if one was a pre-miRNA (i.e. length 1) here they're just identical so drop the second
      duplist.to.keep <- dups[1,]
    } else {
      # they're not simply the same in mir/miR forms, so provide both (i.e. leave duplist.to.keep as dups) but warn
      warning(paste0("Check that both these accessions exist:",paste0(as.character(dups$Accession), collapse=", ")))
    }
  }
  return(duplist.to.keep)
})
duplist.7 <- do.call("rbind", duplist.per.id.7)
good.rows.7 <- rownames(duplist.7)
bad.rows.7 <- setdiff(rownames(dup.poss.7), good.rows.7)
use.rows.7 <- setdiff(rownames(mirgate.mirbase.mapping.7), bad.rows.7)
use.list.7 <- mirgate.mirbase.mapping.7[use.rows.7,]
save(mirgate.output.7, mirgate.mirbase.mapping.7, dup.poss.7, duplist.7, duplist.per.id.7, use.rows.7, use.list.7, file = "mirgate_mapping_7.Rdata")
rm(mirgate.output.7, mirgate.mirbase.mapping.7, duplist.7, duplist.per.id.7) # keep use.list

# 8
mirgate.output.8 <- read.table('mirgate-query-filtered-sorted-names-8.txt', sep='\t')
colnames(mirgate.output.8) <- "miRGate.name"
mirgate.output.8$UpperName <- toupper(mirgate.output.8$miRGate.name)
mirgate.mirbase.mapping.8 <- merge.data.frame(mirgate.output.8, mirbase.aliases)
dup.poss.8 <- do.call("rbind", lapply(mirgate.mirbase.mapping.8[which(base::duplicated(mirgate.mirbase.mapping.8$UpperName)),1], function(upper.name) {
  mirgate.mirbase.mapping.8[which(mirgate.mirbase.mapping.8$UpperName == upper.name),]
}))
duplist.per.id.8 <- lapply(mirgate.mirbase.mapping.8[which(base::duplicated(mirgate.mirbase.mapping.8$UpperName)),1], function(upper.name) {
  dups <- mirgate.mirbase.mapping.8[which(mirgate.mirbase.mapping.8$UpperName == upper.name),]
  duplist.to.keep <- dups[which(as.character(dups$miRGate.name) != as.character(dups$Name)),] # has a pre-miRNA accession, unwanted
  if (nrow(duplist.to.keep == 2)) {
    if (length(unique(as.character(dups$Accession))) == 1) {
      # would be different if one was a pre-miRNA (i.e. length 1) here they're just identical so drop the second
      duplist.to.keep <- dups[1,]
    } else {
      # they're not simply the same in mir/miR forms, so provide both (i.e. leave duplist.to.keep as dups) but warn
      warning(paste0("Check that both these accessions exist:",paste0(as.character(dups$Accession), collapse=", ")))
    }
  }
  return(duplist.to.keep)
})
duplist.8 <- do.call("rbind", duplist.per.id.8)
good.rows.8 <- rownames(duplist.8)
bad.rows.8 <- setdiff(rownames(dup.poss.8), good.rows.8)
use.rows.8 <- setdiff(rownames(mirgate.mirbase.mapping.8), bad.rows.8)
use.list.8 <- mirgate.mirbase.mapping.8[use.rows.8,]
save(mirgate.output.8, mirgate.mirbase.mapping.8, dup.poss.8, duplist.8, duplist.per.id.8, use.rows.8, use.list.8, file = "mirgate_mapping_8.Rdata")
rm(mirgate.output.8, mirgate.mirbase.mapping.8, duplist.8, duplist.per.id.8) # keep use.list

# 9
mirgate.output.9 <- read.table('mirgate-query-filtered-sorted-names-9.txt', sep='\t')
colnames(mirgate.output.9) <- "miRGate.name"
mirgate.output.9$UpperName <- toupper(mirgate.output.9$miRGate.name)
mirgate.mirbase.mapping.9 <- merge.data.frame(mirgate.output.9, mirbase.aliases)
dup.poss.9 <- do.call("rbind", lapply(mirgate.mirbase.mapping.9[which(base::duplicated(mirgate.mirbase.mapping.9$UpperName)),1], function(upper.name) {
  mirgate.mirbase.mapping.9[which(mirgate.mirbase.mapping.9$UpperName == upper.name),]
}))
duplist.per.id.9 <- lapply(mirgate.mirbase.mapping.9[which(base::duplicated(mirgate.mirbase.mapping.9$UpperName)),1], function(upper.name) {
  dups <- mirgate.mirbase.mapping.9[which(mirgate.mirbase.mapping.9$UpperName == upper.name),]
  duplist.to.keep <- dups[which(as.character(dups$miRGate.name) != as.character(dups$Name)),] # has a pre-miRNA accession, unwanted
  if (nrow(duplist.to.keep == 2)) {
    if (length(unique(as.character(dups$Accession))) == 1) {
      # would be different if one was a pre-miRNA (i.e. length 1) here they're just identical so drop the second
      duplist.to.keep <- dups[1,]
    } else {
      # they're not simply the same in mir/miR forms, so provide both (i.e. leave duplist.to.keep as dups) but warn
      warning(paste0("Check that both these accessions exist:",paste0(as.character(dups$Accession), collapse=", ")))
    }
  }
  return(duplist.to.keep)
})
duplist.9 <- do.call("rbind", duplist.per.id.9)
good.rows.9 <- rownames(duplist.9)
bad.rows.9 <- setdiff(rownames(dup.poss.9), good.rows.9)
use.rows.9 <- setdiff(rownames(mirgate.mirbase.mapping.9), bad.rows.9)
use.list.9 <- mirgate.mirbase.mapping.9[use.rows.9,]
save(mirgate.output.9, mirgate.mirbase.mapping.9, dup.poss.9, duplist.9, duplist.per.id.9, use.rows.9, use.list.9, file = "mirgate_mapping_9.Rdata")
rm(mirgate.output.9, mirgate.mirbase.mapping.9, duplist.9, duplist.per.id.9) # keep use.list

# 10
mirgate.output.10 <- read.table('mirgate-query-filtered-sorted-names-10.txt', sep='\t')
colnames(mirgate.output.10) <- "miRGate.name"
mirgate.output.10$UpperName <- toupper(mirgate.output.10$miRGate.name)
mirgate.mirbase.mapping.10 <- merge.data.frame(mirgate.output.10, mirbase.aliases)
dup.poss.10 <- do.call("rbind", lapply(mirgate.mirbase.mapping.10[which(base::duplicated(mirgate.mirbase.mapping.10$UpperName)),1], function(upper.name) {
  mirgate.mirbase.mapping.10[which(mirgate.mirbase.mapping.10$UpperName == upper.name),]
}))
duplist.per.id.10 <- lapply(mirgate.mirbase.mapping.10[which(base::duplicated(mirgate.mirbase.mapping.10$UpperName)),1], function(upper.name) {
  dups <- mirgate.mirbase.mapping.10[which(mirgate.mirbase.mapping.10$UpperName == upper.name),]
  duplist.to.keep <- dups[which(as.character(dups$miRGate.name) != as.character(dups$Name)),] # has a pre-miRNA accession, unwanted
  if (nrow(duplist.to.keep == 2)) {
    if (length(unique(as.character(dups$Accession))) == 1) {
      # would be different if one was a pre-miRNA (i.e. length 1) here they're just identical so drop the second
      duplist.to.keep <- dups[1,]
    } else {
      # they're not simply the same in mir/miR forms, so provide both (i.e. leave duplist.to.keep as dups) but warn
      warning(paste0("Check that both these accessions exist:",paste0(as.character(dups$Accession), collapse=", ")))
    }
  }
  return(duplist.to.keep)
})
duplist.10 <- do.call("rbind", duplist.per.id.10)
good.rows.10 <- rownames(duplist.10)
bad.rows.10 <- setdiff(rownames(dup.poss.10), good.rows.10)
use.rows.10 <- setdiff(rownames(mirgate.mirbase.mapping.10), bad.rows.10)
use.list.10 <- mirgate.mirbase.mapping.10[use.rows.10,]
save(mirgate.output.10, mirgate.mirbase.mapping.10, dup.poss.10, duplist.10, duplist.per.id.10, use.rows.10, use.list.10, file = "mirgate_mapping_10.Rdata")
rm(mirgate.output.10, mirgate.mirbase.mapping.10, duplist.10, duplist.per.id.10) # keep use.list

# 11
mirgate.output.11 <- read.table('mirgate-query-filtered-sorted-names-11.txt', sep='\t')
colnames(mirgate.output.11) <- "miRGate.name"
mirgate.output.11$UpperName <- toupper(mirgate.output.11$miRGate.name)
mirgate.mirbase.mapping.11 <- merge.data.frame(mirgate.output.11, mirbase.aliases)
dup.poss.11 <- do.call("rbind", lapply(mirgate.mirbase.mapping.11[which(base::duplicated(mirgate.mirbase.mapping.11$UpperName)),1], function(upper.name) {
  mirgate.mirbase.mapping.11[which(mirgate.mirbase.mapping.11$UpperName == upper.name),]
}))
duplist.per.id.11 <- lapply(mirgate.mirbase.mapping.11[which(base::duplicated(mirgate.mirbase.mapping.11$UpperName)),1], function(upper.name) {
  dups <- mirgate.mirbase.mapping.11[which(mirgate.mirbase.mapping.11$UpperName == upper.name),]
  duplist.to.keep <- dups[which(as.character(dups$miRGate.name) != as.character(dups$Name)),] # has a pre-miRNA accession, unwanted
  if (nrow(duplist.to.keep == 2)) {
    if (length(unique(as.character(dups$Accession))) == 1) {
      # would be different if one was a pre-miRNA (i.e. length 1) here they're just identical so drop the second
      duplist.to.keep <- dups[1,]
    } else {
      # they're not simply the same in mir/miR forms, so provide both (i.e. leave duplist.to.keep as dups) but warn
      warning(paste0("Check that both these accessions exist:",paste0(as.character(dups$Accession), collapse=", ")))
    }
  }
  return(duplist.to.keep)
})
duplist.11 <- do.call("rbind", duplist.per.id.11)
good.rows.11 <- rownames(duplist.11)
bad.rows.11 <- setdiff(rownames(dup.poss.11), good.rows.11)
use.rows.11 <- setdiff(rownames(mirgate.mirbase.mapping.11), bad.rows.11)
use.list.11 <- mirgate.mirbase.mapping.11[use.rows.11,]
save(mirgate.output.11, mirgate.mirbase.mapping.11, dup.poss.11, duplist.11, duplist.per.id.11, use.rows.11, use.list.11, file = "mirgate_mapping_11.Rdata")
rm(mirgate.output.11, mirgate.mirbase.mapping.11, duplist.11, duplist.per.id.11) # keep use.list

high.conf.set <- read.table('high_conf_miRNA_ID_list.txt', sep='\t')
is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))

# Run from miRBase/20/genomes after decompressing the gff files:
# cat high_conf_miRNA_ID_list.txt | awk '{ print "Derives_from="$0 }' | grep -f - * > high_conf_miRNA_genome_mappings.txt

# cat high_conf_miRNA_genome_mappings.txt | awk '{ split($0,a,"ID="); split(a[2],b,";"); split(b[2],c,"="); split($0,d,"Derives_from="); print d[2]"\t"c[2]; }' > high_conf_miRNA_table.tsv

high.conf.accessions <- read.table('high_conf_miRNA_table.tsv', sep='\t')
colnames(high.conf.accessions) <- c("Pre","Mature")

highConfListings <- function (input.list) {
  # input.list is a subset of mirgate.mirbase.mapping numbered
  confident.list <- lapply(as.character(high.conf.accessions$Mature), function(accession) {
    matched.rows <- input.list[which(input.list$Accession == accession),]
    if (nrow(matched.rows) > 0) {
      return(matched.rows)
    }
  })
  return.list <- do.call("rbind",Filter(Negate(is.NullOb), confident.list))
  return(return.list)
}

write.table(use.list.6, file = "uselist6.tsv", sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

high.conf.listings.6 <- highConfListings(use.list.6)
high.conf.listings.7 <- highConfListings(use.list.7)
high.conf.listings.8 <- highConfListings(use.list.8)
high.conf.listings.9 <- highConfListings(use.list.9)
high.conf.listings.10 <- highConfListings(use.list.10)
high.conf.listings.11 <- highConfListings(use.list.11)

# length(unique(as.character(high.conf.listings.6$Accession)))
# 173
# length(unique(as.character(high.conf.listings.7$Accession)))
# 287
# length(unique(as.character(high.conf.listings.8$Accession)))
# 268
# length(unique(as.character(high.conf.listings.9$Accession)))
# 42
# length(unique(as.character(high.conf.listings.10$Accession)))
# 34
# length(unique(as.character(high.conf.listings.11$Accession)))
# 239

save(high.conf.listings.6, high.conf.listings.7, high.conf.listings.8, high.conf.listings.9,
     high.conf.listings.10, high.conf.listings.11, file="High_confidence_listings.Rdata")

write.table(high.conf.listings.6, file = "high_conf_listings_6.tsv", sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(high.conf.listings.7, file = "high_conf_listings_7.tsv", sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(high.conf.listings.8, file = "high_conf_listings_8.tsv", sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(high.conf.listings.9, file = "high_conf_listings_9.tsv", sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(high.conf.listings.10, file = "high_conf_listings_10.tsv", sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(high.conf.listings.11, file = "high_conf_listings_11.tsv", sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# grep 'mir-9-' mirgate-query-results-*.csv | cut -d\, -f 1 | sort -h | uniq

# mirgate-query-results-6.csv:ENSG00000239857
# mirgate-query-results-7.csv:ENSG00000183647
# mirgate-query-results-7.csv:ENSG00000197951
# mirgate-query-results-8.csv:ENSG00000109991
# mirgate-query-results-8.csv:ENSG00000173575
# mirgate-query-results-11.csv:ENSG00000054219
# mirgate-query-results-11.csv:ENSG00000124785
# mirgate-query-results-11.csv:ENSG00000241399
# mirgate-query-results-11.csv:ENSG00000248672

#   duplist[which(duplist$UpperName == 'HSA-MIR-577'),]
#   duplist[which(duplist$UpperName == 'HSA-MIR-1249'),]
# Confirmed mature used via API query since output lowercases names
#   http://mirgate.bioinfo.cnio.es/ResT/API/human/miRNA_info/hsa-mir-1249/

# library(RColorBrewer)
# library(scales)

showheatmap <- function(n,raw_or_model=c("raw","model")) {
  if (length(raw_or_model)>1) {raw_or_model = "raw"} #default setting would otherwise be both
  checkLoaded(n)
  assign(paste0('poor.fit.',n), eval(as.name(paste0('dsetW.',n)))$waveform$fit$p.value > 0.01)
  assign(paste0('angle.order',n), order(eval(as.name(paste0('dsetW.',n)))$angle[!eval(as.name(paste0('poor.fit.',n))),2]))
  assign(paste0('raw.data.',n), eval(as.name(paste0('dsetW.',n)))$raw[!eval(as.name(paste0('poor.fit.',n))),])
  assign(paste0('model',n), eval(as.name(paste0('dsetW.',n)))$waveform$filtered.data[!eval(as.name(paste0('poor.fit.',n))),])
  assign(paste0('x.',n), 1:ncol(eval(as.name(paste0('raw.data.',n)))))
  assign(paste0('y.',n), 1:nrow(eval(as.name(paste0('raw.data.',n)))))
  assign(paste0('z.raw.data.',n), apply(eval(as.name(paste0('raw.data.',n))), 1, function(x) rescale(x, to=c(0,1))))
  assign(paste0('z.model',n), apply(eval(as.name(paste0('model',n))), 1, function(x) rescale(x, to=c(0,1))))
  hmap <- image(
    eval(as.name(paste0('x.',n))),
    eval(as.name(paste0('y.',n))),
    eval(as.name(paste0('z.raw.data.',n)))[,eval(as.name(paste0('angle.order',n)))],
    col=rev(brewer.pal(11, "RdBu")),
    xlab="Timepoint",
    ylab="Gene ID",
    main=paste0("Raw data (Targeted, ",n," cycles)")
  )
  unloadDset(n)
  return(hmap)
}

png("../heat_maps/903-heat_map_6.png", width=1100, height=1600)
showheatmap(6)
dev.off()

saveheatmap <- function(n, raw_or_model="raw") {
  png(paste0("../heat_maps/903-heat_map_",n,".png"), width=1100, height=1600)
  showheatmap(n, raw_or_model)
  dev.off()
}

savemodelheatmap <- function(n, raw_or_model="model") {
  png(paste0("../heat_maps/903-model-heat_map_",n,".png"), width=1100, height=1600)
  showheatmap(n, raw_or_model)
  dev.off()
}

lapply(6:20, saveheatmap)
lapply(6:20, savemodelheatmap)
