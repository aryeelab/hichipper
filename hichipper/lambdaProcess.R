#!/usr/bin/env Rscript

# File / set of functions that takes restriction enzyme fragment positions
# and bedgraph pileup files from MACS2 and writes modified files

suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(reshape2)))

args <- commandArgs(trailingOnly = TRUE)
resfile <- args[1]
treatmentfile <- args[2]
backgroundfile <- args[3]
outdir <- args[4]
K <- 1000000 # Number of loci to sample from 

# Import Restriction Fragments / Convert to GRanges
if(endsWith(resfile, "gz")){
  resFrags <- data.frame(fread(paste0(input = 'zcat < ', resfile), header = FALSE))
} else {
  resFrags <- data.frame(fread(resfile,  header = FALSE))
}

resSites <- makeGRangesFromDataFrame(setNames(data.frame(
  resFrags[, 1], resFrags[, 2], resFrags[, 2]), c("seqnames", "start", "end")))

# Function that computes a vector of ratios per distance, the total mean,
# and returns the imported background peakset

computeRatioEtc <- function(treatmentfile, backgroundfile){
  
  # Import Treatment / Convert to GRanges
  if(endsWith(treatmentfile, "gz")){
    txt <- data.frame(fread(paste0(input = 'zcat < ', treatmentfile), header = FALSE))
  } else {
    txt <- fread(treatmentfile)
    txt <- txt[txt$V1 %in% seqnames(resSites),]
  }
  txt <- makeGRangesFromDataFrame(txt, seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)
  
  # Import Control / Convert to GRanges
  if(endsWith(backgroundfile, "gz")){
    cont <- data.frame(fread(paste0(input = 'zcat < ', backgroundfile), header = FALSE))
  } else {
    cont <- fread(backgroundfile)
    cont <- cont[cont$V1 %in% seqnames(resSites),]
  }
  cont <- makeGRangesFromDataFrame(cont, seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)
  
  #Sample sites from the genome; make GRanges object
  chroms <- intersect(levels(seqnames(txt)), levels(seqnames(resSites)))
  begins <- tapply(start(txt)[mcols(txt)[,1] != 0 ] , seqnames(txt)[mcols(txt)[,1] != 0 ], min)[chroms]
  ends <- tapply(end(txt)[mcols(txt)[,1] != 0 ] , seqnames(txt)[mcols(txt)[,1] != 0 ], max)[chroms]
  dists <- ends - begins
  nsamples <- round(K * (dists/sum(as.numeric(dists))),digits = 0)
  samplefn <- function (chrom) sample((begins[chrom]:ends[chrom]), nsamples[chrom], replace = FALSE, prob = NULL)
  randGRanges <- makeGRangesFromDataFrame(reshape2::melt(setNames(lapply(chroms, samplefn), chroms)),
                                          seqnames.field = "L1", start.field = "value", end.field = "value")
  
  # Compute vector of distances to cutsite from random peaks; aggregate controls and treatment counts
  nndist <- mcols(distanceToNearest(randGRanges, resSites))[,1]
  txtVals <- mcols(txt)[findOverlaps(randGRanges, txt, select = "arbitrary"),1]
  contVals <- mcols(cont)[findOverlaps(randGRanges, cont, select = "arbitrary"),1]
  
  # Compute nearest neighbor for background
  mid <- (start(cont) + end(cont))/2
  backmid <- as(data.frame(chrom=seqnames(cont), start=mid, end=mid), "GRanges")
  mcols(cont)$dist <- mcols(distanceToNearest(backmid, resSites))[,1]
  
  vals <- tapply(txtVals/contVals, nndist, mean)
  
  # Write adjusted treatment (only filters out regions where no restriction enzyme information is contained)
  outdf <- data.frame(x1 = seqnames(txt), x2 = start(txt), x3 = end(txt), x4 = round(mcols(txt)[,1],2))
  write.table(outdf, file = paste0(outdir, "/adjustedTreatment.bdg.tmp"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  return(list(vals[1:750], mean(txtVals/contVals), cont))
}

hichipall <- computeRatioEtc(treatmentfile, backgroundfile)
hichipratio <- hichipall[[1]]
hichipmean <-  hichipall[[2]]
hichipbg   <-  hichipall[[3]]

# Compute smooth trajectory; scale factor per background set
smoothingSpline <- smooth.spline(as.numeric(names(hichipratio)), unname(hichipratio), spar=0.35)
scaleVec <- c(smoothingSpline$y/hichipmean, rep(1, max(mcols(hichipbg)$dist)- length(smoothingSpline$y) + 1))
names(scaleVec) <- as.character(0:(length(scaleVec)-1))

# Make scale factor; write to new file
dist <- mcols(hichipbg)[,"dist"]
sf <- scaleVec[dist+1] # shift for zero
outdf <- data.frame(x1 = seqnames(hichipbg), x2 = start(hichipbg), x3 = end(hichipbg), x4 = round(mcols(hichipbg)[,1]*sf,2))
globalLambda <- min(mcols(hichipbg)[,1])
outdf[outdf[,4] < globalLambda,4] <- globalLambda
write.table(outdf, file = paste0(outdir, "/adjustedBackground.bdg.tmp"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(data.frame(round(hichipmean,2)),  file = paste0(outdir, "/globalMeanRatio.tmp"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

tdf <- head(hichipratio, 750)
write.table(data.frame(dist = names(tdf), vals = round(unname(tdf),3)), file = paste0(outdir, "/observedRatios.tmp"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

