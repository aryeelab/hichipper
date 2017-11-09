#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(data.table)))

args <- commandArgs(trailingOnly = TRUE)
resfile <- args[1]
peaksfile <- args[2]

pad <- 0 # taken care of in the Python execution

# Import Restriction Fragments / Convert to GRanges
if(endsWith(resfile, "gz")){
  resFrags <- data.frame(fread(paste0(input = 'zcat < ', resfile), header = FALSE))
} else {
  resFrags <- data.frame(fread(resfile,  header = FALSE))
}
resFrags_g <- makeGRangesFromDataFrame(setNames(data.frame(
  resFrags[, 1], resFrags[, 2], resFrags[, 3]), c("seqnames", "start", "end")))

# Import Peaks / Convert to GRanges
if(endsWith(peaksfile, "gz")){
  peaks <- data.frame(fread(paste0(input = 'zcat < ', peaksfile), header = FALSE))
} else {
  peaks <- data.frame(fread(peaksfile,  header = FALSE))
}
peaks$start <- peaks$V2 - pad
peaks$end <- peaks$V3 + pad
peaks <- makeGRangesFromDataFrame(peaks, seqnames.field = "V1", start.field = "start", end.field = "end")

# Overlap the two; create GRanges of extreme fragments
ov <- findOverlaps(resFrags_g, peaks)
mins <- unname(tapply(queryHits(ov), subjectHits(ov), min))
maxs <-  unname(tapply(queryHits(ov), subjectHits(ov), max))

anchorsfin_g <- makeGRangesFromDataFrame(setNames(data.frame(
  resFrags[mins, 1], resFrags[mins, 2], resFrags[maxs, 3]), c("seqnames", "start", "end")))

anchors_g <- reduce(anchorsfin_g)

# Filter anchors that were on terminal restriction fragments
anchors_gfilt <- anchors_g[width(anchors_g) < 50000] 
cat(paste0("Anchors removed due to excessive size (likely at ends of chromosomes): ",
           as.character(length(anchors_g) - length(anchors_gfilt)), " \n"))
write.table(as.data.frame(anchors_gfilt)[,c(1,2,3)], row.names = FALSE, col.names = FALSE, quote = FALSE, 
            sep = "\t", file = paste0(peaksfile, "rf.tmp"))

