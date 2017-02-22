library(GenomicRanges)
library(data.table)

resfile <- "../data/mm9_MboI_resfrag.bed.gz"
chipseqfile <- "../data/mESC_SMC1_chipseq.broadPeak"
pad <- 1000

# Import Restriction Fragments / Convert to GRanges
if(endsWith(resfile, "gz")){
  resFrags <- fread(paste0(input = 'zcat < ', resfile), header = FALSE)
} else {
  resFrags <- fread(resfile,  header = FALSE)
}
resFrags_g <- makeGRangesFromDataFrame(setNames(data.frame(
  resFrags[, 1], resFrags[, 2], resFrags[, 3]), c("seqnames", "start", "end")))

# Import ChIP-Seq / Convert to GRanges
if(endsWith(chipseqfile, "gz")){
  chipseq <- fread(paste0(input = 'zcat < ', chipseqfile), header = FALSE)
} else {
  chipseq <- fread(chipseqfile,  header = FALSE)
}
chipseq$start <- chipseq$V2 - pad
chipseq$End <- chipseq$V3 + pad
chipseq <- makeGRangesFromDataFrame(chipseq, seqnames.field = "V1", start.field = "start", end.field = "end")

# Overlap the two; create GRanges of extreme fragments
ov <- findOverlaps(resFrags_g, chipseq)
mins <- unname(tapply(queryHits(ov), subjectHits(ov), min))
maxs <-  unname(tapply(queryHits(ov), subjectHits(ov), max))

anchorsfin_g <- makeGRangesFromDataFrame(setNames(data.frame(
  resFrags[mins, 1], resFrags[mins, 2], resFrags[maxs, 3]), c("seqnames", "start", "end")))

anchors_g <- reduce(anchorsfin_g)

# Filter anchors that were on terminal restriction fragments
anchors_gfilt <- anchors_g[width(anchors_g) < 50000] 
cat(paste0("Anchors removed due to excessive size: ",
           as.character(length(anchors_g) - length(anchors_gfilt))))
write.table(data.frame(anchors_gfilt)[,c(1,2,3)], row.names = FALSE, col.names = FALSE, quote = FALSE, 
            sep = "\t", file = "../data/hichipper_anchors.bed")

