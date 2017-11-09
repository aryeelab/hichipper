#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

cwd <- args[1]
outdir <- args[2]
samples <- unlist(strsplit(args[c(-1,-2)], split = " "))

suppressMessages(suppressWarnings(library(diffloop)))
suppressMessages(suppressWarnings(library(foreach)))

message("Processing: ", samples)

# Do everything per-sample
tmp <- foreach(sample = samples, .combine="rbind") %do% {
  s <- loopsMake(outdir, samples = paste0(sample, ".filt.intra"), mergegap = 0)
  sampleNames(s) <- sample
  mango_filt_df <- summary(mangoCorrection(s))
  df_filt <- mango_filt_df[,c(1,2,3,4,5,6,7,9)]
  s <- rmchr(s)
  col <- mango_filt_df[,9] < 0.01
  col[col] <- "none"
  col[col=="FALSE"] <- "ns"
  s@rowData$loop.type <- col
  saveRDS(s, file = paste0(outdir, "/", sample, "-HiChIP.rds"))
  write.table(df_filt, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE,
              file = paste0(outdir, "/", sample, ".interactions.all.mango"))
}
