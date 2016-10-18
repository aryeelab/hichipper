#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

cwd <- args[1]
outdir <- args[2]
samples <- unlist(strsplit(args[c(-1,-2)], split = " "))

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(foreach))
suppressMessages(library(gridExtra))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(readr))

pdf <- paste0(outdir, "/", outdir, "_hichipper-qcReport.pdf")

message("Processing: ", samples)
message("Saving QC Report to: ", pdf)

# Creates a dataframe of summary statistics from the individual sample log output files
readstats <- foreach(sample = samples, .combine="rbind") %do% {
  sfilename <- file.path(paste0(outdir, "/", sample, ".stat"))
  rs <- read.table(sfilename, header=FALSE, stringsAsFactors = FALSE, sep = "=")[c(1:7),]
  rs <- cbind(sample=sample, rs)
  colnames(rs) <- c("sample", "metric", "count")
  rs
}

# Get loop lengths for histogram
loop_pets <- foreach(sample = samples, .combine="rbind") %do% {
  sfilename <- file.path(paste0(outdir, "/",sample , ".intra.loop_counts.bedpe"))
  x <- suppressMessages(read_delim(sfilename, " ", col_names = FALSE))
  intra <- x[,1]==x[,4]
  x <- x[intra,]
  loop_length <- rep(as.numeric(((x[, 5]+ x[, 6])/2 - (x[, 2] + x[, 3])/2)[[1]]), as.numeric((x[, 8])[[1]]))
  data.frame(sample = sample, loop_length = pmax(0, loop_length))
}

# Add long-range loop counts to the read stats dataframe.
grouped <- group_by(loop_pets, sample)
summ <- summarise(grouped, long_range=sum(loop_length>=5000))
metrics <- c("Total_PETs", "Mapped_PETs_q30", "Mapped_unique_PETs_q30", "Mapped_unique_intrachromosomal_q30", "Mapped_unique_intrachromosomal_q30_5kb", "Anchor_mapped_PETs_5kb", "Loop_PETs")
readstats$metric <- factor(readstats$metric, levels=metrics)

# Plot read stats
p <- ggplot(readstats, aes(x = sample, y = count, fill=metric))  + 
  geom_bar(stat='identity', position=position_dodge()) + theme_bw() + 
  ggtitle("Counts") + xlab("") + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
p_metrics <- p + scale_y_continuous(labels=comma)

# Plot PET anchor separation distribution
p_hist <- ggplot(loop_pets, aes(loop_length)) + geom_histogram() + scale_x_log10(labels=comma, breaks=10^(3:9)) +
  facet_wrap(~sample, ncol=1, scales = "free_y") + theme_bw()

tab <- suppressMessages(acast(readstats, metric~sample, sum))
tab_percent <- 100*sweep(tab, 2, tab["Total_PETs",], FUN="/")
tab_summary <- rbind(tab[6,]/tab[5,] * 100, tab_percent[c(7,5,2),])
rownames(tab_summary) <- c("% Reads in anchors", "% PETs used in Filtered Loops", "% Long Range Interactions", "% Mapped")

# Output graphics and table to PDF file
pdf(pdf, height=8.5, width=11, onefile = TRUE)
grid.table(format(tab_summary, digits=2, nsmall=2))
p_metrics
plot.new()
grid.table(format(acast(readstats, metric~sample, sum), big.mark=","))
plot.new()
grid.table(format(tab_percent, digits=2, nsmall=2))
p_hist
dev.off()
