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
options(scipen=999)

pdf <- paste0(outdir, "/samples_hichipper-qcReport.pdf")

message("Processing: ", samples)
message("Saving QC Report to: ", pdf)
min_length <- 0
max_length <- 0

# Creates a dataframe of summary statistics from the individual sample log output files
readstats <- foreach(sample = samples, .combine="rbind") %do% {
  sfilename <- file.path(paste0(outdir, "/", sample, ".stat"))
  rs <- read.table(sfilename, header=FALSE, stringsAsFactors = FALSE, sep = "=")
  min_length <- rs[10,2]
  max_length <- rs[11,2]
  rso <- cbind(sample=sample, rs[1:9,])
  colnames(rso) <- c("sample", "metric", "count")
  rso
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

np <- paste0("Mapped_unique_intra_q30_anchor_", round(min_length/1000), "KB-",  round(max_length/1000000), "Mb")
metrics <- c("Total_PETs", "Mapped_PETs_q30", "Mapped_unique_PETs_q30", "Mapped_unique_intrachromosomal_q30", paste0("Mapped_unique_intrachromosomal_q30_>", round(min_length/1000), "KB"),
             np, paste0("Mapped_unique_intra_q30_anchor<", round(min_length/1000), "KB"),  paste0("Mapped_unique_intra_q30_anchor>",round(max_length/1000000), "Mb"), "Number_Peaks")
readstats$metric <- factor(rep(metrics, length(samples)), levels=metrics)

# Plot read stats
cd <- c("#B2182B", "#EF8A62", "#FDDBC7", "#D1E5F0", "#67A9CF", "#2166AC", "#A3A1A1", "gray84","#F7A128")
p_metrics <- ggplot(readstats, aes(x = sample, y = count, fill=metric))  + 
  geom_bar(stat='identity', position=position_dodge()) + theme_bw() + 
  ggtitle("Counts") + xlab("") + scale_fill_manual(values=cd) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))  + scale_y_continuous(labels=comma)

# Plot PET anchor separation distribution
p_hist <- ggplot(loop_pets, aes(loop_length)) + geom_histogram() + scale_x_log10(labels=comma, breaks=10^(3:9)) +
  facet_wrap(~sample, ncol=1, scales = "free_y") + theme_bw()

tab <- suppressMessages(acast(readstats, metric~sample, sum))
tab_percent <- 100*sweep(tab, 2, tab["Total_PETs",], FUN="/")
tab_summary <- rbind(format(tab_percent[c(6,5,2),,drop = FALSE], digits = 2, nsmall = 2), 
                     as.character(tab[1,,drop = FALSE]))
rownames(tab_summary) <- c(paste0("% in Loops"),
                           "% Long Range Interaction",
                           "% Mapped",
                           "Total PETs")
tab_summary <- tab_summary[4:1,]

intraSum <- colSums(tab[c(6,7,8),])
tab_summary2 <- rbind(intraSum,
                      format(tab[6,]/intraSum * 100, digits = 2, nsmall = 2),
                      format(tab[7,]/intraSum * 100, digits = 2, nsmall = 2),
                      format(tab[8,]/intraSum * 100, digits = 2, nsmall = 2))
r <-  paste0(round(min_length/1000), "KB-",  round(max_length/1000000), "Mb")
rownames(tab_summary2) <- c("Intra, anchor-mapped PETs",
                           paste0("% Long Reads in ",r),
                           paste0("% Long Reads < ", round(min_length/1000), "KB"),
                           paste0("% Long Reads > ", round(max_length/1000000), "Mb"))


# Output graphics and table to PDF file
pdf(pdf, height=8.5, width=11, onefile = TRUE)
tt <- ttheme_default(core=list(fg_params=list(hjust=0, x=0.1)))
grid.arrange(
  tableGrob(tab_summary, theme=tt),
  tableGrob(tab_summary2, theme=tt),
  nrow=2)
p_metrics
plot.new()
grid.table(format(acast(readstats, metric~sample, sum), big.mark=","))
plot.new()
grid.table(format(tab_percent[-9,], digits=2, nsmall=2))
p_hist
dev.off()
