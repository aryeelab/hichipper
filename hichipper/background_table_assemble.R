#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

cwd <- args[1]
outdir <- args[2]
sample <- args[3]

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))

message("Assembling table for: ", sample)


hc_names <- c("chr_1", "start_1", "end_1", "chr_2", "start_2", "end_2", "dot", "count")

make_loop_id <- function(df){
  df$id <- paste0(df$chr_1, "_", df$start_1, "_", df$end_1, "_",
                  df$chr_2, "_", df$start_2, "_", df$end_2)
  return(df)
}

# Import 
foreground <- fread(paste0(outdir, "/", sample, ".filt.intra.loop_counts.bedpe"), col.names = hc_names) %>% make_loop_id %>% data.frame
background <- fread(paste0(outdir, "/background_", sample, ".filt.intra.loop_counts.bedpe"), col.names = hc_names) %>% make_loop_id %>% data.frame

# Merge into one dataframe
merge1 <- full_join(foreground, background, by = "id")

# Subset to useful columns
final_df <- merge1[,c("chr_1.x", "start_1.x", "end_1.x", "chr_2.x", "start_2.x", "end_2.x", "count.x", "count.y")]
colnames(final_df) <- c("chr_1", "start_1", "end_1", "chr_2", "start_2", "end_2",paste0(sample,"_count"), "control_count")
final_df[is.na(final_df)] <- 0

# Do a CPM normalization based on total reads (from Caleb's HiC-Pro directory)
#Foxp3_total_pets <- 138755223
#igg_total_pets <- 43818576
#final_df$Foxp3_cpm <- final_df$Foxp3_count/Foxp3_total_pets * 1000000
#final_df$IgG_cpm <- final_df$IgG_count/igg_total_pets * 1000000
#final_df$log2FoldChange <- log2((final_df$Foxp3_cpm + 0.125)/(final_df$IgG_cpm + 0.125))

write.table(final_df, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE,
            file = paste0(outdir, "/", sample, ".backgroundCompare.tsv"))
