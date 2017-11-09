#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
options(warn=-1)

scriptdir <- args[1]
outdir <- args[2]
cwd <- args[3]
version <- args[4]
samples <- unlist(strsplit(args[seq(-1,-4)], split = " "))

cat(getwd())
rmarkdown::render(paste0(outdir,"/qcReport_make.Rmd"), params = list(
  scriptdir = scriptdir, 
  outdir = outdir, 
  cwd = cwd, 
  samples = samples,
  version = version
), quiet = TRUE, output_file=paste0(outdir, ".hichipper.qcreport.html"))

